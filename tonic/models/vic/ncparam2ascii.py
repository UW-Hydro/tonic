#!/usr/bin/env python
"""
Script for writing VIC ascii format parameter files from netcdf format.
"""
description = 'Convert netCDF format VIC parameters to classic ASCII format'
help = description


import numpy as np
from scipy.spatial import cKDTree
from . import grid_params
from tonic.io import read_netcdf

FILL_VALUE = -9999


# -------------------------------------------------------------------- #
def _run(args):

    subset(args.nc_params, upleft=args.upleft, lowright=args.lowright,
           outfiles=args.outfiles, soil_file=args.soil_file,
           snow_file=args.snow_file, veg_file=args.veg_file,
           project=args.project, nijssen2arno=args.nijssen2arno)
    return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def subset(param_file, upleft=False, lowright=False, outfiles=1,
           soil_file=False, snow_file=False,
           veg_file=False, project=None,
           nijssen2arno=False):

    data, attributes = read_netcdf(param_file)

    if nijssen2arno:
        import NIJSSEN2001_to_ARNO
        data = NIJSSEN2001_to_ARNO.convert(data)

    if project:
        print('Project Configuration {0}'.format(project))
        if project == 'RASM':
            outfiles = 1
            cells, yinds, xinds = find_gridcells(data['mask'])
            rasm_soil(data, soil_file)
        else:
            raise ValueError('Unknown project configuration')
        return

    else:
        cells, yinds, xinds = find_gridcells(data['mask'])

    # write snow and veg files
    if veg_file:
        rootzones = data['root_depth'].shape[1]
        veg(data, xinds, yinds, veg_file, rootzones=rootzones, global_lai=True)
    if snow_file:
        snow(data, xinds, yinds, snow_file)

    if (upleft and lowright):
        inds = ((yinds < upleft[0]) and
                (yinds > lowright[0]) and
                (xinds < lowright[1]) and
                (xinds > upleft[1]))
        yinds = yinds[inds]
        xinds = xinds[inds]

    filesize = np.ceil(cells / outfiles)

    for i in range(outfiles):
        start = i * filesize
        end = i * filesize + filesize
        if end > cells:
            end = cells
        if outfiles > 1:
            out_file = '{0}_{1}.txt'.format(soil_file,
                                            str(i).zfill(len(str(outfiles))))
        else:
            out_file = '{0}.txt'.format(soil_file)
        soil(data, xinds[start:end], yinds[start:end], out_file)

    return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def rasm_soil(data, soil_file):
    """Write VIC formatted soil parameter file"""

    message = """
*** ---------------------------------------------------------------------- ***
Notes about RASM soil parameter file generations:
 - To fill in missing grid cells 'mask' variable must be the same as the
   domain file mask.
 - Inactive grid cells will have a dummy line printed for all variables except
   the lons/lats.
 - Any grid cells with nans will be copied from the previous line without nans
   or FILL_VALUEs.
*** --------------------------------------------------------------------- ***\n
    """

    print(message)

    # ---------------------------------------------------------------- #
    c = grid_params.cols(nlayers=3)
    f = grid_params.format(nlayers=3)

    c.soil_param['Nveg'] = np.array([0])
    f.soil_param['Nveg'] = '%1i'

    numcells = data['mask'].size

    arrayshape = (numcells, 1 + np.max([np.max(cols) for v, cols in
                                        c.soil_param.iteritems()]))
    soil_params = np.empty(arrayshape)
    dtypes = ['%1i'] * arrayshape[1]
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # find nearest real grid cell for all grid cells where frozen soil mask is
    # active
    # This is needed because the RASM mask is often different than existing
    # datasets
    print('Finding/filling nearest neighbors for all variables based on '
          'fs_active mask')

    # real mask (from input dataset)
    ry, rx = np.nonzero(data['fs_active'])

    # fill mask (we will fill all of these grid cells with their nearest real
    # value)
    my_mask = np.zeros(data['fs_active'].shape, dtype=int)
    my_mask[ry, rx] = 1

    fy, fx = np.nonzero(my_mask == 0)

    # Find nearest real grid cell
    combined = np.dstack(([data['yc'][ry, rx], data['xc'][ry, rx]]))[0]
    points = list(np.vstack((data['yc'][fy, fx],
                             data['xc'][fy, fx])).transpose())
    mytree = cKDTree(combined)
    dist, inds = mytree.query(points, k=1)

    # loop over all variables and fill in values
    for var in c.soil_param:
        if var not in ['run_cell', 'grid_cell', 'lats', 'lons', 'fs_active',
                       'Nveg']:
            # unmask the variable
            data[var] = np.array(data[var])
            # fill with nearest data
            if data[var].ndim == 2:
                data[var][fy, fx] = data[var][ry[inds], rx[inds]]
            elif data[var].ndim == 3:
                data[var][:, fy, fx] = data[var][:, ry[inds], rx[inds]]
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # This is the domain mask
    yinds, xinds = np.nonzero(data['mask'] == 1)
    ymask, xmask = np.nonzero(data['mask'] == 0)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Fix problematic avg_T values
    print('Finding/filling nearest neighbors for avg_T\n')

    # real mask (from input dataset)
    ry, rx = np.nonzero((data['avg_T'] > -50) & (data['avg_T'] < 99))

    # fill mask (we will fill all of these grid cells with their nearest
    # real value)
    my_mask = np.zeros(data['avg_T'].shape)
    my_mask[ry, rx] = 1

    fy, fx = np.nonzero(my_mask == 0)

    # Find nearest real grid cell
    if len(fy) > 0:
        combined = np.dstack(([data['yc'][ry, rx], data['xc'][ry, rx]]))[0]
        points = list(np.vstack((data['yc'][fy, fx],
                                 data['xc'][fy, fx])).transpose())
        mytree = cKDTree(combined)
        dist, inds = mytree.query(points, k=1)

        data['avg_T'] = np.array(data['avg_T'])
        data['avg_T'][fy, fx] = data['avg_T'][ry[inds], rx[inds]]
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # For rasm, all cols are shifted one to right to make room for nveg in
    # col 0
    i = -1
    for var, cols in c.soil_param.iteritems():
        for col in cols:
            dtypes[col] = f.soil_param[var]

    for (y, x), maskval in np.ndenumerate(data['mask']):
        # advance the row count
        i += 1

        # put real data
        for var in c.soil_param:

            cols = c.soil_param[var]

            if data[var].ndim == 2:
                soil_params[i, cols] = data[var][y, x]
            elif data[var].ndim == 3:
                for j, col in enumerate(cols):
                    soil_params[i, col] = data[var][j, y, x]

        # write the grid cell number
        soil_params[i, 2] = i + 1
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Write phi_s
    soil_params[:, c.soil_param['phi_s']] = -999
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Write nveg
    soil_params[:, 0] = data['Nveg'][:, :].flatten()
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Set all grid cells to run
    soil_params[:, 1] = 1  # run
    soil_params[:, -1] = 1  # run with frozen soils on
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # rewrite the lons/lats columns
    soil_params[:, 3] = data['lats'].flatten()
    soil_params[:, 4] = data['lons'].flatten()

    assert soil_params[-1, 3] == data['lats'][y, x]
    assert soil_params[-1, 4] == data['lons'][y, x]
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # check for nans
    # sometimes RASM has land grid cells that the input file does not.
    # In this case, we will use the previous lnd grid cell
    if np.isnan(soil_params.sum()):
        bad_cells = []
        replacements = []
        for i, cell in enumerate(soil_params):
            if np.isnan(np.sum(cell)):
                bad_cells.append(i)
                j = i
                while True:
                    j -= 1
                    if (FILL_VALUE not in soil_params[j, :]) and \
                            not np.isnan(np.sum(soil_params[j, :])):
                        soil_params[i, 5:] = soil_params[j, 5:]
                        replacements.append(j)
                        break
        print('Fixed {0} bad cells'.format(len(bad_cells)))
        print('Example: {0}-->{1}'.format(bad_cells[0], replacements[0]))
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Print summary of variables
    for var, cols in c.soil_param.iteritems():
        print('{0: <12}--> min: {1:<09.3f}, max: {2:<09.3f}, mean:'
              ' {3:<09.3f}'.format(var,
                                   soil_params[:, cols].min(),
                                   soil_params[:, cols].max(),
                                   soil_params[:, cols].mean()))
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Finish up
    assert soil_params[-1, 3] == data['lats'][y, x]
    assert soil_params[-1, 4] == data['lons'][y, x]
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Write the file
    print('writing soil parameter file: {0}'.format(soil_file))

    np.savetxt(soil_file, soil_params, fmt=dtypes)

    print('done RASM with soil')
    # ---------------------------------------------------------------- #

    return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def soil(data, xinds, yinds, soil_file):
    """Write VIC formatted soil parameter file"""
    c = grid_params.cols(nlayers=3)
    f = grid_params.format(nlayers=3)

    arrayshape = (len(xinds), 1 + np.max([np.max(c.soil_param[var])
                                          for var in c.soil_param]))
    soil_params = np.zeros(arrayshape)
    dtypes = [0] * arrayshape[1]

    for var in c.soil_param:
        if data[var].ndim == 2:
            soil_params[:, c.soil_param[var]] = np.atleast_2d(
                data[var][yinds, xinds]).transpose()
        elif data[var].ndim == 3:
            soil_params[:, c.soil_param[var]] = np.atleast_2d(
                data[var][:, yinds, xinds]).transpose()
        for col in c.soil_param[var]:
            dtypes[col] = f.soil_param[var]

    np.savetxt(soil_file, soil_params, fmt=dtypes)

    print('finished writing soil parameter file: {0}'.format(soil_file))

    return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def snow(data, xinds, yinds, snow_file):
    """Write VIC formatted snowband parameter file"""
    try:
        snow_bands = data['AreaFract'].shape[0]
    except:
        snow_bands = 5

    c = grid_params.cols(snow_bands=snow_bands)
    f = grid_params.format(snow_bands=snow_bands)

    arrayshape = (len(xinds), 1 + np.max([np.max(c.snow_param[var])
                                          for var in c.snow_param]))
    snow_params = np.zeros(arrayshape)
    dtypes = [0] * arrayshape[1]

    for var in c.snow_param:
        if data[var].ndim == 2:
            snow_params[:, c.snow_param[var]] = np.atleast_2d(
                data[var][yinds, xinds]).transpose()
        elif data[var].ndim == 3:
            snow_params[:, c.snow_param[var]] = np.atleast_2d(
                data[var][:, yinds, xinds]).transpose()
        for col in c.snow_param[var]:
            dtypes[col] = f.snow_param[var]

    np.savetxt(snow_file, snow_params, fmt=dtypes)

    print('finished writing snow parameter file: {0}'.format(snow_file))

    return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def veg(data, xinds, yinds, veg_file, rootzones=3, global_lai=True):
    """Write VIC formatted veg parameter file"""

    print('writing veg parameter file: {0}'.format(veg_file))

    # counter for bad grid cells
    count = 0

    f = open(veg_file, 'w')

    for y, x in zip(yinds, xinds):
        gridcell = int(data['gridcell'][y, x])
        n_veg = int(data['Nveg'][y, x])
        cv = data['Cv'][:, y, x]
        veg_class = np.nonzero(cv)[0]

        if not len(veg_class) == n_veg:
            count += 1

        line1 = str(gridcell) + ' ' + str(n_veg) + '\n'
        f.write(line1)
        if n_veg > 0:
            for j in veg_class:
                line2 = [str(j + 1)]
                line2.append(str(cv[j]))
                for k in range(rootzones):
                    line2.append(str(data['root_depth'][j, k, y, x]))
                    line2.append(str(data['root_fract'][j, k, y, x]))
                line2.append('\n')
                f.write(' '.join(line2))
                if global_lai:
                    line3 = []
                    for m in range(12):
                        line3.append(str(data['LAI'][j, m, y, x]))
                    line3.append('\n')
                    f.write(' '.join(line3))
    f.close()

    print('{0} grid cells have unequal veg_classes'.format(count))
    print('finished writing veg parameter file: {0}'.format(veg_file))

    return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def veg_lib():
    raise
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def find_gridcells(mask):
    """ Find number of grid cells in mask"""
    cells = np.sum(mask)

    xinds, yinds = np.nonzero(mask > 0)

    return cells, xinds, yinds
# -------------------------------------------------------------------- #
