#!/usr/bin/env python
"""
Script for writing VIC ascii format parameter files from netcdf format.
"""
import os
import numpy as np
from netCDF4 import Dataset
import grid_params
import argparse
from scipy.spatial import cKDTree

FILL_VALUE = -9999


# -------------------------------------------------------------------- #
def main():

    nc_params, soil_file, UL, LR, outfiles, snow_file, \
        veg_file, project, NIJSSEN2ARNO = process_command_line()

    subset(nc_params, UL=UL, LR=LR, outfiles=outfiles,
           soil_file=soil_file, snow_file=snow_file,
           veg_file=veg_file, project=project,
           NIJSSEN2ARNO=NIJSSEN2ARNO)

    return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def subset(paramNC, UL=False, LR=False, outfiles=1,
           soil_file=False, snow_file=False,
           veg_file=False, project=None,
           NIJSSEN2ARNO=False):

    data, attributes = read_netcdf(paramNC)

    if NIJSSEN2ARNO:
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
        veg(data, xinds, yinds, veg_file, rootzones=rootzones, GLOBAL_LAI=True)
    if snow_file:
        snow(data, xinds, yinds, snow_file)

    if (UL and LR):
        inds = (yinds < UL[0]) & (yinds > LR[0]) & (xinds < LR[1]) \
            & (xinds > UL[1])
        yinds = yinds[inds]
        xinds = xinds[inds]

    filesize = np.ceil(cells/outfiles)

    for i in xrange(outfiles):
        start = i*filesize
        end = i*filesize+filesize
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

    numcells = data['mask'].size

    arrayshape = (numcells, 2+np.max([np.max(c.soil_param[var])
                                      for var in c.soil_param]))
    soil_params = np.empty(arrayshape)
    dtypes = ['%1i']*arrayshape[1]
    # ---------------------------------------------------------------- #

    # # ---------------------------------------------------------------- #
    # # make a dummy line of fill values for unused gridcells
    # dummy = np.zeros(arrayshape[1])+FILL_VALUE
    # dummy[0] = 0
    # dummy[c.soil_param['init_moist']+1] = 9999
    # dummy[c.soil_param['resid_moist']+1] = 0
    # dummy[c.soil_param['quartz']+1] = 0
    # dummy[c.soil_param['depth']+1] = 9999
    # dummy[c.soil_param['avg_T']+1] = 99
    # dummy[c.soil_param['avg_T']+1] = 0
    # dummy[c.soil_param['bulk_density']+1] = 9998
    # dummy[c.soil_param['soil_density']+1] = 9999
    # dummy[c.soil_param['fs_active']+1] = 1

    # print('filling all ocean grid cells with dummy line')
    # # ---------------------------------------------------------------- #

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
    print_flag = 0
    for (y, x), maskval in np.ndenumerate(data['mask']):
        # advance the row count
        i += 1

        # fill with dummy (to be safe)
        # soil_params[i, :] = dummy

        print_flag += 1
        # put real data
        for var in c.soil_param:

            cols = c.soil_param[var]

            if data[var].ndim == 2:
                if print_flag == 1:
                    for col in cols:
                        dtypes[col+1] = f.soil_param[var]

                    print('{0: <12}--> min: {1:<09.3f}, max: {2:<09.3f}, mean:'
                          ' {3:<09.3f}'.format(var,
                                               data[var][yinds, xinds].min(),
                                               data[var][yinds, xinds].max(),
                                               data[var][yinds, xinds].mean()))

                soil_params[i, c.soil_param[var]+1] = data[var][y, x]

            elif data[var].ndim == 3:
                if print_flag == 1:
                    for col in cols:
                        dtypes[col+1] = f.soil_param[var]

                    print('{0: <12}--> min: {1:<09.3f}, max: {2:<09.3f}, mean:'
                          ' {3:<09.3f}'.format(var,
                                               data[var][:, yinds, xinds].min(),
                                               data[var][:, yinds, xinds].max(),
                                               data[var][:, yinds, xinds].mean()))

                for j, col in enumerate(cols):
                    soil_params[i, col+1] = data[var][j, y, x]

                if var == 'phi_s':
                    soil_params[:, c.soil_param['phi_s']+1] = -999

        # write the grid cell number
        soil_params[i, 2] = i+1

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

    # # ---------------------------------------------------------------- #
    # # check for nans
    # # sometimes RASM has land grid cells that the input file does not.
    # # In this case, we will use the previous lnd grid cell
    # for i, cell in enumerate(soil_params):
    #     if np.isnan(np.sum(cell)):
    #         print(cell)
    #         print('fixing nan values in row {}'.format(i))
    #         j = i
    #         while True:
    #             j -= 1
    #             if (FILL_VALUE not in soil_params[j, :]):
    #                 soil_params[i, :] = soil_params[j, :]
    #                 print('replacing with row {}'.format(j))
    #                 break
    # # ---------------------------------------------------------------- #

    print('writing soil parameter file: {0}'.format(soil_file))

    np.savetxt(soil_file, soil_params, fmt=dtypes)

    print('done RASM with soil')

    return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def soil(data, xinds, yinds, soil_file):
    """Write VIC formatted soil parameter file"""
    c = grid_params.cols(nlayers=3)
    f = grid_params.format(nlayers=3)

    arrayshape = (len(xinds), 1+np.max([np.max(c.soil_param[var])
                                        for var in c.soil_param]))
    soil_params = np.zeros(arrayshape)
    dtypes = [0]*arrayshape[1]

    for var in c.soil_param:
        if data[var].ndim == 2:
            soil_params[:, c.soil_param[var]] = np.atleast_2d(data[var][yinds, xinds]).transpose()
        elif data[var].ndim == 3:
            soil_params[:, c.soil_param[var]] = np.atleast_2d(data[var][:, yinds, xinds]).transpose()
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

    arrayshape = (len(xinds), 1+np.max([np.max(c.snow_param[var])
                                        for var in c.snow_param]))
    snow_params = np.zeros(arrayshape)
    dtypes = [0]*arrayshape[1]

    for var in c.snow_param:
        if data[var].ndim == 2:
            snow_params[:, c.snow_param[var]] = np.atleast_2d(data[var][yinds, xinds]).transpose()
        elif data[var].ndim == 3:
            snow_params[:, c.snow_param[var]] = np.atleast_2d(data[var][:, yinds, xinds]).transpose()
        for col in c.snow_param[var]:
            dtypes[col] = f.snow_param[var]

    np.savetxt(snow_file, snow_params, fmt=dtypes)

    print('finished writing snow parameter file: {0}'.format(snow_file))

    return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def veg(data, xinds, yinds, veg_file, rootzones=3, GLOBAL_LAI=True):
    """Write VIC formatted veg parameter file"""

    print('writing veg parameter file: {0}'.format(veg_file))

    # counter for bad grid cells
    count = 0

    f = open(veg_file, 'w')

    for y, x in zip(yinds, xinds):
        gridcell = int(data['gridcell'][y, x])
        Nveg = int(data['Nveg'][y, x])
        Cv = data['Cv'][:, y, x]
        veg_class = np.nonzero(Cv)[0]

        if not len(veg_class) == Nveg:
            count += 1

        line1 = str(gridcell)+' '+str(Nveg)+'\n'
        f.write(line1)
        if Nveg > 0:
            for j in veg_class:
                line2 = [str(j+1)]
                line2.append(str(Cv[j]))
                for k in xrange(rootzones):
                    line2.append(str(data['root_depth'][j, k, y, x]))
                    line2.append(str(data['root_fract'][j, k, y, x]))
                line2.append('\n')
                f.write(' '.join(line2))
                if GLOBAL_LAI:
                    line3 = []
                    for m in xrange(12):
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


# -------------------------------------------------------------------- #
def process_command_line():
    """
    Process command line arguments.
    Optionally may include snow and vegitation parameter files.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("nc_params",
                        type=str,
                        help="Input netCDF VIC parameter file")
    parser.add_argument("--soil_prefix",
                        type=str,
                        help="Output soil param file prefix (default is same "
                             "as nc_params)",
                        default=False)
    parser.add_argument("--veg_prefix",
                        type=int,
                        help="Output veg param file prefix",
                        default=False)
    parser.add_argument("-UL", "--upper_left_corner",
                        type=int,
                        help="Upper left corner for subset",
                        default=False)
    parser.add_argument("-LR", "--lower_right_corner",
                        type=int,
                        help="Lower right corner for subset",
                        default=False)
    parser.add_argument("--outfiles",
                        type=int,
                        help="Number of outfiles",
                        default=1)
    parser.add_argument("--snow_file",
                        type=str,
                        help="Name of output snow file",
                        default=False)
    parser.add_argument("--veg_file",
                        type=str,
                        help="Name of output veg_file",
                        default=False)
    parser.add_argument("--project",
                        type=str,
                        help='Use project configuration options',
                        choices=['RASM'])
    parser.add_argument("--NIJSSEN2ARNO",
                        help='Convert soil parameters from NIJSSEN2001 format '
                             'to ARNO format',
                        action='store_true')

    args = parser.parse_args()

    if args.soil_prefix:
        soil_prefix = args.soil_prefix
    else:
        soil_prefix = os.path.splitext(args.nc_params)[0]

    return args.nc_params, soil_prefix, args.upper_left_corner, \
        args.lower_right_corner, args.outfiles, args.snow_file, \
        args.veg_file, args.project, args.NIJSSEN2ARNO
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Read netCDF Inputs
def read_netcdf(ncFile, var_names=[], coords=False, verbose=True):
    """
    Read data from input netCDF. Will read all variables if none provided.
    Will also return all variable attributes.
    Both variables (data and attributes) are returned as dictionaries named by
    variable
    """

    f = Dataset(ncFile, 'r')
    if var_names == []:
        var_names = f.variables.keys()

    if verbose:
        print('Reading input data var_names:{0} from file: '
              '{1}'.format(var_names, ncFile))

    d = {}
    a = {}

    if coords:
        if isinstance(var_names, str):
            d[var_names] = f.variables[var_names][coords]
            a[var_names] = f.variables[var_names].__dict__
        else:
            for var in var_names:
                d[var] = f.variables[var][coords]
                a[var] = f.variables[var].__dict__
    else:
        if isinstance(var_names, str):
            d[var_names] = f.variables[var_names][:]
            a[var_names] = f.variables[var_names].__dict__
        else:
            for var in var_names:
                d[var] = f.variables[var][:]
                a[var] = f.variables[var].__dict__
    f.close()

    return d, a
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
if __name__ == "__main__":
    main()
# -------------------------------------------------------------------- #
