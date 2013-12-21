#!/usr/bin/env python

import os
import numpy as np
from netCDF4 import Dataset
import grid_params
import argparse

FILL_VALUE = -9999

def main():

    nc_params, soil_file, UL, LR, outfiles, snow_file, veg_file, project = process_command_line()

    subset(nc_params, UL=UL, LR=LR, outfiles=outfiles,
           soil_file=soil_file, snow_file=snow_file, veg_file=veg_file, project=project)

    return

def subset(paramNC, UL=False, LR=False, outfiles=1,
           soil_file=False, snow_file=False, veg_file=False, project=None):

    data, attributes = read_netcdf(paramNC)

    if project:
        print('Project Configuration %s' %project)
        if project == 'RASM':
            outfiles=1
            cells, yinds, xinds = find_gridcells(data['mask'])
            rasm_soil(data, soil_file)
        else:
            raise ValueError('Unknown project configuration')
        return

    else:
        cells, yinds, xinds = find_gridcells(data['mask'])

    if veg_file:
        veg(data, xinds, yinds, veg_file, rootzones=2, GLOBAL_LAI=True)


    if (UL and LR):
        inds = (yinds<UL[0]) & (yinds>LR[0]) & (xinds<LR[1]) & (xinds>UL[1])
        yinds = yinds[inds]
        xinds = xinds[inds]

    filesize = np.ceil(cells/outfiles)

    for i in xrange(outfiles):
        start = i*filesize
        end = i*filesize+filesize
        if end>cells:
            end = cells
        if outfiles>1:
            out_file = '{}_{}.txt'.format(soil_file, str(i).zfill(len(str(outfiles))))
        else:
            out_file = soil_file
        soil(data, xinds[start:end], yinds[start:end], out_file)
        print('finished {}'.format(out_file))
    return

def rasm_soil(data, soil_file):
    """Write VIC formatted soil parameter file"""
    c = grid_params.cols(nlayers=3)
    f = grid_params.format(nlayers=3)

    numcells = data['mask'].size

    arrayshape = (numcells, 2+np.max([np.max(c.soil_param[var]) for var in c.soil_param]))
    soil_params = np.empty(arrayshape)
    dtypes = ['%1i']*arrayshape[1]

    ymask, xmask = np.nonzero(data['mask']==0)

    # make a dummy line of fill values for unused gridcells
    dummy = np.zeros(arrayshape[1])+FILL_VALUE
    dummy[0] = 0
    dummy[c.soil_param['init_moist']+1] = 9999
    dummy[c.soil_param['resid_moist']+1] = 0
    dummy[c.soil_param['quartz']+1] = 0
    dummy[c.soil_param['depth']+1] = 9999
    dummy[c.soil_param['avg_T']+1] = 99
    dummy[c.soil_param['avg_T']+1] = 0
    dummy[c.soil_param['bulk_density']+1] = 9998
    dummy[c.soil_param['soil_density']+1] = 9999
    dummy[c.soil_param['fs_active']+1] = 1

    print('filling all ocean grid cells with dummy line')

    # For rasm, all cols are shifted one to right to make room for nveg in col 0
    i = -1
    for (y, x), maskval in np.ndenumerate(data['mask']):
        # advance the row count
        i += 1

        # fill with dummy (to be safe)
        soil_params[i, :] = dummy

        if maskval > 0:
            # put real data
            for var in c.soil_param:
                if data[var].ndim == 2:

                    soil_params[i, c.soil_param[var]+1] = data[var][y, x]
                elif data[var].ndim == 3:
                    cols = c.soil_param[var]
                    for j, col in enumerate(cols):
                        soil_params[i, col+1] = data[var][j, y, x]
                for col in c.soil_param[var]:
                    dtypes[col+1] = f.soil_param[var]

    # check for nans
    # sometimes RASM has land grid cells that the input file does not.
    # In this case, we will use the previous lnd grid cell
    for i, cell in enumerate(soil_params):
        if np.isnan(np.sum(cell)):
            print('fixing nan values in row {}'.format(i))
            j = i
            while True:
                j -= 1
                if FILL_VALUE not in soil_params[j, :]:
                    soil_params[i, :] = soil_params[j, :]
                    print('replacing with row {}'.format(j))
                    break

    # Write nveg
    soil_params[:, 0] = data['Nveg'][:, :].flatten()

    # Set all grid cells to run
    soil_params[:, 1] = 1

    # rewrite the lons/lats columns
    soil_params[:, 2] = data['gridcell'].flatten()
    soil_params[:, 3] = data['lats'].flatten()
    soil_params[:, 4] = data['lons'].flatten()

    assert soil_params[-1, 3] == data['lats'][y, x]
    assert soil_params[-1, 4] == data['lons'][y, x]

    print('writing soil parameter file: {}'.format(soil_file))

    np.savetxt(soil_file, soil_params, fmt=dtypes)

    print('done RASM with soil')

    return

def soil(data, xinds, yinds, soil_file):
    """Write VIC formatted soil parameter file"""
    c = grid_params.cols(nlayers=3)
    f = grid_params.format(nlayers=3)

    arrayshape = (len(xinds), 1+np.max([np.max(c.soil_param[var]) for var in c.soil_param]))
    soil_params = np.zeros(arrayshape)
    dtypes = [0]*arrayshape[1]

    for var in c.soil_param:
        if data[var].ndim == 2:
            soil_params[:, c.soil_param[var]] = np.atleast_2d(data[var][yinds, xinds]).transpose()
        elif data[var].ndim == 3:
            soil_params[:, c.soil_param[var]] = np.atleast_2d(data[var][:, yinds, xinds]).transpose()
        for col in c.soil_param[var]:
            dtypes[col] = f.soil_param[var]

    print('writing soil parameter file: {}'.format(soil_file))

    np.savetxt(soil_file, soil_params, fmt=dtypes)

    print('done with soil file')

    return

def snow():
    pass

def veg(data, xinds, yinds, veg_file, rootzones=2, GLOBAL_LAI=True):
    """Write VIC formatted veg parameter file"""

    print('writing veg parameter file: {}'.format(veg_file))

    # counter for bad grid cells
    count = 0

    f = open(veg_file, 'w')

    for i in xrange(len(xinds)):
        x,y = xinds[i], yinds[i]
        gridcell = int(data['gridcell'][y,x])
        Nveg = int(data['Nveg'][y,x])
        Cv = data['Cv'][:, y, x]
        veg_class = np.nonzero(Cv)[0]

        if not len(veg_class)==Nveg:
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

    print('{} grid cells have unequal veg_classes'.format(count))
    print('done with veg')

    return

def veg_lib():
    pass

def find_gridcells(mask):
    """ Find number of grid cells in mask"""
    cells = np.sum(mask)

    xinds,yinds = np.nonzero(mask>0)

    return cells, xinds, yinds

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
                        type=int,
                        help="Output soil param file prefix (default is same as nc_params)",
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

    args = parser.parse_args()

    if args.soil_prefix:
        soil_prefix = args.soil_prefix
    else:
        soil_prefix = os.path.splitext(args.nc_params)[0]+'_soil.txt'

    return args.nc_params, soil_prefix, args.upper_left_corner, \
           args.lower_right_corner, args.outfiles, args.snow_file, args.veg_file, args.project

##################################################################################
## Read netCDF Inputs
## Read data from input netCDF.
##################################################################################
def read_netcdf(ncFile, var_names=[], coords=False, verbose=True):
    """
    Read data from input netCDF. Will read all variables if none provided.
    Will also return all variable attributes.
    Both variables (data and attributes) are returned as dictionaries named by variable
    """

    f = Dataset(ncFile,'r')
    if var_names==[]:
        var_names = f.variables.keys()

    if verbose:
        print('Reading input data var_names:{} from file:'.format(var_names, ncFile))

    d={}
    a={}

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

if __name__ == "__main__":
    main()
