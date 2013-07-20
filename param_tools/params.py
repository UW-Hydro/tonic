#!/usr/local/bin/python

import numpy as np

def soil(inFile,Nlayers = 3):
    """
    Load the entire soil file into a dictionary of numpy arrays.
    Also reorders data to match gridcell order of soil file.  
    """
    
    data = np.loadtxt(inFile)
    
    c = cols(Nlayers=Nlayers)
    
    soilDict = {}
    for var in c.soilParam:
        soilDict[var] = np.squeeze(data[:,c.soilParam[var]])

    return soilDict

def snow(snowFile,soilDict,Snow_bands = 5):
    """
    Load the entire snow file into a dictionary of numpy arrays.
    Also reorders data to match gridcell order of soil file.  
    """
    
    data = np.loadtxt(snowFile)
    
    c = cols(Snow_bands=Snow_bands)
    
    snowDict = {}
    for var in c.snowParam:
        snowDict[var] = data[:,c.snowParam[var]]

    target = soilDict['gridcel'].argsort()
    indexes = target[np.searchsorted(soilDict['gridcel'][target], snowDict['cellnum'])]

    for var in snowDict:
        snowDict[var] = np.squeeze(snowDict[var][indexes])

    return snowDict

def veg(vegFile,soilDict,maxRoots = 3,vegClasses = 11,
        cells = False,BLOWING_SNOW = False, LAIndex = False):
    """
    Read the vegetation file from vegFile.  Assumes max length for rootzones and vegclasses
    Also reorders data to match gridcell order of soil file.  
    """
    with open(vegFile) as f:
        lines = f.readlines()

    if not cells:
        cells = len(lines)
        
    gridcel = np.zeros(cells)
    Nveg =  np.zeros(cells)
    Cv = np.zeros((cells,vegClasses))
    root_depth = np.zeros((cells,vegClasses,maxRoots))
    root_fract = np.zeros((cells,vegClasses,maxRoots))
    if BLOWING_SNOW:
        sigma_slope =  np.zeros((cells,vegClasses))
        lag_one = np.zeros((cells,vegClasses))
        fetch = np.zeros((cells,vegClasses))
    if LAIndex:
        lfactor = 2
        LAI = np.zeros((cells,vegClasses,12))
    else:
        lfactor = 1
        
    row, cell = 0,0
    while row < len(lines):
        line = lines[row].strip('\n').split(' ')
        gridcel[cell], Nveg[cell] = np.array(line).astype(int)
        numrows = Nveg[cell]*lfactor+row
        row += 1
        
        while row < numrows:
            lines[row] = lines[row].strip()
            line = lines[row].strip('\n').split(' ')
            temp = np.array(line).astype(float)
            vind = int(temp[0])-1
            Cv[cell,vind] = temp[1]
            if not BLOWING_SNOW:
                rind = (len(temp)-2)/2
            else:
                rind = (len(temp)-5)/2
                sigma_slope[cell,vind,:rind] = temp[-3]
                lag_one[cell,vind,:rind] = temp[-2]
                fetch[cell,vind,:rind] = temp[-1]
            root_depth[cell,vind,:rind] = temp[2::2]
            root_fract[cell,vind,:rind] = temp[3::2]
            row += 1
            if LAIndex:
                lines[row] = lines[row].strip()
                line = lines[row].strip('\n').split(' ')
                LAI[cell,vind,:] = np.array(line).astype(float)
                row += 1
        cell += 1
    vegDict = {}
    vegDict['gridcell'] = gridcel[:cell]
    vegDict['Nveg'] = Nveg[:cell]
    vegDict['Cv'] = Cv[:cell,:]
    vegDict['root_depth'] = root_depth[:cell,:]
    vegDict['root_fract'] = root_fract[:cell,:]
    if BLOWING_SNOW:
        vegDict['sigma_slope'] = sigma_slope[:cell,:]
        vegDict['lag_one'] = lag_one[:cell,:]
        vegDict['fetch'] =fetch[:cell,:]
    if LAIndex:
        vegDict['LAI'] = LAI[:cell,:,:]

    #target = vegDict['gridcell'].argsort()
    #indexes = target[np.searchsorted(soilDict['gridcel'][target], vegDict['gridcell'])]

   # vegDict['Nveg'] = vegDict['Nveg'][inds]
    inds = []
    for sn in soilDict['gridcel']:
        inds.append(np.nonzero(vegDict['gridcell']==sn))

    for var in vegDict:
        vegDict[var] = np.squeeze(vegDict[var][inds])
    return vegDict

def veg_class(inFile,maxcols = 57,skiprows = 1):
    """
    Load the entire vegetation library file into a dictionary of numpy arrays.
    Also reorders data to match gridcell order of soil file.  
    """

    usecols = np.arange(maxcols)
    
    data = np.loadtxt(inFile,usecols = usecols,skiprows = skiprows)
    
    c = cols()
    
    veglDict = {}
    for var in c.vegLib:
        veglDict[var] = np.squeeze(data[:,c.vegLib[var]])

    return veglDict
        
class cols:
    def __init__(self,Nlayers=3,Snow_bands = 5):
        self.soilParam = {'run_cell':[0],
                           'gridcel':[1],
                           'lat':[2],
                           'lon':[3],
                           'infilt':[4],
                           'Ds':[5],
                           'Dsmax':[6],
                           'Ws':[7],
                           'c':[8],
                           'expt':np.arange(9,Nlayers+9),
                           'Ksat':np.arange(Nlayers+9,2*Nlayers+9),
                           'phi_s':np.arange(2*Nlayers+9,3*Nlayers+9),
                           'init_moist':np.arange(3*Nlayers+9,4*Nlayers+9),
                           'elev':[4*Nlayers+9],
                           'depth':np.arange(4*Nlayers+10,5*Nlayers+10),
                           'avg_T':[5*Nlayers+10],
                           'dp':[5*Nlayers+11],
                           'bubble':np.arange(5*Nlayers+12,6*Nlayers+12),
                           'quartz':np.arange(6*Nlayers+12,7*Nlayers+12),
                           'bulk_density':np.arange(7*Nlayers+12,8*Nlayers+12),
                           'soil_density':np.arange(8*Nlayers+12,9*Nlayers+12),
                           'off_gmt':[9*Nlayers+12],
                           'Wcr_FRACT':np.arange(9*Nlayers+13,10*Nlayers+13),
                           'Wpwp_FRACT':np.arange(10*Nlayers+13,11*Nlayers+13),
                           'rough':[11*Nlayers+13],
                           'snow_rough':[11*Nlayers+14],
                           'annual_prec':[11*Nlayers+15],
                           'resid_moist':np.arange(11*Nlayers+16,12*Nlayers+16),
                           'fs_active':[12*Nlayers+16]}

        self.snowParam = {'cellnum':[0],
                          'AreaFract':np.arange(1,Snow_bands+1),
                          'elevation':np.arange(Snow_bands+1,2*Snow_bands+1),
                          'Pfactor':np.arange(2*Snow_bands+1,3*Snow_bands+1)}
        
        self.vegLib = {'Veg_class':[0],
                       'lib_overstory':[1],
                       'lib_rarc':[2],
                       'lib_rmin':[3],
                       'lib_LAI':np.arange(4,16),
                       'lib_albedo':np.arange(16,28),
                       'lib_rough':np.arange(28,40),
                       'lib_displacement':np.arange(40,52),
                       'lib_wind_h':[52],
                       'lib_RGL':[53],
                       'lib_rad_atten':[54],
                       'lib_wind_atten':[55],
                       'lib_trunk_ratio':[56]}
                        #'comment': [57]}

class format:
    def __init__(self,Nlayers=3,Snow_bands = 5):
        self.soilParam = {'run_cell':'%1i',
                           'gridcel':'%1i',
                           'lat':'%1.3f',
                           'lon':'%1.3f',
                           'infilt':'%1.6f',
                           'Ds':'%1.6f',
                           'Dsmax':'%1.6f',
                           'Ws':'%1.6f',
                           'c':'%1.6f',
                           'expt':'%1.6f',
                           'Ksat':'%1.6f',
                           'phi_s':'%1.6f',
                           'init_moist':'%1.6f',
                           'elev':'%1.6f',
                           'depth':'%1.6f',
                           'avg_T':'%1.6f',
                           'dp':'%1.6f',
                           'bubble':'%1.6f',
                           'quartz':'%1.6f',
                           'bulk_density':'%1.6f',
                           'soil_density':'%1.6f',
                           'off_gmt':'%1.6f',
                           'Wcr_FRACT':'%1.6f',
                           'Wpwp_FRACT':'%1.6f',
                           'rough':'%1.6f',
                           'snow_rough':'%1.6f',
                           'annual_prec':'%1.6f',
                           'resid_moist':'%1.6f',
                           'fs_active':'%1i'}

        self.snowParam = {'cellnum':['%1i'],
                          'AreaFract':['%1.6f']*Snow_bands,
                          'elevation':['%1.6f']*Snow_bands,
                          'Pfactor':['%1.6f']*Snow_bands}
        
        self.vegLib = {'Veg_class':['%1i'],
                       'lib_overstory':['%1.6f'],
                       'lib_rarc':['%1.6f'],
                       'lib_rmin':['%1.6f'],
                       'lib_LAI':['%1.6f']*12,
                       'lib_albedo':['%1.6f']*12,
                       'lib_rough':['%1.6f']*12,
                       'lib_displacement':['%1.6f']*12,
                       'lib_wind_h':['%1.6f'],
                       'lib_RGL':['%1.6f'],
                       'lib_rad_atten':['%1.6f'],
                       'lib_wind_atten':['%1.6f'],
                       'lib_trunk_ratio':['%1.6f']}
                        #'comment': [57]}

class description:
    def __init__(self,ORGANIC_FRACT = False,SPATIAL_FROST = False,SPATIAL_SNOW = False,
                 EXCESS_ICE = False, JULY_TAVG_SUPPLIED = False, BLOWING_SNOW = False,GLOBAL_LAI = False):
        self.soilParam = {'run_cell':'1 = Run Grid Cell, 0 = Do Not Run',
                           'gridcel':'Grid cell number',
                           'lat':'Latitude of grid cell',
                           'lon':'Longitude of grid cell',
                           'infilt':'Variable infiltration curve parameter (binfilt)',
                           'Ds':'Fraction of Dsmax where non-linear baseflow begins',
                           'Dsmax':'Maximum velocity of baseflow',
                           'Ws':'Fraction of maximum soil moisture where non-linear baseflow occurs',
                           'c':'Exponent used in baseflow curve, normally set to 2',
                           'expt':'Exponent n (=3+2/lambda) in Campbells eqn for hydraulic conductivity, HBH 5.6 (where lambda = soil pore size distribution parameter).  Values should be > 3.0.',
                           'Ksat':'Saturated hydrologic conductivity',
                           'phi_s':'Soil moisture diffusion parameter',
                           'init_moist':'Initial layer moisture content',
                           'elev':'Average elevation of grid cell',
                           'depth':'Thickness of each soil moisture layer',
                           'avg_T':'Average soil temperature, used as the bottom boundary for soil heat flux solutions',
                           'dp':'Soil thermal damping depth (depth at which soil temperature remains constant through the year, ~4 m)',
                           'bubble':'Bubbling pressure of soil. Values should be > 0.0',
                           'quartz':'Quartz content of soil',
                           'bulk_density':'Bulk density of soil layer',
                           'soil_density':'Soil particle density, normally 2685 kg/m3',
                           'off_gmt':'Time zone offset from GMT. This parameter determines how VIC interprets sub-daily time steps relative to the model start date and time.',
                           'Wcr_FRACT':'Fractional soil moisture content at the critical point (~70% of field capacity) (fraction of maximum moisture)',
                           'Wpwp_FRACT':'Fractional soil moisture content at the wilting point (fraction of maximum moisture)',
                           'rough':'Surface roughness of bare soil',
                           'snow_rough':'Surface roughness of snowpack',
                           'annual_prec':'Average annual precipitation.',
                           'resid_moist':'Soil moisture layer residual moisture.',
                           'fs_active':'If set to 1, then frozen soil algorithm is activated for the grid cell. A 0 indicates that frozen soils are not computed even if soil temperatures fall below 0C.'}

        self.snowParam = {'cellnum':'Grid cell number (should match numbers assigned in soil parameter file)',
                           'AreaFract':'Fraction of grid cell covered by each elevation band. Sum of the fractions must equal 1.',
                           'elevation':'Mean (or median) elevation of elevation band. This is used to compute the change in air temperature from the grid cell mean elevation.',
                           'Pfactor':'Fraction of cell precipitation that falls on each elevation band. Total must equal 1. To ignore effects of elevation on precipitation, set these fractions equal to the area fractions.'}
        
        self.vegLib = {'Veg_class':'Vegetation class identification number (reference index for library table)',
                       'lib_overstory':'Flag to indicate whether or not the current vegetation type has an overstory (TRUE for overstory present [e.g. trees], FALSE for overstory not present [e.g. grass])',
                       'lib_rarc':'Architectural resistance of vegetation type (~2 s/m)',
                       'lib_rmin':'Minimum stomatal resistance of vegetation type (~100 s/m)',
                       'lib_LAI':'Leaf-area index of vegetation type',
                       'lib_albedo':'Shortwave albedo for vegetation type',
                       'lib_rough':'Vegetation roughness length (typically 0.123 * vegetation height)',
                       'lib_displacement':'Vegetation displacement height (typically 0.67 * vegetation height)',
                       'lib_wind_h':'Height at which wind speed is measured.',
                       'lib_RGL':'Minimum incoming shortwave radiation at which there will be transpiration. For trees this is about 30 W/m^2, for crops about 100 W/m^2.',
                       'lib_rad_atten':'Radiation attenuation factor. Normally set to 0.5, though may need to be adjusted for high latitudes.',
                       'lib_wind_atten':'Wind speed attenuation through the overstory. The default value has been 0.5.',
                       'lib_trunk_ratio':'Ratio of total tree height that is trunk (no branches). The default value has been 0.2.',
                       'lib_comment': 'Comment block for vegetation type. Model skips end of line so spaces are valid entrys.'}

        self.vegParam = {'gridcell':'Grid cell number',
                         'Nveg': 'Number of vegetation tiles in the grid cell',
                         'veg_class': 'Vegetation class identification number (reference index to vegetation library)',
                         'Cv':'Fraction of grid cell covered by vegetation tile',
                         'root_depth':'Root zone thickness (sum of depths is total depth of root penetration)',
                         'root_fract':'Fraction of root in the current root zone.',
                         'LAI':'Leaf Area Index, one per month'}

        # if ORGANIC_FRACT:
        #     self.soilParam['organic'] = 'Fraction of soil layer that is organic'
        #     self.soilParam['bul_dens_org'] = 'Bulk density of organic portion of soil'
        #     self.soilParam['soil_dens_org'] = 'Soil particle density of organic portion of soil, normally 1300 kg/m3 '

        # if SPATIAL_FROST:
        #     self.soilParam['frost_slope'] = 'Slope of uniform distribution of soil temperature'

        # if SPATIAL_SNOW:
        #     self.soilParam['max_snow_distrib_slope'] = 'Maximum slope of the snow depth distribution.'

        # if EXCESS_ICE:
        #     self.soilParam['initial_ice_content'] = 'Initial volumetric ice content in soil '

        # if JULY_TAVG_SUPPLIED:
        #     self.soilParam['July_Tavg'] = 'Average July air temperature, used for treeline computations'

        # if BLOWING_SNOW:
        #     self.vegParam['sigma_slope'] = 'Standard deviation of terrain slopes within vegetation tile'
        #     self.vegParam['lag_one'] = 'Standard deviation of terrain slopes within vegetation tile'
        #     self.vegParam['fetch'] = 'Standard deviation of terrain slopes within vegetation tile'

        # if GLOBAL_LAI:
        #     print 'global lai included'
        #     self.vegParam['LAI'] = 'Leaf Area Index, one per month'
                        
class units:
    def __init__(self,ORGANIC_FRACT = False,SPATIAL_FROST = False,SPATIAL_SNOW = False,
                 EXCESS_ICE = False, JULY_TAVG_SUPPLIED = False, BLOWING_SNOW = False,GLOBAL_LAI = False):
        self.soilParam = {'run_cell':'N/A',
                          'gridcel':'N/A',
                          'lat':'degrees',
                          'lon':'degrees',
                          'infilt':'mm/day',
                          'Ds':'fraction',
                          'Dsmax':'mm/day',
                          'Ws':'fraction',
                          'c':'N/A',
                          'expt':'N/A',
                          'Ksat':'mm/day',
                          'phi_s':'mm/mm',
                          'init_moist':'mm',
                          'elev':'m',
                          'depth':'m',
                          'avg_T':'C',
                          'dp':'m',
                          'bubble':'cm',
                          'quartz':'fraction',
                          'bulk_density':'kg/m3',
                          'soil_density':'kg/m3',
                          'off_gmt':'hours',
                          'Wcr_FRACT':'fraction',
                          'Wpwp_FRACT':'fraction',
                          'rough':'m',
                          'snow_rough':'m',
                          'annual_prec':'mm',
                          'resid_moist':'fraction',
                          'fs_active':'binary'}
            
        self.snowParam = {'cellnum':'N/A',
                          'AreaFract':'fraction',
                          'elevation':'m',
                          'Pfactor':'fraction'}
        
        self.vegLib = {'Veg_class':'N/A',
                       'lib_overstory':'N/A',
                       'lib_rarc':'s/m',
                       'lib_rmin':'s/m',
                       'lib_LAI':'N/A',
                       'lib_albedo':'fraction',
                       'lib_rough':'m',
                       'lib_displacement':'m',
                       'lib_wind_h':'m',
                       'lib_RGL':'W/m^2.',
                       'lib_rad_atten':'fraction',
                       'lib_wind_atten':'fraction',
                       'lib_trunk_ratio':'fraction',
                       'lib_comment': 'N/A'}

        self.vegParam = {'gridcell':'N/A',
                         'Nveg': 'N/A',
                         'veg_class': 'N/A',
                         'Cv':'fraction',
                         'root_depth':'m',
                         'root_fract':'fraction',
                         'LAI':'N/A'}
        # if ORGANIC_FRACT:
        #     self.soilParam['organic'] = 'fraction'
        #     self.soilParam['bul_dens_org'] = 'kg/m3'
        #     self.soilParam['soil_dens_org'] = 'kg/m3 '
            
        # if SPATIAL_FROST:
        #     self.soilParam['frost_slope'] = 'C'

        # if SPATIAL_SNOW:
        #     self.soilParam['max_snow_distrib_slope'] = 'm'

        # if EXCESS_ICE:
        #     self.soilParam['initial_ice_content'] = 'N/A'

        # if JULY_TAVG_SUPPLIED:
        #     self.soilParam['July_Tavg'] = 'C'
        
        # if BLOWING_SNOW:
        #     self.vegParam['sigma_slope'] = 'N/A'
        #     self.vegParam['lag_one'] = 'N/A'
        #     self.vegParam['fetch'] = 'm'

        #     if GLOBAL_LAI:
        #         self.vegParam['LAI'] = 'N/A'
            
