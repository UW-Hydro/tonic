#!/usr/local/bin/python

import grid_params

gridFile = "../../../simon/columbia_elev.nc"
soilFile = "/net/juniper/raid/homefc/maca.vic/hb2860_parameter_files/soil/soil.columbia"
snowFile = "/net/juniper/raid/homefc/maca.vic/hb2860_parameter_files/snow/snowbands.0625.PNW.052108"
vegFile = "/net/juniper/raid/homefc/maca.vic/hb2860_parameter_files/veg/veg_param"
vegLFile = "/net/juniper/raid/homefc/maca.vic/hb2860_parameter_files/veg/veglib.LDAS"
ncFile = "test_1"

grid_params.make_grid(gridFile, soilFile, snowFile,vegFile,vegLFile,ncFile)
