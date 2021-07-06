#!/bin/sh

# to limit to region near Calhoun Cut 
# -b 605000 606000 4234000 4236000

python -m stompy.spatial.generate_dem -s ~/src/sfe_dem -o tiles_calhoun -r 2.0 -t 1000 -m -p ~/data/bathy_dwr/gtiff -p /opt/mirrors/ucd-X/Arc_Hydro/CSC_Project/MODELING/1_Hydro_Model_Files/Geometry/Bathymetry_Tiles/NDelta_hydro_2m_v4 -p . -d 2014-04-01 -g ugrid:../dflowfm/CacheSloughComplex_v111-edit19fix-bathy.nc
