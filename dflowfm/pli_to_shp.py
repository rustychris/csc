import os
from stompy.spatial import wkb2shp
import stompy.model.delft.io as dio
from shapely import geometry
import glob

##

shp_dest='gis/model-features.shp'

names=[]
geoms=[]

for fn in glob.glob('*.pli'):
    feats=dio.read_pli(fn)
    for feat in feats: # generally just one per file
        names.append(feat[0])
        geoms.append(geometry.LineString(feat[1]))

wkb2shp.wkb2shp("gis/model-features.shp",geoms,fields=dict(names=names))

# AmericanRiver.pli
# Barker_Pumping_Plant.pli
# DXC.pli
# FlowFMcrs.pli
# Georgiana.pli
# SacramentoRiver.pli
# SRV.pli

##

# Same but for point features
import pandas as pd
df=pd.read_csv("ND_stations.xyn",header=None,names=['x','y','name'],
               sep=r'\s+')

names=[r['name'] for ri,r in df.iterrows()]
for i in range(len(names)):
    if names[i][0] in ['"',"'"] and names[i][0]==names[i][-1]:
        names[i]=names[i].replace(names[i][0],"")

wkb2shp.wkb2shp("gis/point-features.shp",
                [geometry.Point(r.x,r.y) for ri,r in df.iterrows()],
                fields=dict(name=names,
                            type=['obs' for ri,r in df.iterrows()]),
                overwrite=True)


