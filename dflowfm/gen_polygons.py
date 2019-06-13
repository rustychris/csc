# write .pol files for the polygon features in gis/poly-features

from stompy.spatial import wkb2shp
import os
import numpy as np

shp_fn=os.path.join(os.path.dirname(__file__),
                    "gis/poly-features.shp")
polys = wkb2shp.shp2geom(shp_fn)

for feat in polys:
    print(feat['name'])
    with open(feat['name']+".pol",'wt') as fp:
        fp.write(f"* Transcribed from {shp_fn}\n")
        fp.write(f"{feat['name']}\n")
        pnts=np.array( feat['geom'].exterior )
        fp.write(f"{len(pnts)} 2\n")
        for xy in pnts:
            fp.write("%.3f %.3f\n"%(xy[0],xy[1]))

    