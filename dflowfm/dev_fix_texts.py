import run_cal_plotter

import stompy.model.delft.dflow_model as dfm
from stompy.plot import plot_wkb
from stompy import utils
import stompy.grid.unstructured_grid as ugrid

mdupath="runs/grid100_21/flowfm.mdu"
model=dfm.DFlowModel.load(mdupath)
    
his_file=model.his_output()
his_dir=os.path.dirname(his_file)
fig_version="20181116"
plot_dir = os.path.join(his_dir,'figs-%s'%fig_version)
os.path.exists(plot_dir) or os.makedirs(plot_dir)

run_cal_plotter.plot_dir=plot_dir

mr = run_cal_plotter.hcp.DflowfmModelResults([his_file],trim_spinup_days=1.0)

# A bit awkward to mutate stations...
for station in run_cal_plotter.stations:
    station[1]['plot_dir']=plot_dir
    
fig=run_cal_plotter.make_summary_map(model,mr)
    
##

# shift texts to avoid overlaps
ax=fig.axes[0]


txt=texts[0]

txt.get_position() # (x,y)


def adjust_text_position(ax,max_iter=200):
    texts=ax.texts

    bboxes=[]
    for txt in texts:
        ext=txt.get_window_extent()
        bboxes.append( [ [ext.xmin,ext.ymin],
                         [ext.xmax,ext.ymax] ] )
    bboxes=np.array(bboxes)

    # each iteration move overlapping texts by about this
    # much
    dx=2.0
    # come up with pixel offsets for each label.
    pix_offsets=np.zeros( (len(texts),1,2),np.float64 )
    for _ in range(max_iter):
        changed=False

        new_bboxes=bboxes+pix_offsets
        for a in range(len(texts)):
            for b in range(a+1,len(texts)):
                # is there any force between these two labels?
                # check overlap

                int_min=np.maximum(new_bboxes[a,0,:], new_bboxes[b,0,:])
                int_max=np.minimum(new_bboxes[a,1,:], new_bboxes[b,1,:])

                if np.all(int_min<int_max):
                    #print("Collision %s - %s"%(texts[a].get_text(),
                    #                           texts[b].get_text()))

                    # This could probably be faster and less verbose.
                    # separate axis is taken from the overlapping region
                    # and direction
                    
                    # choose the direction that most quickly eliminates the overlap
                    # area. could also just choose the least overlapping direction
                    opt=utils.to_unit( np.array( [int_max[1]-int_min[1],
                                                  int_max[0]-int_min[0]]) )
                    
                    ab=new_bboxes[b].mean(axis=0) - new_bboxes[a].mean(axis=0)
                    if np.dot(opt,ab)<0:
                        opt*=-1

                    pix_offsets[a,0,:] -= dx*opt/2
                    pix_offsets[b,0,:] += dx*opt/2

                    changed=True
        if not changed:
            break

    # Update positions of the texts:
    deltas=np.zeros((len(texts),2),np.float64)
    for i in range(len(texts)):
        txt=texts[i]
        xform=txt.get_transform()
        ixform=xform.inverted()
        p=bboxes[i,0,:]
        pos0=xi.transform_point(p)
        pos_new=xi.transform_point(p+pix_offsets[i,0,:])
        deltas[i]=delta=pos_new-pos0
        txt.set_position( np.array(txt.get_position()) + delta)
    return deltas

