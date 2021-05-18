from stompy.grid import multi_ugrid, unstructured_grid
import glob, os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani

# run name
run = 'sb_rv4_20190601'
# model output directory
out_dir = f"runs\\{run}\\DFM_OUTPUT_flowfm"

paths = glob.glob(os.path.join(out_dir, '*_map.nc'))

mg = multi_ugrid.MultiUgrid(paths, cleanup_dfm=True)

# subset time to animate
#for i in range(16):
#    mg.dss[i] = mg.dss[i].sel(time=slice('2019-07-15', '2019-08-15'))


# plotting final time
"""
scal=mg['mesh2d_NO3'].isel(time=-1).values
plt.figure()
ccoll=mg.grid.plot_cells(values=scal,cmap='jet')
plt.colorbar(ccoll)
plt.axis('equal')
plt.axis('off')
plt.show()

age_conc=mg['mesh2d_NO3'].isel(time=-1).values
conc=mg['mesh2d_ZNit'].isel(time=-1).values

age=age_conc/conc

plt.figure()
valid=conc>0.0001
age[~valid]=0
ccoll=mg.grid.plot_cells(values=age,cmap='jet',mask=valid)
mg.grid.plot_edges(color='k',lw=0.4,alpha=0.1)
#ccoll.set_clim([0,0.05])
plt.colorbar(ccoll)
plt.axis('equal')
plt.show()
"""

####################################
age_conc = mg['mesh2d_NO3']
conc = mg['mesh2d_ZNit']

frames = len(age_conc.sel().values)

fig, ax = plt.subplots()
cb = None

def animate(i):
    ax.clear()

    ac = age_conc.isel(time=i).values
    c = conc.isel(time=i).values

    t = age_conc.isel(time=i).sub_vars[0].time.values  # faffy way to get time...
    print(f'animating t = {t}...')

    age = ac / c
    valid = c > 0.0001
    age[~valid] = 0

    ccoll = mg.grid.plot_cells(values=age, cmap='jet', mask=valid)
    ccoll.set_clim(vmin=0, vmax=i/24)  # set colorbar limits
    ax.add_collection(ccoll)
    lcoll = mg.grid.plot_edges(color='k', lw=0.4, alpha=0.1)
    # ax.add_collection(lcoll)
    plt.axis('equal')
    ax.set(title=f'Time: {t}')
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)

    global cb
    if cb is not None:
        cb.remove()
    cb = plt.colorbar(ccoll, label='Age [days]')
    return ax,

animator = ani.FuncAnimation(fig, animate, interval=1, frames=frames)
animator.save(f'{run}_age.mp4', fps=5, dpi=200)
