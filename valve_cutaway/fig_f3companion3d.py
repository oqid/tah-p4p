import sys; sys.path.insert(0,'.')
from fig_film import *
phases=[16,50,83]
fig=plt.figure(figsize=(rs.WIDTH,4.5))
rows=[crop_common([f"plane_stress_oblique_plane_{k}.png" for k in phases]+[f"plane_dose_oblique_plane_{k}.png" for k in phases])]
ims=rows[0]
asp=ims[0].shape[0]/ims[0].shape[1]
L,R=0.005,0.86; pw=(R-L)/3; ph=pw*rs.WIDTH*asp/4.5
for r in range(2):
    for j,k in enumerate(phases):
        ax=fig.add_axes([L+j*pw,0.95-(r+1)*ph-r*0.02,pw*0.98,ph]); ax.imshow(ims[r*3+j]); ax.axis('off')
        if r==0: ax.set_title(f"{PH[k]:g}% systole · valve {ANG[k]:g}°",fontsize=7.2,pad=2)
y1=0.95-ph; y2=0.95-2*ph-0.02
c1=fig.colorbar(ScalarMappable(LogNorm(1e-2,1e2),'magma'),cax=fig.add_axes([0.885,y1+0.12*ph,0.018,0.76*ph]),extend='both')
c1.set_label("Local scalar stress [Pa]",fontsize=7)
c2=fig.colorbar(ScalarMappable(LogNorm(1e-6,1e-2),'magma_r'),cax=fig.add_axes([0.885,y2+0.12*ph,0.018,0.76*ph]),extend='min')
c2.set_label("Share of linear HI3 dose\nper mm² of slab",fontsize=7)
for c in (c1,c2): c.ax.tick_params(labelsize=6); c.solids.set_alpha(0.72)
for ext in ("png","pdf"): fig.savefig(f"F3_companion_section3d.{ext}",dpi=300)
