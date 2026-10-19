import sys; sys.path.insert(0,'.')
from fig_section import *
from PIL import Image
from render3d import CMAP as CMAP3
from matplotlib.colors import LogNorm
from matplotlib.cm import ScalarMappable
ims={k:np.asarray(Image.open(f"r3d_{k}.png").convert('RGB')) for k in PH}
# common crop box over all renders
nonwhite=np.any(np.stack([np.any(im<245,axis=2) for im in ims.values()]),axis=0)
ys,xs=np.where(nonwhite); y0,y1,x0,x1=ys.min(),ys.max(),xs.min(),xs.max()
fig=plt.figure(figsize=(rs.WIDTH,3.7))
gs=fig.add_gridspec(2,3,hspace=0.10,wspace=0.02,left=0.005,right=0.995,top=0.95,bottom=0.01)
for i,k in enumerate([16,33,50,66,83]):
    ax=fig.add_subplot(gs[i//3,i%3]); ax.imshow(ims[k][y0:y1,x0:x1]); ax.axis('off')
    ax.set_title(TITLE(i,k),fontsize=6.8,loc='left',fontweight='normal',pad=2)
cax=fig.add_axes([0.70,0.36,0.27,0.016])
cb=fig.colorbar(ScalarMappable(LogNorm(0.5,100),CMAP3),cax=cax,orientation='horizontal',extend='both')
cb.set_ticks([0.5,1,10,100]); cb.set_ticklabels(['0.5','1','10','100']); cb.ax.tick_params(labelsize=6,length=2); cb.ax.set_title("Scalar shear stress on path (Pa)",fontsize=6.5,loc='left',pad=3)
fig.text(0.70,0.31,"Rear half of the fluid domain, cut on the plane\nthrough both ports; 700 random inlet-seeded\npaths per phase, clipped to the valve region.",fontsize=6,color=rs.INK,va='top',linespacing=1.35)
for ext in ("png","pdf"): fig.savefig(f"valve_cutaway3d_35cc_120.{ext}",dpi=300)
