import sys; sys.path.insert(0,'.')
from fig_section import *
from PIL import Image
from matplotlib.colors import LogNorm
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
ORDER=[16,33,50,66,83]
def crop_common(files,pad=10):
    ims=[np.asarray(Image.open(f).convert('RGB')) for f in files]
    nw=np.any(np.stack([np.any(im<246,axis=2) for im in ims]),axis=0)
    ys,xs=np.where(nw); y0,y1,x0,x1=max(ys.min()-pad,0),ys.max()+pad,max(xs.min()-pad,0),xs.max()+pad
    return [im[y0:y1,x0:x1] for im in ims]
def film(files,name,legend,phases=ORDER,title_fs=7,strip=0.48):
    ims=crop_common(files)
    asp=ims[0].shape[0]/ims[0].shape[1]
    n=len(phases); W_=rs.WIDTH; pw=W_/n; ih=pw*asp*0.98; tt=0.30
    Ht=tt+ih+strip
    fig=plt.figure(figsize=(W_,Ht))
    for i,(k,im) in enumerate(zip(phases,ims)):
        ax=fig.add_axes([i*pw/W_+0.002,strip/Ht,pw/W_*0.98,ih/Ht]); ax.imshow(im); ax.axis('off')
        ax.set_title(f"{PH[k]:g}% systole\nvalve {ANG[k]:g}°",fontsize=title_fs,loc='center',fontweight='normal',pad=1.5,linespacing=1.15)
    legend(fig,Ht)
    for ext in ("png","pdf"): fig.savefig(f"{name}.{ext}",dpi=300)
    plt.close(fig)
def legend_stress(fig,Ht):
    cax=fig.add_axes([0.30,0.30/Ht,0.40,0.06/Ht])
    cb=fig.colorbar(ScalarMappable(LogNorm(1e-2,1e2),'magma'),cax=cax,orientation='horizontal',extend='both')
    cb.ax.tick_params(labelsize=6,length=2); cb.set_label("Local scalar stress on path [Pa]",fontsize=7,labelpad=1)
def legend_regions(fig,Ht):
    from render_vol import C_VOID
    from render_v2 import C_REC,C_HI,C_CTX
    h=[Patch(fc=C_REC,alpha=0.8),Patch(fc=C_VOID,alpha=0.6),Line2D([],[],color=C_HI,lw=2.2),Line2D([],[],color=C_CTX,lw=1)]
    l=["reversed axial flow","no path reached","scalar stress ≥ 15 Pa","context paths"]
    fig.legend(h,l,loc='center',bbox_to_anchor=(0.5,0.24/Ht),ncol=4,fontsize=6.6,handlelength=1.6,frameon=False,columnspacing=1.6)
if __name__=="__main__":
    film([f"v2_stress_front_valve_{k}.png" for k in ORDER],"cutaway3d_Yup_stress",legend_stress)
    film([f"vol_front_valve_{k}.png" for k in ORDER],"cutaway3d_Yup_regions_front",legend_regions)
    film([f"vol_oblique_down_{k}.png" for k in ORDER],"cutaway3d_Yup_regions_oblique",legend_regions)
