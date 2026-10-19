import sys, pickle; sys.path.insert(0,'.')
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, matplotlib.patheffects as pe
from matplotlib.colors import LinearSegmentedColormap, PowerNorm
from scipy.interpolate import PchipInterpolator
import report_style as rs
rs.apply()
PH={16:16.67,33:33.33,50:50.0,66:66.67,83:83.33}
ANG={16:24.9,33:49.8,50:75.0,66:50.1,83:25.2}         # valve angle from closed, operating_point.py (same as F3)
Q={16:6.78,33:11.67,50:14.51,66:13.94,83:8.86}        # inlet Q, 35 cc 120 BPM (flows_by_phase.csv)
CMAP=LinearSegmentedColormap.from_list("p4p",["#FBF1E8","#F4C7A4","#EE9A6B","#D0668F","#9B4F96","#5B3A8E","#2D1E5C"])
NORM=PowerNorm(0.5,0,2.0)
VOID="#D9DEE2"
def draw(ax,f,k,xl=(-48,26),yl=(-22,37),density=1.6,minW=None,arrows=True):
    gx,gy=f['gx'],f['gy']; m=f['mask']; W=f['W']
    if minW is None: minW=0.02*np.percentile(W[m],90)
    have=m&(W>minW)
    sp=np.where(have,f['speed'],np.nan)
    ax.imshow(np.where(m,1.0,np.nan),extent=(gx[0],gx[-1],gy[0],gy[-1]),origin='lower',cmap=matplotlib.colors.ListedColormap([VOID]),zorder=0,interpolation='nearest')
    im=ax.imshow(sp,extent=(gx[0],gx[-1],gy[0],gy[-1]),origin='lower',cmap=CMAP,norm=NORM,zorder=1,interpolation='bilinear')
    # reversed flow in the outflow tract (downstream of the valve plane): axial velocity < 0
    GX,GY=np.meshgrid(gx,gy)
    rev=have&(GX>-9)&(f['ua']<-0.03)
    from scipy.ndimage import label
    lab,n=label(rev); sizes=np.bincount(lab.ravel())
    rev=np.isin(lab,[i for i in range(1,n+1) if sizes[i]>=12])
    if rev.any():
        ax.contour(gx,gy,rev.astype(float),levels=[0.5],colors='#1F1F1F',linewidths=0.7,linestyles='--',zorder=3)
    if arrows:
        U_=np.where(have,f['ua'],np.nan); V_=np.where(have,f['ub'],np.nan)
        sp_=np.hypot(U_,V_); U_=U_/np.maximum(sp_,1e-9); V_=V_/np.maximum(sp_,1e-9)   # direction only
        st=ax.streamplot(gx,gy,U_,V_,density=density,color='white',linewidth=0.45,arrowsize=0.45,
                         arrowstyle='-|>',minlength=0.08,maxlength=4,broken_streamlines=True,zorder=4)
        st.lines.set_path_effects([pe.withStroke(linewidth=1.1,foreground='#1F1F1F',alpha=0.55)])
    for c in f['wall']: ax.plot(c[:,0],c[:,1],color='#1F1F1F',lw=0.7,zorder=5,solid_capstyle='round')
    for c in f['disc']: ax.fill(c[:,0],c[:,1],color='#3A3A3A',ec='#1F1F1F',lw=0.6,zorder=6)
    ax.set_xlim(*xl); ax.set_ylim(*yl); ax.set_aspect('equal'); ax.axis('off')
    return im
def timeline(ax,cur=None):
    ph=np.array([0]+[PH[k] for k in PH]+[100]); q=np.array([0]+[Q[k] for k in PH]+[0])
    t=np.linspace(0,100,200)
    ax.plot(t,PchipInterpolator(ph,q)(t),color=rs.MUTED,lw=1,ls=':')
    ax.plot([PH[k] for k in PH],[Q[k] for k in PH],'-',color=rs.CASE_STYLE["35 cc 120 BPM"]["color"],lw=1.2)
    for i,k in enumerate(PH):
        ax.plot(PH[k],Q[k],marker=rs.PHASE_MARKERS[i],color=rs.PHASE_COLOURS[i],ms=6 if k==cur else 4,mec='k' if k==cur else 'none',zorder=5)
    ax.set_xlim(0,100); ax.set_ylim(0,16); ax.set_yticks([0,5,10,15]); ax.set_ylabel("Inlet Q (L/min)",color=rs.CASE_STYLE["35 cc 120 BPM"]["color"])
    a2=ax.twinx(); a2.plot([PH[k] for k in PH],[ANG[k] for k in PH],'--',color='#3A3A3A',lw=1)
    a2.set_ylim(0,80); a2.set_ylabel("Disc opening (°) - -"); a2.spines['right'].set_visible(True); a2.grid(False)
    ax.set_xlabel("Systole (%)",labelpad=1); ax.set_xticks([0,16.67,33.33,50,66.67,83.33,100]); ax.set_xticklabels(["0","17","33","50","67","83","100"])
    return a2

TITLE=lambda i,k: f"{PH[k]:g}% systole · valve {ANG[k]:g}°"
def legend_block(fig,im,x,y,w):
    cax=fig.add_axes([x,y,w,0.016]); cb=fig.colorbar(im,cax=cax,orientation='horizontal',extend='max')
    cb.set_ticks([0,0.1,0.5,1,2]); cb.ax.tick_params(labelsize=6,length=2)
    cb.ax.set_title("Path-sampled speed (m/s)",fontsize=6.5,loc='left',pad=3)
    fig.text(x,y-0.05,"white lines: in-plane flow direction\n– – reversed axial flow past the valve\ngrey: no path within ±1 mm of plane",fontsize=6,color=rs.INK,va='top',linespacing=1.35)
def make(F,name,xl,yl,density,figsize,gskw,tl_box,lg):
    fig=plt.figure(figsize=figsize)
    gs=fig.add_gridspec(2,3,**gskw)
    for i,k in enumerate([16,33,50,66,83]):
        ax=fig.add_subplot(gs[i//3,i%3]); im=draw(ax,F[k],k,xl=xl,yl=yl,density=density)
        ax.set_title(TITLE(i,k),fontsize=6.8,loc='left',fontweight='normal',pad=2)
    legend_block(fig,im,*lg)
    for ext in ("png","pdf"): fig.savefig(f"{name}.{ext}",dpi=300)
    plt.close(fig)
if __name__=="__main__":
    F={k:pickle.load(open(f"field{k}.pkl","rb")) for k in PH}
    make(F,"valve_section_full_35cc_120",(-48,26),(-22,37),1.5,(rs.WIDTH,3.95),
         dict(hspace=0.07,wspace=0.03,left=0.005,right=0.995,top=0.955,bottom=0.01),[0.74,0.335,0.19,0.15],(0.70,0.36,0.27))
    make(F,"valve_section_zoom_35cc_120",(-14,26),(-14,11),1.9,(rs.WIDTH,3.3),
         dict(hspace=0.16,wspace=0.03,left=0.005,right=0.995,top=0.95,bottom=0.02),[0.74,0.355,0.19,0.13],(0.70,0.36,0.27))
