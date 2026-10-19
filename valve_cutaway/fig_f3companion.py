import sys, pickle; sys.path.insert(0,'.')
from fig_section import *
from matplotlib.colors import LogNorm
from scipy.ndimage import label
def draw_y(ax,f,field,cmap,norm,line_col='white',line_alpha=0.55,rev=True,density=1.5,xl=(-22,22.5),yl=(-14,27.5)):
    gx,gy=f['gx'],f['gy']; m=f['mask']; W=f['W']
    minW=0.02*np.percentile(W[m],90); have=m&(W>minW)
    ext=(gx[0],gx[-1],gy[0],gy[-1])
    ax.imshow(np.where(m,1.0,np.nan),extent=ext,origin='lower',cmap=matplotlib.colors.ListedColormap([VOID]),zorder=0,interpolation='nearest')
    im=ax.imshow(np.where(have,field,np.nan),extent=ext,origin='lower',cmap=cmap,norm=norm,zorder=1,interpolation='bilinear')
    U_=np.where(have,f['ua'],np.nan); V_=np.where(have,f['ub'],np.nan); s=np.hypot(U_,V_)
    st=ax.streamplot(gx,gy,U_/np.maximum(s,1e-9),V_/np.maximum(s,1e-9),density=density,color=line_col,linewidth=0.4,
                     arrowsize=0.4,minlength=0.08,maxlength=4,zorder=4)
    st.lines.set_alpha(line_alpha); st.arrows.set_alpha(line_alpha)
    if rev:
        GX,GY=np.meshgrid(gx,gy)
        r=have&(f['uax']<-0.01)&(GY>-2)
        lab,n=label(r); sz=np.bincount(lab.ravel()); r=np.isin(lab,[i for i in range(1,n+1) if sz[i]>=12])
        if r.any(): ax.contour(gx,gy,r.astype(float),levels=[0.5],colors=rs.INK,linewidths=0.7,linestyles='--',zorder=3)
    for c in f['wall']: ax.plot(c[:,0],c[:,1],color=rs.INK,lw=0.7,zorder=5)
    for c in f['disc']: ax.fill(c[:,0],c[:,1],color='#3A3A3A',ec=rs.INK,lw=0.6,zorder=6)
    ax.set_xlim(*xl); ax.set_ylim(*yl); ax.set_aspect('equal'); ax.axis('off')
    return im
def scalebar(ax,x0,y0,L=5):
    ax.plot([x0,x0+L],[y0,y0],color=rs.INK,lw=1.2,zorder=7,solid_capstyle='butt')
    ax.text(x0+L/2,y0+0.6,f"{L} mm",ha='center',va='bottom',fontsize=5.8,color=rs.INK,zorder=7)
if __name__=="__main__":
    phases=[16,50,83]
    F={k:pickle.load(open(f"fieldy{k}.pkl","rb")) for k in phases}
    fig,axs=plt.subplots(2,3,figsize=(rs.WIDTH,4.35),gridspec_kw=dict(left=0.005,right=0.84,top=0.94,bottom=0.01,wspace=0.03,hspace=0.10))
    n1=LogNorm(1e-2,1e2); cm2=plt.get_cmap('magma_r').copy(); n2=LogNorm(1e-6,1e-2)
    for j,k in enumerate(phases):
        f=F[k]
        im1=draw_y(axs[0,j],f,10**f['ltau'],'magma',n1)
        axs[0,j].set_title(f"{PH[k]:g}% systole · valve {ANG[k]:g}°",fontsize=7.2,pad=2)
        im2=draw_y(axs[1,j],f,f['dose'],cm2,n2,line_col='#4A4A4A',line_alpha=0.5)
        axs[1,j].text(0.03,0.03,f"mid-plane slab (±1 mm)\nholds {100*f['slab_share']:.0f}% of the dose",transform=axs[1,j].transAxes,
                      fontsize=5.9,color=rs.INK,va='bottom',bbox=dict(fc='white',ec='none',alpha=0.85,pad=1.2),zorder=8)
    scalebar(axs[0,0],12,24)
    axs[0,0].annotate("",xy=(-19,27),xytext=(-19,21),arrowprops=dict(arrowstyle='-|>',color=rs.INK,lw=0.8))
    axs[0,0].text(-18.2,24,"+Y\n(outflow)",fontsize=5.8,color=rs.INK,va='center')
    c1=fig.colorbar(im1,cax=fig.add_axes([0.865,0.575,0.018,0.34]),extend='both'); c1.set_label("Local scalar stress [Pa]",fontsize=7)
    c2=fig.colorbar(im2,cax=fig.add_axes([0.865,0.04,0.018,0.33]),extend='min'); c2.set_label("Share of linear HI3 dose\nper mm² of slab",fontsize=7)
    for c in (c1,c2): c.ax.tick_params(labelsize=6)
    fig.text(0.85,0.525,"lines: in-plane flow\ndirection only\n– – reversed axial flow\n   (> 0.01 m/s)\ngrey: no path ±1 mm",fontsize=5.6,color=rs.INK,va='top',linespacing=1.25)
    for ext in ("png","pdf"): fig.savefig(f"F3_companion_section.{ext}",dpi=300)
