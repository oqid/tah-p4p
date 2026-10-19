import sys, json; sys.path.insert(0,'.')
from fig_film import *
from matplotlib.patches import FancyArrowPatch
from matplotlib.path import Path as MPath
C_BF='#2E6DB4'
SHADE='#8DB8E8'; SHADE_A=0.55
def shaded(f):
    """Multiply-blend a soft tint over the backflow volume's screen footprint (mask_<tag>_<k>.png)."""
    from scipy.ndimage import gaussian_filter
    im=np.asarray(Image.open(f).convert('RGB')).astype(float)/255
    mf=f.replace('v4_','mask_')
    try: m=1-np.asarray(Image.open(mf).convert('L')).astype(float)/255
    except FileNotFoundError: return (im*255).astype(np.uint8)
    m=np.clip(gaussian_filter(m,2.5),0,1)[...,None]*SHADE_A
    tint=np.array(matplotlib.colors.to_rgb(SHADE))
    out=im*(1-m+m*tint)
    return (np.clip(out,0,1)*255).astype(np.uint8)
def crop_box(files,pad=10):
    ims=[shaded(f) for f in files]
    nw=np.any(np.stack([np.any(im<246,axis=2) for im in ims]),axis=0)
    ys,xs=np.where(nw); return ims,(max(ys.min()-pad,0),ys.max()+pad,max(xs.min()-pad,0),xs.max()+pad)
def recirc_glyph(ax,cx,cy,r,sense=1,col='#1F4E8C',lw=0.7):
    """Open circular arrow (≈290°) drawn in image coordinates (y down); sense=+1 is anticlockwise on screen."""
    t=np.linspace(np.radians(60),np.radians(350),60)
    x=cx+r*np.cos(t); y=cy-sense*r*np.sin(t)
    ax.plot(x[:-3],y[:-3],color=col,lw=lw,solid_capstyle='round',zorder=6)
    ax.add_patch(FancyArrowPatch((x[-6],y[-6]),(x[-1],y[-1]),arrowstyle='-|>',mutation_scale=6,color=col,lw=lw,zorder=6,shrinkA=0,shrinkB=0))
def overlay(ax,J,box,backflow=False,glyph=False,labels=None,lab_fs=6.2,glyph_r=None):
    y0,x0=box[0],box[2]
    tr=lambda p:(p[0]-x0,p[1]-y0)
    if backflow and 'backflow_outline' in J:
        L=[np.array(l) for l in J['backflow_outline']]
        lens=[np.sum(np.hypot(*np.diff(l,axis=0).T)) for l in L]; mx=max(lens) if lens else 0
        for ln,le in zip(L,lens):
            if le<0.35*mx: continue
            ax.plot(ln[:,0]-x0,ln[:,1]-y0,color='white',lw=1.6,alpha=0.8,zorder=4)
            ax.plot(ln[:,0]-x0,ln[:,1]-y0,color=C_BF,lw=0.8,ls=(0,(3,2)),zorder=4.5)
    if glyph and 'recirc' in J:
        cx,cy=tr(J['recirc']); r=glyph_r or min(0.7*J.get('recirc_scale',100),60)
        recirc_glyph(ax,cx,cy,r)
    if labels:
        for key,(text,(dx,dy),ha) in labels.items():
            if key not in J: continue
            px,py=tr(J[key])
            ax.annotate(text,xy=(px,py),xytext=(px+dx,py+dy),fontsize=lab_fs,color=rs.INK,ha=ha,va='center',zorder=8,
                        arrowprops=dict(arrowstyle='-',lw=0.55,color=rs.INK,shrinkA=1.5,shrinkB=1.5),
                        bbox=dict(fc='white',ec='none',alpha=0.8,pad=0.6))
            ax.plot(px,py,'o',ms=2.2,color=rs.INK,mec='white',mew=0.5,zorder=9)
def film_v4(tag,name,phases=ORDER,strip=0.5,title_fs=7):
    files=[f"v4_{tag}_{k}.png" for k in phases]; ims,box=crop_box(files)
    y0,y1,x0,x1=box; ims=[im[y0:y1,x0:x1] for im in ims]
    asp=ims[0].shape[0]/ims[0].shape[1]; n=len(phases); W_=rs.WIDTH; pw=W_/n; ih=pw*asp*0.98; tt=0.30; Ht=tt+ih+strip
    fig=plt.figure(figsize=(W_,Ht))
    for i,(k,im) in enumerate(zip(phases,ims)):
        ax=fig.add_axes([i*pw/W_+0.002,strip/Ht,pw/W_*0.98,ih/Ht]); ax.imshow(im); ax.axis('off')
        ax.set_xlim(0,im.shape[1]); ax.set_ylim(im.shape[0],0)
        J=json.load(open(f"v4_{tag}_{k}.json")); overlay(ax,J,box)
        ax.set_title(f"{PH[k]:g}% systole\nvalve {ANG[k]:g}°",fontsize=title_fs,pad=1.5,linespacing=1.15,fontweight='normal')
    cax=fig.add_axes([0.08,0.30/Ht,0.36,0.06/Ht])
    cb=fig.colorbar(ScalarMappable(LogNorm(1e-2,1e2),'magma'),cax=cax,orientation='horizontal',extend='both')
    cb.ax.tick_params(labelsize=6,length=2); cb.set_label("Local scalar stress on path [Pa]",fontsize=7,labelpad=1)
    kx=fig.add_axes([0.52,0.05/Ht,0.46,0.42/Ht]); kx.axis('off'); kx.set_xlim(0,10); kx.set_ylim(0,2)
    kx.add_patch(matplotlib.patches.Rectangle((0.1,1.15),0.8,0.4,fc=SHADE,alpha=SHADE_A,ec='none')); kx.text(1.1,1.35,"backflow region (axial velocity < 0)",fontsize=6.4,va='center')
    recirc_glyph(kx,0.5,0.55,0.28) if False else None
    for ext in ("png","pdf"): fig.savefig(f"{name}.{ext}",dpi=300)
    plt.close(fig)
def key_figure(tag,name,phases=(16,50),labels_for=None):
    files=[f"v4_{tag}_{k}.png" for k in phases]; ims,box=crop_box(files)
    y0,y1,x0,x1=box; ims=[im[y0:y1,x0:x1] for im in ims]
    h,w=ims[0].shape[:2]; m=0.78*w
    fig_w=rs.WIDTH; pan_w=fig_w/2; fig_h=pan_w*h/(w+m)+0.3
    fig=plt.figure(figsize=(fig_w,fig_h))
    for i,(k,im) in enumerate(zip(phases,ims)):
        side='L' if i==0 else 'R'
        ax=fig.add_axes([i*0.5,0,0.5,(fig_h-0.3)/fig_h])
        J=json.load(open(f"v4_{tag}_{k}.json"))
        ax.imshow(im); ax.set_xlim((-m,w) if side=='L' else (0,w+m)); ax.set_ylim(h,0); ax.axis('off')
        overlay(ax,J,box)
        items=[(kk,t) for kk,t in labels_for(k,J).items() if kk in J]
        pts={kk:(J[kk][0]-x0,J[kk][1]-y0) for kk,_ in items}
        items=sorted(items,key=lambda kv:pts[kv[0]][1])
        ys=np.linspace(h*0.18,h*0.86,len(items))
        for (kk,text),yt in zip(items,ys):
            px,py=pts[kk]; xt=-0.06*w if side=='L' else w*1.04
            ax.annotate(text,xy=(px,py),xytext=(xt,yt),fontsize=6.6,color=rs.INK,ha='right' if side=='L' else 'left',va='center',zorder=8,linespacing=1.2,
                        arrowprops=dict(arrowstyle='-',lw=0.5,color=rs.MUTED,shrinkA=3,shrinkB=2.5))
            ax.plot(px,py,'o',ms=2.8,color=rs.INK,mec='white',mew=0.7,zorder=9)
        ax.text(0.5 if True else 0,1.0,f"{PH[k]:g}% systole · valve {ANG[k]:g}°",transform=ax.transAxes,fontsize=7.5,ha='center',va='bottom')
    for ext in ("png","pdf"): fig.savefig(f"{name}.{ext}",dpi=300)
    plt.close(fig)
def LABELS(k,J):
    if k in (16,83,33,66):
        L={'jet':'major-orifice jet\n(high stress)','shear':'shear layer\n(jet edge)','sep':'flow separates\nat disc edge','recirc':'recirculation zone\nbehind disc (backflow)'}
        if k in (33,66): L['jet_minor']='minor-orifice jet'
    else:
        L={'jet':'minor-orifice jet','jet_minor':'major-orifice jet','recirc':'thin separated layer\nat housing step (backflow)'}
    return L
if __name__=="__main__":
    if sys.argv[1]=='key': key_figure(sys.argv[2],sys.argv[3],labels_for=LABELS)
    else: film_v4(sys.argv[1],sys.argv[2])
