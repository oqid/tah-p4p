import sys; sys.path.insert(0,'.')
from section import *
import pyvista as pv
pv.global_theme.allow_empty_mesh=True
from matplotlib.colors import LinearSegmentedColormap
CMAP=LinearSegmentedColormap.from_list("p4p3d",["#7FA3CF","#8C7BC0","#9B4F96","#D0668F","#EE9A6B","#F6C85F"])
RING={38,39,40,42}
PIVOT=STRUT-RING
SHELL=[t for t in range(1,52) if t not in DISC and t not in PIVOT]
AMIN=-0.016
def clip_half(m, keep_back=True):
    m=m.clip(normal=H if keep_back else -H, origin=O, invert=True)
    return m.clip(normal=-A, origin=O+AMIN*A, invert=True)
def path_poly(X,tau,U,off,idx):
    pts=[];cells=[];sc=[];o=0
    for i in idx:
        s,e=off[i],off[i+1]; x=X[s:e]
        if len(x)<3: continue
        pts.append(x); sc.append(tau[s:e]); cells.append(np.r_[len(x),np.arange(o,o+len(x))]); o+=len(x)
    pd=pv.PolyData(np.vstack(pts),lines=np.concatenate(cells)); pd['tau']=np.concatenate(sc); return pd
def render(k,fn,nlines=700,seed=1,window=(1400,1000),cam=None,clim=(0.5,100)):
    fa,ed=geom(k); X,tau,U,off=lines(k)
    rng=np.random.default_rng(seed); idx=rng.choice(len(off)-1,size=min(nlines,len(off)-1),replace=False)
    pl=path_poly(X,tau,U,off,idx)
    pl=clip_half(pl)   # keep paths behind the cut plane, valve region only
    pl['ltau']=np.log10(np.clip(pl['tau'],clim[0],clim[1]))
    p=pv.Plotter(off_screen=True,window_size=window); p.set_background('white')
    shell=clip_half(facemesh(fa,SHELL))
    p.add_mesh(shell,color='#C9D3DA',opacity=0.18,smooth_shading=True,specular=0.2)
    edges=clip_half(edges_mesh(ed))
    p.add_mesh(edges,color='#3A3A3A',line_width=1.2,opacity=0.7)
    # cut outline on the plane
    sl=clip_half(facemesh(fa,[t for t in fa if t not in DISC]).slice(normal=H,origin=O-H*2e-5))
    p.add_mesh(sl,color='#1F1F1F',line_width=2.2)
    p.add_mesh(facemesh(fa,DISC),color='#5E6B75',smooth_shading=True,specular=0.4,ambient=0.25)
    p.add_mesh(clip_half(facemesh(fa,list(PIVOT))),color='#8C979F',smooth_shading=True,ambient=0.25)
    tubes=pl.tube(radius=0.00008,n_sides=8)
    p.add_mesh(tubes,scalars='ltau',cmap=CMAP,clim=np.log10(clim),show_scalar_bar=False,smooth_shading=True,ambient=0.35)
    if cam: p.camera_position=cam
    p.enable_anti_aliasing('ssaa')
    p.screenshot(fn,transparent_background=False); return p.camera_position
if __name__=="__main__":
    # camera: look from the +H side (front), slightly from above and from the outlet side
    focal=O+(5*A-1*E2)/1000
    view=( 1.0*H + 0.50*E2 - 0.45*A); view/=np.linalg.norm(view)
    cam=[tuple(focal+view*0.068), tuple(focal), tuple(E2)]
    k=int(sys.argv[1]) if len(sys.argv)>1 else 50
    render(k,f"r3d_{k}.png",cam=cam)
