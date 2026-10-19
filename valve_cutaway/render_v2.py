import sys; sys.path.insert(0,'.')
from classify import *
import pyvista as pv, matplotlib
pv.global_theme.allow_empty_mesh=True
from matplotlib.colors import LogNorm
Y=np.array([0,1.0,0])
RING={38,39,40,42}; PIVOT=STRUT-RING
SHELL=[t for t in range(1,52) if t not in DISC and t not in PIVOT]
MAGMA=matplotlib.colormaps['magma']
C_REC,C_HI,C_BOTH,C_CTX='#3D7CC9','#D9577F','#7A3E9D','#BAC2C8'
def cam(view,focal,dist):
    v=np.asarray(view,float); v/=np.linalg.norm(v)
    up=Y-(Y@v)*v; up/=np.linalg.norm(up)
    return [tuple(focal+v*dist),tuple(focal),tuple(up)]
def clipper(amin=None,amax=None):
    def f(m):
        m=m.clip(normal=H,origin=O,invert=True)
        if amin is not None: m=m.clip(normal=-A,origin=O+amin*A,invert=True)
        if amax is not None: m=m.clip(normal=A,origin=O+amax*A,invert=True)
        return m
    return f
def polys(X,off,idx,scal=None):
    pts=[];cells=[];sc=[];o=0
    for i in idx:
        s,e=off[i],off[i+1]
        if e-s<3: continue
        pts.append(X[s:e]); cells.append(np.r_[e-s,np.arange(o,o+e-s)]); o+=e-s
        if scal is not None: sc.append(scal[s:e])
    if not pts: return pv.PolyData()
    pd=pv.PolyData(np.vstack(pts),lines=np.concatenate(cells))
    if scal is not None: pd['s']=np.concatenate(sc)
    return pd
def base_scene(k,clip,window,shell_op=0.16):
    fa,ed=geom(k)
    p=pv.Plotter(off_screen=True,window_size=window); p.set_background('white')
    p.add_mesh(clip(facemesh(fa,SHELL)),color='#C9D3DA',opacity=shell_op,smooth_shading=True,specular=0.15)
    p.add_mesh(clip(edges_mesh(ed)),color='#4A4A4A',line_width=1.1,opacity=0.65)
    sl=clip(facemesh(fa,[t for t in fa if t not in DISC]).slice(normal=H,origin=O-H*2e-5))
    p.add_mesh(sl,color='#1F1F1F',line_width=2.0)
    p.add_mesh(facemesh(fa,DISC),color='#56636D',smooth_shading=True,specular=0.4,ambient=0.28)
    p.add_mesh(clip(facemesh(fa,list(PIVOT))),color='#8C979F',smooth_shading=True,ambient=0.25)
    return p
def render_stress(k,fn,camera,clip,n=800,window=(1400,1400),r=0.00008,seed=1):
    c=classify(k); off=c['off']
    idx=np.random.default_rng(seed).choice(len(off)-1,size=min(n,len(off)-1),replace=False)
    pl=clip(polys(c['X'],off,idx,np.log10(np.clip(c['tau'],1e-2,1e2))))
    p=base_scene(k,clip,window)
    p.add_mesh(pl.tube(radius=r,n_sides=8),scalars='s',cmap='magma',clim=(-2,2),show_scalar_bar=False,smooth_shading=True,ambient=0.3)
    p.camera_position=camera; p.enable_anti_aliasing('ssaa'); p.screenshot(fn)
def runs_poly(X,off,flag,idx=None,pad=1):
    """Polylines of consecutive flagged points (padded by one point) on the chosen paths."""
    pts=[];cells=[];o=0
    for i in (range(len(off)-1) if idx is None else idx):
        s,e=off[i],off[i+1]; f=flag[s:e]
        if not f.any(): continue
        d=np.diff(np.r_[0,f.astype(int),0]); st=np.where(d==1)[0]; en=np.where(d==-1)[0]
        for a,b in zip(st,en):
            a=max(a-pad,0); b=min(b+pad,e-s)
            if b-a<2: continue
            pts.append(X[s+a:s+b]); cells.append(np.r_[b-a,np.arange(o,o+b-a)]); o+=b-a
    if not pts: return pv.PolyData()
    return pv.PolyData(np.vstack(pts),lines=np.concatenate(cells))
TAU_HI=15.0
def render_highlight(k,fn,camera,clip,nctx=350,window=(1400,1400),seed=1,ends=True):
    c=classify(k); off=c['off']; n=len(off)-1; X=c['X']
    P,_=to2d(X); pid=np.repeat(np.arange(n),np.diff(off))
    rng=np.random.default_rng(seed)
    ctx=rng.choice(n,size=min(nctx,n),replace=False)
    p=base_scene(k,clip,window,shell_op=0.12)
    p.add_mesh(clip(polys(X,off,ctx)).tube(radius=0.00004,n_sides=6),color=C_CTX,opacity=0.45,smooth_shading=True,ambient=0.5)
    recpt=c['rec'][pid]&(P[:,0]>-3)
    p.add_mesh(clip(runs_poly(X,off,recpt,pad=0)).tube(radius=0.00009,n_sides=8),color=C_REC,smooth_shading=True,ambient=0.3,specular=0.3)
    hipt=c['tau']>=TAU_HI
    p.add_mesh(clip(runs_poly(X,off,hipt)).tube(radius=0.00011,n_sides=8),color=C_HI,smooth_shading=True,ambient=0.3,specular=0.3)
    if ends:
        e=X[off[1:]-1][c['wall']]
        p.add_mesh(clip(pv.PolyData(e)),color='#1F1F1F',point_size=4,render_points_as_spheres=True)
    p.camera_position=camera; p.enable_anti_aliasing('ssaa'); p.screenshot(fn)
VIEWS={
 # name: (view direction from focal, focal offset (A,E2 mm), distance m, amin, amax)
 'front_full':   (H,               (-6,0),  0.150, None, None),
 'front_valve':  (H,               (7,-1),  0.082, -14, None),
 'oblique_down': (H+0.55*A+0.35*E2,(8,0),   0.085, -14, None),
 'oblique_up':   (H-0.45*A+0.45*E2,(5,0),   0.070, -14, None),
}
def view(name):
    v,(fa_,fe),d,amin,amax=VIEWS[name]
    focal=O+(fa_*A+fe*E2)/1000
    return cam(v,focal,d),clipper(None if amin is None else amin/1000,None if amax is None else amax/1000)
if __name__=="__main__":
    k=int(sys.argv[1]); name=sys.argv[2]; mode=sys.argv[3]
    camera,clip=view(name)
    (render_stress if mode=='stress' else render_highlight)(k,f"v2_{mode}_{name}_{k}.png",camera,clip)
