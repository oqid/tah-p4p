import sys; sys.path.insert(0,'.')
from render_v2 import *
C_VOID='#B9A3D6'
def render_volumes(k,fn,camera,amin=-14e-3,window=(1400,1400),nctx=220,seed=1,cut=False,smooth=True):
    fa,ed=geom(k); c=classify(k); off=c['off']; n=len(off)-1
    clipA=lambda m: m.clip(normal=-A,origin=O+amin*A,invert=True)
    clip=(lambda m: clipA(m.clip(normal=H,origin=O,invert=True))) if cut else clipA
    p=pv.Plotter(off_screen=True,window_size=window); p.set_background('white')
    p.add_mesh(clip(facemesh(fa,SHELL)),color='#C9D3DA',opacity=0.10,smooth_shading=True)
    p.add_mesh(clip(edges_mesh(ed)),color='#4A4A4A',line_width=1.1,opacity=0.6)
    p.add_mesh(facemesh(fa,DISC),color='#56636D',smooth_shading=True,specular=0.4,ambient=0.28)
    p.add_mesh(clip(facemesh(fa,list(PIVOT))),color='#8C979F',smooth_shading=True,ambient=0.25)
    g=pv.read(f"vol{k}.vtk")
    vd=g.contour([0.5],scalars='void'); rv=g.contour([0.03],scalars='rev')
    if smooth:
        vd=vd.smooth_taubin(n_iter=30) if vd.n_points else vd; rv=rv.smooth_taubin(n_iter=30) if rv.n_points else rv
    p.add_mesh(clip(vd),color=C_VOID,opacity=0.45,smooth_shading=True,specular=0.2)
    p.add_mesh(clip(rv),color=C_REC,opacity=0.75,smooth_shading=True,specular=0.3)
    ctx=np.random.default_rng(seed).choice(n,size=min(nctx,n),replace=False)
    p.add_mesh(clip(polys(c['X'],off,ctx)).tube(radius=0.00004,n_sides=6),color=C_CTX,opacity=0.5,smooth_shading=True,ambient=0.5)
    P,_=to2d(c['X']); hipt=c['tau']>=TAU_HI
    p.add_mesh(clip(runs_poly(c['X'],off,hipt)).tube(radius=0.0001,n_sides=8),color=C_HI,smooth_shading=True,ambient=0.3,specular=0.3)
    p.camera_position=camera; p.enable_anti_aliasing('ssaa'); p.screenshot(fn)
if __name__=="__main__":
    k=int(sys.argv[1]); name=sys.argv[2]
    camera,_=view(name)
    render_volumes(k,f"vol_{name}_{k}.png",camera,cut=len(sys.argv)>3)
