import sys, pickle; sys.path.insert(0,'.')
from render_v2 import *
from fields_y import RT,UP
HZ=np.cross(Y,H); HZ/=np.linalg.norm(HZ)   # horizontal, perpendicular to Y
def swivel(theta,elev=12,dist=0.085,focal_mm=(7,-1)):
    t=np.radians(theta); e=np.radians(elev)
    Hh=H-(H@Y)*Y; Hh/=np.linalg.norm(Hh)
    v=np.cos(e)*(np.cos(t)*Hh+np.sin(t)*np.cross(Y,Hh))+np.sin(e)*Y
    focal=O+(focal_mm[0]*A+focal_mm[1]*E2)/1000
    return cam(v,focal,dist)
def backflow_layers(k,thr=0.01,minpx=12):
    """Mid-plane reversed-flow region (from fieldy) as a tinted patch + outline, placed just in front of the cut."""
    from scipy.ndimage import label, binary_opening
    f=pickle.load(open(f"fieldy{k}.pkl","rb")); gx,gy=f['gx'],f['gy']
    minW=0.02*np.percentile(f['W'][f['mask']],90); have=f['mask']&(f['W']>minW)
    GX,GY=np.meshgrid(gx,gy)
    r=have&(f['uax']<-thr)&(((GX*0+1)*((O[None,None,:]+(GX[...,None]*RT+GY[...,None]*UP)/1000-O)@A))*1000>-1)
    r=binary_opening(r,iterations=1)
    lab,n=label(r); sz=np.bincount(lab.ravel()); r=np.isin(lab,[i for i in range(1,n+1) if sz[i]>=minpx])
    P3=O+(GX[...,None]*RT+GY[...,None]*UP)/1000+H*4e-5
    g=pv.StructuredGrid(P3[...,0].T,P3[...,1].T,P3[...,2].T); g['r']=gaussian_like(r.astype(float)).T.ravel(order='F')
    fill=g.contour([0.5],scalars='r') if False else g.threshold(0.5,scalars='r')
    edge=g.extract_surface().contour([0.5],scalars='r') if False else None
    # outline: contour of smoothed mask on the structured surface
    surf=g.extract_surface(); line=surf.contour([0.5],scalars='r')
    return surf.threshold(0.5,scalars='r'),line
def gaussian_like(a):
    from scipy.ndimage import gaussian_filter
    return gaussian_filter(a,1.0)
def render_v3(k,fn,camera,n=800,window=(1400,1400),backflow=False,seed=1,amin=-14e-3):
    c=classify(k); off=c['off']; clip=clipper(amin,None)
    idx=np.random.default_rng(seed).choice(len(off)-1,size=min(n,len(off)-1),replace=False)
    pl=clip(polys(c['X'],off,idx,np.log10(np.clip(c['tau'],1e-2,1e2))))
    p=base_scene(k,clip,window)
    p.add_mesh(pl.tube(radius=0.00008,n_sides=8),scalars='s',cmap='magma',clim=(-2,2),show_scalar_bar=False,smooth_shading=True,ambient=0.3)
    if backflow:
        fill,line=backflow_layers(k)
        p.add_mesh(clip(fill.translate(H*1e-5,inplace=False)) if False else fill,color='#7FB3E6',opacity=0.38,ambient=1,diffuse=0)
        p.add_mesh(line,color='#1F4E8C',line_width=3.0,opacity=0.95)
    p.camera_position=camera; p.enable_anti_aliasing('ssaa'); p.screenshot(fn)
    return p
if __name__=="__main__":
    k=int(sys.argv[1]); th=float(sys.argv[2]); bf=len(sys.argv)>3 and sys.argv[3]=='bf'
    render_v3(k,f"v3_{int(th)}{'_bf' if bf else ''}_{k}.png",swivel(th),backflow=bf)

def reversed_glyphs(k,clip,every=3,thr=-0.25,amin_mm=-1):
    """Small cones on path points that move against the outflow axis (backflow), pointing along the local flow."""
    c=classify(k); X=c['X']; off=c['off']; T=tangents(X,off)
    a=((X-O)@A)*1000; m=(a>amin_mm)&((T@A)<thr)
    # keep runs of >=2 points only, thin out
    idx=np.where(m)[0]; idx=idx[(np.isin(idx+1,idx))|(np.isin(idx-1,idx))][::every]
    pts=pv.PolyData(X[idx]); pts['dir']=T[idx]
    pts=clip(pts)
    if pts.n_points==0: return pts
    return pts.glyph(orient='dir',scale=False,factor=0.0011,geom=pv.Cone(radius=0.32,height=1.0,resolution=12))
def render_v3b(k,fn,camera,n=800,window=(1400,1400),seed=1,amin=-14e-3,glyphs=True,outline=False):
    c=classify(k); off=c['off']; clip=clipper(amin,None)
    idx=np.random.default_rng(seed).choice(len(off)-1,size=min(n,len(off)-1),replace=False)
    pl=clip(polys(c['X'],off,idx,np.log10(np.clip(c['tau'],1e-2,1e2))))
    p=base_scene(k,clip,window)
    p.add_mesh(pl.tube(radius=0.00007,n_sides=8),scalars='s',cmap='magma',clim=(-2,2),show_scalar_bar=False,smooth_shading=True,ambient=0.3,opacity=0.9)
    if glyphs:
        gl=reversed_glyphs(k,clip)
        if gl.n_points: p.add_mesh(gl,color='#7FC4F0',smooth_shading=True,ambient=0.45,specular=0.5)
    if outline:
        _,line=backflow_layers(k)
        p.add_mesh(line,color='#2E6DB4',line_width=3.5)
    p.camera_position=camera; p.enable_anti_aliasing('ssaa'); p.screenshot(fn)
    return p
