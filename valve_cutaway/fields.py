import sys; sys.path.insert(0,'.')
from section import *
from scipy.ndimage import gaussian_filter
SLAB=1.0   # mm half-thickness
def domain_mask(fa, gx, gy):
    surf=facemesh(fa,list(fa.keys())).triangulate().clean()
    GX,GY=np.meshgrid(gx,gy)
    P3=O+ (GX.ravel()[:,None]*A + GY.ravel()[:,None]*E2)/1000
    pts=pv.PolyData(P3)
    r=pts.select_enclosed_points(surf, tolerance=1e-6, check_surface=False)
    print(r.array_names); return r[r.array_names[0]].reshape(GX.shape).astype(bool)
def binned(P,vals,w,gx,gy,sigma):
    # normalised convolution on grid
    h=gx[1]-gx[0]
    ix=np.clip(np.round((P[:,0]-gx[0])/h).astype(int),0,len(gx)-1)
    iy=np.clip(np.round((P[:,1]-gy[0])/h).astype(int),0,len(gy)-1)
    W=np.zeros((len(gy),len(gx))); np.add.at(W,(iy,ix),w)
    out=[]
    for v in vals:
        S=np.zeros_like(W); np.add.at(S,(iy,ix),w*v); out.append(gaussian_filter(S,sigma))
    Ws=gaussian_filter(W,sigma)
    return [o/np.maximum(Ws,1e-12) for o in out], Ws
def phase_field(k, h=0.5, sigma=1.4):
    fa,ed=geom(k); X,tau,U,off=lines(k)
    P,d=to2d(X); T=tangents(X,off)
    w=np.exp(-0.5*(d/ (SLAB/1.5))**2); sel=np.abs(d)<SLAB
    gx=np.arange(-48,27+h,h); gy=np.arange(-22,37+h,h)
    V=T*U[:,None]
    (sp,ua,ub,ta),Ws=binned(P[sel],[U[sel],V[sel]@A,V[sel]@E2,tau[sel]],w[sel],gx,gy,sigma)
    m=domain_mask(fa,gx,gy)
    wall=section(fa,[t for t in fa if t not in DISC]); disc=section(fa,DISC)
    return dict(gx=gx,gy=gy,speed=sp,ua=ua,ub=ub,tau=ta,W=Ws,mask=m,wall=wall,disc=disc)
if __name__=="__main__":
    import pickle
    for k in PH:
        f=phase_field(k); pickle.dump(f,open(f"field{k}.pkl","wb")); print(k, f['mask'].mean(), np.nanmax(f['speed']))
