"""Flow features on the mid-plane for annotation: backflow outline, recirculation centre+sense, jet, shear layer, separation tip, sheltered pocket."""
import sys, pickle; sys.path.insert(0,'.')
from render_v2 import *
from fields_y import RT,UP
from scipy.ndimage import label, gaussian_filter
def to3(x,y): return O+(x*RT+y*UP)/1000
def backflow_outline(k,level=0.02):
    g=pv.read(f"vol{k}.vtk")
    s=g.slice(normal=H,origin=O)
    s['rev']=np.nan_to_num(s['rev'],nan=-1)
    ln=s.contour([level],scalars='rev')
    return ln.translate(H*6e-5,inplace=False)
def feats(k):
    f=pickle.load(open(f"fieldy{k}.pkl","rb")); gx,gy=f['gx'],f['gy']; GX,GY=np.meshgrid(gx,gy)
    minW=0.02*np.percentile(f['W'][f['mask']],90); have=f['mask']&(f['W']>minW)
    P3=O[None,None,:]+(GX[...,None]*RT+GY[...,None]*UP)/1000
    a=((P3-O)@A)*1000
    out={}; h=gx[1]-gx[0]
    # backflow region on the plane from the 3D field
    g=pv.read(f"vol{k}.vtk"); probe=pv.PolyData(P3.reshape(-1,3)).sample(g)
    rev=np.nan_to_num(probe['rev'],nan=-1).reshape(GX.shape)
    R=(rev>0.02)&f['mask']&(a>0)
    lab,n=label(R)
    if n:
        sz=np.bincount(lab.ravel()); sz[0]=0; big=lab==sz.argmax()
        if big.sum()>40:
            cy,cx=np.mean(GY[big]),np.mean(GX[big])
            # vorticity sense from smoothed in-plane velocity (speed-weighted)
            u=np.where(have,f['ua']*f['speed'],0); v=np.where(have,f['ub']*f['speed'],0)
            h=gx[1]-gx[0]; w=np.gradient(v,h,axis=1)-np.gradient(u,h,axis=0)
            sense=np.sign(np.nanmean(np.where(big&have,w,np.nan)))
            ext=np.sqrt(big.sum()*h*h/np.pi)
            out['recirc']=dict(p=to3(cx,cy),sense=float(sense),r_mm=float(ext),xy=(cx,cy))
    # jet: fastest point just downstream of the valve
    sp=np.where(have&(a>2)&(a<12),f['speed'],-1); j=np.unravel_index(np.argmax(sp),sp.shape)
    out['jet']=dict(p=to3(GX[j],GY[j]),xy=(GX[j],GY[j]),speed=float(sp[j]))
    # shear layer: between jet and recirculation centre, where speed is halfway
    if 'recirc' in out:
        cx,cy=out['recirc']['xy']; jx,jy=out['jet']['xy']
        t=np.linspace(0,1,200); xs=cx+(jx-cx)*t; ys=cy+(jy-cy)*t
        ii=np.clip(np.round((ys-gy[0])/h).astype(int),0,len(gy)-1); jj=np.clip(np.round((xs-gx[0])/h).astype(int),0,len(gx)-1)
        s=np.nan_to_num(f['speed'][ii,jj]); half=0.5*(s[0]+out['jet']['speed']); m=np.argmax(s>=half)
        out['shear']=dict(p=to3(xs[m],ys[m]),xy=(xs[m],ys[m]))
    # separation: disc tip nearest the recirculation (or the minor side)
    D=np.vstack(f['disc']); 
    ref=np.array(out['recirc']['xy']) if 'recirc' in out else np.array(out['jet']['xy'])*-1
    ends=[D[np.argmin(D@d)] for d in (np.array([1,0.3]),np.array([-1,-0.3]))]
    tip=min(ends,key=lambda e: np.linalg.norm(e-ref))
    out['sep']=dict(p=to3(*tip),xy=tuple(tip))
    # sheltered pocket: largest no-path region downstream of valve
    V=f['mask']&~have&(a>2); lab,n=label(V)
    if n:
        sz=np.bincount(lab.ravel()); sz[0]=0; big=lab==sz.argmax()
        if big.sum()*h*h>15:
            out['pocket']=dict(p=to3(np.mean(GX[big]),np.mean(GY[big])),xy=(np.mean(GX[big]),np.mean(GY[big])),area=float(big.sum()*h*h))
    return out
if __name__=="__main__":
    for k in PH:
        F=feats(k); print(PH[k],{kk:(np.round(v['xy'],1).tolist(), v.get('sense'),round(v.get('r_mm',0),1), round(v.get('area',0))) for kk,v in F.items()})

def largest(mesh):
    if mesh.n_points==0: return None
    c=mesh.connectivity('largest'); return c
def feats3d(k,rear=True):
    g=pv.read(f"vol{k}.vtk"); c=classify(k); X,U,off=c['X'],c['U'],c['off']; T=tangents(X,off)
    side=lambda P: ((P-O)@H)<0 if rear else np.ones(len(P),bool)
    out={}
    rv=g.contour([0.03],scalars='rev'); rv=rv.clip(normal=H,origin=O,invert=True) if rear else rv
    rl=largest(rv)
    if rl is not None and rl.n_points>50:
        cen=np.asarray(rl.points).mean(0)
        b=rl.bounds; size=np.linalg.norm([b[1]-b[0],b[3]-b[2],b[5]-b[4]])*1000
        # rotation sense about the view axis H from path points near the blob
        d=np.linalg.norm(X-cen,axis=1); near=(d<0.004)&side(X)
        Lh=np.cross(X[near]-cen,T[near]*U[near,None])@H
        out['recirc']=dict(p=cen,sense=float(np.sign(Lh.mean())),size_mm=float(size),n=int(near.sum()))
    a=((X-O)@A)*1000; m=side(X)&(a>2)&(a<12)
    thr=np.percentile(U[m],98); jp=X[m&(U>=thr)].mean(0); out['jet']=dict(p=jp,speed=float(U[m].max()))
    e=(X-O)@E2; sj=np.sign((jp-O)@E2); m2=m&(np.sign(e)!=sj)&(np.abs(e)>0.002)
    if m2.sum()>50:
        t2=np.percentile(U[m2],98); jp2=X[m2&(U>=t2)].mean(0)
        out['jet_minor']=dict(p=jp2,speed=float(t2),ratio=float(t2/thr))
    if 'recirc' in out:
        out['shear']=dict(p=0.5*(jp+out['recirc']['p']))
    f=pickle.load(open(f"fieldy{k}.pkl","rb")); D=np.vstack(f['disc'])
    ends=[D[np.argmin(D@d)] for d in (np.array([1,0.3]),np.array([-1,-0.3]))]
    ends3=[to3(*e) for e in ends]
    ref=out['recirc']['p'] if 'recirc' in out else jp
    sep=min(ends3,key=lambda e: np.linalg.norm(e-ref)) if 'recirc' in out else max(ends3,key=lambda e: np.linalg.norm(e-jp))
    out['sep']=dict(p=sep)
    vd=g.contour([0.5],scalars='void'); vd=vd.clip(normal=H,origin=O,invert=True) if rear else vd
    vl=largest(vd)
    if vl is not None and vl.n_points>50:
        out['pocket']=dict(p=np.asarray(vl.points).mean(0))
    return out
if __name__=="__main__" and len(sys.argv)>1:
    for k in PH:
        F=feats3d(k); print(PH[k],{kk:(np.round(((v['p']-O)@np.c_[A,E2,H])*1000,1).tolist(),v.get('sense'),round(v.get('size_mm',0),1)) for kk,v in F.items()})
