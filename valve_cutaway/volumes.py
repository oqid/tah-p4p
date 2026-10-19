import sys; sys.path.insert(0,'.')
from classify import *
import pyvista as pv
from scipy.ndimage import gaussian_filter
def volumes(k,h=0.6e-3,sigma=1.3,amin=-4e-3):
    fa,ed=geom(k); X,tau,U,off=lines(k); T=tangents(X,off)
    surf=facemesh(fa,list(fa)).triangulate().clean()
    lo=surf.points.min(0); hi=surf.points.max(0)
    dims=np.ceil((hi-lo)/h).astype(int)+1
    def idx(P): return np.clip(np.round((P-lo)/h).astype(int),0,dims-1)
    keep=((X-O)@A)>amin-3e-3
    I=idx(X[keep]); ua=(T[keep]@A)*U[keep]
    W=np.zeros(dims); S=np.zeros(dims); Ts=np.zeros(dims)
    np.add.at(W,tuple(I.T),1.0); np.add.at(S,tuple(I.T),ua); np.add.at(Ts,tuple(I.T),np.log10(np.clip(tau[keep],1e-2,1e2)))
    Wg=gaussian_filter(W,sigma); Sg=gaussian_filter(S,sigma)/np.maximum(Wg,1e-9)
    g=pv.ImageData(dimensions=dims,spacing=(h,h,h),origin=lo)
    pts=g.points
    ins=pv.PolyData(pts).select_enclosed_points(surf,tolerance=1e-6,check_surface=False)['SelectedPoints'].astype(bool)
    a=(pts-O)@A
    region=ins&(a>amin)
    Wf=Wg.ravel(order='F'); Sf=Sg.ravel(order='F')
    thr=0.03*np.percentile(Wf[region&(Wf>0)],90)
    rev=np.where(region&(Wf>thr),-Sf,-1.0)     # >0 where axial flow reversed
    void=np.where(region,(Wf<thr).astype(float),0.0)
    g['rev']=rev; g['void']=gaussian_filter(void.reshape(dims,order='F'),1.0).ravel(order='F')
    return g
if __name__=="__main__":
    import pickle
    for k in PH:
        g=volumes(k); g.save(f"vol{k}.vtk")
        rv=g.contour([0.03],scalars='rev'); vd=g.contour([0.5],scalars='void')
        print(PH[k],"reversed vol mm3",round(rv.volume*1e9 if rv.n_points else 0,1) , "void-iso pts",vd.n_points)
