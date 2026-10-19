"""Mid-plane fields in a Y-up basis (outflow up, as in the CAD drawings), incl. GW HI3 dose."""
import sys; sys.path.insert(0,'.')
from fields import *
GW=dict(A=3.62e-5,alpha=2.416,beta=0.785)
UP=np.array([0,1.0,0]); UP=UP-(UP@H)*H; UP/=np.linalg.norm(UP)
RT=np.cross(UP,H); RT*= np.sign(RT[0]) if RT[0]!=0 else 1
def to2y(P3): Q=P3-O; return np.c_[Q@RT,Q@UP]*1000,(Q@H)*1000
def node_dose(X,U,off):
    """Linearised GW HI3 dose per node, as explore_core.inner (dt from the incoming segment, velocity floor)."""
    n=len(off)-1; lens=np.diff(off); pid=np.repeat(np.arange(n),lens)
    first=np.zeros(len(X),bool); first[off[:-1]]=True
    sl=np.r_[0,np.linalg.norm(np.diff(X,axis=0),axis=1)]; sl[first]=0
    sv=np.r_[0,0.5*(U[1:]+U[:-1])]; sv[first]=0
    pos=U>0; mv=np.bincount(pid,np.where(pos,U,0),n)/np.maximum(np.bincount(pid,pos.astype(float),n),1)
    floor=np.maximum(0.01*mv,0.001)[pid]; sv=np.maximum(sv,floor)
    dt=np.where(first,0,sl/sv)
    return dt
def phase_field_y(k,h=0.5,sigma=1.4,xl=(-22,22),yl=(-14,30)):
    fa,ed=geom(k); X,tau,U,off=lines(k); T=tangents(X,off)
    dt=node_dose(X,U,off)
    inner=GW['A']**(1/GW['beta'])*np.clip(tau,0,None)**(GW['alpha']/GW['beta'])*dt
    P,d=to2y(X); sel=np.abs(d)<SLAB; w=np.exp(-0.5*(d/(SLAB/1.5))**2)
    gx=np.arange(xl[0],xl[1]+h,h); gy=np.arange(yl[0],yl[1]+h,h)
    V=T*U[:,None]
    (sp,ux,uy,ua,lt),Ws=binned(P[sel],[U[sel],V[sel]@RT,V[sel]@UP,V[sel]@A,np.log10(np.clip(tau[sel],1e-3,None))],w[sel],gx,gy,sigma)
    # dose: plain sums (top-hat slab), smoothed, as share of the phase's total dose per mm^2
    ix=np.clip(np.round((P[sel,0]-gx[0])/h).astype(int),0,len(gx)-1); iy=np.clip(np.round((P[sel,1]-gy[0])/h).astype(int),0,len(gy)-1)
    D=np.zeros((len(gy),len(gx))); np.add.at(D,(iy,ix),inner[sel]); D=gaussian_filter(D,sigma)/(h*h)/inner.sum()
    # mask in this basis
    surf=facemesh(fa,list(fa)).triangulate().clean()
    GX,GY=np.meshgrid(gx,gy); P3=O+(GX.ravel()[:,None]*RT+GY.ravel()[:,None]*UP)/1000
    m=pv.PolyData(P3).select_enclosed_points(surf,tolerance=1e-6,check_surface=False)['SelectedPoints'].reshape(GX.shape).astype(bool)
    def sec(ids):
        s=facemesh(fa,ids).slice(normal=H,origin=O)
        if s.n_points==0: return []
        s=s.strip(); return [to2y(np.asarray(s.get_cell(i).points))[0] for i in range(s.n_cells)]
    return dict(gx=gx,gy=gy,speed=sp,ua=ux,ub=uy,uax=ua,ltau=lt,dose=D,W=Ws,mask=m,
                wall=sec([t for t in fa if t not in DISC]),disc=sec(DISC),slab_share=inner[sel].sum()/inner.sum())
if __name__=="__main__":
    import pickle
    for k in PH:
        f=phase_field_y(k); pickle.dump(f,open(f"fieldy{k}.pkl","wb"))
        print(PH[k],"slab dose share %.3f"%f['slab_share'], "max dose/mm2 %.2e"%np.nanmax(f['dose']))
