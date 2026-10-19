import sys; sys.path.insert(0,'.')
from section import *
OUT=np.array([20.24,30.52,13.84])/1000
VALVE_ANGLE={16:24.9,33:49.8,50:75.0,66:50.1,83:25.2}   # operating_point.py (as in F3)
def classify(k):
    X,tau,U,off=lines(k); T=tangents(X,off); P,d=to2d(X)
    n=len(off)-1; rec=np.zeros(n,bool); peak=np.zeros(n); wall=np.zeros(n,bool)
    for i,(s,e) in enumerate(zip(off[:-1],off[1:])):
        a=P[s:e,0]; ta=T[s:e]@A
        r=(a>1)&(ta<-0.2)
        if r.sum()>=2:
            seg=np.diff(np.r_[0,r.astype(int),0]); st=np.where(seg==1)[0]; en=np.where(seg==-1)[0]
            rec[i]=(en-st).max()>=2
        peak[i]=tau[s:e].max()
        wall[i]=((X[e-1]-OUT)@A)<-0.002
    hi=peak>=np.percentile(peak,95)
    return dict(X=X,tau=tau,U=U,off=off,rec=rec,hi=hi,wall=wall,peak=peak)
if __name__=="__main__":
    for k in PH:
        c=classify(k); print(PH[k],"recirc",c['rec'].sum(),"hi",c['hi'].sum(),"both",(c['rec']&c['hi']).sum(),"wall-ended",c['wall'].sum(), "P95 peak",np.percentile(c['peak'],95).round(1))
