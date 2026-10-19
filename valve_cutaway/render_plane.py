import sys, pickle; sys.path.insert(0,'.')
from render_v2 import *
from fields_y import RT,UP
import matplotlib.pyplot as plt
def plane_mesh(f,vals):
    gx,gy=f['gx'],f['gy']; GX,GY=np.meshgrid(gx,gy)
    P3=O+(GX[...,None]*RT+GY[...,None]*UP)/1000
    g=pv.StructuredGrid(P3[...,0].T,P3[...,1].T,P3[...,2].T)
    minW=0.02*np.percentile(f['W'][f['mask']],90); have=f['mask']&(f['W']>minW)
    g['v']=np.where(have,vals,np.nan).T.ravel(order='F'); g['m']=f['mask'].T.ravel(order='F').astype(float)
    return g
def stream_lines3d(f,density=1.4):
    gx,gy=f['gx'],f['gy']; minW=0.02*np.percentile(f['W'][f['mask']],90); have=f['mask']&(f['W']>minW)
    U_=np.where(have,f['ua'],np.nan); V_=np.where(have,f['ub'],np.nan); s=np.hypot(U_,V_)
    fig,ax=plt.subplots(); st=ax.streamplot(gx,gy,U_/np.maximum(s,1e-9),V_/np.maximum(s,1e-9),density=density,minlength=0.08,maxlength=4); plt.close(fig)
    segs=st.lines.get_segments()
    pts=[];cells=[];o=0
    for sg in segs:
        sg=np.asarray(sg); P3=O+(sg[:,0:1]*RT+sg[:,1:2]*UP)/1000+H*3e-5
        pts.append(P3); cells.append(np.r_[len(P3),np.arange(o,o+len(P3))]); o+=len(P3)
    return pv.PolyData(np.vstack(pts),lines=np.concatenate(cells))
def render_plane(k,fn,camera,which='stress',window=(1400,1400)):
    fa,ed=geom(k); f=pickle.load(open(f"fieldy{k}.pkl","rb"))
    clip=clipper(-16e-3,None); clipA=lambda m: m.clip(normal=-A,origin=O-16e-3*A,invert=True)
    p=pv.Plotter(off_screen=True,window_size=window); p.set_background('white')
    p.add_mesh(clip(facemesh(fa,SHELL)),color='#C9D3DA',opacity=0.20,smooth_shading=True)
    p.add_mesh(clip(edges_mesh(ed)),color='#4A4A4A',line_width=1.1,opacity=0.6)
    g=plane_mesh(f,10**f['ltau'] if which=='stress' else f['dose'])
    inside=g.threshold(0.5,scalars='m')
    p.add_mesh(clipA(inside),color=VOID_3D,ambient=1.0,diffuse=0.0,opacity=0.30)
    data=g.threshold(value=(-1e30,1e30),scalars='v').translate(H*1.5e-5,inplace=False)   # drops NaN
    if which=='stress':
        p.add_mesh(clipA(data),scalars='v',cmap='magma',log_scale=True,clim=(1e-2,1e2),show_scalar_bar=False,ambient=1.0,diffuse=0.0,nan_opacity=0,opacity=OP)
        col='white'
    else:
        p.add_mesh(clipA(data),scalars='v',cmap='magma_r',log_scale=True,clim=(1e-6,1e-2),show_scalar_bar=False,ambient=1.0,diffuse=0.0,nan_opacity=0,opacity=OP)
        col='#4A4A4A'
    p.add_mesh(clipA(stream_lines3d(f)),color=col,line_width=1.8,opacity=0.85)
    p.add_mesh(facemesh(fa,DISC),color='#56636D',smooth_shading=True,specular=0.4,ambient=0.28)
    p.add_mesh(clip(facemesh(fa,list(PIVOT))),color='#8C979F',smooth_shading=True,ambient=0.25)
    p.camera_position=camera; p.enable_anti_aliasing('ssaa'); p.screenshot(fn)
VOID_3D='#D9DEE2'
OP=0.72
if __name__=="__main__":
    k=int(sys.argv[1]); which=sys.argv[2]; vname=sys.argv[3] if len(sys.argv)>3 else 'oblique_plane'
    VIEWS['oblique_plane']=(H+0.45*A-0.55*E2,(7,0),0.085,-16,None)
    camera,_=view(vname)
    render_plane(k,f"plane_{which}_{vname}_{k}.png",camera,which)
