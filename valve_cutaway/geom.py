import gmsh, numpy as np, sys
def tess(f, size=0.0006):
    gmsh.initialize(); gmsh.option.setNumber("General.Terminal",0)
    gmsh.option.setString("Geometry.OCCTargetUnit","M")
    gmsh.model.occ.importShapes(f); gmsh.model.occ.synchronize()
    gmsh.option.setNumber("Mesh.MeshSizeMax",0.0016); gmsh.option.setNumber("Mesh.MeshSizeMin",0.0004)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature",16)
    gmsh.model.mesh.generate(2)
    faces={}
    for d,t in gmsh.model.getEntities(2):
        nt,nc,_=gmsh.model.mesh.getNodes(2,t,includeBoundary=True)
        et,etags,en=gmsh.model.mesh.getElements(2,t)
        tri=np.array(en[0]).reshape(-1,3)
        m={n:i for i,n in enumerate(nt)}
        tri=np.vectorize(m.get)(tri)
        faces[t]=dict(V=np.array(nc).reshape(-1,3),F=tri,type=gmsh.model.getType(2,t),area=gmsh.model.occ.getMass(2,t))
    edges=[]
    for d,t in gmsh.model.getEntities(1):
        nt,nc,_=gmsh.model.mesh.getNodes(1,t,includeBoundary=True)
        et,etags,en=gmsh.model.mesh.getElements(1,t)
        if not len(en): continue
        m={n:i for i,n in enumerate(nt)}; V=np.array(nc).reshape(-1,3)
        seg=np.vectorize(m.get)(np.array(en[0]).reshape(-1,2))
        up=[s for s in gmsh.model.getAdjacencies(1,t)[0]]
        edges.append(dict(V=V,S=seg,faces=up))
    gmsh.finalize()
    return faces,edges
if __name__=="__main__":
    import pickle
    names={16:"16.67",33:"33",50:"50",66:"66.6",83:"83.3"}
    for k,n in names.items():
        import os
        stp=os.environ.get("STP_DIR", r"E:\35cc Retry\3D 35cc STP Files")
        fa,ed=tess(os.path.join(stp,f"3D 35cc {n}% systole.stp"))
        pickle.dump((fa,ed),open(f"geom{k}.pkl","wb"))
        print(k,len(fa),sum(len(v['F']) for v in fa.values()), [ (t,v['type'],round(v['area']*1e6,1)) for t,v in fa.items() if 150<v['area']*1e6<260 and v['type']=='Plane'])
