"""Optional STEP surface tessellation for streamline figures (requires Gmsh)."""
from dataclasses import dataclass
from pathlib import Path
from uuid import uuid4

import numpy as np


@dataclass
class FluidDomainMesh:
    vertices: np.ndarray
    triangles: np.ndarray
    border_segments: np.ndarray


def camera_border_collection(vertices, triangles, cad_segments, **style):
    """Combine CAD curves with tessellated apparent contours on every draw.

    Adjacent faces form a silhouette when their opposite corners project to
    the same side of their shared edge. This also tolerates reversed CAD face
    winding. Contours remain visible through the translucent shell.
    """
    from mpl_toolkits.mplot3d.art3d import Line3DCollection
    from mpl_toolkits.mplot3d.proj3d import proj_transform

    adjacency = {}
    for a, b, c in triangles:
        for u, v, opposite in ((a, b, c), (b, c, a), (c, a, b)):
            adjacency.setdefault(tuple(sorted((u, v))), []).append(opposite)
    pairs = [(u, v, corners[0], corners[1])
             for (u, v), corners in adjacency.items() if len(corners) == 2]
    pairs = np.asarray(pairs, dtype=int).reshape(-1, 4)

    class CameraBorders(Line3DCollection):
        def do_3d_projection(self):
            if len(pairs):
                x, y, _ = proj_transform(*vertices.T, self.axes.get_proj())
                projected = np.column_stack((x, y))
                a, b, c, d = (projected[pairs[:, i]] for i in range(4))
                edge = b - a
                side_c = edge[:, 0] * (c - a)[:, 1] - edge[:, 1] * (c - a)[:, 0]
                side_d = edge[:, 0] * (d - a)[:, 1] - edge[:, 1] * (d - a)[:, 0]
                contour = vertices[pairs[side_c * side_d > 0, :2]]
                self.set_segments(np.concatenate((cad_segments, contour), axis=0))
            return super().do_3d_projection()

    return CameraBorders(cad_segments, **style)


def domain_surface_collection(vertices, triangles, shading="flat", **style):
    """Camera-lit, two-sided matte or satin shading for a translucent shell."""
    from matplotlib.colors import to_rgb
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    if shading not in {"flat", "matte", "glossy"}:
        raise ValueError("Domain shading must be flat, matte, or glossy.")
    faces = vertices[triangles]
    normals = np.cross(faces[:, 1] - faces[:, 0], faces[:, 2] - faces[:, 0])
    normals /= np.maximum(np.linalg.norm(normals, axis=1, keepdims=True), 1e-30)
    base = np.array(to_rgb("#9caeb8"))

    class LitSurface(Poly3DCollection):
        def do_3d_projection(self):
            if shading != "flat":
                elev, azim = np.deg2rad([self.axes.elev, self.axes.azim])
                view = np.array([np.cos(elev) * np.cos(azim),
                                 np.cos(elev) * np.sin(azim), np.sin(elev)])
                right = np.array([-np.sin(azim), np.cos(azim), 0.])
                up = np.cross(view, right)
                light = view - 0.6 * right + 0.8 * up
                light /= np.linalg.norm(light)
                facing = normals * np.where(normals @ view >= 0, 1., -1.)[:, None]
                diffuse = np.maximum(facing @ light, 0)
                colors = base * (0.3 + 0.7 * diffuse[:, None])
                if shading == "glossy":
                    halfway = light + view
                    halfway /= np.linalg.norm(halfway)
                    shine = 0.25 * np.maximum(facing @ halfway, 0) ** 24
                    colors += shine[:, None]
                self.set_facecolor(np.clip(colors, 0, 1))
            return super().do_3d_projection()

    return LitSurface(faces, facecolors=base, edgecolors="none", **style)


def load_fluid_domain(path, scale=1.0):
    """Read STEP/STP into metres, preserving CAD curves as outline segments."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Fluid-domain file not found: {path}")
    if path.suffix.lower() not in {".step", ".stp"}:
        raise ValueError("Fluid domain must be a .step or .stp file.")
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("Domain scale must be finite and positive.")
    try:
        import gmsh
    except (ImportError, OSError) as exc:
        raise RuntimeError("STEP overlays require Gmsh: python -m pip install -r requirements-domain.txt") from exc

    owns_session = not gmsh.isInitialized()
    if owns_session:
        gmsh.initialize()
    previous_model = gmsh.model.getCurrent()
    options = {"General.Terminal": 0, "Mesh.ElementOrder": 1,
               "Mesh.RecombineAll": 0, "Mesh.MeshSizeFromCurvature": 12,
               "Mesh.MeshSizeMin": 0, "Mesh.MeshSizeMax": 1}
    old_options = {key: gmsh.option.getNumber(key) for key in options}
    old_unit = gmsh.option.getString("Geometry.OCCTargetUnit")
    added_model = False
    try:
        gmsh.model.add("fluid_domain_" + uuid4().hex)
        added_model = True
        gmsh.option.setString("Geometry.OCCTargetUnit", "M")
        gmsh.model.occ.importShapes(str(path.resolve()), highestDimOnly=False)
        gmsh.model.occ.synchronize()
        bounds = np.array(gmsh.model.getBoundingBox(-1, -1))
        extent = float(np.max(bounds[3:] - bounds[:3]))
        if not np.isfinite(extent) or extent <= 0:
            raise ValueError("Fluid-domain model has no usable geometry.")
        options["Mesh.MeshSizeMin"] = extent / 100
        options["Mesh.MeshSizeMax"] = extent / 30
        for key, value in options.items():
            gmsh.option.setNumber(key, value)
        gmsh.model.mesh.generate(2)
        node_tags, coordinates, _ = gmsh.model.mesh.getNodes()
        vertices = np.asarray(coordinates).reshape(-1, 3) * scale
        order = np.argsort(node_tags)
        sorted_tags = np.asarray(node_tags)[order]

        def indices(tags, width):
            return order[np.searchsorted(sorted_tags, np.asarray(tags))].reshape(-1, width)

        _, triangle_nodes = gmsh.model.mesh.getElementsByType(2)
        triangles = indices(triangle_nodes, 3)
        if not len(triangles):
            raise ValueError("Fluid-domain model contains no meshable surfaces.")
        _, curve_nodes = gmsh.model.mesh.getElementsByType(1)
        segments = vertices[indices(curve_nodes, 2)]
        return FluidDomainMesh(vertices, triangles, segments)
    finally:
        if owns_session:
            gmsh.finalize()
        else:
            if added_model:
                gmsh.model.remove()
            for key, value in old_options.items():
                gmsh.option.setNumber(key, value)
            gmsh.option.setString("Geometry.OCCTargetUnit", old_unit)
            if previous_model:
                gmsh.model.setCurrent(previous_model)
