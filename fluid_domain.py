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
