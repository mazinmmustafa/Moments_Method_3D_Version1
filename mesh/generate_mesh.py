from engine import mesh
import os

is_view = False

# export .brep files from FreeCAD

# filename_FreeCAD = "mesh/FreeCAD/sphere.brep"
filename_FreeCAD = "mesh/FreeCAD/test_shape.brep"

filename_Mesh = "mesh/shape.vtk"

mesh.call_gmsh(filename_FreeCAD, 3, 0.1)
if is_view:
    os.system(f"gmsh {filename_Mesh}")
    pass
elements = mesh.read_mesh(filename_Mesh)
mesh.write_mesh(elements)