from engine import mesh
import os

is_view = False
is_simple= False

# export .brep files from FreeCAD
# call gmsh to obtain a mesh
# define physicals groups and save .geo
# obtain physical groups numbers

# simple geometries without physical groups
# filename_FreeCAD = "mesh/FreeCAD/sphere.brep"
# filename_FreeCAD = "mesh/FreeCAD/test_shape.brep"
# filename_FreeCAD = "mesh/FreeCAD/test_shape_labels-GND.brep"
# calling gmsh
# mesh.call_gmsh(filename_FreeCAD, 3, 2.0)

filename_Mesh = "mesh/shape.vtk"

if is_view:
    os.system(f"gmsh {filename_Mesh}")
    pass

elements = mesh.read_mesh(filename_Mesh)
mesh.write_mesh(elements)