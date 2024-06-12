import os
import sys 
import numpy as np

tol = 1.0E-6

class vertex:
    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.x = np.round(x, -int(np.log10(tol)))
        self.y = np.round(y, -int(np.log10(tol)))
        self.z = np.round(z, -int(np.log10(tol)))
        pass
    def __repr__(self):
        str = "({:21.14E}, {:21.14E}, {:21.14E})".format(self.x, self.y, self.z)
        return str

class element:
    def __init__(self):
        self.vertices = []
        pass
    def add_vertex(self, v):
        self.vertices.append(v)
        pass
    def __repr__(self):
        str = f"element {len(self.vertices)-1}d:"
        for vertex in self.vertices:
            str+="\n"+vertex.__repr__()
            continue
        return str

def read_mesh(filename):
    try:
        with open(filename) as my_file:
            text = my_file.read()
        pass
    except:
        print("error: unable to open mesh file!")
        exit(0)
    print(f"reading {filename}..., done!")
    text = text.splitlines()
    vertices = []
    elements = []
    i = 0
    i_break1 = 4
    i_break2 = 6
    N_vertices = 0
    N_elements = 0
    flag1 = True
    flag2 = False
    for line in text:
        if i==i_break1:
            N_vertices = int(line.split()[1])
            print(f"found {N_vertices} vertices..., done!")
            pass
        if flag1 and i>i_break1:
            values = line.split()
            new_vertex = vertex(float(values[0]), float(values[1]), float(values[2]))
            vertices.append(new_vertex)
            pass
        if i==N_vertices+i_break1:
            flag1 = False
        if i==N_vertices+i_break2:
            N_elements = int(line.split()[1])
            print(f"found {N_elements} elements..., done!")
            flag2 = True
            pass
        if flag2 and i>N_vertices+i_break2:
            values = line.split()
            new_element = element()
            N_dimensions = int(values[0])
            for ii in range(1, N_dimensions+1):
                x = vertices[int(values[ii])].x
                y = vertices[int(values[ii])].y
                z = vertices[int(values[ii])].z
                new_vertex = vertex(x, y, z)
                new_element.add_vertex(new_vertex)
                continue
            elements.append(new_element)
            pass
        if flag2 and i==N_vertices+N_elements+i_break2:
            break
        i+=1
        continue
    return elements

def call_gmsh(filename, opt, res):
    # opt: 2 or 3 for triangles or tetraherons
    output = "mesh/shape.vtk"
    try:
        assert(opt==2 or opt==3)
        pass
    except:
        print("error: invalid meshing option")
        exit(1)
    try:
        assert(res>0.0)
        pass
    except:
        print("error: invalid meshing resolution")
        exit(1)
    os.system(f"gmsh {filename} -clmax {res} -{opt} -format vtk -save_all -o {output}")

def write_mesh(elements):
    count_0d = 0
    count_1d = 0
    count_2d = 0
    count_3d = 0
    try:
        file_0d = open("mesh/mesh/mesh_0d.txt", "w")
        file_1d = open("mesh/mesh/mesh_1d.txt", "w")
        file_2d = open("mesh/mesh/mesh_2d.txt", "w")
        file_3d = open("mesh/mesh/mesh_3d.txt", "w")
        file_data = open("mesh/mesh/mesh_data.txt", "w")
        for element in elements:
            # 0d
            if len(element.vertices)==1:
                file_0d.write("{:21.14E}, ".format(element.vertices[0].x))
                file_0d.write("{:21.14E}, ".format(element.vertices[0].y))
                file_0d.write("{:21.14E}\n".format(element.vertices[0].z))
                file_0d.write("\n")
                count_0d+=1
                pass
            # 1d
            if len(element.vertices)==2:
                file_1d.write("{:21.14E}, ".format(element.vertices[0].x))
                file_1d.write("{:21.14E}, ".format(element.vertices[0].y))
                file_1d.write("{:21.14E}\n".format(element.vertices[0].z))
                file_1d.write("{:21.14E}, ".format(element.vertices[1].x))
                file_1d.write("{:21.14E}, ".format(element.vertices[1].y))
                file_1d.write("{:21.14E}\n".format(element.vertices[1].z))
                file_1d.write("\n")
                count_1d+=1
                pass
            # 2d
            if len(element.vertices)==3:
                file_2d.write("{:21.14E}, ".format(element.vertices[0].x))
                file_2d.write("{:21.14E}, ".format(element.vertices[0].y))
                file_2d.write("{:21.14E}\n".format(element.vertices[0].z))
                file_2d.write("{:21.14E}, ".format(element.vertices[1].x))
                file_2d.write("{:21.14E}, ".format(element.vertices[1].y))
                file_2d.write("{:21.14E}\n".format(element.vertices[1].z))
                file_2d.write("{:21.14E}, ".format(element.vertices[2].x))
                file_2d.write("{:21.14E}, ".format(element.vertices[2].y))
                file_2d.write("{:21.14E}\n".format(element.vertices[2].z))
                file_2d.write("\n")
                count_2d+=1
                pass
            # 3d
            if len(element.vertices)==4:
                file_3d.write("{:21.14E}, ".format(element.vertices[0].x))
                file_3d.write("{:21.14E}, ".format(element.vertices[0].y))
                file_3d.write("{:21.14E}\n".format(element.vertices[0].z))
                file_3d.write("{:21.14E}, ".format(element.vertices[1].x))
                file_3d.write("{:21.14E}, ".format(element.vertices[1].y))
                file_3d.write("{:21.14E}\n".format(element.vertices[1].z))
                file_3d.write("{:21.14E}, ".format(element.vertices[2].x))
                file_3d.write("{:21.14E}, ".format(element.vertices[2].y))
                file_3d.write("{:21.14E}\n".format(element.vertices[2].z))
                file_3d.write("{:21.14E}, ".format(element.vertices[3].x))
                file_3d.write("{:21.14E}, ".format(element.vertices[3].y))
                file_3d.write("{:21.14E}\n".format(element.vertices[3].z))
                file_3d.write("\n")
                count_3d+=1
                pass
            # 
            continue
        file_0d.close()
        file_1d.close()
        file_2d.close()
        file_3d.close()
        file_data.write("{:d}\n".format(count_0d))
        file_data.write("{:d}\n".format(count_1d))
        file_data.write("{:d}\n".format(count_2d))
        file_data.write("{:d}\n".format(count_3d))
        file_data.close()
        print("meshing completed successfully")
        print(f"found {count_0d} 0d elements")
        print(f"found {count_1d} 1d elements")
        print(f"found {count_2d} 2d elements")
        print(f"found {count_3d} 3d elements")
        pass
    except:
        print("error: unable to open file!")
        exit(0)
    pass