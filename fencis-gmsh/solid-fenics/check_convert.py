
from fenics import *

xml_file = "mesh.xml"
mesh = Mesh(xml_file)

U = FunctionSpace(mesh, 'P', 1)
u = Function(U)

vtkfile = File('test.pvd')
vtkfile << (u)

cd = MeshFunction('size_t', mesh, "mesh_physical_region.xml");
fd = MeshFunction('size_t', mesh, "mesh_facet_region.xml");
