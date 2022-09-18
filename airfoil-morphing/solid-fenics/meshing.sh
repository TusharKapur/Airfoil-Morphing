#!/bin/sh
gmsh mesh.geo -2 mesh.msh -format msh2
dolfin-convert mesh.msh mesh.xml
python3 check_convert.py
