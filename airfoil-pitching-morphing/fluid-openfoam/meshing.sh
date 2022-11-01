#!/bin/sh
gmsh mesh.geo -3 mesh.msh -format msh2
gmshToFoam mesh.msh
changeDictionary
checkMesh
