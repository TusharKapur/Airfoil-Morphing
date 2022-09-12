Geometry.OldNewReg=0;
General.ExpertMode=1;
size_beam = 1.5E-00;
size_wall = 1.5E-00;
Point(1) = {0, 1, 0, size_beam};
Point(2) = {0, 0, 0, size_beam};
Point(3) = {4, 0, 0, size_beam};
Point(4) = {4, 1, 0, size_beam};

Point(5) = {-20, 10, 0, size_beam};
Point(6) = {-20, -10, 0, size_beam};
Point(7) = {30, -10, 0, size_beam};
Point(8) = {30, 10, 0, size_beam};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8};

Plane Surface(1) = {1};
//+
Extrude {0, 0, 1} {
  Surface{1}; Layers {1}; Recombine;
}
//+
Physical Surface("inlet", 1) = {7};
//+
Physical Surface("outlet", 2) = {9};
//+
Physical Surface("flap", 4) = {3, 6, 5, 4};
//+
Physical Surface("frontAndBack", 5) = {1, 11};
//+
Physical Surface("upperWall", 6) = {10};
//+
Physical Surface("lowerWall", 7) = {8};
//+
Physical Volume("domain", 8) = {1};
