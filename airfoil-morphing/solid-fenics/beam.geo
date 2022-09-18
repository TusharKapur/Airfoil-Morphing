Geometry.OldNewReg=0;
General.ExpertMode=1;
size_beam = 1.5E-01;
Point(1) = {0, 1, 0, size_beam};
Point(2) = {0, 0, 0, size_beam};
Point(3) = {4, 0, 0, size_beam};
Point(4) = {4, 1, 0, size_beam};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Physical Curve("left", 1) = {1};
//+
Physical Curve("up", 2) = {4};
//+
Physical Curve("down", 3) = {2};
//+
Physical Curve("right", 4) = {3};
//+
Physical Surface("beam", 5) = {1};
