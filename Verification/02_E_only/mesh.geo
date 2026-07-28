Point(1) = {0, 0, 0, 0.005};
Point(2) = {1, 0, 0, 0.005};
Point(3) = {1, 1, 0, 0.005};
Point(4) = {0, 1, 0, 0.005};

Line(1) = {4, 3};
Line(2) = {3, 2};
Line(3) = {2, 1};
Line(4) = {1, 4};

Curve Loop(1) = {-3, -2, -1, -4};

Plane Surface(1) = {1};
Physical Curve("wall", 5) = {4, 1, 2, 3};
Physical Surface("plasma", 6) = {1};
