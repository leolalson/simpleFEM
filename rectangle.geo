cl__1 = 0.005;
meshsize = 0.04;
p1 = 0.3;
p2 = 0.120;
Point(1) = {0, 0, 0, cl__1};
Point(2) = {p1, 0, 0, cl__1};
Point(3) = {p1, p2, 0, cl__1};
Point(4) = {0, p2, 0, cl__1};
Line(5) = {1, 2};
Transfinite Line {5} = (p1/meshsize) Using Progression 1;
Line(6) = {2, 3};
Transfinite Line {6} = (p2/meshsize) Using Progression 1;
Line(7) = {3, 4};
Transfinite Line {7} = (p1/meshsize) Using Progression 1;
Line(8) = {4, 1};
Transfinite Line {8} = (p2/meshsize) Using Progression 1;
Curve Loop(10) = {5, 6, 7, 8};
Plane Surface(10) = {10};
Transfinite Surface {10};//+
Physical Surface(5) = {10};
//+
Physical Curve(1) = {5};
//+
Physical Curve(2) = {6};
//+
Physical Curve(3) = {7};//+
Physical Curve(4) = {8};