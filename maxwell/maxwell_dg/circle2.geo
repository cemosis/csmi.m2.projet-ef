h = 0.1;
r=1;
r2=0.3;
Point(1) = {0,0,0,h};
Point(2) = {r,0,0,h};
Point(3) = {0,r,0,h};
Point(4) = {-r,0,0,h};
Point(5) = {0,-r,0,h};

Point(22) = {r2,0,0,h};
Point(33) = {0,r2,0,h};
Point(44) = {-r2,0,0,h};
Point(55) = {0,-r2,0,h};



Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};

Circle(11) = {22,1,33};
Circle(22) = {33,1,44};
Circle(33) = {44,1,55};
Circle(44) = {55,1,22};


Line Loop(5) = {2, 3, 4, 1};

Line Loop(6) = {22, 33, 44, 11};


Plane Surface(1) = {5,-6};
//Plane Surface(2) = {6};

Physical Line("Dirichlet") = {2,3,4,5}
Physical Line("Metal") = {22,33,44,55}
