// From GMSH tutorials

lc = 0.005;

Point(1) = {-0.025,  0.025, -0.06, lc};
Point(2) = {-0.025, -0.025, -0.06, lc};
Point(3) = { 0.025, -0.025, -0.06, lc};
Point(4) = { 0.025,  0.025, -0.06, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};

Plane Surface(1) = {1};

// Extrude towards the z direction
// Save extruded elements in v[]
v[] = Extrude {0, 0, 0.05} { Surface{1}; };

Physical Volume(111) = {v[1]};
Physical Volume('tets') = {111};
Physical Volume(1) = {v[1]};


Mesh.Algorithm3D = 1;
Mesh.OptimizeThreshold = 0.4;
Mesh.MeshSizeMin = 0.002;
Mesh.MeshSizeMax = 0.0035;
// Mesh.Smoothing = 100;

// Needs MMG:
// RefineMesh;

Mesh 3;
