// SetFactory("OpenCASCADE");
length = 0.1;
inner_rad = 0.01/2;
outer_rad = 0.015/2;
DefineConstant[Nquadrant = {3, Min 1, Step 1, Name "Nquadrant"}];
DefineConstant[Nradial = {2, Min 1, Step 1, Name "Nradial"}];
DefineConstant[Nlength = {30, Min 1, Step 1, Name "Nlength"}];

For i In {0:3}
    theta = i * Pi / 2;
    Point(2*i+1) = {inner_rad*Cos(theta), inner_rad*Sin(theta), -length/2};
    Point(2*i+2) = {outer_rad*Cos(theta), outer_rad*Sin(theta), -length/2};
EndFor
center = newp; Point(center) = {0, 0, -length/2};
For i In {0:3}
    Circle(3*i+1) = {2*i+1, center, 2*((i+1)%4)+1}; // inner circle
    Circle(3*i+2) = {2*i+2, center, 2*((i+1)%4)+2}; // outer circle
    Line(3*i+3) = {2*i+1, 2*i+2};
EndFor
For i In {0:3}
    Transfinite Curve {3*i+1, 3*i+2} = Nquadrant + 1;
    Transfinite Curve {3*i+3} = Nradial + 1;
    Curve Loop(i+1) = {3*i+3, 3*i+2, -(3*((i+1)%4)+3), -(3*i+1)};
    Plane Surface(i+1) = {i+1};
EndFor
cyl[] = Extrude {0, 0, length} {Surface {1:4}; Layers{Nlength+1}; Recombine;};
Transfinite Surface {1:4};
Transfinite Surface {cyl[0], cyl[6], cyl[12], cyl[18]};
Printf("New cyl: %g %g %g %g %g %g %g %g", cyl[0], cyl[1], cyl[2], cyl[3], cyl[4], cyl[5], cyl[6], cyl[7]);

Transfinite Volume {:};
Recombine Surface {:};

bottom = 998;
top = 999;
Physical Surface(bottom) = {1:4};
Physical Surface(top) = {cyl[0], cyl[6], cyl[12], cyl[18]};
Physical Volume(1) = {1:4};