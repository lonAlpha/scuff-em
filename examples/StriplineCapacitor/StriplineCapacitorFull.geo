//////////////////////////////////////////////////
// GMSH geometry file for a stripline capacitor
// 
// Homer Reid 2/15/2017
//////////////////////////////////////////////////

//////////////////////////////////////////////////
// user-definable constants //////////////////////
//////////////////////////////////////////////////
DefineConstant[ W = 1.0  ];  // trace width
DefineConstant[ L = 10.0 ];  // trace length
DefineConstant[ T = 1.0  ];  // dielectric thickness
DefineConstant[ B = 10.0 ];  // PCB edge length

DefineConstant[ N = 2    ];  // triangles per unit length

DefineConstant[ UseSymmetry = 0 ];

//////////////////////////////////////////////////
// derived constants /////////////////////////////
//////////////////////////////////////////////////
NW   = Ceil[N*W];
NT   = Ceil[N*T];
NBW2 = Ceil[N*0.5*(B-W)];
NB   = Ceil[N*B];
NB2  = Ceil[0.5*N*B];
NL   = Ceil[N*L];

//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////
Point(1) = { -0.5*W, -L/2, T };
Point(2) = {  0.5*W, -L/2, T };
Point(3) = {  0.5*B, -L/2, T };
Point(4) = {  0.5*B, -L/2, 0 };
Point(5) = { -0.5*B, -L/2, 0 };
Point(6) = { -0.5*B, -L/2, T };

// upper trace
Line(1) = {1, 2};
Transfinite Line{1} = NW+1;

// top and side dielectric-air interfaces
Line(2)  = {2, 3};
Transfinite Line{2} = NBW2+1;
Line(3)  = {3, 4};
Transfinite Line{3} = NT+1;
Line(20) = {1,6}; 
Transfinite Line{20} = NBW2+1;
Line(21) = {6,5}; 
Transfinite Line{21} = NT+1;

// ground plane
Line(4) = {4, 5};
Transfinite Line{4} = NB+1;

Extrude(0,L,0)  { Line{1};   Layers{NL}; }
Extrude(0,L,0)  { Line{2,3,20,21}; Layers{NL}; }
Extrude(0,L,0)  { Line{4};   Layers{NL}; }

Extrude(-B,0,0)  { Line{3,30}; Layers{NB}; }
//Translate(0,L,0) { Duplicata{Surface{24};}}

Physical Surface(1) = { 25     };        // trace
Physical Surface(2) = { 29, 33, 37, 41, 49, 53}; //dielectric-air interface
Physical Surface(3) = { 45     };        // ground plane
