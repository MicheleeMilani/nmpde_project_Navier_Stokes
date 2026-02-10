// turek_lite.geo (Versione alleggerita del 30%)

// 1. Parametri Geometrici
L = 2.2;
H = 0.41;
cx = 0.2;
cy = 0.2;
r = 0.05;

// --- PARAMETRI DI DENSITÀ (ALLEGGERITI) ---
// Più alti sono questi numeri, meno elementi avrai.
lc_far = 0.08;       // Prima era 0.05 -> Sfondo più leggero
lc_cylinder = 0.003; // Prima era 0.002 -> Meno punti sul cerchio
lc_wake = 0.006;     // Prima era 0.004 -> Scia meno fitta

// --- GEOMETRIA ---
Point(1) = {0, 0, 0, lc_far};
Point(2) = {L, 0, 0, lc_far};
Point(3) = {L, H, 0, lc_far};
Point(4) = {0, H, 0, lc_far};

Point(5) = {cx, cy, 0, lc_cylinder};
Point(6) = {cx-r, cy, 0, lc_cylinder};
Point(7) = {cx, cy-r, 0, lc_cylinder};
Point(8) = {cx+r, cy, 0, lc_cylinder};
Point(9) = {cx, cy+r, 0, lc_cylinder};

Line(1) = {1, 2}; Line(2) = {2, 3}; Line(3) = {3, 4}; Line(4) = {4, 1};
Circle(5) = {6, 5, 7}; Circle(6) = {7, 5, 8}; Circle(7) = {8, 5, 9}; Circle(8) = {9, 5, 6};
Curve Loop(1) = {1, 2, 3, 4}; Curve Loop(2) = {5, 6, 7, 8};
Plane Surface(1) = {1, 2};

// IDs per Deal.II
Physical Curve(0) = {4};       // Inlet
Physical Curve(1) = {2};       // Outlet
Physical Curve(2) = {1, 3};    // Walls
Physical Curve(3) = {5, 6, 7, 8}; // Cylinder
Physical Surface(0) = {1};

// =====================================================
// CAMPI DI RAFFINAMENTO (FIELDS)
// =====================================================

// --- 1. Gradiente attorno al CILINDRO ---
Field[1] = Distance;
Field[1].CurvesList = {5, 6, 7, 8};
Field[1].Sampling = 100;

Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = lc_cylinder; 
Field[2].LcMax = lc_far;      
Field[2].DistMin = 0.05;      // Zona densa vicinissima al cilindro
Field[2].DistMax = 0.25;      // Transizione fino a questa distanza (ridotto leggermente)

// --- 2. Scatola per la SCIA (Wake) ---
Field[3] = Box;
Field[3].XMin = 0.15; 
Field[3].XMax = 2.2;
Field[3].YMin = 0.12;
Field[3].YMax = 0.28;
Field[3].VIn  = lc_wake; // Qui usiamo il valore alleggerito (0.006)
Field[3].VOut = lc_far;  

// Mantengo la transizione morbida
Field[3].Thickness = 0.2; 

// --- 3. Combinazione (Min) ---
Field[4] = Min;
Field[4].FieldsList = {2, 3};

Background Field = 4;

// Disattiviamo l'interpolazione dai bordi
Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.CharacteristicLengthFromPoints = 0;
Mesh.CharacteristicLengthFromCurvature = 0;
