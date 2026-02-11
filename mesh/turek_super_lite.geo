// turek_smooth_1k.geo (Target: ~1000 elementi con transizione dolce)

// 1. Parametri Geometrici
L = 2.2;
H = 0.41;
cx = 0.2;
cy = 0.2;
r = 0.05;

// --- NUOVI PARAMETRI DI DENSITÀ ---
lc_far = 0.12;       // Ridotto (era 0.18): evita triangoli troppo grandi ai bordi
lc_cylinder = 0.02;  // Aumentato (era 0.012): risparmiamo elementi qui
lc_wake = 0.045;     // Scia: via di mezzo per mantenere la scia fluida

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
Physical Curve(0) = {4};         
Physical Curve(1) = {2};         
Physical Curve(2) = {1, 3};      
Physical Curve(3) = {5, 6, 7, 8}; 
Physical Surface(0) = {1};

// =====================================================
// CAMPI DI RAFFINAMENTO (MASSIMA GRADUALITÀ)
// =====================================================

Field[1] = Distance;
Field[1].CurvesList = {5, 6, 7, 8};

Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = lc_cylinder; 
Field[2].LcMax = lc_far;      
Field[2].DistMin = 0.05;      
Field[2].DistMax = 0.4;       // Transizione spalmata per evitare l'effetto "esplosione"

Field[3] = Box;
Field[3].XMin = 0.1; 
Field[3].XMax = 2.2;
Field[3].YMin = 0.05;         // Allargata quasi a tutto il canale
Field[3].YMax = 0.35;         // Allargata quasi a tutto il canale
Field[3].VIn  = lc_wake;
Field[3].VOut = lc_far;
Field[3].Thickness = 0.2;     

Field[4] = Min;
Field[4].FieldsList = {2, 3};

Background Field = 4;

// Cruciale per evitare che i punti agli angoli comandino la mesh
Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.CharacteristicLengthFromPoints = 0;
Mesh.CharacteristicLengthFromCurvature = 0;

Mesh.Algorithm = 6;