// turek_smooth.geo

// 1. Parametri Geometrici
L = 2.2;
H = 0.41;
cx = 0.2;
cy = 0.2;
r = 0.05;

// Dimensioni base
lc_far = 0.05;      // Dimensione celle lontano (sfondo)
lc_cylinder = 0.002; // Dimensione MOLTO fine sul bordo del cilindro
lc_wake = 0.004;     // Dimensione fine dentro la scia

// --- GEOMETRIA (Punti e Linee) ---
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
// STRATEGIA DI RAFFINAMENTO GRADUALE
// =====================================================

// --- STEP 1: Gradiente attorno al CILINDRO ---
// Calcoliamo la distanza dalle curve del cilindro (5,6,7,8)
Field[1] = Distance;
Field[1].CurvesList = {5, 6, 7, 8};
Field[1].Sampling = 100;

// Creiamo un gradiente basato sulla distanza
Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = lc_cylinder; // Dimensione a contatto col cilindro
Field[2].LcMax = lc_far;      // Dimensione lontano
Field[2].DistMin = 0.05;      // Fino a questa distanza, mantieni lc_cylinder
Field[2].DistMax = 0.2;       // Da 0.05 a 0.2, cresci gradualmente fino a lc_far

// --- STEP 2: Scatola per la SCIA (Wake) ---
Field[3] = Box;
Field[3].XMin = 0.15; 
Field[3].XMax = 2.2;
Field[3].YMin = 0.12;
Field[3].YMax = 0.28;
Field[3].VIn  = lc_wake; // Dimensione dentro la scia
Field[3].VOut = lc_far;  // Dimensione fuori dalla scia

// *** IL TRUCCO PER LA GRADUALITÀ ***
// Thickness definisce quanto è "larga" la zona di sfumatura tra VIn e VOut.
// Prima avevi 0.08 (troppo poco). Mettiamo 0.15 o 0.2 per sfumare dolcemente.
Field[3].Thickness = 0.2; 

// --- STEP 3: Combinazione Finale (Min) ---
// Prendi il minimo tra il gradiente del cilindro (Field 2) e la scatola (Field 3)
Field[4] = Min;
Field[4].FieldsList = {2, 3};

// Impostiamo il campo risultante come background
Background Field = 4;

// Disattiviamo l'interpolazione dai bordi per far comandare i Campi
Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.CharacteristicLengthFromPoints = 0;
Mesh.CharacteristicLengthFromCurvature = 0;
