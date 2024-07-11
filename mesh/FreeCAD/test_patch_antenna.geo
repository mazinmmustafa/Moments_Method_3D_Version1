Merge "test_patch_antenna.brep";
//+
MeshSize {27, 28, 29, 30, 31, 32, 21, 22, 23, 24, 25, 26} = 4.0; // Patch
//+
MeshSize {6, 8, 5, 7, 2, 4, 1, 3} = 6.0; // Substrate Corners
//+
MeshSize {15, 16, 17, 18, 19, 20, 9, 10, 11, 12, 13, 14} = 6.0; // GND Patch
//+
MeshSize {20, 32, 21, 9} = 4.0; // Port
//+
Physical Surface("Patch", 1) = {9};
//+
Physical Surface("GND", 2) = {8, 4};
//+
Physical Surface("Port", 3) = {10};
//+
Physical Volume("Substrate", 4) = {1};
