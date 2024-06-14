Merge "test_patch_antenna.brep";
//+
MeshSize {27, 28, 29, 30, 31, 32, 21, 22, 23, 24, 25, 26} = 0.6;
//+
MeshSize {6, 8, 5, 7, 2, 4, 1, 3} = 2.0;
//+
MeshSize {15, 16, 17, 18, 19, 20, 9, 10, 11, 12, 13, 14} = 2.0;
//+
MeshSize {20, 32, 21, 9} = 0.6;
//+
Physical Surface("Patch", 1) = {9};
//+
Physical Surface("GND", 2) = {8, 4};
//+
Physical Surface("Port", 3) = {10};
//+
Physical Volume("Substrate", 4) = {1};
