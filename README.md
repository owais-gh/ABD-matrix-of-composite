# ABD-matrix-of-composite
Calculating ABD matrix for single material composite with equal thickness layers, works for any number of layers and orientation angle


Possible error messages and their fixes

1.
{Index exceeds array bounds.
Error in ABD_Matrix_calculation_for_composite_materials (line 49)
    c=cosd(theta(k));} 
=Number of layers(n)>Number of layer orientations entered

2.
Missing values for 1 or more layes=Number of layers(n)<Number of layer orientations entered
