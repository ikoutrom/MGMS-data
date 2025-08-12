This set of scripts runs the analyses for the 2D mesh without a mass on the top, subjected to a harmonic horizontal force at the top.

There are two separate master scripts:

- VMS_2DharmonicPsingleMesh.m: this runs the case with a SINGLE MESH (no ACM or selective mass scaling used)
- VMS_2DharmonicP.m: this runs the case with selective mass scaling and use of ACM

We also have an auxiliary routine for the element-level computations in Q4 elastic elements (this are called by the above scripts as necessary).

Please see the script files for further annotations and information.