This set of scripts runs the analyses for the 2D mesh with a mass on the top, subjected to a slowly increasing displacement at the top and then left to conduct free vibration.

There are three separate master scripts:

- VMS_2DsingleMeshMtop.m: this runs the case with a SINGLE MESH (no ACM or selective mass scaling used)
- VMS_2DMtop.m: this runs the case where the ACM is conforming to the actual mesh
- VMS_2DnoncomfMtop.m: this runs the case where the ACM is non-conforming to the actual mesh

We also have two auxiliary routines for the element-level computations in Q4 elastic elements (these are called by the above scripts as necessary).

Please see the script files for further annotations and information.