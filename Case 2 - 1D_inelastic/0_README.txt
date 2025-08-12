This set of scripts runs the analyses for the 1D mesh with monotonically increasing tip displacement and INELASTIC (elastoplastic or viscoplastic) material behavior.

There are three separate master scripts:

- VMS_1DsingleMeshYield.m: this runs the case with a SINGLE MESH (no ACM or selective mass scaling used)
- VMS_1DYield.m: this runs the case of ELASTOPLASTIC material, with selective mass scaling and use of ACM
- VMS_1DYieldVP.m: this runs the case of VISCOPLASTIC SOFTENING material, with selective mass scaling and use of ACM

We also have two auxiliary routines for the element-level computations: truss1D.m (elastoplastic material) and truss1Dvp.m (viscoplastic material)

Please see the script files for further annotations and information.