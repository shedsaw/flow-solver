//============================================================
// 
//  Defined_Variables.h
//  
//  List of defined variables for compile time.
//
//  Written by - Shane Sawyer
//
//=============================================================

#ifndef DEF_VAR_H
#define DEF_VAR_H

#define PERCENTILE_CONVERGENCE 0          // Rather than look at all points, check the convergence of points below the percentile.
#define PERCENTILE 99                     // The upper limit on the percentile check. Looks at all values below this threshold.

#define LIMITER_NEIGHBOR_AVERAGE 0        // When firing up the limiter, attempt to smooth the results by averaging the value of all nearest neighbors.

#define DEBUG 0
#define DEBUG_GRAD 0
#define DEBUG_HESS 0
#define DEBUG_DERIVS 0

#define ANNULUS 0
#define NACA 0
#define CYLINDER 0
#define CYLINDER_HACK 0                   // Logic flag to enable a hack for the cylinder cases where part of the grid
                                          // will be ran at 2nd order and the remaining portion will be processed as higher order.

#define MMS 0                             // Do method of manufactured solutions to verify code.

#define MMS_EXP 0

#define NEQ 4                             // Number of equations to solve.
#define NDIM 2                            // Number of spatial dimensions.
#define NUM_ELEM_TYPE 3                   // Number of allowable element types.
                                          //  Intentions:
                                          //    2D - 3 elements:
                                          //           edge,
                                          //           triangle,
                                          //           quad.
                                          //    -> Edge segements are handled implicitly with the generic mesh
                                          //       format of Dr. Karman and do not need connectivity in the cell to
                                          //       node array.  7/3/07 - Made an executive decision to handle them.
#define NUM_VAR 4                         // Number of dependent variables.
                                          //  4 for the 2D Navier-Stokes equations.

#define MAX_NUM_VERT 4                    // Maximum number of vertices - the most nodes an element can have. 4 for quads.

#define NUM_MOM 9                         // Number of control volume moments. 5 for 3rd order, 9 for 4th order.
#define NUM_HESS 5                        // Number of entries stored in the hessian array (includes gradients).
#define NUM_GRAD 2                        // Number of spatial derivatives.

#define HESS_FULL_SYSTEM 0                // Solve the full or reduced LSQR system for the hessian.
#define GEOM_WEIGHT 1                     // The exponent applied to the distance weighting term in the hessian computation.

#define SVDFAC 1                          // Flag to determine if the least squares is to be solved by the SVD factorization.
#define QRFAC 0                           // ---- QR factorization.
#define GSFAC 0                           // ---- Gram-Schmidt orthonormalization.

#define PV_EXT 1                          // Flag for determining if the reconstruction for third order uses the CV averges found in grid->Q
                                          // or if I use the point value from the least squares problem.

#define FULL_GAUSS_FLUX 0                 // Do Gaussian integration for the flux for all orders of reconstruction.

#define RECON_PRIM 0                      // Use the primitive variables to construct the least squares problem.

#define EDGE_BASED 0                      // Turn the subedge solver model into a pseudo-edge based solver.

#define Edge    0
#define Tri     1
#define Quad    2

#define INTERIOR 0
#define BOUNDARY 1
#define GHOST 2
#define PHANTOM 3

#define freenull(p)   { if ( p!=0) free(p); (p) = NULL; }

#define CHECKPT {printf("Checkpoint: %s, line %d\n",__FILE__,__LINE__);  fflush(stdout);}

#define MIN(x,y) ((x) <= (y) ? (x) : (y))
#define MAX(x,y) ((x) >= (y) ? (x) : (y))

#define sign(x) (x < 0.0) ? -1.0 : 1.0

#endif
