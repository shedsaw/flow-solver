//=============================================================
// 
//  compute_derivatives..C
//  
//  Functions to calculate the derivatives of the conserved variables.
//  This code is derived from the work of Ollivier-Gooch et al.
//
//  Written by - Shane Sawyer
//
//=============================================================

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "grid.h"
#include "mesh_utils.h"
#include "Defined_Variables.h"
#include "linear_algebra.h"
#include "hessian.h"
#include "liblapack_interface.h"
#include "grid_io.h"
#include "curved_boundaries.h"
#include "reconstruct.h"
#include "compute_derivatives.h"


//=============================================================
// 
//  Compute_Derivatives()
//
//  Calculates the derivatives of the dependent variables using a
//  least squares technique for k-exact reconstruction.
//  
//  GRID *grid                        // The grid.
//  PARAMS params                     // Parameters for the run.
//  int order_desired                 // The order to which the least squares problem is applied.
//
//=============================================================

void Compute_Derivatives ( GRID *grid, PARAMS params, int order_desired )
{
  int i,j,k;                             // Loop counters.
  int v;                                 // Vector loop.
  int ind;                               // Component loop.
  int nodeL, nodeR;                      // Edge nodes.
  int num_neighbors;                     // Number of neighbors around a node.
  int e;                                 // Element type/loop counter.
  int var, row;                          // More loop counters.
  int node;                              // Neighboring node.
  double dx[NDIM];                       // dx vector.
  double p;                              // Exponent for weighting.
  double dist;                           // Distance from node.
  double w;                              // Weight term.
  double temp;                           // Intermediate calculation.
  int count;                             // Count up offending derivatives.
  int mtail;                             // How much to add to the end of num_neighbors ( 0 or 1 ).
  int bytes_allocated = 0;               // Total number of bytes allocated by the QR factorization.
  int factorize = 0;                     // Whether or not to factorize, 0 by default set to 1 if memory is not alloced.
  double innerprod;                      // Vector inner product.
  int start_index, end_index;

  // Variables used for the SVD routine.
  double *A = NULL;
  double *U = NULL;
  double *sol = NULL;
  double *temp1 = NULL, *temp2 = NULL, *temp3 = NULL;
  double *aa=NULL;
  double *s=NULL;
  double *u=NULL;
  double *vt=NULL;
  double *work=NULL;
  int lda,lds,ldu,ldvt,lwork,info;
  int n;
  int m;
  int MAX_NEIGHBORS, MAX_LWORK;
  int Adim, Udim, tempdim, aadim, udim, vtdim, sdim, workdim;
  char jobu = 'A';
  char jobvt = 'A';
  int nbytes = 0;

  // Debug variables.
  int imax;
  double *qprime = NULL;
  double term;
  double norm;
  char **derivs = NULL;
  int nodes[MAX_NUM_VERT];               // Element vertices.
  int num_vert;                          // Number of vertices for an element.
  int left_id;                           // ID of left edge.
  int right_id;                          // ID of right edge.
  double xc[NDIM];                       // Vectors.
  double xmidL[NDIM];
  double xmidR[NDIM];
  double v_left[NDIM];
  double v_right[NDIM];
  double v_middle[NDIM];
  double bary_coord[7][3];               // The barycentric coordinates.
  double ct_gp[7][2];                    // Gauss integration points in (r,s) space (xi,eta).
  double ct_nodes[7][2];                 // The x,y locations of the 6 nodes defining the p2 triangle (curved side).
  double N1,N2,N3,N4,N5,N6;              // Shape functions.
  double N1r,N1s,N2r,N2s,N3r,N3s;        // Shape function derivatives.
  double N4r,N4s,N5r,N5s,N6r,N6s;
  double xr,xs,yr,ys,jacobian;           // Transformation jacobian.
  double wg[7];                          // The associated weights;
  double b_r = (6. - sqrt(15.))/21.;     // Precompute some the coordinates.
  double b_t = (6. + sqrt(15.))/21.;
  double b_s = (9. + 2.*sqrt(15.))/21.;
  double b_u = (9. - 2.*sqrt(15.))/21.;
  double wA = (155. - sqrt(15.))/1200.;
  double wB = (155. + sqrt(15.))/1200.;

  double XI,ETA;                         // Triangle coordinates.
  double tri_int;                        // Triangle integration.

  double gp[7][3];                       // Gauss points for integration.
  double dA;                             // Differential area (for a triangle).
  double ext_val[7];                     // Extrapolation values at Gauss points.
  double max_diff = 0.;                  // Maximum difference between extrapolated value
                                         // and the analytical value.
  double *ext_cv_avg = NULL;             // CV average from the extrapolation.
  double total_integral;

  int b,seg,bc;
  int inf_subedge, inf_bedge;
  int iedge, gelem;
  int real_node, ghost_node;
  int ct_flag;

  double norm_1_s, norm_2_s, norm_i_s;
  double norm_1_b, norm_2_b, norm_i_b;
  double xmid,ymid;
  double x1,x2,x3;
  double y1,y2,y3;
  double xi,yi;
  double t1,t2,t3;
  double recon_val;
  double XL[2],XR[2],GP[6];
  double xL,xR,yL,yR;

  double QL[NUM_VAR], QR[NUM_VAR];
  
  char buff[100];

  // Pointers.
  int nn = grid->nn;
  int *nnsn2 = grid->nnsn2;
  int *nsn2 = grid->nsn2;
  int *gnnsn1 = grid->gnnsn1;
  int *gnnsn = grid->gnnsn;
  int *gnsn  = grid->gnsn;
  double *x = grid->x;
  double *q = grid->nQ;
  double *hess = grid->hess;
  double *moments = grid->Moments;
  double *point_value = grid->point_value;
  double *cv_avg = NULL;

  factorize = 0;
  p = (double)GEOM_WEIGHT;

  if ( order_desired < 1 || order_desired > 4 )
    {
      printf("FATAL ERROR: The desired order requested to Compute_Derivatives is invalid. <%d>\n",order_desired);
      fflush(stdout);
      exit(1);
    }
  
  if (DEBUG_DERIVS)
    {
      derivs = (char**)malloc(NUM_MOM*sizeof(char*));
      if ( derivs == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'derivs'.\n"); exit(1); }
      
      for ( i=0; i < NUM_MOM; i++ )
	{
	  derivs[i] = NULL;
	  derivs[i] = (char*)malloc(10*sizeof(char));
	  
	  if ( derivs[i] == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'derivs[%d]'.\n",i); exit(1); }
	}

      sprintf((derivs[0]), "Fx");
      sprintf((derivs[1]), "Fy");
      sprintf((derivs[2]), "Fxx");
      sprintf((derivs[3]), "Fyy");
      sprintf((derivs[4]), "Fxy");
      sprintf((derivs[5]), "Fxxx");
      sprintf((derivs[6]), "Fyyy");
      sprintf((derivs[7]), "Fxxy");
      sprintf((derivs[8]), "Fxyy");
      factorize = 1;
    }

  
  // Do some logic checking to make sure the compile time parameters are kosher.
  
  i = SVDFAC + QRFAC + GSFAC;
  
  if ( i == 0 )
    {
      printf("FATAL ERROR: NO METHOD WAS SELECTED TO SOLVE THE LEAST SQAURES PROBLEM.\n");
      exit(0);
    }
  
  if ( i > 1 )
    {
      printf("CRITICAL ERROR: MORE THAN ONE METHOD WAS SELECTED TO SOLVE THE LEAST SQUARES PROBLEM OR THE INPUT IS INCORRECT.\n");
      exit(0);
    }
  
  // Allocate some memory that we need.
  sol = (double*)malloc((NUM_MOM+1)*sizeof(double));
  if ( sol == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'sol'.\n"); exit(1); }
  
  cv_avg = (double*)malloc((nn+1)*sizeof(double));
  if ( cv_avg == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'cv_avg'.\n"); exit(1); }

  ext_cv_avg = (double*)malloc((nn+1)*sizeof(double));
  if ( ext_cv_avg == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'ext_cv_avg'.\n"); exit(1); }
  
  if ( grid->hess == NULL )
    {
      grid->hess = (double*)malloc((grid->nn + 1)*NUM_MOM*NUM_VAR*sizeof(double));
      if ( grid->hess == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'hess'.\n"); exit(1); }
    }
  
  if ( grid->point_value == NULL )
    {
      grid->point_value = (double*)malloc((grid->nn + 1)*NUM_VAR*sizeof(double));
      if ( grid->point_value == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'point_value'.\n"); exit(1); }
    }

  if ( grid->SVD_U == NULL )
    {
      grid->SVD_U = (double**)malloc( (nn + 1)*sizeof(double*));
      if ( grid->SVD_U == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'SVD_U'.\n"); exit(1); }

      for ( i=0; i <= nn; i++ )
	{
	  grid->SVD_U[i] = NULL;
	}

      factorize = 1;
    }

  if ( grid->SVD_V == NULL )
    {
      grid->SVD_V = (double**)malloc( (nn + 1)*sizeof(double*));
      if ( grid->SVD_V == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'SVD_V'.\n"); exit(1); }

      for ( i=0; i <= nn; i++ )
	{
	  grid->SVD_V[i] = NULL;
	}
    }
  
  if ( grid->SVD_S == NULL )
    {
      grid->SVD_S = (double**)malloc( (nn + 1)*sizeof(double*));
      if ( grid->SVD_S == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'SVD_S'.\n"); exit(1); }

      for ( i=0; i <= nn; i++ )
	{
	  grid->SVD_S[i] = NULL;
	}
    }

  if ( grid->SVD_size == NULL )
    {
      grid->SVD_size = (double*)malloc( (nn + 1)*3*sizeof(double*));
      if ( grid->SVD_size == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'SVD_size'.\n"); exit(1); }
    }

  // Do we need to construct and decompose the matrix?
  if ( order_desired != grid->SVDinv_order )
    {
      factorize = 1;
    }
  
  if ( factorize )
    {
      printf("Compute_derivatives report: SVDinv_order = %d. The desired order is %d.\n",grid->SVDinv_order, order_desired);
      fflush(stdout);
    }
  
  // Reset pointers.
  hess = grid->hess;
  point_value = grid->point_value;

  // Clean out all memory that will be reset in this function.
  for ( i=0; i < (NUM_MOM+1); i++ )
    sol[i] = 0.;
  
  for ( i=0; i <= grid->nn; i++ )
    cv_avg[i] = 0.;
  
  for ( i=0; i <= grid->nn; i++ )
    ext_cv_avg[i] = 0.;
  
  for ( i=0; i <= grid->nn; i++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  point_value[i*NUM_VAR+j] = 0.;
	}
    }

  for ( i=0; i <= grid->nn; i++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  for ( k=0; k < NUM_MOM; k++ )
	    hess[i*NUM_MOM*NUM_VAR + j*NUM_MOM + k] = 0.; 
	}
    }
  
  // Get the area averaged values of the test function.
  if (DEBUG_DERIVS)
    {
      Compute_CV_Averages(grid,cv_avg);
      
      // Copy the cv averages over to the Q array.
      for ( i=1; i <= grid->nn; i++ )
	{
	  grid->Q[i*NUM_VAR+0] = cv_avg[i];
	  grid->nQ[i*NUM_VAR+0] = cv_avg[i];
	}
      
      if ( grid->phi == NULL )
	{
	  grid->phi = (double*)malloc( (grid->nn+1)*2*4*sizeof(double) );
	  if ( grid->phi == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'phi'.\n"); exit(1); }
	}
      
      for ( i=1; i <= grid->nn; i++ )
	{
	  for ( j=0; j < 8; j++ )
	    {
	      grid->phi[i*8+j] = 1.;
	    }
	}
      
      grid->citer = 1;
    }

  // Find the maximum degree of neighbors in the mesh and allocate temp memory for dot products.
  MAX_NEIGHBORS = 0;
  for ( i=1; i <= nn; i++ )
    {
      num_neighbors = gnnsn[ gnnsn1[i+1] ] - gnnsn[ ( gnnsn1[i] + 1 ) ];
      MAX_NEIGHBORS = MAX(MAX_NEIGHBORS,num_neighbors);
    }
  
  MAX_NEIGHBORS++;  // Account for the node itself.

  if (DEBUG)
    printf("  Maximum number of 2nd degree neighbors in the mesh is %d.\n",MAX_NEIGHBORS-1);

  tempdim = MAX_NEIGHBORS;

  temp1 = (double*)malloc(tempdim*sizeof(double));
  if ( temp1 == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'temp1'.\n"); exit(1); }
  
  temp2 = (double*)malloc(tempdim*sizeof(double));
  if ( temp2 == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'temp2'.\n"); exit(1); }
  
  temp3 = (double*)malloc(tempdim*sizeof(double));
  if ( temp3 == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'temp3'.\n"); exit(1); }

  for ( i=0; i < tempdim; i++ )
    {
      temp1[i] = 0.;  temp2[i] = 0.;  temp3[i] = 0.;
    }

  Udim = MAX_NEIGHBORS;

  U = (double*)malloc(Udim*sizeof(double));
  if ( U == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'U'.\n"); exit(1); }

  for ( i=0; i < Udim; i++ )
    U[i] = 0.;

  // if we need to recompute (or compute the first pass) the svd, allocate memory for lapack call, find the storage requirements,
  // and allocate space for the pseudoinverse.
  if ( SVDFAC && factorize )
    {
      if ( grid->SVD_U != NULL )
	{
	  for ( i=1; i <= nn; i++ )
	    {
	      /*
	      if ( grid->SVD_U[i] != NULL )
		freenull( grid->SVD_U[i] );
	      
	      if ( grid->SVD_V[i] != NULL )
		freenull( grid->SVD_V[i] );
	      
	      if ( grid->SVD_S[i] != NULL )
		freenull( grid->SVD_S[i] );
	      */
	      
	      freenull( grid->SVD_U[i] );
	      freenull( grid->SVD_V[i] );
	      freenull( grid->SVD_S[i] );
	    }
	}
      
      // First we want to allocate memory used in the SVD routine.
      n = NUM_MOM+1;
      i = MIN(MAX_NEIGHBORS,n) * 3 + MAX(MAX_NEIGHBORS,n);
      j = MIN(MAX_NEIGHBORS,n) * 5;
      MAX_LWORK = MAX(i,j) * 2;

      Adim = n*MAX_NEIGHBORS;
      aadim = MAX_NEIGHBORS*n;
      udim = MAX_NEIGHBORS*MAX_NEIGHBORS;
      vtdim = n*n;
      sdim = MAX_NEIGHBORS;
      workdim = MAX_LWORK;
  
      A = (double*)malloc(Adim*sizeof(double));
      if ( A == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'A'.\n"); exit(1); }

      aa = (double*)malloc(aadim*sizeof(double));
      if ( aa == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'aa'.\n"); exit(1); }

      u = (double*)malloc(udim*sizeof(double));
      if ( u == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'u'.\n"); exit(1); }

      vt = (double*)malloc(vtdim*sizeof(double));
      if ( vt == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'vt'.\n"); exit(1); }

      s = (double*)malloc(sdim*sizeof(double));
      if ( s == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 's'.\n"); exit(1); }

      work = (double*)malloc(workdim*sizeof(double));
      if ( work == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'work'.\n"); exit(1); }

      // Clean the memory.
      for ( i=0; i < Adim; i++ )
	A[i] = 0.;

      for ( i=0; i < aadim; i++ )
	aa[i] = 0.;
      
      for ( i=0; i < udim; i++ )
	u[i] = 0.;
      
      for ( i=0; i < vtdim; i++ )
	vt[i] = 0.;
      
      for ( i=0; i < sdim; i++ )
	s[i] = 0.;
      
      for ( i=0; i < workdim; i++ )
	work[i] = 0.;

      // now we need to allocate space for the pseudoinverse.
      for ( i=1; i <= nn; i++ )
	{
	  num_neighbors = gnnsn[ gnnsn1[i+1] ] - gnnsn[ ( gnnsn1[i] + 1 ) ];
	  
	  // how many columns in the system?
	  switch ( order_desired )
	    {
	    case 1:
	      n = 2;
	      break;
	    case 2:    // second order
	      n = 2;
	      break;
	    case 3:    // third order
	      n = 5;
	      break;
	    case 4:    // fourth order
	      n = 9;
	      break;
	    default:
	      printf("FATAL ERROR: order_desired is < %d >?\n",order_desired);
	      fflush(stdout);
	      exit(1);
	    }

	  if ( HESS_FULL_SYSTEM )
	    {
	      n++;
	      num_neighbors++;
	    }

	  // Allocate the space.
	  m = num_neighbors;

	  // U is mxm, V is nxn, S is MIN(m,n);
	  grid->SVD_size[i*3+0] = m;
	  grid->SVD_size[i*3+1] = n;
	  grid->SVD_size[i*3+2] = MIN(m,n);

	  grid->SVD_U[i] = (double*)malloc((m*m)*sizeof(double));
	  if ( grid->SVD_U[i] == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'SVD_U[%d]'.\n",i); exit(1); }

	  grid->SVD_V[i] = (double*)malloc((n*n)*sizeof(double));
	  if ( grid->SVD_V[i] == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'SVD_V[%d]'.\n",i); exit(1); }

	  grid->SVD_S[i] = (double*)malloc((MIN(m,n))*sizeof(double));
	  if ( grid->SVD_S[i] == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'SVD_S[%d]'.\n",i); exit(1); }
	  
	  nbytes += (m*m) + (n*n) + MIN(m,n);

	}

      printf("Finished allocating memory for the SVD in Compute_Derivatives(). A total of %d doubles were allocated, or a total of %d bytes.\n",nbytes,nbytes*8);
    }


  // Now we need to build the linear system and decompose it.
  if ( SVDFAC && factorize )
    {
      if ( HESS_FULL_SYSTEM )
	{
	  if ( order_desired == 1 || order_desired == 2 )  // 2nd
	    {
	      // know that n is 2+1
	      n = 3;
	      
	      for ( i=1; i <= nn; i++ )
		{
		  num_neighbors = gnnsn[ gnnsn1[i+1] ] - gnnsn[ gnnsn1[i] ];
		  
		  // Put in the mean constraint.
		  A[0 + 0*num_neighbors] = 1.0;
		  A[0 + 1*num_neighbors] = moments[i*NUM_MOM+0]; // x
		  A[0 + 2*num_neighbors] = moments[i*NUM_MOM+1]; // y
		  
		  start_index = gnnsn[ (gnnsn1[i]+1) ];
		  end_index   = gnnsn[ gnnsn1[i+1] ];
		  
		  for ( j=start_index; j < end_index; j++ )
		    {
		      node = gnsn[j];
		      
		      // Get the geometric weight term for current node.
		      dx[0] = x[node*NDIM+0] - x[i*NDIM+0];
		      dx[1] = x[node*NDIM+1] - x[i*NDIM+1];
		      
		      dist = sqrt(dx[0]*dx[0] + dx[1]*dx[1]);
		      
		      w = 1.0 / pow(dist,p);
		      
		      row = j - start_index + 1;
		      
		      A[row+0*num_neighbors] = w*1.0;
		      
		      A[row+1*num_neighbors] = w*(x[node*NDIM+0] - x[i*NDIM+0] +
						  moments[node*NUM_MOM+0]);
		      
		      A[row+2*num_neighbors] = w*(x[node*NDIM+1] - x[i*NDIM+1] +
						  moments[node*NUM_MOM+1]);
		    }
		  
		  // Set some variables.
		  m = num_neighbors;
		  lwork = MAX(3*MIN(m,n)+MAX(m,n),5*MIN(m,n));
		  lds = MIN(m,n);
		  lda=m;
		  ldu=m;
		  ldvt=n;
		  
		  // Copy the matrix over.
		  for ( k=0; k < m*n; k++ )
		    aa[k] = A[k];
		  
		  // Get the SVD
		  dgesvd_( &jobu, &jobvt, &m, &n, aa, &lda, s, u, &ldu, vt, &ldvt, work, &MAX_LWORK, &info );

		  // now store the pseudoinverse.
		  for ( j=0; j < (m*m); j++ )
		    grid->SVD_U[i][j] = u[j];

		  for ( j=0; j < (n*n); j++ )
		    grid->SVD_V[i][j] = vt[j];

		  for ( j=0; j < lds; j++ )
		    grid->SVD_S[i][j] = s[j];
		}
	    }

	  else if ( order_desired == 3 )  // 3rd
	    {
	      // know that n is 5+1
	      n = 6;

	      for ( i=1; i <= nn; i++ )
		{
		  if ( grid->node_order[i] < 3 )
		    continue;
		  
		  num_neighbors = gnnsn[ gnnsn1[i+1] ] - gnnsn[ gnnsn1[i] ];

		  // Put in the mean constraint.
		  A[0 + 0*num_neighbors] = 1.0;
		  A[0 + 1*num_neighbors] = moments[i*NUM_MOM+0]; // x
		  A[0 + 2*num_neighbors] = moments[i*NUM_MOM+1]; // y
		  A[0 + 3*num_neighbors] = 0.5*moments[i*NUM_MOM+2]; // xx
		  A[0 + 4*num_neighbors] = 0.5*moments[i*NUM_MOM+3]; // yy
		  A[0 + 5*num_neighbors] = moments[i*NUM_MOM+4]; // xy

		  start_index = gnnsn[ (gnnsn1[i]+1) ];
		  end_index   = gnnsn[ gnnsn1[i+1] ];

		  for ( j=start_index; j < end_index; j++ )
		    {
		      node = gnsn[j];

		      // Get the geometric weight term for current node.
		      dx[0] = x[node*NDIM+0] - x[i*NDIM+0];
		      dx[1] = x[node*NDIM+1] - x[i*NDIM+1];

		      dist = sqrt(dx[0]*dx[0] + dx[1]*dx[1]);

		      w = 1.0 / pow(dist,p);
		  
		      row = j - start_index + 1;
	      
		      A[row+0*num_neighbors] = w*1.0;
		      
		      A[row+1*num_neighbors] = w*(x[node*NDIM+0] - x[i*NDIM+0] +
						  moments[node*NUM_MOM+0]);

		      A[row+2*num_neighbors] = w*(x[node*NDIM+1] - x[i*NDIM+1] +
						  moments[node*NUM_MOM+1]);

		      A[row+3*num_neighbors] = 0.5*w*( (x[node*NDIM+0] - x[i*NDIM+0])*
						       (x[node*NDIM+0] - x[i*NDIM+0]) +
						       2.0*moments[node*NUM_MOM+0]*(x[node*NDIM+0] - x[i*NDIM+0]) +
						       moments[node*NUM_MOM+2] );

		      A[row+4*num_neighbors] = 0.5*w*( (x[node*NDIM+1] - x[i*NDIM+1])*
						       (x[node*NDIM+1] - x[i*NDIM+1]) +
						       2.0*moments[node*NUM_MOM+1]*(x[node*NDIM+1] - x[i*NDIM+1]) +
						       moments[node*NUM_MOM+3] );

		      A[row+5*num_neighbors] = w*(moments[node*NUM_MOM+4] +
						  moments[node*NUM_MOM+0]*(x[node*NDIM+1] - x[i*NDIM+1]) +
						  moments[node*NUM_MOM+1]*(x[node*NDIM+0] - x[i*NDIM+0]) +
						  (x[node*NDIM+1] - x[i*NDIM+1])*(x[node*NDIM+0] - x[i*NDIM+0]));
		    }

		  // Set some variables.
		  m = num_neighbors;
		  lwork = MAX(3*MIN(m,n)+MAX(m,n),5*MIN(m,n));
		  lds = MIN(m,n);
		  lda=m;
		  ldu=m;
		  ldvt=n;
		  
		  // Copy the matrix over.
		  for ( k=0; k < m*n; k++ )
		    aa[k] = A[k];
		  
		  // Get the SVD
		  dgesvd_( &jobu, &jobvt, &m, &n, aa, &lda, s, u, &ldu, vt, &ldvt, work, &MAX_LWORK, &info );
		  
		  // now store the pseudoinverse.
		  for ( j=0; j < (m*m); j++ )
		    grid->SVD_U[i][j] = u[j];

		  for ( j=0; j < (n*n); j++ )
		    grid->SVD_V[i][j] = vt[j];

		  for ( j=0; j < lds; j++ )
		    grid->SVD_S[i][j] = s[j];
		}
	    }

	  else  // 4th
	    {
	      // know that n is 9 + 1
	      n = 10;
	      
	      for ( i=1; i <= nn; i++ )
		{
		  if ( grid->node_order[i] < 3 )
		    continue;

		  num_neighbors = gnnsn[ gnnsn1[i+1] ] - gnnsn[ gnnsn1[i] ];
		  
		  // Put in the mean constraint.
		  A[0 + 0*num_neighbors] = 1.0;
		  A[0 + 1*num_neighbors] = moments[i*NUM_MOM+0]; // x
		  A[0 + 2*num_neighbors] = moments[i*NUM_MOM+1]; // y
		  A[0 + 3*num_neighbors] = 0.5*moments[i*NUM_MOM+2]; // xx
		  A[0 + 4*num_neighbors] = 0.5*moments[i*NUM_MOM+3]; // yy
		  A[0 + 5*num_neighbors] = moments[i*NUM_MOM+4]; // xy
		  A[0 + 6*num_neighbors] = (1./6.)*moments[i*NUM_MOM+5]; // xxx
		  A[0 + 7*num_neighbors] = (1./6.)*moments[i*NUM_MOM+6]; // yyy
		  A[0 + 8*num_neighbors] = 0.5*moments[i*NUM_MOM+7]; // xxy
		  A[0 + 9*num_neighbors] = 0.5*moments[i*NUM_MOM+8]; // xyy

		  start_index = gnnsn[ (gnnsn1[i]+1) ];
		  end_index   = gnnsn[ gnnsn1[i+1] ];

		  for ( j=start_index; j < end_index; j++ )
		    {
		      node = gnsn[j];

		      // Get the geometric weight term for current node.
		      dx[0] = x[node*NDIM+0] - x[i*NDIM+0];
		      dx[1] = x[node*NDIM+1] - x[i*NDIM+1];

		      dist = sqrt(dx[0]*dx[0] + dx[1]*dx[1]);

		      w = 1.0 / pow(dist,p);
		  
		      row = j - start_index + 1;
	      
		      A[row+0*num_neighbors] = w*1.0;
		      
		      A[row+1*num_neighbors] = w*(x[node*NDIM+0] - x[i*NDIM+0] +
						  moments[node*NUM_MOM+0]);

		      A[row+2*num_neighbors] = w*(x[node*NDIM+1] - x[i*NDIM+1] +
						  moments[node*NUM_MOM+1]);

		      A[row+3*num_neighbors] = 0.5*w*( (x[node*NDIM+0] - x[i*NDIM+0])*
						       (x[node*NDIM+0] - x[i*NDIM+0]) +
						       2.0*moments[node*NUM_MOM+0]*(x[node*NDIM+0] - x[i*NDIM+0]) +
						       moments[node*NUM_MOM+2] );

		      A[row+4*num_neighbors] = 0.5*w*( (x[node*NDIM+1] - x[i*NDIM+1])*
						       (x[node*NDIM+1] - x[i*NDIM+1]) +
						       2.0*moments[node*NUM_MOM+1]*(x[node*NDIM+1] - x[i*NDIM+1]) +
						       moments[node*NUM_MOM+3] );

		      A[row+5*num_neighbors] = w*(moments[node*NUM_MOM+4] +
						  moments[node*NUM_MOM+0]*(x[node*NDIM+1] - x[i*NDIM+1]) +
						  moments[node*NUM_MOM+1]*(x[node*NDIM+0] - x[i*NDIM+0]) +
						  (x[node*NDIM+1] - x[i*NDIM+1])*(x[node*NDIM+0] - x[i*NDIM+0]));

		      A[row+6*num_neighbors] = (1./6.)*w*(moments[node*NUM_MOM+5] +
							  3.0*moments[node*NUM_MOM+2]*(x[node*NDIM+0] - x[i*NDIM+0]) +
							  3.0*moments[node*NUM_MOM+0]*(x[node*NDIM+0] - x[i*NDIM+0])*(x[node*NDIM+0] - x[i*NDIM+0]) +
							  (x[node*NDIM+0] - x[i*NDIM+0])*(x[node*NDIM+0] - x[i*NDIM+0])*(x[node*NDIM+0] - x[i*NDIM+0]));

		      
		      A[row+7*num_neighbors] = (1./6.)*w*(moments[node*NUM_MOM+6] +
							  3.0*moments[node*NUM_MOM+3]*(x[node*NDIM+1] - x[i*NDIM+1]) +
							  3.0*moments[node*NUM_MOM+1]*(x[node*NDIM+1] - x[i*NDIM+1])*(x[node*NDIM+1] - x[i*NDIM+1]) +
							  (x[node*NDIM+1] - x[i*NDIM+1])*(x[node*NDIM+1] - x[i*NDIM+1])*(x[node*NDIM+1] - x[i*NDIM+1]));

		      A[row+8*num_neighbors] = 0.5*w*(moments[node*NUM_MOM+7] +
						      2.0*moments[node*NUM_MOM+4]*(x[node*NDIM+0] - x[i*NDIM+0]) +
						      2.0*moments[node*NUM_MOM+0]*(x[node*NDIM+0] - x[i*NDIM+0])*(x[node*NDIM+1] - x[i*NDIM+1]) +
						      moments[node*NUM_MOM+1]*(x[node*NDIM+0] - x[i*NDIM+0])*(x[node*NDIM+0] - x[i*NDIM+0]) +
						      moments[node*NUM_MOM+2]*(x[node*NDIM+1] - x[i*NDIM+1]) +
						      (x[node*NDIM+0] - x[i*NDIM+0])*(x[node*NDIM+0] - x[i*NDIM+0])*(x[node*NDIM+1] - x[i*NDIM+1]) );

		      A[row+9*num_neighbors] = 0.5*w*(moments[node*NUM_MOM+8] +
						      2.0*moments[node*NUM_MOM+4]*(x[node*NDIM+1] - x[i*NDIM+1]) +
						      2.0*moments[node*NUM_MOM+1]*(x[node*NDIM+0] - x[i*NDIM+0])*(x[node*NDIM+1] - x[i*NDIM+1]) +
						      moments[node*NUM_MOM+0]*(x[node*NDIM+1] - x[i*NDIM+1])*(x[node*NDIM+1] - x[i*NDIM+1]) +
						      moments[node*NUM_MOM+3]*(x[node*NDIM+0] - x[i*NDIM+0]) +
						      (x[node*NDIM+0] - x[i*NDIM+0])*(x[node*NDIM+1] - x[i*NDIM+1])*(x[node*NDIM+1] - x[i*NDIM+1]) );
						      

		    }

		  // Set some variables.
		  m = num_neighbors;
		  lwork = MAX(3*MIN(m,n)+MAX(m,n),5*MIN(m,n));
		  lds = MIN(m,n);
		  lda=m;
		  ldu=m;
		  ldvt=n;
		  
		  // Copy the matrix over.
		  for ( k=0; k < m*n; k++ )
		    aa[k] = A[k];
		  
		  // Get the SVD
		  dgesvd_( &jobu, &jobvt, &m, &n, aa, &lda, s, u, &ldu, vt, &ldvt, work, &MAX_LWORK, &info );
		  
		  // now store the pseudoinverse.
		  for ( j=0; j < (m*m); j++ )
		    grid->SVD_U[i][j] = u[j];

		  for ( j=0; j < (n*n); j++ )
		    grid->SVD_V[i][j] = vt[j];

		  for ( j=0; j < lds; j++ )
		    grid->SVD_S[i][j] = s[j];
		}
	    }
	}

      // reduced system
      else
	{
	  if ( order_desired == 1 || order_desired == 2 )  // 2nd
	    {
	      // know that n is 2
	      n = 2;
	      
	      for ( i=1; i <= nn; i++ )
		{
		  num_neighbors = gnnsn[ gnnsn1[i+1] ] - gnnsn[ ( gnnsn1[i] + 1 ) ];
		  
		  start_index = gnnsn[ (gnnsn1[i]+1) ];
		  end_index   = gnnsn[ gnnsn1[i+1] ];
		  
		  for ( j=start_index; j < end_index; j++ )
		    {
		      node = gnsn[j];
		      
		      // Get the geometric weight term for current node.
		      dx[0] = x[node*NDIM+0] - x[i*NDIM+0];
		      dx[1] = x[node*NDIM+1] - x[i*NDIM+1];
		      
		      dist = sqrt(dx[0]*dx[0] + dx[1]*dx[1]);
		      
		      w = 1.0 / pow(dist,p);
		      
		      row = j - start_index;
		      
		      A[row+0*num_neighbors] = w*(x[node*NDIM+0] - x[i*NDIM+0] +
						  moments[node*NUM_MOM+0] -
						  moments[i*NUM_MOM+0] );
		      
		      A[row+1*num_neighbors] = w*(x[node*NDIM+1] - x[i*NDIM+1] +
						  moments[node*NUM_MOM+1] -
						  moments[i*NUM_MOM+1] );
		    }
		  
		  // Set some variables.
		  m = num_neighbors;
		  lwork = MAX(3*MIN(m,n)+MAX(m,n),5*MIN(m,n));
		  lds = MIN(m,n);
		  lda=m;
		  ldu=m;
		  ldvt=n;
		  
		  // Copy the matrix over.
		  for ( k=0; k < m*n; k++ )
		    aa[k] = A[k];
		  
		  // Get the SVD
		  dgesvd_( &jobu, &jobvt, &m, &n, aa, &lda, s, u, &ldu, vt, &ldvt, work, &MAX_LWORK, &info );

		  // now store the pseudoinverse.
		  for ( j=0; j < (m*m); j++ )
		    grid->SVD_U[i][j] = u[j];

		  for ( j=0; j < (n*n); j++ )
		    grid->SVD_V[i][j] = vt[j];

		  for ( j=0; j < lds; j++ )
		    grid->SVD_S[i][j] = s[j];
		}
	    }

	  else if ( order_desired == 3 )  // 3rd
	    {
	      // know that n is 5
	      n = 5;

	      for ( i=1; i <= nn; i++ )
		{
		  if ( grid->node_order[i] < 3 )
		    continue;

		  num_neighbors = gnnsn[ gnnsn1[i+1] ] - gnnsn[ ( gnnsn1[i] + 1 ) ];

		  start_index = gnnsn[ (gnnsn1[i]+1) ];
		  end_index   = gnnsn[ gnnsn1[i+1] ];

		  for ( j=start_index; j < end_index; j++ )
		    {
		      node = gnsn[j];
		      
		      // Get the geometric weight term for current node.
		      dx[0] = x[node*NDIM+0] - x[i*NDIM+0];
		      dx[1] = x[node*NDIM+1] - x[i*NDIM+1];
		      
		      dist = sqrt(dx[0]*dx[0] + dx[1]*dx[1]);
		      
		      w = 1.0 / pow(dist,p);
		      
		      row = j - start_index;
		      
		      A[row+0*num_neighbors] = w*(x[node*NDIM+0] - x[i*NDIM+0] +
						  moments[node*NUM_MOM+0] -
						  moments[i*NUM_MOM+0]);

		      A[row+1*num_neighbors] = w*(x[node*NDIM+1] - x[i*NDIM+1] +
						  moments[node*NUM_MOM+1] -
						  moments[i*NUM_MOM+1]);

		      A[row+2*num_neighbors] = 0.5*w*( (x[node*NDIM+0] - x[i*NDIM+0])*
						       (x[node*NDIM+0] - x[i*NDIM+0]) +
						       2.0*moments[node*NUM_MOM+0]*(x[node*NDIM+0] - x[i*NDIM+0]) +
						       moments[node*NUM_MOM+2] -
						       moments[i*NUM_MOM+2]);

		      A[row+3*num_neighbors] = 0.5*w*( (x[node*NDIM+1] - x[i*NDIM+1])*
						       (x[node*NDIM+1] - x[i*NDIM+1]) +
						       2.0*moments[node*NUM_MOM+1]*(x[node*NDIM+1] - x[i*NDIM+1]) +
						       moments[node*NUM_MOM+3] -
						       moments[i*NUM_MOM+3]);

		      A[row+4*num_neighbors] = w*(moments[node*NUM_MOM+4] +
						  moments[node*NUM_MOM+0]*(x[node*NDIM+1] - x[i*NDIM+1]) +
						  moments[node*NUM_MOM+1]*(x[node*NDIM+0] - x[i*NDIM+0]) +
						  (x[node*NDIM+1] - x[i*NDIM+1])*(x[node*NDIM+0] - x[i*NDIM+0]) -
						  moments[i*NUM_MOM+4]);
		    }

		  // Set some variables.
		  m = num_neighbors;
		  lwork = MAX(3*MIN(m,n)+MAX(m,n),5*MIN(m,n));
		  lds = MIN(m,n);
		  lda=m;
		  ldu=m;
		  ldvt=n;
		  
		  // Copy the matrix over.
		  for ( k=0; k < m*n; k++ )
		    aa[k] = A[k];
		  
		  // Get the SVD
		  dgesvd_( &jobu, &jobvt, &m, &n, aa, &lda, s, u, &ldu, vt, &ldvt, work, &MAX_LWORK, &info );
		  
		  // now store the pseudoinverse.
		  for ( j=0; j < (m*m); j++ )
		    grid->SVD_U[i][j] = u[j];

		  for ( j=0; j < (n*n); j++ )
		    grid->SVD_V[i][j] = vt[j];

		  for ( j=0; j < lds; j++ )
		    grid->SVD_S[i][j] = s[j];
		}
	    }

	  else  // 4th
	    {
	      // know that n is 9
	      n = 9;
	      
	      for ( i=1; i <= nn; i++ )
		{
		  if ( grid->node_order[i] < 3 )
		    continue;

		  num_neighbors = gnnsn[ gnnsn1[i+1] ] - gnnsn[ ( gnnsn1[i] + 1 ) ];
		 
		  start_index = gnnsn[ (gnnsn1[i]+1) ];
		  end_index   = gnnsn[ gnnsn1[i+1] ];

		  for ( j=start_index; j < end_index; j++ )
		    {
		      node = gnsn[j];

		      // Get the geometric weight term for current node.
		      dx[0] = x[node*NDIM+0] - x[i*NDIM+0];
		      dx[1] = x[node*NDIM+1] - x[i*NDIM+1];

		      dist = sqrt(dx[0]*dx[0] + dx[1]*dx[1]);

		      w = 1.0 / pow(dist,p);
		  
		      row = j - start_index;
	      
		      A[row+0*num_neighbors] = w*(x[node*NDIM+0] - x[i*NDIM+0] +
						  moments[node*NUM_MOM+0] -
						  moments[i*NUM_MOM+0]);

		      A[row+1*num_neighbors] = w*(x[node*NDIM+1] - x[i*NDIM+1] +
						  moments[node*NUM_MOM+1] -
						  moments[i*NUM_MOM+1]);

		      A[row+2*num_neighbors] = 0.5*w*( (x[node*NDIM+0] - x[i*NDIM+0])*
						       (x[node*NDIM+0] - x[i*NDIM+0]) +
						       2.0*moments[node*NUM_MOM+0]*(x[node*NDIM+0] - x[i*NDIM+0]) +
						       moments[node*NUM_MOM+2] -
						       moments[i*NUM_MOM+2]);

		      A[row+3*num_neighbors] = 0.5*w*( (x[node*NDIM+1] - x[i*NDIM+1])*
						       (x[node*NDIM+1] - x[i*NDIM+1]) +
						       2.0*moments[node*NUM_MOM+1]*(x[node*NDIM+1] - x[i*NDIM+1]) +
						       moments[node*NUM_MOM+3] -
						       moments[i*NUM_MOM+3]);

		      A[row+4*num_neighbors] = w*(moments[node*NUM_MOM+4] +
						  moments[node*NUM_MOM+0]*(x[node*NDIM+1] - x[i*NDIM+1]) +
						  moments[node*NUM_MOM+1]*(x[node*NDIM+0] - x[i*NDIM+0]) +
						  (x[node*NDIM+1] - x[i*NDIM+1])*(x[node*NDIM+0] - x[i*NDIM+0]) -
						  moments[i*NUM_MOM+4]);

		      A[row+5*num_neighbors] = (1./6.)*w*(moments[node*NUM_MOM+5] +
							  3.0*moments[node*NUM_MOM+2]*(x[node*NDIM+0] - x[i*NDIM+0]) +
							  3.0*moments[node*NUM_MOM+0]*(x[node*NDIM+0] - x[i*NDIM+0])*(x[node*NDIM+0] - x[i*NDIM+0]) +
							  (x[node*NDIM+0] - x[i*NDIM+0])*(x[node*NDIM+0] - x[i*NDIM+0])*(x[node*NDIM+0] - x[i*NDIM+0]) -
							  moments[i*NUM_MOM+5]);

		      
		      A[row+6*num_neighbors] = (1./6.)*w*(moments[node*NUM_MOM+6] +
							  3.0*moments[node*NUM_MOM+3]*(x[node*NDIM+1] - x[i*NDIM+1]) +
							  3.0*moments[node*NUM_MOM+1]*(x[node*NDIM+1] - x[i*NDIM+1])*(x[node*NDIM+1] - x[i*NDIM+1]) +
							  (x[node*NDIM+1] - x[i*NDIM+1])*(x[node*NDIM+1] - x[i*NDIM+1])*(x[node*NDIM+1] - x[i*NDIM+1]) -
							  moments[i*NUM_MOM+6]);

		      A[row+7*num_neighbors] = 0.5*w*(moments[node*NUM_MOM+7] +
						      2.0*moments[node*NUM_MOM+4]*(x[node*NDIM+0] - x[i*NDIM+0]) +
						      2.0*moments[node*NUM_MOM+0]*(x[node*NDIM+0] - x[i*NDIM+0])*(x[node*NDIM+1] - x[i*NDIM+1]) +
						      moments[node*NUM_MOM+1]*(x[node*NDIM+0] - x[i*NDIM+0])*(x[node*NDIM+0] - x[i*NDIM+0]) +
						      moments[node*NUM_MOM+2]*(x[node*NDIM+1] - x[i*NDIM+1]) +
						      (x[node*NDIM+0] - x[i*NDIM+0])*(x[node*NDIM+0] - x[i*NDIM+0])*(x[node*NDIM+1] - x[i*NDIM+1]) -
						      moments[i*NUM_MOM+7]);
		      
		      A[row+8*num_neighbors] = 0.5*w*(moments[node*NUM_MOM+8] +
						      2.0*moments[node*NUM_MOM+4]*(x[node*NDIM+1] - x[i*NDIM+1]) +
						      2.0*moments[node*NUM_MOM+1]*(x[node*NDIM+0] - x[i*NDIM+0])*(x[node*NDIM+1] - x[i*NDIM+1]) +
						      moments[node*NUM_MOM+0]*(x[node*NDIM+1] - x[i*NDIM+1])*(x[node*NDIM+1] - x[i*NDIM+1]) +
						      moments[node*NUM_MOM+3]*(x[node*NDIM+0] - x[i*NDIM+0]) +
						      (x[node*NDIM+0] - x[i*NDIM+0])*(x[node*NDIM+1] - x[i*NDIM+1])*(x[node*NDIM+1] - x[i*NDIM+1]) -
						      moments[i*NUM_MOM+8]);
		      
		      
		    }

		  if ( i == 0 )
		    {
		      for ( j=0; j < num_neighbors; j++ )
			{
			  printf("%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",
				 A[j + 0*num_neighbors], A[j + 1*num_neighbors], A[j + 2*num_neighbors], A[j + 3*num_neighbors],
				 A[j + 4*num_neighbors], A[j + 5*num_neighbors], A[j + 6*num_neighbors], A[j + 7*num_neighbors],
				 A[j + 8*num_neighbors] );
			}
		    }
		  
		  // Set some variable.
		  m = num_neighbors;
		  lwork = MAX(3*MIN(m,n)+MAX(m,n),5*MIN(m,n));
		  lds = MIN(m,n);
		  lda=m;
		  ldu=m;
		  ldvt=n;

		  m = grid->SVD_size[i*3+0];
		  n = grid->SVD_size[i*3+1];
		  lds =  grid->SVD_size[i*3+2];
		  
		  // Copy the matrix over.
		  for ( k=0; k < m*n; k++ )
		    aa[k] = A[k];
		  
		  // Get the SVD
		  dgesvd_( &jobu, &jobvt, &m, &n, aa, &lda, s, u, &ldu, vt, &ldvt, work, &MAX_LWORK, &info );
		  
		  // now store the pseudoinverse.
		  for ( j=0; j < (m*m); j++ )
		    grid->SVD_U[i][j] = u[j];

		  for ( j=0; j < (n*n); j++ )
		    grid->SVD_V[i][j] = vt[j];

		  for ( j=0; j < lds; j++ )
		    grid->SVD_S[i][j] = s[j];
		}
	    }
	}

      grid->SVDinv_order = order_desired;

    }
  else if ( QRFAC && factorize )
    {
      // do nothing now.
    }
  else   // Gram-Schmidt orthonormalization.
    {
      // do nothing now.
    }

  // At this point, all the matrix decompositions have been performed either above or in a previous function call, so now I need to
  // use the saved pseudoinverse to get the least squares solution.

  if ( SVDFAC )
    {
      if ( HESS_FULL_SYSTEM )
	{
	  if ( order_desired == 1 || order_desired == 2 )  // 2nd
	    {
	      for ( i=1; i <= nn; i++ )
		{
		  // set some variables.
		  m = grid->SVD_size[i*3+0];
		  n = grid->SVD_size[i*3+1];
		  lds = grid->SVD_size[i*3+2];

		  // loop over all the variables and compute the repsective derivatives.
		  for( var=0; var < NUM_VAR; var++)
		    {
		      U[0] = q[i*NUM_VAR+var];
		      
		      if (DEBUG_DERIVS)
			U[0] = cv_avg[i];
		      
		      start_index = gnnsn[ (gnnsn1[i]+1) ];
		      end_index   = gnnsn[ gnnsn1[i+1] ];

		      for ( j=start_index; j < end_index; j++ )
			{
			  node = gnsn[j];

			  // Get the geometric weight term for current node.
			  dx[0] = x[node*NDIM+0] - x[i*NDIM+0];
			  dx[1] = x[node*NDIM+1] - x[i*NDIM+1];
			  
			  dist = sqrt(dx[0]*dx[0] + dx[1]*dx[1]);
		      
			  w = 1.0 / pow(dist,p);
			  
			  row = j - start_index + 1;

			  U[row] = w*q[node*NUM_VAR+var];
			  
			  if (DEBUG_DERIVS)
			    U[row] = w*cv_avg[node];
			  
			}

		      // Perform the linear algebra for the system Ax=b -> (SVD)x=b
		      for ( j=0; j < m; j++ )
			{
			  for( k=0; k < m; k++ )
			    {
			      temp1[k]=grid->SVD_U[i][j*m+k];
			    }
			  temp2[j]=dotProductGeneral(U,temp1,m);
			}
		      
		      for ( j=0; j < tempdim; j++ )
			temp1[j] = 0.;
	      
		      for ( j=0; j < lds; j++ )
			{
			  temp1[j]=temp2[j] / ( grid->SVD_S[i][j] );
			}
	      
		      for ( j=0; j < tempdim; j++ )
			temp2[j] = 0.;
		      for ( j=0; j < tempdim; j++ )
			temp3[j] = 0.;

	      
		      for ( j=0; j < n; j++ )
			{
			  for( k=0; k < n; k++ )
			    {
			      temp3[k]= grid->SVD_V[i][j*n+k];
			    }
			  sol[j] = dotProductGeneral(temp3,temp1,n);
			}
		      
		      hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 0] = sol[1];
		      hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 1] = sol[2];
		      
		      point_value[i*NUM_VAR+var] = sol[0];
		    }
		}
	    }
	  else if ( order_desired == 3 )
	    {
	      for ( i=1; i <= nn; i++ )
		{
		  if ( grid->node_order[i] < 3 )
		    continue;

		  // set some variables.
		  m = grid->SVD_size[i*3+0];
		  n = grid->SVD_size[i*3+1];
		  lds = grid->SVD_size[i*3+2];

		  // loop over all the variables and compute the repsective derivatives.
		  for( var=0; var < NUM_VAR; var++)
		    {
		      U[0] = q[i*NUM_VAR+var];
		      
		      if (DEBUG_DERIVS)
			U[0] = cv_avg[i];
		      
		      start_index = gnnsn[ (gnnsn1[i]+1) ];
		      end_index   = gnnsn[ gnnsn1[i+1] ];

		      for ( j=start_index; j < end_index; j++ )
			{
			  node = gnsn[j];

			  // Get the geometric weight term for current node.
			  dx[0] = x[node*NDIM+0] - x[i*NDIM+0];
			  dx[1] = x[node*NDIM+1] - x[i*NDIM+1];
			  
			  dist = sqrt(dx[0]*dx[0] + dx[1]*dx[1]);
		      
			  w = 1.0 / pow(dist,p);
			  
			  row = j - start_index + 1;

			  U[row] = w*q[node*NUM_VAR+var];
			  
			  if (DEBUG_DERIVS)
			    U[row] = w*cv_avg[node];
			  
			}

		      // Perform the linear algebra for the system Ax=b -> (SVD)x=b
		      for ( j=0; j < m; j++ )
			{
			  for( k=0; k < m; k++ )
			    {
			      temp1[k]=grid->SVD_U[i][j*m+k];
			    }
			  temp2[j]=dotProductGeneral(U,temp1,m);
			}
		      
		      for ( j=0; j < tempdim; j++ )
			temp1[j] = 0.;
	      
		      for ( j=0; j < lds; j++ )
			{
			  temp1[j]=temp2[j] / ( grid->SVD_S[i][j] );
			}
	      
		      for ( j=0; j < tempdim; j++ )
			temp2[j] = 0.;
		      for ( j=0; j < tempdim; j++ )
			temp3[j] = 0.;

	      
		      for ( j=0; j < n; j++ )
			{
			  for( k=0; k < n; k++ )
			    {
			      temp3[k]= grid->SVD_V[i][j*n+k];
			    }
			  sol[j] = dotProductGeneral(temp3,temp1,n);
			}
		      
		      hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 0] = sol[1];
		      hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 1] = sol[2];
		      hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 2] = sol[3];
		      hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 3] = sol[4];
		      hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 4] = sol[5];
		      
		      point_value[i*NUM_VAR+var] = sol[0];
		    }
		}
	    }
	  else // 4th
	    {
	      for ( i=1; i <= nn; i++ )
		{
		  if ( grid->node_order[i] < 3 )
		    continue;

		  // set some variables.
		  m = grid->SVD_size[i*3+0];
		  n = grid->SVD_size[i*3+1];
		  lds = grid->SVD_size[i*3+2];

		  // loop over all the variables and compute the repsective derivatives.
		  for( var=0; var < NUM_VAR; var++)
		    {
		      U[0] = q[i*NUM_VAR+var];
		      
		      if (DEBUG_DERIVS)
			U[0] = cv_avg[i];
		      
		      start_index = gnnsn[ (gnnsn1[i]+1) ];
		      end_index   = gnnsn[ gnnsn1[i+1] ];

		      for ( j=start_index; j < end_index; j++ )
			{
			  node = gnsn[j];

			  // Get the geometric weight term for current node.
			  dx[0] = x[node*NDIM+0] - x[i*NDIM+0];
			  dx[1] = x[node*NDIM+1] - x[i*NDIM+1];
			  
			  dist = sqrt(dx[0]*dx[0] + dx[1]*dx[1]);
		      
			  w = 1.0 / pow(dist,p);
			  
			  row = j - start_index + 1;

			  U[row] = w*q[node*NUM_VAR+var];
			  
			  if (DEBUG_DERIVS)
			    U[row] = w*cv_avg[node];
			  
			}

		      // Perform the linear algebra for the system Ax=b -> (SVD)x=b
		      for ( j=0; j < m; j++ )
			{
			  for( k=0; k < m; k++ )
			    {
			      temp1[k]=grid->SVD_U[i][j*m+k];
			    }
			  temp2[j]=dotProductGeneral(U,temp1,m);
			}
		      
		      for ( j=0; j < tempdim; j++ )
			temp1[j] = 0.;
	      
		      for ( j=0; j < lds; j++ )
			{
			  temp1[j]=temp2[j] / ( grid->SVD_S[i][j] );
			}
	      
		      for ( j=0; j < tempdim; j++ )
			temp2[j] = 0.;
		      for ( j=0; j < tempdim; j++ )
			temp3[j] = 0.;

	      
		      for ( j=0; j < n; j++ )
			{
			  for( k=0; k < n; k++ )
			    {
			      temp3[k]= grid->SVD_V[i][j*n+k];
			    }
			  sol[j] = dotProductGeneral(temp3,temp1,n);
			}
		      
		      hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 0] = sol[1];
		      hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 1] = sol[2];
		      hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 2] = sol[3];
		      hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 3] = sol[4];
		      hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 4] = sol[5];
		      hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 5] = sol[6];
		      hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 6] = sol[7];
		      hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 7] = sol[8];
		      hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 8] = sol[9];
		      
		      point_value[i*NUM_VAR+var] = sol[0];
		    }
		}
	    }
	}
      else // reduced system solve.
	{
	  if ( order_desired == 1 || order_desired == 2 )  // 2nd
	    {
	      for ( i=1; i <= nn; i++ )
		{
		  // set some variables.
		  m = grid->SVD_size[i*3+0];
		  n = grid->SVD_size[i*3+1];
		  lds = grid->SVD_size[i*3+2];

		  // loop over all the variables and compute the repsective derivatives.
		  for( var=0; var < NUM_VAR; var++)
		    {
		      start_index = gnnsn[ (gnnsn1[i]+1) ];
		      end_index   = gnnsn[ gnnsn1[i+1] ];
		      
		      for ( j=start_index; j < end_index; j++ )
			{
			  node = gnsn[j];
			  
			  // Get the geometric weight term for current node.
			  dx[0] = x[node*NDIM+0] - x[i*NDIM+0];
			  dx[1] = x[node*NDIM+1] - x[i*NDIM+1];
			  
			  dist = sqrt(dx[0]*dx[0] + dx[1]*dx[1]);
			  
			  w = 1.0 / pow(dist,p);
			  
			  row = j - start_index;

			  U[row] = w*(q[node*NUM_VAR+var] - q[i*NUM_VAR+var]);
			  
			  if (DEBUG_DERIVS)
			    U[row] = w*(cv_avg[node]-cv_avg[i]);
			  
			}

		      // Perform the linear algebra for the system Ax=b -> (SVD)x=b
		      for ( j=0; j < m; j++ )
			{
			  for( k=0; k < m; k++ )
			    {
			      temp1[k]=grid->SVD_U[i][j*m+k];
			    }
			  temp2[j]=dotProductGeneral(U,temp1,m);
			}
		      
		      for ( j=0; j < tempdim; j++ )
			temp1[j] = 0.;
	      
		      for ( j=0; j < lds; j++ )
			{
			  temp1[j]=temp2[j] / ( grid->SVD_S[i][j] );
			}
	      
		      for ( j=0; j < tempdim; j++ )
			temp2[j] = 0.;
		      for ( j=0; j < tempdim; j++ )
			temp3[j] = 0.;

	      
		      for ( j=0; j < n; j++ )
			{
			  for( k=0; k < n; k++ )
			    {
			      temp3[k]= grid->SVD_V[i][j*n+k];
			    }
			  sol[j] = dotProductGeneral(temp3,temp1,n);
			}
		      
		      hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 0] = sol[0];
		      hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 1] = sol[1];
		      
		      point_value[i*NUM_VAR+var] = q[i*NUM_VAR+var] - 
			                           moments[i*NUM_MOM+0]*sol[0] -
			                           moments[i*NUM_MOM+1]*sol[1];

		      if ( DEBUG_DERIVS )
			{
			  point_value[i*NUM_VAR+var] = cv_avg[i] - 
			                               moments[i*NUM_MOM+0]*sol[0] -
			                               moments[i*NUM_MOM+1]*sol[1];
			}
		    }
		}
	    }
	  else if ( order_desired == 3 )
	    {
	      for ( i=1; i <= nn; i++ )
		{
		  if ( grid->node_order[i] < 3 )
		    continue;

		  // set some variables.
		  m = grid->SVD_size[i*3+0];
		  n = grid->SVD_size[i*3+1];
		  lds = grid->SVD_size[i*3+2];

		  // loop over all the variables and compute the repsective derivatives.
		  for( var=0; var < NUM_VAR; var++)
		    {
		      start_index = gnnsn[ (gnnsn1[i]+1) ];
		      end_index   = gnnsn[ gnnsn1[i+1] ];
		      
		      for ( j=start_index; j < end_index; j++ )
			{
			  node = gnsn[j];
			  
			  // Get the geometric weight term for current node.
			  dx[0] = x[node*NDIM+0] - x[i*NDIM+0];
			  dx[1] = x[node*NDIM+1] - x[i*NDIM+1];
			  
			  dist = sqrt(dx[0]*dx[0] + dx[1]*dx[1]);
			  
			  w = 1.0 / pow(dist,p);
			  
			  row = j - start_index;

			  U[row] = w*(q[node*NUM_VAR+var] - q[i*NUM_VAR+var]);
			  
			  if (DEBUG_DERIVS)
			    U[row] = w*(cv_avg[node]-cv_avg[i]);
			  
			}

		      // Perform the linear algebra for the system Ax=b -> (SVD)x=b
		      for ( j=0; j < m; j++ )
			{
			  for( k=0; k < m; k++ )
			    {
			      temp1[k]=grid->SVD_U[i][j*m+k];
			    }
			  temp2[j]=dotProductGeneral(U,temp1,m);
			}
		      
		      for ( j=0; j < tempdim; j++ )
			temp1[j] = 0.;
	      
		      for ( j=0; j < lds; j++ )
			{
			  temp1[j]=temp2[j] / ( grid->SVD_S[i][j] );
			}
	      
		      for ( j=0; j < tempdim; j++ )
			temp2[j] = 0.;
		      for ( j=0; j < tempdim; j++ )
			temp3[j] = 0.;

	      
		      for ( j=0; j < n; j++ )
			{
			  for( k=0; k < n; k++ )
			    {
			      temp3[k]= grid->SVD_V[i][j*n+k];
			    }
			  sol[j] = dotProductGeneral(temp3,temp1,n);
			}
		      
		      hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 0] = sol[0];
		      hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 1] = sol[1];
		      hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 2] = sol[2];
		      hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 3] = sol[3];
		      hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 4] = sol[4];
		      
		      point_value[i*NUM_VAR+var] = q[i*NUM_VAR+var] - 
		                                   moments[i*NUM_MOM+0]*sol[0] -
		                                   moments[i*NUM_MOM+1]*sol[1] -
		                                   moments[i*NUM_MOM+2]*sol[2]*0.5 -
		                                   moments[i*NUM_MOM+3]*sol[3]*0.5 -
		                                   moments[i*NUM_MOM+4]*sol[4];

		      if ( DEBUG_DERIVS )
			{
			  point_value[i*NUM_VAR+var] = cv_avg[i] -
			                               moments[i*NUM_MOM+0]*sol[0] -
		                                       moments[i*NUM_MOM+1]*sol[1] -
		                                       moments[i*NUM_MOM+2]*sol[2]*0.5 -
		                                       moments[i*NUM_MOM+3]*sol[3]*0.5 -
		                                       moments[i*NUM_MOM+4]*sol[4];
			}
		    }
		}
	    }
	  else // 4th
	    {
	      for ( i=1; i <= nn; i++ )
		{
		  if ( grid->node_order[i] < 3 )
		    continue;

		  // set some variables.
		  m = grid->SVD_size[i*3+0];
		  n = grid->SVD_size[i*3+1];
		  lds = grid->SVD_size[i*3+2];

		  // loop over all the variables and compute the repsective derivatives.
		  for( var=0; var < NUM_VAR; var++)
		    {
		      start_index = gnnsn[ (gnnsn1[i]+1) ];
		      end_index   = gnnsn[ gnnsn1[i+1] ];
		      
		      for ( j=start_index; j < end_index; j++ )
			{
			  node = gnsn[j];
			  
			  // Get the geometric weight term for current node.
			  dx[0] = x[node*NDIM+0] - x[i*NDIM+0];
			  dx[1] = x[node*NDIM+1] - x[i*NDIM+1];
			  
			  dist = sqrt(dx[0]*dx[0] + dx[1]*dx[1]);
			  
			  w = 1.0 / pow(dist,p);
			  
			  row = j - start_index;

			  U[row] = w*(q[node*NUM_VAR+var] - q[i*NUM_VAR+var]);
			  
			  if (DEBUG_DERIVS)
			    U[row] = w*(cv_avg[node]-cv_avg[i]);
			  
			}

		      if ( i == 0 && var == 0 )
			{
			  for ( j=0; j < m; j++ )
			    {
			      printf("%.15e\n", U[j]);
			    }
			}

		      // Perform the linear algebra for the system Ax=b -> (SVD)x=b
		      for ( j=0; j < m; j++ )
			{
			  for( k=0; k < m; k++ )
			    {
			      temp1[k]=grid->SVD_U[i][j*m+k];
			    }
			  temp2[j]=dotProductGeneral(U,temp1,m);
			}
		      
		      for ( j=0; j < tempdim; j++ )
			temp1[j] = 0.;
	      
		      for ( j=0; j < lds; j++ )
			{
			  temp1[j]=temp2[j] / ( grid->SVD_S[i][j] );
			}
	      
		      for ( j=0; j < tempdim; j++ )
			temp2[j] = 0.;
		      for ( j=0; j < tempdim; j++ )
			temp3[j] = 0.;

	      
		      for ( j=0; j < n; j++ )
			{
			  for( k=0; k < n; k++ )
			    {
			      temp3[k]= grid->SVD_V[i][j*n+k];
			    }
			  sol[j] = dotProductGeneral(temp3,temp1,n);
			}
		      
		      hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 0] = sol[0];
		      hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 1] = sol[1];
		      hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 2] = sol[2];
		      hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 3] = sol[3];
		      hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 4] = sol[4];
		      hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 5] = sol[5];
		      hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 6] = sol[6];
		      hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 7] = sol[7];
		      hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 8] = sol[8];
		      
		      point_value[i*NUM_VAR+var] = q[i*NUM_VAR+var] - 
		                                   moments[i*NUM_MOM+0]*sol[0] -
		                                   moments[i*NUM_MOM+1]*sol[1] -
		                                   moments[i*NUM_MOM+2]*sol[2]*0.5 -
		                                   moments[i*NUM_MOM+3]*sol[3]*0.5 -
		                                   moments[i*NUM_MOM+4]*sol[4] -
			                           moments[i*NUM_MOM+5]*sol[5]*(1./6.) -
			                           moments[i*NUM_MOM+6]*sol[6]*(1./6.) -
			                           moments[i*NUM_MOM+7]*sol[7]*0.5 -
			                           moments[i*NUM_MOM+8]*sol[8]*0.5;

		      if ( DEBUG_DERIVS )
			{
			  point_value[i*NUM_VAR+var] = cv_avg[i] -
			                               moments[i*NUM_MOM+0]*sol[0] -
		                                       moments[i*NUM_MOM+1]*sol[1] -
		                                       moments[i*NUM_MOM+2]*sol[2]*0.5 -
		                                       moments[i*NUM_MOM+3]*sol[3]*0.5 -
		                                       moments[i*NUM_MOM+4]*sol[4] -
			                               moments[i*NUM_MOM+5]*sol[5]*(1./6.) -
			                               moments[i*NUM_MOM+6]*sol[6]*(1./6.) -
			                               moments[i*NUM_MOM+7]*sol[7]*0.5 -
			                               moments[i*NUM_MOM+8]*sol[8]*0.5;
			}

		      if ( i == 0 && var == 0 )
			{
			  printf("solution\n");
			  for ( j=0; j < n; j++ )
			    printf("%.15e\n",sol[j]);
			  
			  printf("point value = %.15e\n",point_value[i*NUM_VAR]);
			}
		    }
		}
	      
	    }
	}
    }	  
  else
    {
      printf("FATAL ERROR: NO METHOD WAS SELECTED FOR SOLVING THE LEAST-SQUARES PROBLEM.\n");
      exit(1);
    }

  // debug check for the annulus case.
  if ( ANNULUS && 0 )
    {
      Test_Annulus_Derivatives( grid , params );
      
      printf("\n\n");
      printf("FINISHED TESTING THE ANNULUS.\n");
      printf("exiting...\n");
      fflush(stdout);
      exit(1);
    }
  
  
  if ( DEBUG_DERIVS )
    {
      FILE *fp_pv = NULL;
      fp_pv = fopen("point_values.dat","w");
      
      for ( i=1; i <= grid->nn; i++ )
	{
	  fprintf(fp_pv,"%.15e\n",point_value[i*NUM_VAR+0]);
	}
      fclose(fp_pv);
    }
  else if (0)
    {
      FILE *fp_pv1 = NULL;
      FILE *fp_pv2 = NULL;
      FILE *fp_pv3 = NULL;
      FILE *fp_pv4 = NULL;
      
      fp_pv1 = fopen("derivs_pv1.dat","w");
      fp_pv2 = fopen("derivs_pv2.dat","w");
      fp_pv3 = fopen("derivs_pv3.dat","w");
      fp_pv4 = fopen("derivs_pv4.dat","w");
      
      for ( i=1; i <= grid->nn; i++ )
	{
	  fprintf(fp_pv1,"%.15e\n",point_value[i*NUM_VAR+0]);
	  fprintf(fp_pv2,"%.15e\n",point_value[i*NUM_VAR+1]);
	  fprintf(fp_pv3,"%.15e\n",point_value[i*NUM_VAR+2]);
	  fprintf(fp_pv4,"%.15e\n",point_value[i*NUM_VAR+3]);
	}
      fclose(fp_pv1);
      fclose(fp_pv2);
      fclose(fp_pv3);
      fclose(fp_pv4);
    }
  
  if ( DEBUG_DERIVS )
    {
      FILE *fp_derivs = NULL;
      fp_derivs = fopen("derivs_values.dat","w");
      
      for ( i=1; i <= grid->nn; i++ )
	{
	  if ( order_desired == 1 || order_desired == 2 )
	    k = 2;

	  if ( order_desired == 3 )
	    k = 5;

	  if ( order_desired == 4 )
	    k = 9;

	  for ( j=0; j < k; j++ )
	    {
	      fprintf(fp_derivs,"%.15e\n",hess[i*NUM_VAR*NUM_MOM + 0*NUM_MOM + j]);
	    }
	}
      fclose(fp_derivs);
    }
  else if (0)
    {
      FILE *fp_h1 = NULL;
      FILE *fp_h2 = NULL;
      FILE *fp_h3 = NULL;
      FILE *fp_h4 = NULL;
      
      fp_h1 = fopen("derivs_q1.dat","w");
      fp_h2 = fopen("derivs_q2.dat","w");
      fp_h3 = fopen("derivs_q3.dat","w");
      fp_h4 = fopen("derivs_q4.dat","w");
      
      for ( i=1; i <= grid->nn; i++ )
	{
	  if ( order_desired == 1 || order_desired == 2 )
	    k = 2;

	  if ( order_desired == 3 )
	    k = 5;

	  if ( order_desired == 4 )
	    k = 9;
	  
	  for ( j=0; j < k; j++ )
	    {
	      fprintf(fp_h1,"%.15e\n",hess[i*NUM_VAR*NUM_MOM + 0*NUM_MOM + j]);
	      fprintf(fp_h2,"%.15e\n",hess[i*NUM_VAR*NUM_MOM + 1*NUM_MOM + j]);
	      fprintf(fp_h3,"%.15e\n",hess[i*NUM_VAR*NUM_MOM + 2*NUM_MOM + j]);
	      fprintf(fp_h4,"%.15e\n",hess[i*NUM_VAR*NUM_MOM + 3*NUM_MOM + j]);
	    }
	}
      fclose(fp_h1);
      fclose(fp_h2);
      fclose(fp_h3);
      fclose(fp_h4);
    }

  if ( (MMS && RECON_PRIM) && 0 )
    {
      double pi = M_PI;
      
      norm = 0.;
      temp = 0.;
      imax = 0;
      term = 0.;
      
      // compute the error in predicting the point value of the node.
      printf("Density\n");
      for ( i=1; i <= grid->nn; i++ )
	{
	  x1 = grid->x[i*NDIM+0];
	  y1 = grid->x[i*NDIM+1];
	  
	  term = fabs( point_value[i*NUM_VAR+0] - ( 1.0 + sin(pi*x1)*sin(pi*y1) ) );
	  norm += term;
	}
      printf("L1 Norm of point value difference = %.15E\n",norm);

      norm = 0.;
      for ( i=1; i <= grid->nn; i++ )
	{
	  x1 = grid->x[i*NDIM+0];
	  y1 = grid->x[i*NDIM+1];

	  term = fabs( point_value[i*NUM_VAR+0] - ( 1.0 + sin(pi*x1)*sin(pi*y1) ) );
	  norm += (term*term);
	}
      printf("L2 Norm of point value difference = %.15E\n",norm);
      
      x1 = grid->x[1*NDIM+0];
      y1 = grid->x[1*NDIM+1];
      norm = fabs( point_value[1*NUM_VAR+0] - ( 1.0 + sin(pi*x1)*sin(pi*y1) ) );
      imax = 1;
      
      for ( i=2; i <= grid->nn; i++ )
	{
	  x1 = grid->x[i*NDIM+0];
	  y1 = grid->x[i*NDIM+1];
	  
	  term = fabs( point_value[i*NUM_VAR+0] - ( 1.0 + sin(pi*x1)*sin(pi*y1) ) );
	  if ( term > norm )
	    {
	      norm = term;
	      imax = i;
	    }
	}
      printf("Infinity norm of the point values difference = %.15E\n",norm);
      printf("This occured at node %d: %f   %f\n",imax,(float)x[imax*NDIM+0],(float)x[imax*NDIM+1]);
      
      printf("X velocity: u\n");
      for ( i=1; i <= grid->nn; i++ )
	{
	  x1 = grid->x[i*NDIM+0];
	  y1 = grid->x[i*NDIM+1];

	  term = fabs( point_value[i*NUM_VAR+1] - ( 0.5 + 0.5*sin(pi*x1)*cos(2.*pi*y1) ) );
	  norm += term;
	}
      printf("L1 Norm of point value difference = %.15E\n",norm);

      norm = 0.;
      for ( i=1; i <= grid->nn; i++ )
	{
	  x1 = grid->x[i*NDIM+0];
	  y1 = grid->x[i*NDIM+1];

	  term = fabs( point_value[i*NUM_VAR+1] - ( 0.5 + 0.5*sin(pi*x1)*cos(2.*pi*y1) ) );
	  norm += (term*term);
	}
      printf("L2 Norm of point value difference = %.15E\n",norm);

      x1 = grid->x[1*NDIM+0];
      y1 = grid->x[1*NDIM+1];
      norm = fabs( point_value[1*NUM_VAR+1] - ( 0.5 + 0.5*sin(pi*x1)*cos(2.*pi*y1) ) );
      imax = 1;

      for ( i=2; i <= grid->nn; i++ )
	{
	  x1 = grid->x[i*NDIM+0];
	  y1 = grid->x[i*NDIM+1];

	  term = fabs( point_value[i*NUM_VAR+1] - ( 0.5 + 0.5*sin(pi*x1)*cos(2.*pi*y1) ) );
	  if ( term > norm )
	    {
	      norm = term;
	      imax = i;
	    }
	}
      printf("Infinity norm of the point values difference = %.15E\n",norm);
      printf("This occured at node %d: %f   %f\n",imax,(float)x[imax*NDIM+0],(float)x[imax*NDIM+1]);


      printf("Y velocity: v\n");
      for ( i=1; i <= grid->nn; i++ )
	{
	  x1 = grid->x[i*NDIM+0];
	  y1 = grid->x[i*NDIM+1];

	  term = fabs( point_value[i*NUM_VAR+2] - ( 0.5 + 0.5*cos(2.*pi*x1)*sin(pi*y1) ) );
	  norm += term;
	}
      printf("L1 Norm of point value difference = %.15E\n",norm);

      norm = 0.;
      for ( i=1; i <= grid->nn; i++ )
	{
	  x1 = grid->x[i*NDIM+0];
	  y1 = grid->x[i*NDIM+1];

	  term = fabs( point_value[i*NUM_VAR+2] - ( 0.5 + 0.5*cos(2.*pi*x1)*sin(pi*y1) ) );
	  norm += (term*term);
	}
      printf("L2 Norm of point value difference = %.15E\n",norm);

      x1 = grid->x[1*NDIM+0];
      y1 = grid->x[1*NDIM+1];
      norm = fabs( point_value[1*NUM_VAR+2] - ( 0.5 + 0.5*cos(2.*pi*x1)*sin(pi*y1) ) );
      imax = 1;

      for ( i=2; i <= grid->nn; i++ )
	{
	  x1 = grid->x[i*NDIM+0];
	  y1 = grid->x[i*NDIM+1];

	  term = fabs( point_value[i*NUM_VAR+2] - ( 0.5 + 0.5*cos(2.*pi*x1)*sin(pi*y1) ) );
	  if ( term > norm )
	    {
	      norm = term;
	      imax = i;
	    }
	}
      printf("Infinity norm of the point values difference = %.15E\n",norm);
      printf("This occured at node %d: %f   %f\n",imax,(float)x[imax*NDIM+0],(float)x[imax*NDIM+1]);

      printf("Pressure:\n");
      for ( i=1; i <= grid->nn; i++ )
	{
	  x1 = grid->x[i*NDIM+0];
	  y1 = grid->x[i*NDIM+1];

	  term = fabs( point_value[i*NUM_VAR+3] - ( 1.0/1.4 + 0.05*cos(2.*pi*x1)*cos(2.*pi*y1) ) );
	  norm += term;
	}
      printf("L1 Norm of point value difference = %.15E\n",norm);

      norm = 0.;
      for ( i=1; i <= grid->nn; i++ )
	{
	  x1 = grid->x[i*NDIM+0];
	  y1 = grid->x[i*NDIM+1];

	  term = fabs( point_value[i*NUM_VAR+3] - ( 1.0/1.4 + 0.05*cos(2.*pi*x1)*cos(2.*pi*y1) ) );
	  norm += (term*term);
	}
      printf("L2 Norm of point value difference = %.15E\n",norm);

      x1 = grid->x[1*NDIM+0];
      y1 = grid->x[1*NDIM+1];
      norm = fabs( point_value[1*NUM_VAR+3] - ( 1.0/1.4 + 0.05*cos(2.*pi*x1)*cos(2.*pi*y1) ) );
      imax = 1;

      for ( i=2; i <= grid->nn; i++ )
	{
	  x1 = grid->x[i*NDIM+0];
	  y1 = grid->x[i*NDIM+1];

	  term = fabs( point_value[i*NUM_VAR+3] - ( 1.0/1.4 + 0.05*cos(2.*pi*x1)*cos(2.*pi*y1) ) );
	  if ( term > norm )
	    {
	      norm = term;
	      imax = i;
	    }
	}
      printf("Infinity norm of the point values difference = %.15E\n",norm);
      printf("This occured at node %d: %f   %f\n",imax,(float)x[imax*NDIM+0],(float)x[imax*NDIM+1]);
      
    }

  //#if DEBUG_DERIVS
  if (0)
    {
      printf("DERIVS DEBUG:\n");
      for ( i=1; i <= grid->nn; i++ )
	{
	  printf("Node %d:\n",i);
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      printf("  Var %d: qx = %f\n",j,(float)grid->hess[i*NUM_MOM*NUM_VAR + j*NUM_MOM + 0]);
	      printf("          qy = %f\n",(float)grid->hess[i*NUM_MOM*NUM_VAR + j*NUM_MOM + 1]);

	      if ( order_desired > 2 )
		{
		  printf("         qxx = %f\n",(float)grid->hess[i*NUM_MOM*NUM_VAR + j*NUM_MOM + 2]);
		  printf("         qyy = %f\n",(float)grid->hess[i*NUM_MOM*NUM_VAR + j*NUM_MOM + 3]);
		  printf("         qxy = %f\n",(float)grid->hess[i*NUM_MOM*NUM_VAR + j*NUM_MOM + 4]);
		}

	      if ( order_desired > 3 )
		{
		  printf("         qxxx = %f\n",(float)grid->hess[i*NUM_MOM*NUM_VAR + j*NUM_MOM + 5]);
		  printf("         qyyy = %f\n",(float)grid->hess[i*NUM_MOM*NUM_VAR + j*NUM_MOM + 6]);
		  printf("         qxxy = %f\n",(float)grid->hess[i*NUM_MOM*NUM_VAR + j*NUM_MOM + 7]);
		  printf("         qxyy = %f\n",(float)grid->hess[i*NUM_MOM*NUM_VAR + j*NUM_MOM + 8]);
		}
	    }
	}
      fflush(stdout);
    }


  // Now we may do some analysis on the result for debug cases.
  if ( DEBUG_DERIVS )
    {
      norm = 0.;  term = 0.;

      // Test the point value convergence.
      for ( i=1; i <= grid->nn; i++ )
	{
	  term = fabs( point_value[i*NUM_VAR+0] - Test_Function(&(grid->x[i*NDIM])) );
	  norm += term;
	}
      printf("L1 Norm of point value difference = %.15E\n",norm);
      
      norm = 0.;
      for ( i=1; i <= grid->nn; i++ )
	{
	  term = fabs( point_value[i*NUM_VAR+0] - Test_Function(&(grid->x[i*NDIM])) );
	  norm += (term*term);
	}
      printf("L2 Norm of point value difference = %.15E\n",norm);
      
      norm = fabs( point_value[1*NUM_VAR+0] - Test_Function(&(grid->x[1*NDIM])) );
      imax = 1;
      
      for ( i=2; i <= grid->nn; i++ )
	{
	  term = fabs( point_value[i*NUM_VAR+0] - Test_Function(&(grid->x[i*NDIM])) );
	  if ( term > norm )
	    {
	      norm = term;
	      imax = i;
	    }
	}
      printf("Infinity norm of the point values difference = %.15E\n",norm);
      printf("This occured at node %d: %f   %f\n",imax,(float)x[imax*NDIM+0],(float)x[imax*NDIM+1]);
      
      // Compare what the Compute_Derivative routine returned against
      // what I'm expecting it to be. Print out any large differences.
      // Print out the norm.

      qprime = (double*)malloc((nn+1)*NUM_MOM*sizeof(double));
      if ( qprime == NULL ) { printf("MEMORY ERROR: Could not allocate 'qprime'.\n"); exit(0); }
      
      // Set the analytical values.
      for ( i=1; i <= nn; i++ )
	{
	  Test_Function_Derivatives(&(x[i*NDIM]),&(qprime[i*NUM_MOM]));
	}
      
      for ( i=0; i < NUM_MOM; i++ )
	{
	  printf("Checking %s derivatives...  ",derivs[i]);
	  count=0;
	  for( j=1; j <= nn; j++ )
	    {
	      if ( fabs( qprime[j*NUM_MOM+i] - hess[j*NUM_VAR*NUM_MOM+0*NUM_MOM+i] ) > 1.E-10 )
		{
		  count++;
		}
	    }
	  printf("Found %d problems!\n",count);
	}
      
      // Find the l1 norm of the differences.
      for ( i=0; i < NUM_MOM; i++ )
	{
	  printf("L1 norm of %s derivative difference = ",derivs[i]);
	  norm = 0.;
	  for( j=1; j <= nn; j++ )
	    {
	      term = fabs( qprime[j*NUM_MOM+i] - hess[j*NUM_VAR*NUM_MOM+0*NUM_MOM+i] );
	      norm += term;
	    }
	  printf("%.15E\n",norm);
	}

      // Find the l2 norm of the differences.
      for ( i=0; i < NUM_MOM; i++ )
	{
	  printf("L2 norm of %s derivative difference = ",derivs[i]);
	  norm = 0.;
	  for( j=1; j <= nn; j++ )
	    {
	      term = fabs( qprime[j*NUM_MOM+i] - hess[j*NUM_VAR*NUM_MOM+0*NUM_MOM+i] );
	      norm += (term*term);
	    }
	  printf("%.15E\n",sqrt(norm));
	}
      
      // Find the infinity norm of the differences.
      for( i=0; i < NUM_MOM; i++ )
	{
	  printf("Infinty norm of %s derivative difference = ",derivs[i]);
	  norm = fabs( qprime[1*NUM_MOM+i] - hess[1*NUM_VAR*NUM_MOM+0*NUM_MOM+i] );
	  imax = 1;
	  for( j=1; j <= nn; j++ )
	    {
	      term = fabs( qprime[j*NUM_MOM+i] - hess[j*NUM_VAR*NUM_MOM+0*NUM_MOM+i] );
	      if ( term > norm )
		{
		  norm = term;
		  imax = j;
		}
	    }
	  printf("%.15E\n",norm);
	  printf("This occured at node %d: %f   %f\n",imax,(float)x[imax*NDIM+0],(float)x[imax*NDIM+1]);
	}

      // Find the l1 norm of the derivative vector.
      norm = 0.;
      for( j=1; j <= nn; j++ )
	{
	  for ( i=0; i < NUM_MOM; i++ )
	    {
	      term = fabs( qprime[j*NUM_MOM+i] - hess[j*NUM_VAR*NUM_MOM+0*NUM_MOM+i] );
	      norm += term;
	    }
	}
      printf("L1 Norm of Derivative Difference Vector = %.15E\n",norm);
      
      // Find the l2 norm of the derivative vector.
      norm = 0.;
      for( j=1; j <= nn; j++ )
	{
	  for ( i=0; i < NUM_MOM; i++ )
	    {
	      term = fabs( qprime[j*NUM_MOM+i] - hess[j*NUM_VAR*NUM_MOM+0*NUM_MOM+i] );
	      norm += (term*term);
	    }
	}
      printf("L2 Norm of Derivative Difference Vector = %.15E\n",sqrt(norm));
      
      // Find the infinity norm of the derivative vector.
      norm = fabs( qprime[1*NUM_MOM+0] - hess[1*NUM_VAR*NUM_MOM+0*NUM_MOM+0] );
      imax = 1;
      k = 0;
      for ( i=1; i < NUM_MOM; i++ )
	{
	  term = fabs( qprime[1*NUM_MOM+i] - hess[1*NUM_VAR*NUM_MOM+0*NUM_MOM+i] );
	  if ( term > norm )
	    {
	      norm = term;
	      k = i;
	    }
	}
      // The above finds the first max in the first node.
      
      // Now search through the rest.
      for ( j=2; j <= nn; j++ )
	{
	  for ( i=0; i < NUM_MOM; i++ )
	    {
	      term = fabs( qprime[j*NUM_MOM+i] - hess[j*NUM_VAR*NUM_MOM+0*NUM_MOM+i] );
	      if ( term > norm )
		{
		  norm = term;
		  imax = j;
		  k = i;
		}
	    }
	}
      printf("Infinity norm of the derivative difference vector = %.15E\n",norm);
      printf("This occured at node %d: %f   %f\n",imax,(float)x[imax*NDIM+0],(float)x[imax*NDIM+1]);
      printf("In the %s derivative.\n",derivs[k]);

      printf("\n");
      printf("Reconstruction test...\n");
      
      // Loop over all the subedges and reconstruct to the Gauss points from both directions
      // and look at the error created.
      
      norm_1_s = 0.;
      norm_2_s = 0.;
      norm_i_s = 0.;
      
      t1 = sqrt(15.0)/5.0;
      t2 = 0.;
      t3 = sqrt(15.0)/(-5.0);

      double norm_rm = 0.0;
      FILE *fp1 = NULL;
      fp1 = fopen("debug_gauss_points_extrap_4th.dat","w");

      
      for ( i=1; i <= grid->nsubedges; i++ )
	{
	  iedge = grid->subedges[i*2];
	  
	  nodeL = grid->edges[iedge*2 + 0];
	  nodeR = grid->edges[iedge*2 + 1];
	  
	  xmid = grid->xm_subedges[i*NDIM+0];
	  ymid = grid->xm_subedges[i*NDIM+1];
	  
	  gelem = grid->subedges[i*2+1];
	  xc[0] = grid->el_cent[gelem*2];
	  xc[1] = grid->el_cent[gelem*2+1];

	  x1 = (xmid+xc[0])*0.5 + (xc[0]-xmid)*0.5*t1;
	  y1 = (ymid+xc[1])*0.5 + (xc[1]-ymid)*0.5*t1;
      
	  x2 = (xmid+xc[0])*0.5 + (xc[0]-xmid)*0.5*t2;
	  y2 = (ymid+xc[1])*0.5 + (xc[1]-ymid)*0.5*t2;
      
	  x3 = (xmid+xc[0])*0.5 + (xc[0]-xmid)*0.5*t3;
	  y3 = (ymid+xc[1])*0.5 + (xc[1]-ymid)*0.5*t3;

	  // LEFT SIDE FIRST.
	  
	  // First gp.
	  dx[0] = x1 - x[nodeL*NDIM+0];
	  dx[1] = y1 - x[nodeL*NDIM+1];
	  gp[0][0] = x1;  gp[0][1] = y1;

	  Reconstruct_Gauss ( i, nodeL, nodeR, gp[0], QL, QR, grid, params, 1 );

	  recon_val = QL[0];

	  fprintf(fp1,"%.15e\n",recon_val);

	  temp = Test_Function ( gp[0] );

	  temp = fabs( recon_val - temp );
	  
	  norm_1_s += temp;
	  norm_2_s += (temp*temp);
	  if ( temp > norm_i_s )
	    {
	      norm_i_s = temp;
	      inf_subedge = i;
	    }
	  norm_rm += fabs(recon_val - QL[0]);

	  // Second gp.
	  dx[0] = x2 - x[nodeL*NDIM+0];
	  dx[1] = y2 - x[nodeL*NDIM+1];
	  gp[0][0] = x2;  gp[0][1] = y2;
	  
	  Reconstruct_Gauss ( i, nodeL, nodeR, gp[0], QL, QR, grid, params, 1 );
	  recon_val = QL[0];

	  fprintf(fp1,"%.15e\n",recon_val);
	  
	  temp = Test_Function ( gp[0] );

	  temp = fabs( recon_val - temp );

	  norm_1_s += temp;
	  norm_2_s += (temp*temp);
	  if ( temp > norm_i_s )
	    {
	      norm_i_s = temp;
	      inf_subedge = i;
	    }
	  norm_rm += fabs(recon_val - QL[0]);

	  // Third gp.
	  dx[0] = x3 - x[nodeL*NDIM+0];
	  dx[1] = y3 - x[nodeL*NDIM+1];
	  gp[0][0] = x3;  gp[0][1] = y3;

	  Reconstruct_Gauss ( i, nodeL, nodeR, gp[0], QL, QR, grid, params, 1 );
	  recon_val = QL[0];

	  fprintf(fp1,"%.15e\n",recon_val);
	  
	  temp = Test_Function ( gp[0] );
	  
	  temp = fabs( recon_val - temp );
	  
	  norm_1_s += temp;
	  norm_2_s += (temp*temp);
	  if ( temp > norm_i_s )
	    {
	      norm_i_s = temp;
	      inf_subedge = i;
	    }
	  norm_rm += fabs(recon_val - QL[0]);
	  
	  // RIGHT NODE.
	  
	  // First gp.
	  dx[0] = x1 - x[nodeR*NDIM+0];
	  dx[1] = y1 - x[nodeR*NDIM+1];
	  gp[0][0] = x1;  gp[0][1] = y1;
	  
	  Reconstruct_Gauss ( i, nodeL, nodeR, gp[0], QL, QR, grid, params, 1 );
	  recon_val = QR[0];
	  
	  fprintf(fp1,"%.15e\n",recon_val);

	  temp = Test_Function ( gp[0] );

	  temp = fabs( recon_val - temp );
	  
	  norm_1_s += temp;
	  norm_2_s += (temp*temp);
	  if ( temp > norm_i_s )
	    {
	      norm_i_s = temp;
	      inf_subedge = i;
	    }
	  norm_rm += fabs(recon_val - QR[0]);
	  
	  // Second gp.
	  dx[0] = x2 - x[nodeR*NDIM+0];
	  dx[1] = y2 - x[nodeR*NDIM+1];
	  gp[0][0] = x2;  gp[0][1] = y2;

	  Reconstruct_Gauss ( i, nodeL, nodeR, gp[0], QL, QR, grid, params, 1 );
	  recon_val = QR[0];

	  fprintf(fp1,"%.15e\n",recon_val);
	  
	  temp = Test_Function ( gp[0] );

	  temp = fabs( recon_val - temp );

	  norm_1_s += temp;
	  norm_2_s += (temp*temp);
	  if ( temp > norm_i_s )
	    {
	      norm_i_s = temp;
	      inf_subedge = i;
	    }
	  norm_rm += fabs(recon_val - QR[0]);

	  // Third gp.
	  dx[0] = x3 - x[nodeR*NDIM+0];
	  dx[1] = y3 - x[nodeR*NDIM+1];
	  gp[0][0] = x3;  gp[0][1] = y3;

	  Reconstruct_Gauss ( i, nodeL, nodeR, gp[0], QL, QR, grid, params, 1 );
	  recon_val = QR[0];

	  fprintf(fp1,"%.15e\n",recon_val);

	  temp = Test_Function ( gp[0] );

	  temp = fabs( recon_val - temp );

	  norm_1_s += temp;
	  norm_2_s += (temp*temp);
	  if ( temp > norm_i_s )
	    {
	      norm_i_s = temp;
	      inf_subedge = i;
	    }
	  norm_rm += fabs(recon_val - QR[0]);
	}

      fclose(fp1);
      // Loop over all the boundary edges and reconstruct to the Gauss points and look at the error created.
      
      norm_1_b = 0.;
      norm_2_b = 0.;
      norm_i_b = 0.;
      
      fp1 = fopen("debug_gauss_points_extrap_boundary_4th.dat","w");

      for ( i=1; i <= grid->nbedges; i++ )
	{
	  // Grab the data from the structure.
	  real_node = grid->bedges[i*5+0];
	  ghost_node = grid->bedges[i*5+1] + grid->nn;
	  
	  b = grid->bedges[i*5+3];
	  bc = grid->bbc[b];

	  // Get the segment.
	  seg = grid->bedges[i*5+4];
	  
	  // Now get the nodes attached to the edge.
	  nodeL = grid->bs[b][seg][0];
	  nodeR = grid->bs[b][seg][1];
	  
	  // Get the node coordinates.
	  xi = grid->x[real_node*NDIM+0];
	  yi = grid->x[real_node*NDIM+1];
	  
	  // Get the edge midpoint.
	  xmid = 0.5*( grid->x[nodeL*NDIM+0] + grid->x[nodeR*NDIM+0] );
	  ymid = 0.5*( grid->x[nodeL*NDIM+1] + grid->x[nodeR*NDIM+1] );
	  
	  if ( bc <= 10 ) // Linear boundary.
	    {
	      if ( real_node == nodeL )
		{
		  xL = xi; xR = xmid;
		  yL = yi; yR = ymid;
		}
	      else
		{
		  xL = xmid; xR = xi;
		  yL = ymid; yR = yi;
		}
	      
	      t1 = sqrt(15.0)/5.0;
	      t2 = 0.;
	      t3 = sqrt(15.0)/(-5.0);

	      x1 = (xL+xR)*0.5 + (xR-xL)*0.5*t1;
	      y1 = (yL+yR)*0.5 + (yR-yL)*0.5*t1;
      
	      x2 = (xL+xR)*0.5 + (xR-xL)*0.5*t2;
	      y2 = (yL+yR)*0.5 + (yR-yL)*0.5*t2;

	      x3 = (xL+xR)*0.5 + (xR-xL)*0.5*t3;
	      y3 = (yL+yR)*0.5 + (yR-yL)*0.5*t3;

	      // First gp.
	      dx[0] = x1 - x[real_node*NDIM+0];
	      dx[1] = y1 - x[real_node*NDIM+1];
	      gp[0][0] = x1;  gp[0][1] = y1;

	      Reconstruct_Gauss_Boundary ( real_node, gp[0], QL, grid, params);

	      recon_val = QL[0];

	      fprintf(fp1,"%.15e\n",recon_val);
	      
	      temp = Test_Function ( gp[0] );
	      
	      temp = fabs( recon_val - temp );
	  
	      norm_1_b += temp;
	      norm_2_b += (temp*temp);
	      if ( temp > norm_i_b )
		{
		  norm_i_b = temp;
		  inf_bedge = i;
		}
	      norm_rm += fabs(recon_val - QL[0]);
	      
	      // Second gp.
	      dx[0] = x2 - x[real_node*NDIM+0];
	      dx[1] = y2 - x[real_node*NDIM+1];
	      gp[0][0] = x2;  gp[0][1] = y2;

	      Reconstruct_Gauss_Boundary ( real_node, gp[0], QL, grid, params);

	      recon_val = QL[0];

	      fprintf(fp1,"%.15e\n",recon_val);
	      
	      temp = Test_Function ( gp[0] );

	      temp = fabs( recon_val - temp );

	      norm_1_b += temp;
	      norm_2_b += (temp*temp);
	      if ( temp > norm_i_b )
		{
		  norm_i_b = temp;
		  inf_bedge = i;
		}
	      norm_rm += fabs(recon_val - QL[0]);

	      // Third gp.
	      dx[0] = x3 - x[real_node*NDIM+0];
	      dx[1] = y3 - x[real_node*NDIM+1];
	      gp[0][0] = x3;  gp[0][1] = y3;

	      Reconstruct_Gauss_Boundary ( real_node, gp[0], QL, grid, params);

	      recon_val = QL[0];

	      fprintf(fp1,"%.15e\n",recon_val);
	      
	      temp = Test_Function ( gp[0] );
	      
	      temp = fabs( recon_val - temp );

	      norm_1_b += temp;
	      norm_2_b += (temp*temp);
	      if ( temp > norm_i_b )
		{
		  norm_i_b = temp;
		  inf_bedge = i;
		}
	      norm_rm += fabs(recon_val - QL[0]);
	    }
	  else // curved boundaries.
	    {
	      if ( real_node == nodeL )
		{
		  XL[0] = grid->x[real_node*NDIM+0];
		  XL[1] = grid->x[real_node*NDIM+1];
		  XR[0] = grid->x[ghost_node*NDIM+0];
		  XR[1] = grid->x[ghost_node*NDIM+1];
		}
	      else
		{
		  XR[0] = grid->x[real_node*NDIM+0];
		  XR[1] = grid->x[real_node*NDIM+1];
		  XL[0] = grid->x[ghost_node*NDIM+0];
		  XL[1] = grid->x[ghost_node*NDIM+1];
		}
	      
	      // Get the Gauss points.
	      curved_boundary_gauss_points ( grid, bc, XL, XR, GP );
	      
	      // First gp.
	      dx[0] = GP[0] - x[real_node*NDIM+0];
	      dx[1] = GP[1] - x[real_node*NDIM+1];
	      gp[0][0] = GP[0];  gp[0][1] = GP[1];

	      Reconstruct_Gauss_Boundary ( real_node, gp[0], QL, grid, params);
	      
	      recon_val = QL[0];

	      temp = Test_Function ( gp[0] );
	      
	      temp = fabs( recon_val - temp );
	  
	      norm_1_b += temp;
	      norm_2_b += (temp*temp);
	      if ( temp > norm_i_b )
		{
		  norm_i_b = temp;
		  inf_bedge = i;
		}

	      // Second gp.
	      dx[0] = GP[2] - x[real_node*NDIM+0];
	      dx[1] = GP[3] - x[real_node*NDIM+1];
	      gp[0][0] = GP[2];  gp[0][1] = GP[3];

	      Reconstruct_Gauss_Boundary ( real_node, gp[0], QL, grid, params);
	      
	      recon_val = QL[0];
	      
	      temp = Test_Function ( gp[0] );
	      
	      temp = fabs( recon_val - temp );
	  
	      norm_1_b += temp;
	      norm_2_b += (temp*temp);
	      if ( temp > norm_i_b )
		{
		  norm_i_b = temp;
		  inf_bedge = i;
		}

	      // Third gp.
	      dx[0] = GP[4] - x[real_node*NDIM+0];
	      dx[1] = GP[5] - x[real_node*NDIM+1];
	      gp[0][0] = GP[4];  gp[0][1] = GP[5];

	      Reconstruct_Gauss_Boundary ( real_node, gp[0], QL, grid, params);
	      
	      recon_val = QL[0];
	      
	      temp = Test_Function ( gp[0] );
	      
	      temp = fabs( recon_val - temp );
	  
	      norm_1_b += temp;
	      norm_2_b += (temp*temp);
	      if ( temp > norm_i_b )
		{
		  norm_i_b = temp;
		  inf_bedge = i;
		}
	    }

	}

      fclose(fp1);

      // Now I can print out all this stuff to the screen.

      //printf("L1 norm of difference between reconstructions is %.15E\n",norm_rm);

      iedge = grid->subedges[inf_subedge*2];
      nodeL = grid->edges[iedge*2 + 0];
      nodeR = grid->edges[iedge*2 + 1];
      printf("L1 norm of the interior edges reconstruction error is %.15E\n",norm_1_s);
      printf("L2 norm of the interior edges reconstruction error is %.15E\n",sqrt(norm_2_s));
      printf("The infinity norm of the interior edges reconstruction error is %.15E\n",norm_i_s);
      printf("  This occured on subedge %d: nodeL = %d  nodeR = %d\n",inf_subedge,nodeL,nodeR);
      printf("  nodeL : < %f , %f >\n",(float)grid->x[nodeL*NDIM+0],(float)grid->x[nodeL*NDIM+1]);
      printf("  nodeR : < %f , %f >\n",(float)grid->x[nodeR*NDIM+0],(float)grid->x[nodeR*NDIM+1]);
      printf("\n");
      real_node = grid->bedges[inf_bedge*5+0];
      b = grid->bedges[inf_bedge*5+3];
      seg = grid->bedges[inf_bedge*5+4];
      printf("L1 norm of the boundary edges reconstruction error is %.15E\n",norm_1_b);
      printf("L2 norm of the boundary edges reconstruction error is %.15E\n",sqrt(norm_2_b));
      printf("The infinity norm of the boundary edges reconstruction error is %.15E\n",norm_i_b);
      printf("  This occured on bedge %d: real_node = %d , boundary = %d , segment = %d\n",inf_bedge,real_node,b,seg);
      printf("  real_node : < %f , %f >\n",(float)grid->x[real_node*NDIM+0],(float)grid->x[real_node*NDIM+1]);
      printf("\n");
      printf("L1 norm for all edges reconstruction error is %.15E\n",norm_1_s + norm_1_b);
      printf("L2 norm for all edges reconstruction error is %.15E\n",sqrt(norm_2_s + norm_2_b));
      if ( norm_i_s > norm_i_b )
	{
	  printf("The infinity norm for all edges reconstruction error is %.15E\n",norm_i_s);
	  printf("  This occured on subedge %d: nodeL = %d  nodeR = %d\n",inf_subedge,nodeL,nodeR);
	  printf("  nodeL : < %f , %f >\n",(float)grid->x[nodeL*NDIM+0],(float)grid->x[nodeL*NDIM+1]);
	  printf("  nodeR : < %f , %f >\n",(float)grid->x[nodeR*NDIM+0],(float)grid->x[nodeR*NDIM+1]);
	  printf("\n\n");
	}
      else
	{
	  printf("The infinity norm for all edges reconstruction error is %.15E\n",norm_i_b);
	  printf("  This occured on bedge %d: real_node = %d , boundary = %d , segment = %d\n",inf_bedge,real_node,b,seg);
	  printf("  real_node : < %f , %f >\n",(float)grid->x[real_node*NDIM+0],(float)grid->x[real_node*NDIM+1]);
	  printf("\n\n");
	}

      
      // Now we calculate the area-averaged value of the test function using Gaussian quadrature
      // from the extrapolated values. Track the largest extrapolation error and report the norms
      // of the differences.

      // Set up the barycentric coordinates and their respective weights.
      bary_coord[0][0] = 1./3.; bary_coord[0][1] = 1./3.; bary_coord[0][2] = 1./3.;
      bary_coord[1][0] = b_r;   bary_coord[1][1] = b_r;   bary_coord[1][2] = b_s;
      bary_coord[2][0] = b_r;   bary_coord[2][1] = b_s;   bary_coord[2][2] = b_r;
      bary_coord[3][0] = b_s;   bary_coord[3][1] = b_r;   bary_coord[3][2] = b_r;
      bary_coord[4][0] = b_t;   bary_coord[4][1] = b_t;   bary_coord[4][2] = b_u;
      bary_coord[5][0] = b_t;   bary_coord[5][1] = b_u;   bary_coord[5][2] = b_t;
      bary_coord[6][0] = b_u;   bary_coord[6][1] = b_t;   bary_coord[6][2] = b_t;

      // Lets also map these guys onto the reference triangle for integration on the triangles with curved sides.
      for ( j=0; j < 7; j++ )
	{
	  ct_gp[j][0] = (bary_coord[j][0] * 0. +
			 bary_coord[j][1] * 1. +
			 bary_coord[j][2] * 0. );

	  ct_gp[j][1] = (bary_coord[j][0] * 0. +
			 bary_coord[j][1] * 0. +
			 bary_coord[j][2] * 1. );
			
	}  // The weights can be reused since everything is consistent in the ordering.
      
      wg[0] = 9./40.;
      wg[1] = wA;
      wg[2] = wA;
      wg[3] = wA;
      wg[4] = wB;
      wg[5] = wB;
      wg[6] = wB;
      double tri_contrib_1, tri_contrib_2;
      for ( j=1; j <= nn; j++ )
	ext_cv_avg[j] = 0.;

      // Start the process by looping over all the elements and looping through the internal median dual pieces.
      for ( e=Tri; e <= Quad; e++ )                 // Element types loop.
	{
	  // Find the number of nodes and edges for this element type.
	  num_vert = NumberNodesForElement(e);

	  for ( j=1; j <= grid->num_elem[e]; j++ )              // Element loop.
	    {
	      // Retrieve the element node ids. This is for convenience.
	      for ( k=0; k < num_vert; k++ )              // Element node extraction loop.
		{
		  nodes[k] = grid->c2n[e][j*num_vert+k];
		}                                         // End element node extraction loop.
	      
	      // Find the element centroid.
	      xc[0] = 0.; xc[1] = 0.;
	      for ( k=0; k < num_vert; k++ )
		{
		  for ( i=0; i < NDIM; i++ )
		    {
		      xc[i] += grid->x[NDIM*nodes[k]+i];
		    }
		}
	      xc[0] /= (num_vert*1.);
	      xc[1] /= (num_vert*1.);

	      // Loop through each node in the element.
	      for ( k=0; k < num_vert; k++ )              // Element node loop.
		{
		  node = nodes[k];
		  
		  // Get the two neighboring vertices.
		  left_id = (k+1)%num_vert;
		  right_id = (k+num_vert-1)%num_vert;
		  
		  left_id = grid->c2n[e][j*num_vert+left_id];
		  right_id = grid->c2n[e][j*num_vert+right_id];
		  
		  // Get the edge midpoints.
		  for ( i=0; i < NDIM; i++ )
		    {
		      xmidL[i] = grid->x[NDIM*node+i] + grid->x[NDIM*left_id+i];
		      xmidR[i] = grid->x[NDIM*node+i] + grid->x[NDIM*right_id+i];
		    }

		  for ( i=0; i < NDIM; i++ )
		    {
		      xmidL[i] *= 0.5;
		      xmidR[i] *= 0.5;
		    }

		  // Make the vectors.
		  for ( i=0; i < NDIM; i++ )
		    {
		      v_left[i]   = xmidL[i] - grid->x[NDIM*node+i];
		      v_right[i]  = xmidR[i] - grid->x[NDIM*node+i];
		      v_middle[i] = xc[i]    - grid->x[NDIM*node+i];
		    }

		  // Process the first triangle.

		  ct_flag = 0;

		  // If the node and its neighbor are both boundary nodes, we need to gather some information.
		  if ( grid->node_state[node] == BOUNDARY && grid->node_state[left_id] == BOUNDARY )
		    {
		      // Retrieve the relevant information.
		      
		      // First find the correct bedge. It will have both nodes.
		      for ( i=1; i <= grid->nbedges; i++ )
			{
			  if ( node == grid->bedges[i*5] )  // candidate bedge.
			    {
			      b = grid->bedges[i*5+3];
			      seg = grid->bedges[i*5+4];

			      if ( (  node == grid->bs[b][seg][0]   ||  node == grid->bs[b][seg][1]   ) &&
				   ( left_id == grid->bs[b][seg][0] || left_id == grid->bs[b][seg][1] )   )  // we have a winner.
				{
				  real_node = grid->bedges[i*5+0];
				  ghost_node = grid->bedges[i*5+1] + grid->nn;
				  b = grid->bedges[i*5+3];
				  bc = grid->bbc[b];
				  seg = grid->bedges[i*5+4];

				  if ( bc > 10 )
				    ct_flag = 1;
			      
				  break;
				}
			    }
			}
		    }

		  // Is this node a boundary node with a boundary neighbor AND is it a curved boundary -> special treatment.
		  //if ( bc > 10 && ( grid->node_state[node] == BOUNDARY && grid->node_state[left_id] == BOUNDARY ) )
		  if ( ct_flag )
		    {
		      // Now we need to set up the six nodes used by the shape functions.
		      
		      // Node 1 is the element centroid.
		      ct_nodes[1][0] = xc[0];
		      ct_nodes[1][1] = xc[1];

		      // Node 2 is the node we are working on.
		      ct_nodes[2][0] = grid->x[NDIM*node+0];
		      ct_nodes[2][1] = grid->x[NDIM*node+1];

		      // Third node is the ghost node.
		      ct_nodes[3][0] = grid->x[NDIM*ghost_node+0];
		      ct_nodes[3][1] = grid->x[NDIM*ghost_node+1];

		      // Node 4 is the midpoint of 1,2.
		      ct_nodes[4][0] = 0.5*(ct_nodes[1][0] + ct_nodes[2][0]);
		      ct_nodes[4][1] = 0.5*(ct_nodes[1][1] + ct_nodes[2][1]);

		      // Node 5 is the midpoint of the curved boundary between node and ghost node.
		      curved_boundary_midpoint(grid, bc, &(ct_nodes[2][0]), &(ct_nodes[3][0]), &(ct_nodes[5][0]) );

		      // Node 6 is the midpoint of 1,3.
		      ct_nodes[6][0] = 0.5*(ct_nodes[1][0] + ct_nodes[3][0]);
		      ct_nodes[6][1] = 0.5*(ct_nodes[1][1] + ct_nodes[3][1]);

		      // First we need to make sure that the concavity constraint is met. Under a linear transformation
		      // to the reference element, node 5 must be greater than 0.25 in both r and s coordinates.
		      // See Flaherty's notes online or Claes Johnson's book.

		      XI = ( (ct_nodes[1][0]*ct_nodes[3][1] - ct_nodes[1][1]*ct_nodes[3][0]) -
			     (ct_nodes[5][0]*ct_nodes[3][1] - ct_nodes[5][1]*ct_nodes[3][0]) +
			     (ct_nodes[5][0]*ct_nodes[1][1] - ct_nodes[5][1]*ct_nodes[1][0])  ) /
			   ( (ct_nodes[1][0]*ct_nodes[3][1] - ct_nodes[1][1]*ct_nodes[3][0]) -
			     (ct_nodes[2][0]*ct_nodes[3][1] - ct_nodes[2][1]*ct_nodes[3][0]) +
			     (ct_nodes[2][0]*ct_nodes[1][1] - ct_nodes[2][1]*ct_nodes[1][0])  );

		      ETA = ( (ct_nodes[1][0]*ct_nodes[2][1] - ct_nodes[1][1]*ct_nodes[2][0]) -
			      (ct_nodes[5][0]*ct_nodes[2][1] - ct_nodes[5][1]*ct_nodes[2][0]) +
			      (ct_nodes[5][0]*ct_nodes[1][1] - ct_nodes[5][1]*ct_nodes[1][0])  ) /
			    ( (ct_nodes[1][0]*ct_nodes[2][1] - ct_nodes[1][1]*ct_nodes[2][0]) -
			      (ct_nodes[3][0]*ct_nodes[2][1] - ct_nodes[3][1]*ct_nodes[2][0]) +
			      (ct_nodes[3][0]*ct_nodes[1][1] - ct_nodes[1][0]*ct_nodes[3][1])  );

		      if ( ETA < 0.25 || XI < 0.25 )  // The transformation will not be invertible.
			{
			  printf("CRITICAL ERROR: Triangle 1 of element %d of type %d is too concave (hessian.C)!\n",j,e);
			  printf("  xi = %.15e     eta = %.15e\n",XI,ETA);
			  fflush(stdout);
			  exit(1);
			}

		      tri_int = 0.;

		      // Now I can loop over the Gaussian integration nodes and apply the quadratic transformation.
		      for ( i=0; i < 7; i++ )
			{
			  XI = ct_gp[i][0];
			  ETA = ct_gp[i][1];
			  
			  // Compute basis functions.
			  N1 = 1. - 3.*XI - 3.*ETA + 4.*XI*ETA + 2.*XI*XI + 2.*ETA*ETA;
			  N2 = XI*(2.*XI - 1.);
			  N3 = ETA*(2.*ETA - 1.);
			  N4 = 4.*XI*(1. - XI - ETA);
			  N5 = 4.*XI*ETA;
			  N6 = 4.*ETA*(1. - XI - ETA);

			  // Compute their derivatives.
			  N1r = -3. + 4.*ETA + 4.*XI;
			  N1s = -3. + 4.*XI + 4.*ETA;
			  N2r = 4.*XI - 1.;
			  N2s = 0.;
			  N3r = 0.;
			  N3s = 4.*ETA - 1.;
			  N4r = 4.*(1. - 2.*XI - ETA);
			  N4s = -4.*XI;
			  N5r = 4.*ETA;
			  N5s = 4.*XI;
			  N6r = -4.*ETA;
			  N6s = 4.*(1. - 2.*ETA - XI);

			  // Compute the derivatives of x,y wrt to the transformation.
			  xr = N1r*ct_nodes[1][0] + N2r*ct_nodes[2][0] + N3r*ct_nodes[3][0] + N4r*ct_nodes[4][0]
			     + N5r*ct_nodes[5][0] + N6r*ct_nodes[6][0];
			  xs = N1s*ct_nodes[1][0] + N2s*ct_nodes[2][0] + N3s*ct_nodes[3][0] + N4s*ct_nodes[4][0]
			     + N5s*ct_nodes[5][0] + N6s*ct_nodes[6][0];
			  yr = N1r*ct_nodes[1][1] + N2r*ct_nodes[2][1] + N3r*ct_nodes[3][1] + N4r*ct_nodes[4][1]
			     + N5r*ct_nodes[5][1] + N6r*ct_nodes[6][1];
			  ys = N1s*ct_nodes[1][1] + N2s*ct_nodes[2][1] + N3s*ct_nodes[3][1] + N4s*ct_nodes[4][1]
			     + N5s*ct_nodes[5][1] + N6s*ct_nodes[6][1];

			  jacobian = xr*ys - xs*yr;

			  // Now I need to extrapolate to this point. First we get the real coordinates of the
			  // Gaussian point.
			  gp[0][0] = N1*ct_nodes[1][0] + N2*ct_nodes[2][0] + N3*ct_nodes[3][0] +
                                     N4*ct_nodes[4][0] + N5*ct_nodes[5][0] + N6*ct_nodes[6][0];

			  gp[0][1] = N1*ct_nodes[1][1] + N2*ct_nodes[2][1] + N3*ct_nodes[3][1] +
                                     N4*ct_nodes[4][1] + N5*ct_nodes[5][1] + N6*ct_nodes[6][1];

			  dx[0] = gp[0][0] - x[node*NDIM+0];
			  dx[1] = gp[0][1] - x[node*NDIM+1];

			  // Get the extapolated value at the Gauss point.
			  Reconstruct_Gauss_Boundary ( node, gp[0], QL, grid, params);  // one-side reconstruction. It really doesn't care if its on a boundary.

			  ext_val[0] = QL[0];

			  temp = Test_Function(gp[0]);
			  
			  temp = fabs( ext_val[0] - temp );
			  
			  // Compare to the max.
			  max_diff = MAX(max_diff,temp);

			  // Now we accumulate to the triangle.
			  tri_int += ( wg[i] * ext_val[0] * jacobian * 0.5 );
			}
		      tri_contrib_1 = tri_int;
		      // Add to the running total for the node.
		      ext_cv_avg[node] += tri_int;
      
		    }
		  else
		    {
		      // Get the area.
		      dA = 0.5*( v_left[0]*v_middle[1] - v_left[1]*v_middle[0] );

		      // Initialize the Gauss points.
		      for ( v=0; v < 7; v++ )
			{
			  for ( ind=0; ind < NDIM; ind++ )
			    {
			      gp[v][ind] = (bary_coord[v][0] * x[node*NDIM+ind] +
					    bary_coord[v][1] * xmidL[ind] +
					    bary_coord[v][2] * xc[ind]);
			    }
			  gp[v][2] = 0.;   // Always 0 since 2D vectors.
			}

		      // Now we need to extrapolate to the Gauss points using the Reconstruction
		      // function.
		      for ( i=0; i < 7; i++ )
			{
			  dx[0] = gp[i][0] - x[node*NDIM+0];
			  dx[1] = gp[i][1] - x[node*NDIM+1];
			  
			  // Get the extapolated value at the Gauss point.
			  Reconstruct_Gauss_Boundary ( node, gp[i], QL, grid, params);  // one-side reconstruction. It really doesn't care if its on a boundary.

			  ext_val[i] = QL[0];
			  
			  temp = Test_Function(gp[i]);
			  
			  temp = fabs( ext_val[i] - temp );
			  
			  // Compare to the max.
			  max_diff = MAX(max_diff,temp);
			}
		      
		      // Now we can compute the area integral of the triangle.
		      temp = 0.;
		      for ( i=0; i < 7; i++ )
			{
			  temp += ( wg[i]*ext_val[i] );
			}
		      
		      temp *= dA;
		      
		      // Accumulate the integral to the node.
		      ext_cv_avg[node] += temp;
		    }
		  
		  
		  // Repeat the process for the second triangle.
		  
		  ct_flag = 0;
		  
		  // If the node and its neighbor are both boundary nodes, we need to gather some information.
		  if ( grid->node_state[node] == BOUNDARY && grid->node_state[right_id] == BOUNDARY )
		    {
		      // Retrieve the relevant information.
		      
		      // First find the correct bedge. It will have both nodes.
		      for ( i=1; i <= grid->nbedges; i++ )
			{
			  if ( node == grid->bedges[i*5] )  // candidate bedge.
			    {
			      b = grid->bedges[i*5+3];
			      seg = grid->bedges[i*5+4];
			      
			      if ( (   node == grid->bs[b][seg][0]   ||   node == grid->bs[b][seg][1]   ) &&
				   ( right_id == grid->bs[b][seg][0] || right_id == grid->bs[b][seg][1] )   )  // we have a winner.
				{
				  real_node = grid->bedges[i*5+0];
				  ghost_node = grid->bedges[i*5+1] + grid->nn;
				  b = grid->bedges[i*5+3];
				  bc = grid->bbc[b];
				  seg = grid->bedges[i*5+4];

				  if ( bc > 10 )
				    ct_flag = 1;
			      
				  break;
				}
			    }
			}
		    }
		  

		  // Is this node a boundary node with a boundary neighbor AND is it a curved boundary -> special treatment.
		  //if ( bc > 10 && ( grid->node_state[node] == BOUNDARY && grid->node_state[right_id] == BOUNDARY ) )
		  if ( ct_flag )
		    {
		      // Now we need to set up the six nodes used by the shape functions.
		      
		      // Node 1 is the element centroid.
		      ct_nodes[1][0] = xc[0];
		      ct_nodes[1][1] = xc[1];

		      // Node 2 is the ghost node.
		      ct_nodes[2][0] = grid->x[NDIM*ghost_node+0];
		      ct_nodes[2][1] = grid->x[NDIM*ghost_node+1];

		      // Third node is the node we are working on.
		      ct_nodes[3][0] = grid->x[NDIM*node+0];
		      ct_nodes[3][1] = grid->x[NDIM*node+1];

		      // Node 4 is the midpoint of 1,2.
		      ct_nodes[4][0] = 0.5*(ct_nodes[1][0] + ct_nodes[2][0]);
		      ct_nodes[4][1] = 0.5*(ct_nodes[1][1] + ct_nodes[2][1]);

		      // Node 5 is the midpoint of the curved boundary between node and ghost node.
		      curved_boundary_midpoint(grid, bc, &(ct_nodes[2][0]), &(ct_nodes[3][0]), &(ct_nodes[5][0]) );

		      // Node 6 is the midpoint of 1,3.
		      ct_nodes[6][0] = 0.5*(ct_nodes[1][0] + ct_nodes[3][0]);
		      ct_nodes[6][1] = 0.5*(ct_nodes[1][1] + ct_nodes[3][1]);

		      // First we need to make sure that the concavity constraint is met. Under a linear transformation
		      // to the reference element, node 5 must be greater than 0.25 in both r and s coordinates.
		      // See Flaherty's notes online or Claes Johnson's book.

		      XI = ( (ct_nodes[1][0]*ct_nodes[3][1] - ct_nodes[1][1]*ct_nodes[3][0]) -
			     (ct_nodes[5][0]*ct_nodes[3][1] - ct_nodes[5][1]*ct_nodes[3][0]) +
			     (ct_nodes[5][0]*ct_nodes[1][1] - ct_nodes[5][1]*ct_nodes[1][0])  ) /
			   ( (ct_nodes[1][0]*ct_nodes[3][1] - ct_nodes[1][1]*ct_nodes[3][0]) -
			     (ct_nodes[2][0]*ct_nodes[3][1] - ct_nodes[2][1]*ct_nodes[3][0]) +
			     (ct_nodes[2][0]*ct_nodes[1][1] - ct_nodes[2][1]*ct_nodes[1][0])  );

		      ETA = ( (ct_nodes[1][0]*ct_nodes[2][1] - ct_nodes[1][1]*ct_nodes[2][0]) -
			      (ct_nodes[5][0]*ct_nodes[2][1] - ct_nodes[5][1]*ct_nodes[2][0]) +
			      (ct_nodes[5][0]*ct_nodes[1][1] - ct_nodes[5][1]*ct_nodes[1][0])  ) /
			    ( (ct_nodes[1][0]*ct_nodes[2][1] - ct_nodes[1][1]*ct_nodes[2][0]) -
			      (ct_nodes[3][0]*ct_nodes[2][1] - ct_nodes[3][1]*ct_nodes[2][0]) +
			      (ct_nodes[3][0]*ct_nodes[1][1] - ct_nodes[1][0]*ct_nodes[3][1])  );

		      if ( ETA < 0.25 || XI < 0.25 )  // The transformation will not be invertible.
			{
			  printf("CRITICAL ERROR: Triangle 2 of element %d of type %d is too concave (hessian.C)!\n",j,e);
			  printf("  xi = %.15e     eta = %.15e\n",XI,ETA);
			  fflush(stdout);
			  exit(1);
			}

		      tri_int = 0.;

		      // Now I can loop over the Gaussian integration nodes and apply the quadratic transformation.
		      for ( i=0; i < 7; i++ )
			{
			  XI = ct_gp[i][0];
			  ETA = ct_gp[i][1];
			  
			  // Compute basis functions.
			  N1 = 1. - 3.*XI - 3.*ETA + 4.*XI*ETA + 2.*XI*XI + 2.*ETA*ETA;
			  N2 = XI*(2.*XI - 1.);
			  N3 = ETA*(2.*ETA - 1.);
			  N4 = 4.*XI*(1. - XI - ETA);
			  N5 = 4.*XI*ETA;
			  N6 = 4.*ETA*(1. - XI - ETA);

			  // Compute their derivatives.
			  N1r = -3. + 4.*ETA + 4.*XI;
			  N1s = -3. + 4.*XI + 4.*ETA;
			  N2r = 4.*XI - 1.;
			  N2s = 0.;
			  N3r = 0.;
			  N3s = 4.*ETA - 1.;
			  N4r = 4.*(1. - 2.*XI - ETA);
			  N4s = -4.*XI;
			  N5r = 4.*ETA;
			  N5s = 4.*XI;
			  N6r = -4.*ETA;
			  N6s = 4.*(1. - 2.*ETA - XI);

			  // Compute the derivatives of x,y wrt to the transformation.
			  xr = N1r*ct_nodes[1][0] + N2r*ct_nodes[2][0] + N3r*ct_nodes[3][0] + N4r*ct_nodes[4][0]
			     + N5r*ct_nodes[5][0] + N6r*ct_nodes[6][0];
			  xs = N1s*ct_nodes[1][0] + N2s*ct_nodes[2][0] + N3s*ct_nodes[3][0] + N4s*ct_nodes[4][0]
			     + N5s*ct_nodes[5][0] + N6s*ct_nodes[6][0];
			  yr = N1r*ct_nodes[1][1] + N2r*ct_nodes[2][1] + N3r*ct_nodes[3][1] + N4r*ct_nodes[4][1]
			     + N5r*ct_nodes[5][1] + N6r*ct_nodes[6][1];
			  ys = N1s*ct_nodes[1][1] + N2s*ct_nodes[2][1] + N3s*ct_nodes[3][1] + N4s*ct_nodes[4][1]
			     + N5s*ct_nodes[5][1] + N6s*ct_nodes[6][1];

			  jacobian = xr*ys - xs*yr;

			  // Now I need to extrapolate to this point. First we get the real coordinates of the
			  // Gaussian point.
			  gp[0][0] = N1*ct_nodes[1][0] + N2*ct_nodes[2][0] + N3*ct_nodes[3][0] +
                                     N4*ct_nodes[4][0] + N5*ct_nodes[5][0] + N6*ct_nodes[6][0];

			  gp[0][1] = N1*ct_nodes[1][1] + N2*ct_nodes[2][1] + N3*ct_nodes[3][1] +
                                     N4*ct_nodes[4][1] + N5*ct_nodes[5][1] + N6*ct_nodes[6][1];

			  dx[0] = gp[0][0] - x[node*NDIM+0];
			  dx[1] = gp[0][1] - x[node*NDIM+1];

			  // Get the extapolated value at the Gauss point.
			  Reconstruct_Gauss_Boundary ( node, gp[0], QL, grid, params);  // one-side reconstruction. It really doesn't care if its on a boundary.

			  ext_val[0] = QL[0];

			  temp = Test_Function(gp[0]);
			  
			  temp = fabs( ext_val[0] - temp );
			  
			  // Compare to the max.
			  max_diff = MAX(max_diff,temp);

			  // Now we accumulate to the triangle.
			  tri_int += ( wg[i] * ext_val[0] * jacobian * 0.5 );
			}
		      tri_contrib_2 = tri_int;
		      // Add to the running total for the node.
		      ext_cv_avg[node] += tri_int;
		    }
		  else
		    {
		      // Get the area.
		      dA = 0.5*( v_middle[0]*v_right[1] - v_middle[1]*v_right[0] );
		      
		      // Initialize the Gauss points.
		      for ( v=0; v < 7; v++ )
			{
			  for ( ind=0; ind < NDIM; ind++ )
			    {
			      gp[v][ind] = (bary_coord[v][0] * x[node*NDIM+ind] +
					    bary_coord[v][1] * xmidR[ind] +
					    bary_coord[v][2] * xc[ind]);
			    }
			  gp[v][2] = 0.;   // Always 0 since 2D vectors.
			}
		      
		      // Now we need to extrapolate to the Gauss points using the Reconstruction
		      // function.
		      for ( i=0; i < 7; i++ )
			{
			  dx[0] = gp[i][0] - x[node*NDIM+0];
			  dx[1] = gp[i][1] - x[node*NDIM+1];
			  
			  // Get the extapolated value at the Gauss point.
			  Reconstruct_Gauss_Boundary ( node, gp[i], QL, grid, params);  // one-side reconstruction. It really doesn't care if its on a boundary.

			  ext_val[i] = QL[0];

			  temp = Test_Function(gp[i]);
			  
			  temp = fabs( ext_val[i] - temp );
			  
			  // Compare to the max.
			  max_diff = MAX(max_diff,temp);
			}
		      
		      // Now we can the area integral of the triangle.
		      temp = 0.;
		      for ( i=0; i < 7; i++ )
			{
			  temp += ( wg[i]*ext_val[i] );
			}
		      
		      temp *= dA;
		      
		      // Accumulate the integral to the node.
		      ext_cv_avg[node] += temp;
		    }
		}                                       // End element node loop.
	      
	    }                                           // End element loop.

	}                                               // End element type loop.

      // Now we need to get the total integral values not just area averaged.
      for ( i=1; i <= nn; i++ )
	{
	  cv_avg[i] *= grid->cv_area[i];
	}

      //for ( i=1; i <= nn; i++ )
	//printf("ext_cv_avg[%d] = %.15E\n",i,ext_cv_avg[i]);

      // Report the largest difference.
      printf("EXTRAPOLATION TEST: The largest error at the Gauss points is: %.15e\n",max_diff);

      // Compute the norms of the difference.
      norm = 0.;
      for ( i=1; i <= nn; i++ )
	{
	  norm += fabs( cv_avg[i] - ext_cv_avg[i] );
	}
      printf("EXTRAPOLATION TEST: L1 norm of difference: %.15e\n",norm);

      norm = 0.;
      for ( i=1; i <= nn; i++ )
	{
	  temp = cv_avg[i] - ext_cv_avg[i];
	  norm += (temp*temp);
	}
      printf("EXTRAPOLATION TEST: L2 norm of difference: %.15e\n",sqrt(norm));

      norm = 0.;
      for ( i=1; i <= nn; i++ )
	{
	  temp = fabs( cv_avg[i] - ext_cv_avg[i] );
	  if ( temp > norm )
	    {
	      norm = temp;
	      j = i;
	    }
	}
      printf("EXTRAPOLATION TEST: L_inf norm of difference: %.15e   at node %d\n",norm,j);
      
      // Compute the total integral of the function over the domain.
      total_integral = 0.;
      for ( i=1; i <= nn; i++ )
	{
	  total_integral += ( cv_avg[i] );
	}

      printf("Total integral of the test function over the domain using divergence and Gaussian Quadrature is %.15E\n",
	     total_integral);
      
      temp = total_integral;
      total_integral = 0.;
      for ( i=1; i <= nn; i++ )
	{
	  total_integral += ( ext_cv_avg[i] );
	}
      
      printf("Total integral of the test function over the domain using extrapolation and Gaussian Quadrature is %.15E\n",
	     total_integral);

      printf("The difference in the integral values is %.15E\n",fabs(total_integral - temp));

      for ( i=1; i <= nn; i++ )
	{
	  grid->Q[i*NUM_VAR] = fabs(cv_avg[i] - ext_cv_avg[i]);
	}

      // Now write a solution file with this.
      sprintf(buff,"tecplot_diff_cv_ext.dat");
      write_tecplot_solutionC ( buff, grid );
      
    }

  if ( DEBUG_DERIVS )
    {
      for ( i=0; i < NUM_MOM; i++ )
	freenull(derivs[i]);

      freenull(derivs);

      freenull(qprime);
    }

  freenull(U);
  freenull(sol);
  freenull(temp1);
  freenull(temp2);
  freenull(temp3);
  freenull(cv_avg);
  freenull(ext_cv_avg);
  
  if ( factorize )
    {
      freenull(A);
      freenull(aa);
      freenull(s);
      freenull(u);
      freenull(vt);
      freenull(work);
    }

  //scanf("%c",&jobu);

  return;
}
