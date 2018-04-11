//=============================================================
// 
//  hessian.C
//  
//  Functions to calculate the Hessian of the conserved variables.
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


//=============================================================
// 
//  Compute_Hessian()
//
//  Calculates the Hessian of the dependent variables using a
//  least squares technique for k-exact reconstruction.
//  
//  GRID *grid                        // The grid.
//
//=============================================================

void Compute_Hessian ( GRID *grid, PARAMS params )
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
  double *x = grid->x;
  double *q = grid->nQ;
  double *hess = grid->hess;
  double *moments = grid->Moments;
  double *point_value = grid->point_value;
  double *cv_avg = NULL;

  if (DEBUG_HESS)
    {
      derivs = (char**)malloc(NUM_MOM*sizeof(char*));
      if ( derivs == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'derivs'.\n"); exit(1); }

      for ( i=0; i < 5; i++ )
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

      //for ( i=1; i <= grid->nn; i++ )
      //{
      //grid->Q[i*NUM_VAR+0] = Test_Function(&(grid->x[i*NDIM]));
      //grid->Q[i*NUM_VAR+1] = 4.*grid->x[i*NDIM] + 5.*grid->x[i*NDIM+1];
      //grid->Q[i*NUM_VAR+2] = 7.*grid->x[i*NDIM] + 8.*grid->x[i*NDIM+1];
      //grid->Q[i*NUM_VAR+3] = 10.*grid->x[i*NDIM] + 11.*grid->x[i*NDIM+1];
      //}
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
  sol = (double*)malloc(6*sizeof(double));
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

  if (!SVDFAC && 0) // Allocate memory for the QR factorization if we are using it.
    {
      if ( grid->Qfac == NULL )
	{
	  grid->Qfac = (double**)malloc((nn+1)*sizeof(double*));
	  if ( grid->Qfac == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'Qfac'.\n"); exit(1); }
	}
      if ( grid->Rfac == NULL )
	{
	  grid->Rfac = (double**)malloc((nn+1)*sizeof(double*));
	  if ( grid->Rfac == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'Rfac'.\n"); exit(1); }
	}

      for ( i=0; i <= grid->nn; i++ )
	{
	  grid->Qfac[i] = 0;
	  grid->Rfac[i] = 0;
	}
    }

  // Reset pointers.
  hess = grid->hess;
  point_value = grid->point_value;

  // Clean out all memory that will be reset in this function.
  for ( i=0; i < 6; i++ )
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
  if (DEBUG_HESS)
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

  // Find maximum neighbor degree.
  MAX_NEIGHBORS = 0;
  for ( i=1; i <= nn; i++ )
    {
      num_neighbors = nnsn2[i+1]-nnsn2[i];
      MAX_NEIGHBORS = MAX(MAX_NEIGHBORS,num_neighbors);
    }

  MAX_NEIGHBORS++;  // Account for the node itself.

  if (DEBUG)
    printf("  Maximum number of 2nd degree neighbors in the mesh is %d.\n",MAX_NEIGHBORS-1);
  
  // We can now begin to solve the system. I'm including the full system as well as the system after
  // the first row has been reduced using Gaussian elimination. I also want to leave it general enough
  // so I'll include the a geometric weight term in the calculations.

  // The system will currently be solved every time-step using the SVD (called from lapack).
  

  // First we want to allocate memory used in the SVD routine.
  n = 6;
  i = MIN(MAX_NEIGHBORS,n) * 3 + MAX(MAX_NEIGHBORS,n);
  j = MIN(MAX_NEIGHBORS,n) * 5;
  MAX_LWORK = MAX(i,j) * 2;

  Adim = n*MAX_NEIGHBORS;
  Udim = MAX_NEIGHBORS;
  tempdim = MAX_NEIGHBORS;
  aadim = MAX_NEIGHBORS*n;
  udim = MAX_NEIGHBORS*MAX_NEIGHBORS;
  vtdim = n*n;
  sdim = MAX_NEIGHBORS;
  workdim = MAX_LWORK;
  
  A = (double*)malloc(Adim*sizeof(double));
  if ( A == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'A'.\n"); exit(1); }

  U = (double*)malloc(Udim*sizeof(double));
  if ( U == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'U'.\n"); exit(1); }

  temp1 = (double*)malloc(tempdim*sizeof(double));
  if ( temp1 == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'temp1'.\n"); exit(1); }

  temp2 = (double*)malloc(tempdim*sizeof(double));
  if ( temp2 == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'temp2'.\n"); exit(1); }

  temp3 = (double*)malloc(tempdim*sizeof(double));
  if ( temp3 == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'temp3'.\n"); exit(1); }

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

  for ( i=0; i < Udim; i++ )
    U[i] = 0.;

  for ( i=0; i < tempdim; i++ )
    {
      temp1[i] = 0.;  temp2[i] = 0.;  temp3[i] = 0.;
    }

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

  p = (double)GEOM_WEIGHT;

  /*
  double Ad[24];
  double Qd[36];
  double Bd[6];
  double Xd[4];

  Ad[0] =  1.0;
  Ad[1] =  3.0;
  Ad[2] =  2.0;
  Ad[3] =  4.0;
  
  Ad[4] =  6.0;
  Ad[5] =  5.0;
  Ad[6] =  7.0;
  Ad[7] =  8.0;
      
  Ad[8] =   9.0;
  Ad[9] =   8.0;
  Ad[10] =  7.0;
  Ad[11] =  6.0;
  
  Ad[12] =  4.0;
  Ad[13] =  6.0;
  Ad[14] =  7.0;
  Ad[15] =  8.0;

  Ad[16] =  1.0;
  Ad[17] =  0.0;
  Ad[18] =  0.0;
  Ad[19] =  0.0;
      
  Ad[20] =  0.0;
  Ad[21] =  0.0;
  Ad[22] =  0.0;
  Ad[23] =  1.0;

  Bd[0] = 5.;
  Bd[1] = 5.;
  Bd[2] = 4.;
  Bd[3] = 2.;
  Bd[4] = 5.;
  Bd[5] = 1.;

  QR_House(Ad,Qd,6,4);
  QR_QtB ( Qd, Bd, Xd, 6 );
  QR_Backsub ( Ad, Xd, 4 );
  exit(1);
  */

  if ( DEBUG_HESS )
    {
      FILE *fp_cva = NULL;
      fp_cva = fopen("cv_averages.dat","w");

      for ( i=1; i <= grid->nn; i++ )
	{
	  fprintf(fp_cva,"%.15e\n",cv_avg[i]);
	}
      fclose(fp_cva);
    }
  else if ( 0 )
    {
      FILE *fp_q1 = NULL;
      FILE *fp_q2 = NULL;
      FILE *fp_q3 = NULL;
      FILE *fp_q4 = NULL;
      
      fp_q1 = fopen("hess_qv1.dat","w");
      fp_q2 = fopen("hess_qv2.dat","w");
      fp_q3 = fopen("hess_qv3.dat","w");
      fp_q4 = fopen("hess_qv4.dat","w");
      
      for ( i=1; i <= grid->nn; i++ )
	{
	  fprintf(fp_q1,"%.15e\n",q[i*NUM_VAR+0]);
	  fprintf(fp_q2,"%.15e\n",q[i*NUM_VAR+1]);
	  fprintf(fp_q3,"%.15e\n",q[i*NUM_VAR+2]);
	  fprintf(fp_q4,"%.15e\n",q[i*NUM_VAR+3]);
	}
      fclose(fp_q1);
      fclose(fp_q2);
      fclose(fp_q3);
      fclose(fp_q4);
    }

  if ( SVDFAC )
    {
      if ( HESS_FULL_SYSTEM )
	{
	  for ( i=1; i <= nn; i++ )
	    {
	      if ( grid->node_order[i] < 3 )
		continue;

	      num_neighbors = nnsn2[i+1] - nnsn2[i];
	      n = 6;

	      num_neighbors++;  // Include the node itself.

	      // Put in the mean constraint.
	      A[0 + 0*num_neighbors] = 1.0;
	      A[0 + 1*num_neighbors] = moments[i*NUM_MOM+0]; // x
	      A[0 + 2*num_neighbors] = moments[i*NUM_MOM+1]; // y
	      A[0 + 3*num_neighbors] = 0.5*moments[i*NUM_MOM+2]; // xx
	      A[0 + 4*num_neighbors] = 0.5*moments[i*NUM_MOM+3]; // yy
	      A[0 + 5*num_neighbors] = moments[i*NUM_MOM+4]; // xy

	      for ( j=nnsn2[i]; j < nnsn2[i+1]; j++ )
		{
		  node = nsn2[j];

		  // Get the geometric weight term for current node.
		  dx[0] = x[node*NDIM+0] - x[i*NDIM+0];
		  dx[1] = x[node*NDIM+1] - x[i*NDIM+1];

		  dist = sqrt(dx[0]*dx[0] + dx[1]*dx[1]);

		  w = 1.0 / pow(dist,p);
		  
		  row = j - nnsn2[i] + 1;
	      
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

	      // Finished contruction of LHS.

	      if ( i == 0 )
		{
		  for ( j=0; j < num_neighbors; j++ )
		    {
		      printf("%.15e %.15e %.15e %.15e %.15e %.15e\n", A[j + 0*num_neighbors], A[j + 1*num_neighbors], A[j + 2*num_neighbors],
			     A[j + 3*num_neighbors], A[j + 4*num_neighbors], A[j + 5*num_neighbors] );
		    }
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

	      // Now presumably we have the SVD of the matrix which is constant for each variable.
	      // So, we may now loop over the number of variables to get the gradiant/hessian of
	      // each of them.
	      for( var=0; var < NUM_VAR; var++)
		{
		  U[0] = q[i*NUM_VAR+var];

                  if (DEBUG_HESS)
		    U[0] = cv_avg[i];

		  for ( j=nnsn2[i]; j < nnsn2[i+1]; j++ )
		    {
		      node = nsn2[j];

		      // Get the geometric weight term for current node.
		      dx[0] = x[node*NDIM+0] - x[i*NDIM+0];
		      dx[1] = x[node*NDIM+1] - x[i*NDIM+1];

		      dist = sqrt(dx[0]*dx[0] + dx[1]*dx[1]);
		      
		      w = 1.0 / pow(dist,p);
		      
		      row = j - nnsn2[i] + 1;
		      
		      U[row] = w*q[node*NUM_VAR+var];

		      if (DEBUG_HESS)
			U[row] = w*cv_avg[node];
		    }

		  if ( i == 0 && var == 0 )
		    {
		      for ( j=0; j < num_neighbors; j++ )
			{
			  printf("%.15e\n", U[j]);
			}
		    }

		  // Perform the linear algebra for the system Ax=b -> (SVD)x=b
		  for ( j=0; j < m; j++ )
		    {
		      for( k=0; k < m; k++ )
			{
			  temp1[k]=u[j*m+k];
			}
		      temp2[j]=dotProductGeneral(U,temp1,m);
		    }
	      
		  for ( j=0; j < tempdim; j++ )
		    temp1[j] = 0.;
	      
		  for ( j=0; j < lds; j++ )
		    {
		      temp1[j]=temp2[j]/s[j];
		    }
	      
		  for ( j=0; j < tempdim; j++ )
		    temp2[j] = 0.;
		  for ( j=0; j < tempdim; j++ )
		    temp3[j] = 0.;

	      
		  for ( j=0; j < n; j++ )
		    {
		      for( k=0; k < n; k++ )
			{
			  temp3[k]=vt[j*n+k];
			}
		      sol[j] = dotProductGeneral(temp3,temp1,n);
		    }
	      
		  hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 0] = sol[1];
		  hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 1] = sol[2];
		  hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 2] = sol[3];
		  hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 3] = sol[4];
		  hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 4] = sol[5];
	      
		  point_value[i*NUM_VAR+var] = sol[0];
		  
		  if ( DEBUG_HESS && var == 0 )
		    {
		      if ( (fabs( x[i*NDIM+0] - 0.5 ) < 1.0e-10) && (fabs( x[i*NDIM+1] ) < 1.0e-10 ) ) // the node I'll study.
			{
			  printf("Node %d is at 0.5,0.0 :\n",i);
			  
			  for ( j=0; j < lds; j++ )
			    {
			      printf("sigma[%d] = %.15E\n",j,s[j]);
			    }
			  
			  printf("point values = %.15E\n",sol[0]);
			  printf("fx = %.15E\n",sol[1]);
			  printf("fy = %.15E\n",sol[2]);
			  printf("fxx = %.15E\n",sol[3]);
			  printf("fyy = %.15E\n",sol[4]);
			  printf("fxy = %.15E\n",sol[5]);
			}		      
		    }		  
		}
	    }
	}
      else
	{
	  for ( i=1; i <= nn; i++ )
	    {
	      if ( grid->node_order[i] < 3 )
		continue;
	      
	      num_neighbors = nnsn2[i+1] - nnsn2[i];
	      n = 5;

	      for ( j=nnsn2[i]; j < nnsn2[i+1]; j++ )
		{
		  node = nsn2[j];
		  
		  // Get the geometric weight term for current node.
		  dx[0] = x[node*NDIM+0] - x[i*NDIM+0];
		  dx[1] = x[node*NDIM+1] - x[i*NDIM+1];
		  
		  dist = sqrt(dx[0]*dx[0] + dx[1]*dx[1]);
		  
		  w = 1.0 / pow(dist,p);

		  row = j - nnsn2[i];
		  
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

	      // Finished contruction of LHS.
	      
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

	      // Now presumably we have the SVD of the matrix which is constant for each variable.
	      // So, we may now loop over the number of variables to get the gradiant/hessian of
	      // each of them.
	      for( var=0; var < NUM_VAR; var++)
		{
		  for ( j=nnsn2[i]; j < nnsn2[i+1]; j++ )
		    {
		      node = nsn2[j];

		      // Get the geometric weight term for current node.
		      dx[0] = x[node*NDIM+0] - x[i*NDIM+0];
		      dx[1] = x[node*NDIM+1] - x[i*NDIM+1];

		      dist = sqrt(dx[0]*dx[0] + dx[1]*dx[1]);
		      
		      w = 1.0 / pow(dist,p);

		      row = j - nnsn2[i];

		      U[row] = w*(q[node*NUM_VAR+var] - q[i*NUM_VAR+var]);

                      if (DEBUG_HESS)
			U[row] = w*(cv_avg[node]-cv_avg[i]);
		    }

		  // Perform the linear algebra for the system Ax=b -> (SVD)x=b
		  for ( j=0; j < m; j++ )
		    {
		      for( k=0; k < m; k++ )
			{
			  temp1[k]=u[j*m+k];
			}
		      temp2[j]=dotProductGeneral(U,temp1,m);
		    }
		  
		  for ( j=0; j < tempdim; j++ )
		    temp1[j] = 0.;
	      
		  for ( j=0; j < lds; j++ )
		    {
		      temp1[j]=temp2[j]/s[j];
		    }
	      
		  for ( j=0; j < tempdim; j++ )
		    temp2[j] = 0.;
		  for ( j=0; j < tempdim; j++ )
		    temp3[j] = 0.;

	      
		  for ( j=0; j < n; j++ )
		    {
		      for( k=0; k < n; k++ )
			{
			  temp3[k]=vt[j*n+k];
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

		  if ( DEBUG_HESS)
		    {
		      point_value[i*NUM_VAR+var] = cv_avg[i] -
			                           moments[i*NUM_MOM+0]*sol[0] -
		                                   moments[i*NUM_MOM+1]*sol[1] -
		                                   moments[i*NUM_MOM+2]*sol[2]*0.5 -
		                                   moments[i*NUM_MOM+3]*sol[3]*0.5 -
		                                   moments[i*NUM_MOM+4]*sol[4];
		    }

		  if ( DEBUG_HESS && var == 0 )
		    {
		      if ( (fabs( x[i*NDIM+0] - 0.5 ) < 1.0e-10) && (fabs( x[i*NDIM+1] ) < 1.0e-10 ) ) // the node I'll study.
			{
			  printf("Node %d is at 0.5,0.0 :\n",i);
			  
			  for ( j=0; j < lds; j++ )
			    {
			      printf("sigma[%d] = %.15E\n",j,s[j]);
			    }

			  printf("point values = %.15E\n",point_value[i*NUM_VAR+var]);
			  printf("fx = %.15E\n",sol[0]);
			  printf("fy = %.15E\n",sol[1]);
			  printf("fxx = %.15E\n",sol[2]);
			  printf("fyy = %.15E\n",sol[3]);
			  printf("fxy = %.15E\n",sol[4]);
			}		      
		    }
		}
	    }
	}
    }
  else if ( QRFAC ) // Do QR factorization.
    {
      if ( HESS_FULL_SYSTEM )  // Get the size of the matrices (at least n).
	{
	  n = 6;
	  mtail = 1;
	}
      else
	{
	  n = 5;
	  mtail = 0;
	}
      
      // Do the memory allocation once.
      if ( grid->Qfac == NULL )
	{
	  printf("  Allocating space for the QR factorization.\n");

	  grid->Qfac = (double**)malloc((nn+1)*sizeof(double*));
	  if ( grid->Qfac == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'Qfac'.\n"); exit(1); }
	  
	  grid->Rfac = (double**)malloc((nn+1)*sizeof(double*));
	  if ( grid->Rfac == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'Rfac'.\n"); exit(1); }

	  // Loop over the nodes and figure out how big the matrices are.
	  
	  for ( i=1; i <= grid->nn; i++ )
	    {
	      m = nnsn2[i+1] - nnsn2[i] + mtail;
	      
	      // Now allocate for Q (m by m) and R (m by n).
	      
	      grid->Qfac[i] = (double*)malloc((m*m)*sizeof(double));
	      if ( grid->Qfac[i] == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'Qfac[%d]'.\n",i); exit(1); }

	      grid->Rfac[i] = (double*)malloc((m*n)*sizeof(double));
	      if ( grid->Rfac[i] == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'Rfac[%d]'.\n",i); exit(1); }
	      
	      bytes_allocated += ( m*m + m*n );
	    }

	  printf("  Done allocating memory. Total memory allocated is %d bytes.\n",bytes_allocated*((int)sizeof(double)));
	  fflush(stdout);
	  factorize = 1;
	}

      // Do the QR factorization if necessary.
      if ( factorize )
	{
	  if ( HESS_FULL_SYSTEM )
	    {
	      for ( i=1; i <= nn; i++ )
		{
		  m = nnsn2[i+1] - nnsn2[i] + 1;
		  n = 6;

		  // Put in the mean constraint.
		  grid->Rfac[i][0*n + 0] = 1.0;
		  grid->Rfac[i][0*n + 1] = moments[i*NUM_MOM+0]; // x
		  grid->Rfac[i][0*n + 2] = moments[i*NUM_MOM+1]; // y
		  grid->Rfac[i][0*n + 3] = 0.5*moments[i*NUM_MOM+2]; // xx
		  grid->Rfac[i][0*n + 4] = 0.5*moments[i*NUM_MOM+3]; // yy
		  grid->Rfac[i][0*n + 5] = moments[i*NUM_MOM+4]; // xy
		  
		  for ( j=nnsn2[i]; j < nnsn2[i+1]; j++ )
		    {
		      node = nsn2[j];

		      // Get the geometric weight term for current node.
		      dx[0] = x[node*NDIM+0] - x[i*NDIM+0];
		      dx[1] = x[node*NDIM+1] - x[i*NDIM+1];

		      dist = sqrt(dx[0]*dx[0] + dx[1]*dx[1]);
		      
		      w = 1.0 / pow(dist,p);
		  
		      row = j - nnsn2[i] + 1;
		      
		      grid->Rfac[i][row*n+0] = w*1.0;
		      
		      grid->Rfac[i][row*n+1] = w*(x[node*NDIM+0] - x[i*NDIM+0] +
						  moments[node*NUM_MOM+0]);
		      
		      grid->Rfac[i][row*n+2] = w*(x[node*NDIM+1] - x[i*NDIM+1] +
						  moments[node*NUM_MOM+1]);

		      grid->Rfac[i][row*n+3] = 0.5*w*( (x[node*NDIM+0] - x[i*NDIM+0])*
						       (x[node*NDIM+0] - x[i*NDIM+0]) +
						       2.0*moments[node*NUM_MOM+0]*(x[node*NDIM+0] - x[i*NDIM+0]) +
						       moments[node*NUM_MOM+2] );

		      grid->Rfac[i][row*n+4] = 0.5*w*( (x[node*NDIM+1] - x[i*NDIM+1])*
						       (x[node*NDIM+1] - x[i*NDIM+1]) +
						       2.0*moments[node*NUM_MOM+1]*(x[node*NDIM+1] - x[i*NDIM+1]) +
						       moments[node*NUM_MOM+3] );

		      grid->Rfac[i][row*n+5] = w*(moments[node*NUM_MOM+4] +
						  moments[node*NUM_MOM+0]*(x[node*NDIM+1] - x[i*NDIM+1]) +
						  moments[node*NUM_MOM+1]*(x[node*NDIM+0] - x[i*NDIM+0]) +
						  (x[node*NDIM+1] - x[i*NDIM+1])*(x[node*NDIM+0] - x[i*NDIM+0]));
		    }

		  // The system is built so just have to factor it.
		  QR_House ( grid->Rfac[i], grid->Qfac[i], m, n );

		}
	    }
	  else  // Mean constrained is removed by Gaussian elemination.
	    {
	      for ( i=1; i <= nn; i++ )
		{
		  m = nnsn2[i+1] - nnsn2[i];
		  n = 5;
		  
		  for ( j=nnsn2[i]; j < nnsn2[i+1]; j++ )
		    {
		      node = nsn2[j];
		      
		      // Get the geometric weight term for current node.
		      dx[0] = x[node*NDIM+0] - x[i*NDIM+0];
		      dx[1] = x[node*NDIM+1] - x[i*NDIM+1];
		      
		      dist = sqrt(dx[0]*dx[0] + dx[1]*dx[1]);
		      
		      w = 1.0 / pow(dist,p);
		      
		      row = j - nnsn2[i];
		      
		      grid->Rfac[i][row*n+0] = w*(x[node*NDIM+0] - x[i*NDIM+0] +
						  moments[node*NUM_MOM+0] -
						  moments[i*NUM_MOM+0]);
		      
		      grid->Rfac[i][row*n+1] = w*(x[node*NDIM+1] - x[i*NDIM+1] +
						  moments[node*NUM_MOM+1] -
						  moments[i*NUM_MOM+1]);
		      
		      grid->Rfac[i][row*n+2] = 0.5*w*( (x[node*NDIM+0] - x[i*NDIM+0])*
						       (x[node*NDIM+0] - x[i*NDIM+0]) +
						       2.0*moments[node*NUM_MOM+0]*(x[node*NDIM+0] - x[i*NDIM+0]) +
						       moments[node*NUM_MOM+2] -
						       moments[i*NUM_MOM+2]);
		      
		      grid->Rfac[i][row*n+3] = 0.5*w*( (x[node*NDIM+1] - x[i*NDIM+1])*
						       (x[node*NDIM+1] - x[i*NDIM+1]) +
						       2.0*moments[node*NUM_MOM+1]*(x[node*NDIM+1] - x[i*NDIM+1]) +
						       moments[node*NUM_MOM+3] -
						       moments[i*NUM_MOM+3]);
		      
		      grid->Rfac[i][row*n+4] = w*(moments[node*NUM_MOM+4] +
						  moments[node*NUM_MOM+0]*(x[node*NDIM+1] - x[i*NDIM+1]) +
						  moments[node*NUM_MOM+1]*(x[node*NDIM+0] - x[i*NDIM+0]) +
						  (x[node*NDIM+1] - x[i*NDIM+1])*(x[node*NDIM+0] - x[i*NDIM+0]) -
						  moments[i*NUM_MOM+4]);
		    }
		  
		  // The system is built so just have to factor it.
		  QR_House ( grid->Rfac[i], grid->Qfac[i], m, n );
		  
		}
	    }
	  
	}

      // Loop through the nodes, build the right hand side, and back-substitute.
      if ( HESS_FULL_SYSTEM )
	{
	  for ( i=1; i <= grid->nn; i++ )
	    {
	      if ( grid->node_order[i] < 3 )
		continue;

	      m = nnsn2[i+1] - nnsn2[i] + 1;
	      n = 6;
	      
	      for( var=0; var < NUM_VAR; var++)
		{
		  U[0] = q[i*NUM_VAR+var];
		  
                  if (DEBUG_HESS)
		    U[0] = cv_avg[i];

		  for ( j=nnsn2[i]; j < nnsn2[i+1]; j++ )
		    {
		      node = nsn2[j];
		  
		      // Get the geometric weight term for current node.
		      dx[0] = x[node*NDIM+0] - x[i*NDIM+0];
		      dx[1] = x[node*NDIM+1] - x[i*NDIM+1];
		  
		      dist = sqrt(dx[0]*dx[0] + dx[1]*dx[1]);
		      
		      w = 1.0 / pow(dist,p);
		      
		      row = j - nnsn2[i] + 1;
		  
		      U[row] = w*q[node*NUM_VAR+var];

		      if (DEBUG_HESS)
			U[row] = w*cv_avg[node];
		    }

		  // Do the multiplication Q^t*b
		  QR_QtB ( grid->Qfac[i], U, temp1, m );
		  QR_Backsub ( grid->Rfac[i], temp1, n );
		  
		  hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 0] = temp1[1];
		  hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 1] = temp1[2];
		  hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 2] = temp1[3];
		  hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 3] = temp1[4];
		  hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 4] = temp1[5];
		  
		  point_value[i*NUM_VAR+var] = temp1[0];

		  if ( DEBUG_HESS && var == 0 )
		    {
		      if ( (fabs( x[i*NDIM+0] - 0.5 ) < 1.0e-10) && (fabs( x[i*NDIM+1] ) < 1.0e-10 ) ) // the node I'll study.
			{
			  printf("Node %d is at 0.5,0.0 :\n",i);
			  printf("point values = %.15E\n",temp1[0]);
			  printf("fx = %.15E\n",temp1[1]);
			  printf("fy = %.15E\n",temp1[2]);
			  printf("fxx = %.15E\n",temp1[3]);
			  printf("fyy = %.15E\n",temp1[4]);
			  printf("fxy = %.15E\n",temp1[5]);
			}
		    }

		}
	    }
	}
      else
	{
	  for ( i=1; i <= grid->nn; i++ )
	    {
	      if ( grid->node_order[i] < 3 )
		continue;

	      m = nnsn2[i+1] - nnsn2[i];
	      n = 5;

	      for( var=0; var < NUM_VAR; var++)
		{
		  for ( j=nnsn2[i]; j < nnsn2[i+1]; j++ )
		    {
		      node = nsn2[j];

		      // Get the geometric weight term for current node.
		      dx[0] = x[node*NDIM+0] - x[i*NDIM+0];
		      dx[1] = x[node*NDIM+1] - x[i*NDIM+1];
		      
		      dist = sqrt(dx[0]*dx[0] + dx[1]*dx[1]);
		      
		      w = 1.0 / pow(dist,p);

		      row = j - nnsn2[i];

		      U[row] = w*(q[node*NUM_VAR+var] - q[i*NUM_VAR+var]);

                      if (DEBUG_HESS)
			U[row] = w*(cv_avg[node]-cv_avg[i]);
		    }

		  // Do the multiplication Q^t*b
		  QR_QtB ( grid->Qfac[i], U, temp1, m );
		  QR_Backsub ( grid->Rfac[i], temp1, n );
	      
		  hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 0] = temp1[0];
		  hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 1] = temp1[1];
		  hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 2] = temp1[2];
		  hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 3] = temp1[3];
		  hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 4] = temp1[4];
	      
		  point_value[i*NUM_VAR+var] = q[i*NUM_VAR+var] - 
		                               moments[i*NUM_MOM+0]*temp1[0] -
		                               moments[i*NUM_MOM+1]*temp1[1] -
		                               moments[i*NUM_MOM+2]*temp1[2]*0.5 -
		                               moments[i*NUM_MOM+3]*temp1[3]*0.5 -
		                               moments[i*NUM_MOM+4]*temp1[4];

		  if ( DEBUG_HESS)
		    {
		      point_value[i*NUM_VAR+var] = cv_avg[i] -
			                           moments[i*NUM_MOM+0]*temp1[0] -
		                                   moments[i*NUM_MOM+1]*temp1[1] -
		                                   moments[i*NUM_MOM+2]*temp1[2]*0.5 -
		                                   moments[i*NUM_MOM+3]*temp1[3]*0.5 -
		                                   moments[i*NUM_MOM+4]*temp1[4];
		    }

		  if ( DEBUG_HESS && var == 0 )
		    {
		      if ( (fabs( x[i*NDIM+0] - 0.5 ) < 1.0e-10) && (fabs( x[i*NDIM+1] ) < 1.0e-10 ) ) // the node I'll study.
			{
			  printf("Node %d is at 0.5,0.0 :\n",i);
			  printf("point values = %.15E\n",point_value[i*NUM_VAR+0]);
			  printf("fx = %.15E\n",temp1[0]);
			  printf("fy = %.15E\n",temp1[1]);
			  printf("fxx = %.15E\n",temp1[2]);
			  printf("fyy = %.15E\n",temp1[3]);
			  printf("fxy = %.15E\n",temp1[4]);
			}
		    }
		}
	    }
	}

    }
  else if ( GSFAC)
    {
      if ( HESS_FULL_SYSTEM )
	{
	  for ( i=1; i <= nn; i++ )
	    {
	      if ( grid->node_order[i] < 3 )
		continue;

	      num_neighbors = nnsn2[i+1] - nnsn2[i];
	      n = 6;
	  
	      num_neighbors++;  // Include the node itself.
	      m = num_neighbors;

	      // Put in the mean constraint.
	      A[0 + 0*num_neighbors] = 1.0;
	      A[0 + 1*num_neighbors] = moments[i*NUM_MOM+0]; // x
	      A[0 + 2*num_neighbors] = moments[i*NUM_MOM+1]; // y
	      A[0 + 3*num_neighbors] = 0.5*moments[i*NUM_MOM+2]; // xx
	      A[0 + 4*num_neighbors] = 0.5*moments[i*NUM_MOM+3]; // yy
	      A[0 + 5*num_neighbors] = moments[i*NUM_MOM+4]; // xy

	      for ( j=nnsn2[i]; j < nnsn2[i+1]; j++ )
		{
		  node = nsn2[j];
		  
		  // Get the geometric weight term for current node.
		  dx[0] = x[node*NDIM+0] - x[i*NDIM+0];
		  dx[1] = x[node*NDIM+1] - x[i*NDIM+1];

		  dist = sqrt(dx[0]*dx[0] + dx[1]*dx[1]);
		  
		  w = 1.0 / pow(dist,p);
		  
		  row = j - nnsn2[i] + 1;
	      
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

	      // Finished contruction of LHS.
	      

	      // Now we start the Gram-Schmidt Orthonormalization algorithm. The matrix 'aa' will be used here to
	      // contain the orthogonal basis.
	      
	      // Set the first basis vector.
	      for ( j=0; j < m; j++ )
		{
		  aa[j] = A[j];
		}
	      
	      normalizeVectorGeneral(&(aa[0]),m);

	      // Construct the remaining basis vectors.
	      for ( j=1; j < n; j++ )
		{
		  // Copy over the current column vector from A.
		  for ( k=0; k < m; k++ )
		    {
		      aa[j*m + k] = A[j*m + k];
		    }

		  // Subtract off components from the previous vectors.
		  for ( k=0; k < j; k++ )
		    {
		      innerprod = dotProductGeneral( &(aa[k*m]) , &(A[j*m]) , m );

		      for ( row=0; row < m; row++ )
			{
			  aa[j*m + row] -= innerprod*aa[k*m+row];
			}
		    }

		  // Normalize the vector.
		  normalizeVectorGeneral( &(aa[j*m]),m );
		}

	      // We now have constructed the orthogonal matrix Q. Now work on the
	      // right hand side.
	      for( var=0; var < NUM_VAR; var++)
		{
		  U[0] = q[i*NUM_VAR+var];
		  
                  if (DEBUG_HESS)
		    U[0] = cv_avg[i];

		  for ( j=nnsn2[i]; j < nnsn2[i+1]; j++ )
		    {
		      node = nsn2[j];

		      // Get the geometric weight term for current node.
		      dx[0] = x[node*NDIM+0] - x[i*NDIM+0];
		      dx[1] = x[node*NDIM+1] - x[i*NDIM+1];

		      dist = sqrt(dx[0]*dx[0] + dx[1]*dx[1]);
		      
		      w = 1.0 / pow(dist,p);

		      row = j - nnsn2[i] + 1;

		      U[row] = w*q[node*NUM_VAR+var];

		      if (DEBUG_HESS)
			U[row] = w*cv_avg[node];
		    }

		  // Now do the linear algebra Rx = Q'b.
		  for ( j=0; j < n; j++ )
		    {
		      vt[j] = dotProductGeneral( &(aa[j*m]), U, m );
		    }

		  // And finally do back substitution to get the result.
		  for ( j=n-1; j >= 0; j-- )
		    {
		      sol[j] = vt[j];

		      for ( k=j+1; k < n; k++ )
			{
			  innerprod = dotProductGeneral( &(aa[j*m]), &(A[k*m]), m );

			  sol[j] -= sol[k]*innerprod;
			}

		      innerprod = dotProductGeneral( &(aa[j*m]), &(A[j*m]), m );

		      sol[j] = sol[j] / innerprod;
		    }

		  // Store the solution.
		  hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 0] = sol[1];
		  hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 1] = sol[2];
		  hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 2] = sol[3];
		  hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 3] = sol[4];
		  hess[i*NUM_VAR*NUM_MOM + var*NUM_MOM + 4] = sol[5];
	      
		  point_value[i*NUM_VAR+var] = sol[0];

		  if ( DEBUG_HESS && var == 0 )
		    {
		      if ( (fabs( x[i*NDIM+0] - 0.5 ) < 1.0e-10) && (fabs( x[i*NDIM+1] ) < 1.0e-10 ) ) // the node I'll study.
			{
			  printf("Node %d is at 0.5,0.0 :\n",i);
			  printf("point values = %.15E\n",sol[0]);
			  printf("fx = %.15E\n",sol[1]);
			  printf("fy = %.15E\n",sol[2]);
			  printf("fxx = %.15E\n",sol[3]);
			  printf("fyy = %.15E\n",sol[4]);
			  printf("fxy = %.15E\n",sol[5]);
			}
		    }
		}
	    }
	}
      else // Reduced system.
	{
	  for ( i=1; i <= nn; i++ )
	    {
	      if ( grid->node_order[i] < 3 )
		continue; 
	      
	      num_neighbors = nnsn2[i+1] - nnsn2[i];
	      n = 5;
	      m = num_neighbors;
	      
	      for ( j=nnsn2[i]; j < nnsn2[i+1]; j++ )
		{
		  node = nsn2[j];
		  
		  // Get the geometric weight term for current node.
		  dx[0] = x[node*NDIM+0] - x[i*NDIM+0];
		  dx[1] = x[node*NDIM+1] - x[i*NDIM+1];
		  
		  dist = sqrt(dx[0]*dx[0] + dx[1]*dx[1]);
		  
		  w = 1.0 / pow(dist,p);

		  row = j - nnsn2[i];
		  
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

	      // Finished contruction of LHS.

	      // Now we start the Gram-Schmidt Orthonormalization algorithm. The matrix 'aa' will be used here to
	      // contain the orthogonal basis.
	      
	      // Set the first basis vector.
	      for ( j=0; j < m; j++ )
		{
		  aa[j] = A[j];
		}
	      
	      normalizeVectorGeneral(&(aa[0]),m);

	      // Construct the remaining basis vectors.
	      for ( j=1; j < n; j++ )
		{
		  // Copy over the current column vector from A.
		  for ( k=0; k < m; k++ )
		    {
		      aa[j*m + k] = A[j*m + k];
		    }

		  // Subtract off components from the previous vectors.
		  for ( k=0; k < j; k++ )
		    {
		      innerprod = dotProductGeneral( &(aa[k*m]) , &(A[j*m]) , m );

		      for ( row=0; row < m; row++ )
			{
			  aa[j*m + row] -= innerprod*aa[k*m+row];
			}
		    }

		  // Normalize the vector.
		  normalizeVectorGeneral( &(aa[j*m]),m );
		}

	      // Work on the right hand side now.
	      for( var=0; var < NUM_VAR; var++)
		{
		  for ( j=nnsn2[i]; j < nnsn2[i+1]; j++ )
		    {
		      node = nsn2[j];

		      // Get the geometric weight term for current node.
		      dx[0] = x[node*NDIM+0] - x[i*NDIM+0];
		      dx[1] = x[node*NDIM+1] - x[i*NDIM+1];

		      dist = sqrt(dx[0]*dx[0] + dx[1]*dx[1]);
		      
		      w = 1.0 / pow(dist,p);

		      row = j - nnsn2[i];

		      U[row] = w*(q[node*NUM_VAR+var] - q[i*NUM_VAR+var]);

		      if (DEBUG_HESS)
			U[row] = w*(cv_avg[node]-cv_avg[i]);
		    }

		  // Now do the linear algebra Rx = Q'b.
		  for ( j=0; j < n; j++ )
		    {
		      vt[j] = dotProductGeneral( &(aa[j*m]), U, m );
		    }

		  // And finally do back substitution to get the result.
		  for ( j=n-1; j >= 0; j-- )
		    {
		      sol[j] = vt[j];

		      for ( k=j+1; k < n; k++ )
			{
			  innerprod = dotProductGeneral( &(aa[j*m]), &(A[k*m]), m );

			  sol[j] -= sol[k]*innerprod;
			}

		      innerprod = dotProductGeneral( &(aa[j*m]), &(A[j*m]), m );

		      sol[j] = sol[j] / innerprod;
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

		  if ( DEBUG_HESS)
		    {
		      point_value[i*NUM_VAR+var] = cv_avg[i] -
			                           moments[i*NUM_MOM+0]*sol[0] -
		                                   moments[i*NUM_MOM+1]*sol[1] -
		                                   moments[i*NUM_MOM+2]*sol[2]*0.5 -
		                                   moments[i*NUM_MOM+3]*sol[3]*0.5 -
		                                   moments[i*NUM_MOM+4]*sol[4];
		    }


		  if ( DEBUG_HESS && var == 0 )
		    {
		      if ( (fabs( x[i*NDIM+0] - 0.5 ) < 1.0e-10) && (fabs( x[i*NDIM+1] ) < 1.0e-10 ) ) // the node I'll study.
			{
			  printf("Node %d is at 0.5,0.0 :\n",i);
			  printf("point values = %.15E\n",point_value[i*NUM_VAR+0]);
			  printf("fx = %.15E\n",sol[0]);
			  printf("fy = %.15E\n",sol[1]);
			  printf("fxx = %.15E\n",sol[2]);
			  printf("fyy = %.15E\n",sol[3]);
			  printf("fxy = %.15E\n",sol[4]);
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

  if ( DEBUG_HESS )
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
      
      fp_pv1 = fopen("hess_pv1.dat","w");
      fp_pv2 = fopen("hess_pv2.dat","w");
      fp_pv3 = fopen("hess_pv3.dat","w");
      fp_pv4 = fopen("hess_pv4.dat","w");
      
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

  if ( DEBUG_HESS )
    {
      FILE *fp_hess = NULL;
      fp_hess = fopen("hess_values.dat","w");

      for ( i=1; i <= grid->nn; i++ )
	{
	  for ( j=0; j < 5; j++ )
	    {
	      fprintf(fp_hess,"%.15e\n",hess[i*NUM_VAR*NUM_MOM + 0*NUM_MOM + j]);
	    }
	}
      fclose(fp_hess);
    }
  else if (0)
    {
      FILE *fp_h1 = NULL;
      FILE *fp_h2 = NULL;
      FILE *fp_h3 = NULL;
      FILE *fp_h4 = NULL;
      
      fp_h1 = fopen("hess_q1.dat","w");
      fp_h2 = fopen("hess_q2.dat","w");
      fp_h3 = fopen("hess_q3.dat","w");
      fp_h4 = fopen("hess_q4.dat","w");
      
      for ( i=1; i <= grid->nn; i++ )
	{
	  for ( j=0; j < 5; j++ )
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
  
  //#if DEBUG_HESS
  if (0)
    {
      printf("HESS DEBUG:\n");
      for ( i=1; i <= grid->nn; i++ )
	{
	  printf("Node %d:\n",i);
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      printf("  Var %d: qx = %f\n",j,(float)grid->hess[i*NUM_MOM*NUM_VAR + j*NUM_MOM + 0]);
	      printf("          qy = %f\n",(float)grid->hess[i*NUM_MOM*NUM_VAR + j*NUM_MOM + 1]);
	      printf("         qxx = %f\n",(float)grid->hess[i*NUM_MOM*NUM_VAR + j*NUM_MOM + 2]);
	      printf("         qyy = %f\n",(float)grid->hess[i*NUM_MOM*NUM_VAR + j*NUM_MOM + 3]);
	      printf("         qxy = %f\n",(float)grid->hess[i*NUM_MOM*NUM_VAR + j*NUM_MOM + 4]);
	    }
	}
      fflush(stdout);
    }

  // Now we may do some analysis on the result for debug cases.
  if ( DEBUG_HESS )
    {
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

      // Compare what the Hessian routine returned against
      // what I'm expecting it to be. Print out any large differences.
      // Print out the norm.
      qprime = (double*)malloc((nn+1)*NUM_MOM*sizeof(double));
      if ( qprime == NULL ) { printf("MEMORY ERROR: Could not allocate 'qprime'.\n"); exit(0); }

      // Set the analytical values.
      for ( i=1; i <= nn; i++ )
	{
	  Test_Function_Derivatives(&(x[i*NDIM]),&(qprime[i*NUM_MOM]));
	}

      for ( i=0; i < 5; i++ )
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
      for ( i=0; i < 5; i++ )
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
      for ( i=0; i < 5; i++ )
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
      for( i=0; i < 5; i++ )
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
	  for ( i=0; i < 5; i++ )
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
	  for ( i=0; i < 5; i++ )
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
      for ( i=1; i < 5; i++ )
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
	  for ( i=0; i < 5; i++ )
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
      fp1 = fopen("debug_gauss_points_extrap.dat","w");

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
	  //fprintf(fp1,"%.15e\n",QL[0]);

	  recon_val = (point_value[nodeL*NUM_VAR] +
		       hess[nodeL*NUM_VAR*NUM_MOM+0*NUM_MOM+0]*dx[0] +
		       hess[nodeL*NUM_VAR*NUM_MOM+0*NUM_MOM+1]*dx[1] +
		       hess[nodeL*NUM_VAR*NUM_MOM+0*NUM_MOM+2]*dx[0]*dx[0]*0.5 +
		       hess[nodeL*NUM_VAR*NUM_MOM+0*NUM_MOM+3]*dx[1]*dx[1]*0.5 +
		       hess[nodeL*NUM_VAR*NUM_MOM+0*NUM_MOM+4]*dx[0]*dx[1]);

	  fprintf(fp1,"%.15e\n",recon_val);

	  temp = Test_Function ( gp[0] );

	  temp = fabs( recon_val - temp );
	  
	  norm_1_s += fabs(temp);
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
	  //fprintf(fp1,"%.15e\n",QL[0]);

	  recon_val = (point_value[nodeL*NUM_VAR] +
		       hess[nodeL*NUM_VAR*NUM_MOM+0*NUM_MOM+0]*dx[0] +
		       hess[nodeL*NUM_VAR*NUM_MOM+0*NUM_MOM+1]*dx[1] +
		       hess[nodeL*NUM_VAR*NUM_MOM+0*NUM_MOM+2]*dx[0]*dx[0]*0.5 +
		       hess[nodeL*NUM_VAR*NUM_MOM+0*NUM_MOM+3]*dx[1]*dx[1]*0.5 +
		       hess[nodeL*NUM_VAR*NUM_MOM+0*NUM_MOM+4]*dx[0]*dx[1]);

	  fprintf(fp1,"%.15e\n",recon_val);
	  
	  temp = Test_Function ( gp[0] );

	  temp = fabs( recon_val - temp );

	  norm_1_s += fabs(temp);
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
	  //fprintf(fp1,"%.15e\n",QL[0]);

	  recon_val = (point_value[nodeL*NUM_VAR] +
		       hess[nodeL*NUM_VAR*NUM_MOM+0*NUM_MOM+0]*dx[0] +
		       hess[nodeL*NUM_VAR*NUM_MOM+0*NUM_MOM+1]*dx[1] +
		       hess[nodeL*NUM_VAR*NUM_MOM+0*NUM_MOM+2]*dx[0]*dx[0]*0.5 +
		       hess[nodeL*NUM_VAR*NUM_MOM+0*NUM_MOM+3]*dx[1]*dx[1]*0.5 +
		       hess[nodeL*NUM_VAR*NUM_MOM+0*NUM_MOM+4]*dx[0]*dx[1]);

	  fprintf(fp1,"%.15e\n",recon_val);
 
	  temp = Test_Function ( gp[0] );

	  temp = fabs( recon_val - temp );

	  norm_1_s += fabs(temp);
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
	  //fprintf(fp1,"%.15e\n",QR[0]);
 
	  recon_val = (point_value[nodeR*NUM_VAR] +
		       hess[nodeR*NUM_VAR*NUM_MOM+0*NUM_MOM+0]*dx[0] +
		       hess[nodeR*NUM_VAR*NUM_MOM+0*NUM_MOM+1]*dx[1] +
		       hess[nodeR*NUM_VAR*NUM_MOM+0*NUM_MOM+2]*dx[0]*dx[0]*0.5 +
		       hess[nodeR*NUM_VAR*NUM_MOM+0*NUM_MOM+3]*dx[1]*dx[1]*0.5 +
		       hess[nodeR*NUM_VAR*NUM_MOM+0*NUM_MOM+4]*dx[0]*dx[1]);

	  fprintf(fp1,"%.15e\n",recon_val);

	  temp = Test_Function ( gp[0] );

	  temp = fabs( recon_val - temp );
	  
	  norm_1_s += fabs(temp);
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
	  //fprintf(fp1,"%.15e\n",QR[0]);

	  recon_val = (point_value[nodeR*NUM_VAR] +
		       hess[nodeR*NUM_VAR*NUM_MOM+0*NUM_MOM+0]*dx[0] +
		       hess[nodeR*NUM_VAR*NUM_MOM+0*NUM_MOM+1]*dx[1] +
		       hess[nodeR*NUM_VAR*NUM_MOM+0*NUM_MOM+2]*dx[0]*dx[0]*0.5 +
		       hess[nodeR*NUM_VAR*NUM_MOM+0*NUM_MOM+3]*dx[1]*dx[1]*0.5 +
		       hess[nodeR*NUM_VAR*NUM_MOM+0*NUM_MOM+4]*dx[0]*dx[1]);

	  fprintf(fp1,"%.15e\n",recon_val);
	  
	  temp = Test_Function ( gp[0] );

	  temp = fabs( recon_val - temp );

	  norm_1_s += fabs(temp);
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
	  //fprintf(fp1,"%.15e\n",QR[0]);
 
	  recon_val = (point_value[nodeR*NUM_VAR] +
		       hess[nodeR*NUM_VAR*NUM_MOM+0*NUM_MOM+0]*dx[0] +
		       hess[nodeR*NUM_VAR*NUM_MOM+0*NUM_MOM+1]*dx[1] +
		       hess[nodeR*NUM_VAR*NUM_MOM+0*NUM_MOM+2]*dx[0]*dx[0]*0.5 +
		       hess[nodeR*NUM_VAR*NUM_MOM+0*NUM_MOM+3]*dx[1]*dx[1]*0.5 +
		       hess[nodeR*NUM_VAR*NUM_MOM+0*NUM_MOM+4]*dx[0]*dx[1]);

	  fprintf(fp1,"%.15e\n",recon_val);

	  temp = Test_Function ( gp[0] );

	  temp = fabs( recon_val - temp );

	  norm_1_s += fabs(temp);
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

      fp1 = fopen("debug_gauss_points_extrap_boundary.dat","w");

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

	      recon_val = (point_value[real_node*NUM_VAR] +
			   hess[real_node*NUM_VAR*NUM_MOM+0*NUM_MOM+0]*dx[0] +
			   hess[real_node*NUM_VAR*NUM_MOM+0*NUM_MOM+1]*dx[1] +
			   hess[real_node*NUM_VAR*NUM_MOM+0*NUM_MOM+2]*dx[0]*dx[0]*0.5 +
			   hess[real_node*NUM_VAR*NUM_MOM+0*NUM_MOM+3]*dx[1]*dx[1]*0.5 +
			   hess[real_node*NUM_VAR*NUM_MOM+0*NUM_MOM+4]*dx[0]*dx[1]);

	      fprintf(fp1,"%.15e\n",recon_val);
	      
	      temp = Test_Function ( gp[0] );
	      
	      temp = fabs( recon_val - temp );
	  
	      norm_1_b += fabs(temp);
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

	      recon_val = (point_value[real_node*NUM_VAR] +
			   hess[real_node*NUM_VAR*NUM_MOM+0*NUM_MOM+0]*dx[0] +
			   hess[real_node*NUM_VAR*NUM_MOM+0*NUM_MOM+1]*dx[1] +
			   hess[real_node*NUM_VAR*NUM_MOM+0*NUM_MOM+2]*dx[0]*dx[0]*0.5 +
			   hess[real_node*NUM_VAR*NUM_MOM+0*NUM_MOM+3]*dx[1]*dx[1]*0.5 +
			   hess[real_node*NUM_VAR*NUM_MOM+0*NUM_MOM+4]*dx[0]*dx[1]);

	      fprintf(fp1,"%.15e\n",recon_val);
	      
	      temp = Test_Function ( gp[0] );

	      temp = fabs( recon_val - temp );

	      norm_1_b += fabs(temp);
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

	      recon_val = (point_value[real_node*NUM_VAR] +
			   hess[real_node*NUM_VAR*NUM_MOM+0*NUM_MOM+0]*dx[0] +
			   hess[real_node*NUM_VAR*NUM_MOM+0*NUM_MOM+1]*dx[1] +
			   hess[real_node*NUM_VAR*NUM_MOM+0*NUM_MOM+2]*dx[0]*dx[0]*0.5 +
			   hess[real_node*NUM_VAR*NUM_MOM+0*NUM_MOM+3]*dx[1]*dx[1]*0.5 +
			   hess[real_node*NUM_VAR*NUM_MOM+0*NUM_MOM+4]*dx[0]*dx[1]);

	      fprintf(fp1,"%.15e\n",recon_val);
	      
	      temp = Test_Function ( gp[0] );
	      
	      temp = fabs( recon_val - temp );

	      norm_1_b += fabs(temp);
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

	      recon_val = (point_value[real_node*NUM_VAR] +
			   hess[real_node*NUM_VAR*NUM_MOM+0*NUM_MOM+0]*dx[0] +
			   hess[real_node*NUM_VAR*NUM_MOM+0*NUM_MOM+1]*dx[1] +
			   hess[real_node*NUM_VAR*NUM_MOM+0*NUM_MOM+2]*dx[0]*dx[0]*0.5 +
			   hess[real_node*NUM_VAR*NUM_MOM+0*NUM_MOM+3]*dx[1]*dx[1]*0.5 +
			   hess[real_node*NUM_VAR*NUM_MOM+0*NUM_MOM+4]*dx[0]*dx[1]);
	      
	      temp = Test_Function ( gp[0] );
	      
	      temp = fabs( recon_val - temp );
	  
	      norm_1_b += fabs(temp);
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

	      recon_val = (point_value[real_node*NUM_VAR] +
			   hess[real_node*NUM_VAR*NUM_MOM+0*NUM_MOM+0]*dx[0] +
			   hess[real_node*NUM_VAR*NUM_MOM+0*NUM_MOM+1]*dx[1] +
			   hess[real_node*NUM_VAR*NUM_MOM+0*NUM_MOM+2]*dx[0]*dx[0]*0.5 +
			   hess[real_node*NUM_VAR*NUM_MOM+0*NUM_MOM+3]*dx[1]*dx[1]*0.5 +
			   hess[real_node*NUM_VAR*NUM_MOM+0*NUM_MOM+4]*dx[0]*dx[1]);
	      
	      temp = Test_Function ( gp[0] );
	      
	      temp = fabs( recon_val - temp );
	  
	      norm_1_b += fabs(temp);
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

	      recon_val = (point_value[real_node*NUM_VAR] +
			   hess[real_node*NUM_VAR*NUM_MOM+0*NUM_MOM+0]*dx[0] +
			   hess[real_node*NUM_VAR*NUM_MOM+0*NUM_MOM+1]*dx[1] +
			   hess[real_node*NUM_VAR*NUM_MOM+0*NUM_MOM+2]*dx[0]*dx[0]*0.5 +
			   hess[real_node*NUM_VAR*NUM_MOM+0*NUM_MOM+3]*dx[1]*dx[1]*0.5 +
			   hess[real_node*NUM_VAR*NUM_MOM+0*NUM_MOM+4]*dx[0]*dx[1]);
	      
	      temp = Test_Function ( gp[0] );
	      
	      temp = fabs( recon_val - temp );
	  
	      norm_1_b += fabs(temp);
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

      printf("L1 norm of difference between reconstructions is %.15E\n",norm_rm);

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
			  ext_val[0] = (point_value[node*NUM_VAR] +
					hess[node*NUM_VAR*NUM_MOM+0*NUM_MOM+0]*dx[0] +
					hess[node*NUM_VAR*NUM_MOM+0*NUM_MOM+1]*dx[1] +
					hess[node*NUM_VAR*NUM_MOM+0*NUM_MOM+2]*dx[0]*dx[0]*0.5 +
					hess[node*NUM_VAR*NUM_MOM+0*NUM_MOM+3]*dx[1]*dx[1]*0.5 +
					hess[node*NUM_VAR*NUM_MOM+0*NUM_MOM+4]*dx[0]*dx[1]);

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
			  ext_val[i] = (point_value[node*NUM_VAR] +
					hess[node*NUM_VAR*NUM_MOM+0*NUM_MOM+0]*dx[0] +
					hess[node*NUM_VAR*NUM_MOM+0*NUM_MOM+1]*dx[1] +
					hess[node*NUM_VAR*NUM_MOM+0*NUM_MOM+2]*dx[0]*dx[0]*0.5 +
					hess[node*NUM_VAR*NUM_MOM+0*NUM_MOM+3]*dx[1]*dx[1]*0.5 +
					hess[node*NUM_VAR*NUM_MOM+0*NUM_MOM+4]*dx[0]*dx[1]);

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
			  ext_val[0] = (point_value[node*NUM_VAR] +
					hess[node*NUM_VAR*NUM_MOM+0*NUM_MOM+0]*dx[0] +
					hess[node*NUM_VAR*NUM_MOM+0*NUM_MOM+1]*dx[1] +
					hess[node*NUM_VAR*NUM_MOM+0*NUM_MOM+2]*dx[0]*dx[0]*0.5 +
					hess[node*NUM_VAR*NUM_MOM+0*NUM_MOM+3]*dx[1]*dx[1]*0.5 +
					hess[node*NUM_VAR*NUM_MOM+0*NUM_MOM+4]*dx[0]*dx[1]);

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
			  ext_val[i] = (point_value[node*NUM_VAR] +
					hess[node*NUM_VAR*NUM_MOM+0*NUM_MOM+0]*dx[0] +
					hess[node*NUM_VAR*NUM_MOM+0*NUM_MOM+1]*dx[1] +
					hess[node*NUM_VAR*NUM_MOM+0*NUM_MOM+2]*dx[0]*dx[0]*0.5 +
					hess[node*NUM_VAR*NUM_MOM+0*NUM_MOM+3]*dx[1]*dx[1]*0.5 +
					hess[node*NUM_VAR*NUM_MOM+0*NUM_MOM+4]*dx[0]*dx[1]);

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

  if ( DEBUG_HESS )
    {
      for ( i=0; i < NUM_MOM; i++ )
	freenull(derivs[i]);

      freenull(derivs);
    }

  freenull(A);
  freenull(U);
  freenull(sol);
  freenull(temp1);
  freenull(temp2);
  freenull(temp3);
  freenull(aa);
  freenull(s);
  freenull(u);
  freenull(vt);
  freenull(work);
  freenull(qprime);
  freenull(cv_avg);
  freenull(ext_cv_avg);
  
  return;
}




//=============================================================
// 
//  Test_Function()
//
//  Test a simple scalar function.
//
//  double *x;                       // Coordinates.
//
//=============================================================

double Test_Function ( double *x )
{
  double f;                         // Function value.
  double pi = M_PI;
  double r, r2;

  if ( ANNULUS )
    {
      double rho,u,v,P,E;
      double rhoi,Mi,U,Ui,Ri;
      rhoi = 1.0;
      Mi = 2.0;
      Ui = 2.0;
      Ri = 2.0;
      
      r2 = (x[0]*x[0]) + (x[1]*x[1]);
      r = sqrt(r2);
      
      rho = pow( ( 1.0 + ((1.4 - 1.0)/2.0)*Mi*Mi*( 1.0 - (Ri*Ri)/r2)) , (1.0/(1.4-1.0)) );
      U   = (Ui*Ri)/r;
      u   = (x[1]*U)/r;
      v   = (-x[0]*U)/r;
      P   = ( pow( rho , 1.4) )/1.4;
      
      E   = P / (1.4 - 1.) + 0.5*rho*(u*u + v*v);

      // conserved: rho,rho*u,rho*v,E
      // primitive: rho,u,v,P

      return (E);
    }

  if ( 0 )
    {
      double rho,u,v,P,E;

      rho = 1.0 + sin(pi*x[0])*sin(pi*x[1]);
      u   = 0.25 + 0.25*sin(pi*x[0])*cos(2.*pi*x[1]);
      v   = 0.25 + 0.25*cos(2.*pi*x[0])*sin(pi*x[1]);
      P   = 1.0/1.4 + 0.05*cos(2.*pi*x[0])*cos(2.*pi*x[1]);

      E   = P/(1.4-1.0) + 0.5*rho*(u*u + v*v);

      // conserved: rho,rho*u,rho*v,E
      // primitive: rho,u,v,P

      return (rho);
    }

  // f = x
  //f = x[0];

  // f = y
  //f = x[1];

  // f = x^2
  //f = x[0]*x[0];

  // f = x^3
  //f = x[0]*x[0]*x[0];

  // f = x + 5y + 3x^2 + 4y^2 + 10xy
  //f = x[0] + 5.*x[1] + 3.*x[0]*x[0] + 4.*x[1]*x[1] + 10.*x[0]*x[1];

  // f = x + 2y + 3x^2 + 4y^2 + 7xy + 1.2yx^2 + 3.7xy^2 + 9x^3 + 12y^3
  //f = x[0] + 2.*x[1] + 3.*x[0]*x[0] + 4.*x[1]*x[1] + 7.*x[0]*x[1] + 1.2*x[0]*x[0]*x[1] +
  //    3.7*x[0]*x[1]*x[1] + 9.*x[0]*x[0]*x[0] + 12.*x[1]*x[1]*x[1];


  // Test functions used in my thesis.
  // f = x^2 + y^2 + xy + x + y
  //f = x[0]*x[0] + x[1]*x[1] + x[0]*x[1] + x[0] + x[1];

  // f = 3x^3 + 5xy^2
  //f = 3.0*x[0]*x[0]*x[0] + 5.0*x[0]*x[1]*x[1];

  // f = sin(pi x) cos(pi y )
  //f = sin(pi*x[0])*cos(pi*x[1]);

  // f = e^(-r^2)
  r2 = (x[0]-0.5)*(x[0]-0.5) + (x[1]-0.5)*(x[1]-0.5);
  r = sqrt ( r2 );
  f = exp(-r2);

  return f;
}


//=============================================================
// 
//  Test_Function_Divergent()
//
//  Integral of a simple scalar function wrt x.
//
//  double *x;                       // Coordinates.
//
//=============================================================

double Test_Function_Divergent ( double *x )
{
  double f;                         // Function value.
  double pi = M_PI;

  //rho
  //f = x[0] - (1.0/pi)*(cos(pi*x[0])*sin(pi*x[1]));

  // u
  //f = 0.25*x[0] - (0.25/pi)*(cos(pi*x[0])*cos(2.*pi*x[1]));

  // v
  //f = 0.25*x[0] + (0.125/pi)*sin(2.*pi*x[0])*sin(pi*x[1]);

  // P
  //f = (1.0/1.4)*x[0] + (0.025/pi)*sin(2.*pi*x[0])*cos(2.*pi*x[1]);

  // f = x
  //f = 0.5*x[0]*x[0];

  // f = y
  //f = x[0]*x[1];

  // f = x^2
  //f = (1./3.)*x[0]*x[0]*x[0];

  // f = x^3
  //f = 0.25*x[0]*x[0]*x[0]*x[0];

  // f = x + 5y + 3x^2 + 4y^2 + 10xy
  //f = 0.5*x[0]*x[0] + 5.*x[1]*x[0] + x[0]*x[0]*x[0] + 4.*x[0]*x[1]*x[1] + 5.*x[0]*x[0]*x[1];

  // f = x + 2y + 3x^2 + 4y^2 + 7xy + 1.2yx^2 + 3.7xy^2 + 9x^3 + 12y^3
  //f = 0.5*x[0]*x[0] + 2.*x[0]*x[1] + x[0]*x[0]*x[0] + 4.*x[1]*x[1]*x[0] + 3.5*x[0]*x[0]*x[1] + 0.4*x[0]*x[0]*x[0]*x[1] +
  //  (3.7/2.)*x[0]*x[0]*x[1]*x[1] + (9./4.)*x[0]*x[0]*x[0]*x[0] + 12.*x[0]*x[1]*x[1]*x[1];

  // Test functions used in my thesis.
  
  // f =  1/3 x^3 + xy^2 + 1/2 x^2 y + 1/2 x^2 + xy
  //f = (1./3.)*x[0]*x[0]*x[0] + x[0]*x[1]*x[1] + 0.5*x[0]*x[0]*x[1] + 0.5*x[0]*x[0] + x[0]*x[1];

  // f = 3/4 x^4 + 5/2 x^2 y^2
  //f = 0.75*x[0]*x[0]*x[0]*x[0] + 2.5*x[0]*x[0]*x[1]*x[1];

  // f = -pi^-1 * cos(pi*x) cos(pi*y)
  //f = (-1./pi)*cos(pi*x[0])*cos(pi*x[1]);

  // f = e^(-r^2) Does not have a closed form solution. -->> Error function.

  return f;
}

//=============================================================
// 
//  Test_Function_Derivatives()
//
//  Integral of a simple scalar function wrt x.
//
//  double *x;                       // Coordinates.
//  double *derivs;                  // Function derivatives.
//
//=============================================================

void Test_Function_Derivatives ( double *x, double *derivs )
{
  double pi = M_PI,r,r2,f;

  // derivs[0] = fx
  // derivs[1] = fy
  // derivs[2] = fxx
  // derivs[3] = fyy
  // derivs[4] = fxy
  // derivs[5] = fxxx
  // derivs[6] = fyyy
  // derivs[7] = fxxy
  // derivs[8] = fxyy

  // rho
  //derivs[0] = pi*cos(pi*x[0])*sin(pi*x[1]);
  //derivs[1] = pi*sin(pi*x[0])*cos(pi*x[1]);
  //derivs[2] = -pi*pi*sin(pi*x[0])*sin(pi*x[1]);
  //derivs[3] = -pi*pi*sin(pi*x[0])*sin(pi*x[1]);
  //derivs[4] = pi*pi*cos(pi*x[0])*cos(pi*x[1]);
  //derivs[5] = -pi*pi*pi*cos(pi*x[0])*sin(pi*x[1]);
  //derivs[6] = -pi*pi*pi*sin(pi*x[0])*sin(pi*x[1]);
  //derivs[7] = -pi*pi*pi*sin(pi*x[0])*sin(pi*x[1]);
  //derivs[8] = -pi*pi*pi*cos(pi*x[0])*sin(pi*x[1]);

  // u
  //derivs[0] = 0.25*pi*cos(pi*x[0])*cos(2.*pi*x[1]);
  //derivs[1] = -0.5*pi*sin(pi*x[0])*sin(2.*pi*x[1]);
  //derivs[2] = -0.25*pi*pi*sin(pi*x[0])*cos(2.*pi*x[1]);
  //derivs[3] = -pi*pi*sin(pi*x[0])*cos(2.*pi*x[1]);
  //derivs[4] = -0.5*pi*pi*cos(pi*x[0])*sin(2.*pi*x[1]);
  //derivs[5] = -0.25*pi*pi*pi*cos(pi*x[0])*cos(2.*pi*x[1]);
  //derivs[6] = 2.*pi*pi*pi*sin(pi*x[0])*sin(2.*pi*x[1]);
  //derivs[7] = 0.5*pi*pi*pi*sin(pi*x[0])*sin(2.*pi*x[1]);
  //derivs[8] = -pi*pi*pi*cos(pi*x[0])*cos(2.*pi*x[1]);

  // v
  //derivs[0] = -0.5*pi*sin(2.*pi*x[0])*sin(pi*x[1]);
  //derivs[1] = 0.25*pi*cos(2.*pi*x[0])*cos(pi*x[1]);
  //derivs[2] = -pi*pi*cos(2.*pi*x[0])*sin(pi*x[1]);
  //derivs[3] = -0.25*pi*pi*cos(2.*pi*x[0])*sin(pi*x[1]);
  //derivs[4] = -0.5*pi*pi*sin(2.*pi*x[0])*cos(pi*x[1]);
  //derivs[5] = 2.*pi*pi*pi*sin(2.*pi*x[0])*sin(pi*x[1]);
  //derivs[6] = -0.25*pi*pi*pi*cos(2.*pi*x[0])*cos(pi*x[1]);
  //derivs[7] = -pi*pi*pi*cos(2.*pi*x[0])*cos(pi*x[1]);
  //derivs[8] = 0.5*pi*pi*pi*sin(2.*pi*x[0])*sin(pi*x[1]);

  // P
  //derivs[0] = -0.1*pi*sin(2.*pi*x[0])*cos(2.*pi*x[1]);
  //derivs[1] = -0.1*pi*cos(2.*pi*x[0])*sin(2.*pi*x[1]);
  //derivs[2] = -0.2*pi*pi*cos(2.*pi*x[0])*cos(2.*pi*x[1]);
  //derivs[3] = -0.2*pi*pi*cos(2.*pi*x[0])*cos(2.*pi*x[1]);
  //derivs[4] = 0.2*pi*pi*sin(2.*pi*x[0])*sin(2.*pi*x[1]);
  //derivs[5] = 0.4*pi*pi*pi*sin(2.*pi*x[0])*cos(2.*pi*x[1]);
  //derivs[6] = 0.4*pi*pi*pi*cos(2.*pi*x[0])*sin(2.*pi*x[1]);
  //derivs[7] = 0.4*pi*pi*pi*cos(2.*pi*x[0])*sin(2.*pi*x[1]);
  //derivs[8] = 0.4*pi*pi*pi*sin(2.*pi*x[0])*cos(2.*pi*x[1]);


  // f = x
  //derivs[0] = 1.0;
  //derivs[1] = 0.;
  //derivs[2] = 0.;
  //derivs[3] = 0.;
  //derivs[4] = 0.;
  //derivs[5] = 0.;
  //derivs[6] = 0.;
  //derivs[7] = 0.;
  //derivs[8] = 0.;
  
  // f = y
  //derivs[0] = 0.;
  //derivs[1] = 1.0;
  //derivs[2] = 0.;
  //derivs[3] = 0.;
  //derivs[4] = 0.;
  //derivs[5] = 0.;
  //derivs[6] = 0.;
  //derivs[7] = 0.;
  //derivs[8] = 0.;

  // f = x^2
  //derivs[0] = 2.*x[0];
  //derivs[1] = 0.;
  //derivs[2] = 2.;
  //derivs[3] = 0.;
  //derivs[4] = 0.;
  //derivs[5] = 0.;
  //derivs[6] = 0.;
  //derivs[7] = 0.;
  //derivs[8] = 0.;

  // f = x^3
  //derivs[0] = 3.*x[0]*x[0];
  //derivs[1] = 0.;
  //derivs[2] = 6.*x[0];
  //derivs[3] = 0.;
  //derivs[4] = 0.;
  //derivs[5] = 6.0;
  //derivs[6] = 0.;
  //derivs[7] = 0.;
  //derivs[8] = 0.;

  // f = x + 5y + 3x^2 + 4y^2 + 10xy
  //derivs[0] = 1. + 6.*x[0] + 10.*x[1];
  //derivs[1] = 5. + 8.*x[1] + 10.*x[0];
  //derivs[2] = 6.;
  //derivs[3] = 8.;
  //derivs[4] = 10.;
  //derivs[5] = 0.;
  //derivs[6] = 0.;
  //derivs[7] = 0.;
  //derivs[8] = 0.;

  // f = x + 2y + 3x^2 + 4y^2 + 7xy + 1.2yx^2 + 3.7xy^2 + 9x^3 + 12y^3
  //derivs[0] = 1. + 6.*x[0] + 7.*x[1] + 2.4*x[0]*x[1] + 3.7*x[1]*x[1] + 27.*x[0]*x[0];
  //derivs[1] = 2. + 8.*x[1] + 7.*x[0] + 1.2*x[0]*x[0] + 7.4*x[0]*x[1] + 36.*x[1]*x[1];
  //derivs[2] = 6. + 2.4*x[1] + 54.*x[0];
  //derivs[3] = 8. + 7.4*x[0] + 72.*x[1];
  //derivs[4] = 7. + 2.4*x[0] + 7.4*x[1];
  //derivs[5] = 54.0;
  //derivs[6] = 72.0;
  //derivs[7] = 2.4;
  //derivs[8] = 7.4;
  

  // Test functions from my thesis.

  // f = x^2 + y^2 + xy + x + y
  //derivs[0] = 2.*x[0] + x[1] + 1.;
  //derivs[1] = 2.*x[1] + x[0] + 1.;
  //derivs[2] = 2.;
  //derivs[3] = 2.;
  //derivs[4] = 1.;
  //derivs[5] = 0.;
  //derivs[6] = 0.;
  //derivs[7] = 0.;
  //derivs[8] = 0.;

  // f = 3x^3 + 5xy^2
  //derivs[0] = 9.*x[0]*x[0] + 5.*x[1]*x[1];
  //derivs[1] = 10.*x[0]*x[1];
  //derivs[2] = 18.*x[0];
  //derivs[3] = 10.*x[0];
  //derivs[4] = 10.*x[1];
  //derivs[5] = 18.;
  //derivs[6] = 0.;
  //derivs[7] = 0.;
  //derivs[8] = 10.;

  // f = sin(pi x) cos(pi y )
  //derivs[0] = pi*cos(pi*x[0])*cos(pi*x[1]);
  //derivs[1] = -pi*sin(pi*x[0])*sin(pi*x[1]);
  //derivs[2] = -pi*pi*sin(pi*x[0])*cos(pi*x[1]);
  //derivs[3] = -pi*pi*sin(pi*x[0])*cos(pi*x[1]);
  //derivs[4] = -pi*pi*cos(pi*x[0])*sin(pi*x[1]);
  //derivs[5] = -pi*pi*pi*cos(pi*x[0])*cos(pi*x[1]);
  //derivs[6] = pi*pi*pi*sin(pi*x[0])*sin(pi*x[1]);
  //derivs[7] = pi*pi*pi*sin(pi*x[0])*sin(pi*x[1]);
  //derivs[8] = -pi*pi*pi*cos(pi*x[0])*cos(pi*x[1]);

  // f = e^(-r^2)
  r2 = (x[0]-0.5)*(x[0]-0.5) + (x[1]-0.5)*(x[1]-0.5);
  r = sqrt ( r2 );
  f = exp(-r2);
  derivs[0] = (-2.*x[0] + 1.) * f;
  derivs[1] = (-2.*x[1] + 1.) * f;
  derivs[2] = -2. * f + (-2.*x[0] + 1.)*(-2.*x[0] + 1.)*f;
  derivs[3] = -2. * f + (-2.*x[1] + 1.)*(-2.*x[1] + 1.)*f;
  derivs[4] = (-2.*x[0]+1.)*(-2.*x[1]+1.)*f;
  derivs[5] = -6.*(-2.*x[0] + 1.)*f + (-2.*x[0] + 1.)*(-2.*x[0] + 1.)*(-2.*x[0] + 1.)*f;
  derivs[6] = -6.*(-2.*x[1] + 1.)*f + (-2.*x[1] + 1.)*(-2.*x[1] + 1.)*(-2.*x[1] + 1.)*f;
  derivs[7] = -2.*(-2.*x[1] + 1.)*f + (-2.*x[0] + 1.)*(-2.*x[0] + 1.)*(-2.*x[1] + 1.)*f;
  derivs[8] = -2.*(-2.*x[0] + 1.)*f + (-2.*x[0] + 1.)*(-2.*x[1] + 1.)*(-2.*x[1] + 1.)*f; 


  // LEFT OVER STUFF -----------------------------------------------------------------
  //derivs[0] = 2.5;
  //derivs[1] = 12.0;
  //derivs[2] = 0.;
  //derivs[3] = 0.;
  //derivs[4] = 0.;

  //derivs[0] = 20.*x[0];
  //derivs[1] = 24.*x[1];
  //derivs[2] = 20.;
  //derivs[3] = 24.;
  //derivs[4] = 0.;
  
  return;
}

//=============================================================
// 
//  Get_Annulus_Derivatives
//
//  Computes the derivatives of the analytical solution to the
//  annulus problem for all variables.
//
//  double *X;                       // Coordinates.
//  double *derivs;                  // Function derivatives.
//
//=============================================================

void Get_Annulus_Derivatives ( double *X, double *derivs )
{
  double gamma, gm1, oogm1;
  double Ri, r, Ri2, r2, Mi, Ui, U;
  double rhoi, rho, u, v, rhou, rhov;
  double P, E;
  double x,y;
  
  double drdx,drdy,dr2dxx,dr2dyy,dr2dxy,dr3dxxx,dr3dyyy,dr3dxxy,dr3dxyy;
  double dudx,dudy,du2dxx,du2dyy,du2dxy,du3dxxx,du3dyyy,du3dxxy,du3dxyy;
  double dvdx,dvdy,dv2dxx,dv2dyy,dv2dxy,dv3dxxx,dv3dyyy,dv3dxxy,dv3dxyy;
  double dPdx,dPdy,dP2dxx,dP2dyy,dP2dxy,dP3dxxx,dP3dyyy,dP3dxxy,dP3dxyy;
  
  double drudx,drudy,dru2dxx,dru2dyy,dru2dxy,dru3dxxx,dru3dyyy,dru3dxxy,dru3dxyy;
  double drvdx,drvdy,drv2dxx,drv2dyy,drv2dxy,drv3dxxx,drv3dyyy,drv3dxxy,drv3dxyy;
  double dEdx,dEdy,dE2dxx,dE2dyy,dE2dxy,dE3dxxx,dE3dyyy,dE3dxxy,dE3dxyy;

  double term, term2;

  // Set the variables.
  gamma = 1.4;
  gm1 = gamma - 1.0;
  oogm1 = 1.0 / ( gamma - 1.0 );

  Ri = 2.0;
  Ri2 = Ri*Ri;
  Mi = 2.0;
  rhoi = 1.0;
  Ui = Mi * ( pow( rhoi , (gm1*0.5) ));

  x = X[0];
  y = X[1];

  r2 = x*x + y*y;
  r = sqrt(r2);

  rho = pow( ( 1.0 + gm1*0.5*Mi*Mi*(1.0 - (Ri2/r2))) , (1.0/gm1) );
  U = (Ui*Ri)/r;
  u = (y*U)/r;
  v = -1.0 * (x*U)/r;
  P = ( pow( rho , gamma ) )/gamma;

  E = P / gm1 + 0.5*rho*(u*u + v*v);
  
  // set up the value of a common term that appears frequently in the derivatives.
  term = ( 1.0 + gm1*0.5*Mi*Mi*(1.0 - (Ri2/r2)));

  term2 = ( (y * y * Ui * Ui * Ri * Ri)/(pow(r2,2.)) + (x * x * Ui * Ui * Ri * Ri)/(pow(r2,2.)) );

  // now I reset what Pressure is to simplify the derivatives below. This term appears frequently.
  P = ( pow( rho , gamma ) );

  if ( RECON_PRIM )
    {
      // rho derivatives.
      drdx = ( 2.0 * rho * gm1 * 0.5 * Mi*Mi * Ri2 * x ) / ( gm1 * r2*r2 * term );

      drdy = ( 2.0 * rho * gm1 * 0.5 * Mi*Mi * Ri2 * y ) / ( gm1 * r2*r2 * term );

      dr2dxx =   (  4.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x*x )    / ( (pow(gm1,2.)) * (pow(r2,4.)) * (pow(term,2.)) )
	       - (  8.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * x*x )    / ( (pow(gm1,1.)) * (pow(r2,3.)) * (pow(term,1.)) )
               + (  2.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) )          / ( (pow(gm1,1.)) * (pow(r2,2.)) * (pow(term,1.)) )
               - (  4.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x*x )    / ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,2.)) );

      dr2dyy =   (  4.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y*y )    / ( (pow(gm1,2.)) * (pow(r2,4.)) * (pow(term,2.)) )
	       - (  8.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * y*y )    / ( (pow(gm1,1.)) * (pow(r2,3.)) * (pow(term,1.)) )
               + (  2.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) )          / ( (pow(gm1,1.)) * (pow(r2,2.)) * (pow(term,1.)) )
               - (  4.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y*y )    / ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,2.)) );

      dr2dxy =   (  4.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y*x )    / ( (pow(gm1,2.)) * (pow(r2,4.)) * (pow(term,2.)) )
	       - (  8.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * x*y )    / ( (pow(gm1,1.)) * (pow(r2,3.)) * (pow(term,1.)) )
       	       - (  4.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x*y )    / ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,2.)) );

      dr3dxxx =   (  8.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * x*x*x ) / ( (pow(gm1,3.)) * (pow(r2,6.)) * (pow(term,3.)) )
                - ( 48.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x*x*x ) / ( (pow(gm1,2.)) * (pow(r2,5.)) * (pow(term,2.)) )
	        + ( 12.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x )     / ( (pow(gm1,2.)) * (pow(r2,4.)) * (pow(term,2.)) )
                - ( 24.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * x*x*x ) / ( (pow(gm1,2.)) * (pow(r2,6.)) * (pow(term,3.)) )
	        + ( 48.0 * rho *       gm1*0.5       * (pow(Mi,2.)) * (pow(Ri,2.)) * x*x*x ) / (      gm1      * (pow(r2,4.)) * term )
	        - ( 24.0 * rho *       gm1*0.5       * (pow(Mi,2.)) * (pow(Ri,2.)) * x )     / (      gm1      * (pow(r2,3.)) * term )
	        + ( 48.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x*x*x ) / (      gm1      * (pow(r2,5.)) * (pow(term,2.)) )
	        - ( 12.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x )     / (      gm1      * (pow(r2,4.)) * (pow(term,2.)) )
        	+ ( 16.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * x*x*x ) / (      gm1      * (pow(r2,6.)) * (pow(term,3.)) );

      dr3dyyy =   (  8.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * y*y*y ) / ( (pow(gm1,3.)) * (pow(r2,6.)) * (pow(term,3.)) )
                - ( 48.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y*y*y ) / ( (pow(gm1,2.)) * (pow(r2,5.)) * (pow(term,2.)) )
	        + ( 12.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y )     / ( (pow(gm1,2.)) * (pow(r2,4.)) * (pow(term,2.)) )
                - ( 24.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * y*y*y ) / ( (pow(gm1,2.)) * (pow(r2,6.)) * (pow(term,3.)) )
	        + ( 48.0 * rho *       gm1*0.5       * (pow(Mi,2.)) * (pow(Ri,2.)) * y*y*y ) / (      gm1      * (pow(r2,4.)) * term )
	        - ( 24.0 * rho *       gm1*0.5       * (pow(Mi,2.)) * (pow(Ri,2.)) * y )     / (      gm1      * (pow(r2,3.)) * term )
	        + ( 48.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y*y*y ) / (      gm1      * (pow(r2,5.)) * (pow(term,2.)) )
	        - ( 12.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y )     / (      gm1      * (pow(r2,4.)) * (pow(term,2.)) )
        	+ ( 16.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * y*y*y ) / (      gm1      * (pow(r2,6.)) * (pow(term,3.)) );

      dr3dxxy =   (  8.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * y*x*x ) / ( (pow(gm1,3.)) * (pow(r2,6.)) * (pow(term,3.)) )
                - ( 48.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x*x*y ) / ( (pow(gm1,2.)) * (pow(r2,5.)) * (pow(term,2.)) )
	        - ( 24.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * x*x*y ) / ( (pow(gm1,2.)) * (pow(r2,6.)) * (pow(term,3.)) )
                + ( 48.0 * rho *       gm1*0.5       * (pow(Mi,2.)) * (pow(Ri,2.)) * x*x*y ) / (      gm1      * (pow(r2,4.)) * term )
         	+ ( 48.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x*x*y ) / (      gm1      * (pow(r2,5.)) * (pow(term,2.)) )
	        + (  4.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y )     / ( (pow(gm1,2.)) * (pow(r2,4.)) * (pow(term,2.)) )
	        - (  8.0 * rho *       gm1*0.5       * (pow(Mi,2.)) * (pow(Ri,2.)) * y )     / (      gm1      * (pow(r2,3.)) * term )
	        - (  4.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y )     / (      gm1      * (pow(r2,4.)) * (pow(term,2.)) )
        	+ ( 16.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * x*x*y ) / (      gm1      * (pow(r2,6.)) * (pow(term,3.)) );

      
      dr3dxyy =   (  8.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * y*y*x ) / ( (pow(gm1,3.)) * (pow(r2,6.)) * (pow(term,3.)) )
                - ( 48.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y*y*x ) / ( (pow(gm1,2.)) * (pow(r2,5.)) * (pow(term,2.)) )
	        + (  4.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x )     / ( (pow(gm1,2.)) * (pow(r2,4.)) * (pow(term,2.)) )
                - ( 24.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * y*y*x ) / ( (pow(gm1,2.)) * (pow(r2,6.)) * (pow(term,3.)) )
         	+ ( 48.0 * rho *       gm1*0.5       * (pow(Mi,2.)) * (pow(Ri,2.)) * x*y*y ) / (      gm1      * (pow(r2,4.)) * term )
	        + ( 48.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x*y*y ) / (      gm1      * (pow(r2,5.)) * (pow(term,2.)) )
	        - (  8.0 * rho *       gm1*0.5       * (pow(Mi,2.)) * (pow(Ri,2.)) * x )     / (      gm1      * (pow(r2,3.)) * term )
	        + ( 16.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * x*y*y ) / (      gm1      * (pow(r2,6.)) * (pow(term,3.)) )
        	- (  4.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x )     / (      gm1      * (pow(r2,4.)) * (pow(term,2.)) );

      // u derivatives.
      dudx = ( -2.0 * y * Ui * Ri * x ) / ( pow(r2,2.) );
      dudy = ( Ui * Ri ) / ( pow(r2,1.) ) - ( 2.0 * y * y * Ui * Ri ) / ( pow(r2,2.) );

      du2dxx = ( 8.0 * y * Ui * Ri * x * x ) / ( pow(r2,3.) ) - ( 2.0 * y * Ui * Ri ) / ( pow(r2,2.) );
      du2dyy = (-6.0 * y * Ui * Ri ) / ( pow(r2,2.) ) + ( 8.0 * y * y * y * Ui * Ri ) / ( pow(r2,3.) );
      du2dxy = (-2.0 * Ui * Ri * x ) / ( pow(r2,2.) ) + ( 8.0 * y * y * Ui * Ri * x ) / ( pow(r2,3.) );

      du3dxxx = (-48.0 * y * Ui * Ri * x * x * x ) / ( pow(r2,4.) ) + ( 24.0 * y * Ui * Ri * x ) / ( pow(r2,3.) );
      du3dyyy = (-6.0 * Ui * Ri ) / ( pow(r2,2.) ) + ( 48.0 * y * y * Ui * Ri ) / ( pow(r2,3.) ) - ( 48.0 * y * y * y * y * Ui * Ri ) / ( pow(r2,4.) );
      du3dxxy = ( 8.0 * Ui * Ri * x * x ) / ( pow(r2,3.) ) - ( 48.0 * y * y * Ui * Ri * x * x ) / ( pow(r2,4.) ) - ( 2.0 * Ui * Ri ) / ( pow(r2,2.) ) + ( 8.0 * y * y * Ui * Ri ) / ( pow(r2,3.) );
      du3dxyy = ( 24.0 * y * Ui * Ri * x ) / ( pow(r2,3.) ) - ( 48.0 * y * y * y * Ui * Ri * x ) / ( pow(r2,4.) );

      // v derivatives.
      dvdx = -1.*( Ui * Ri ) / ( pow(r2,1.) ) + ( 2.0 * x * x * Ui * Ri ) / ( pow(r2,2.) );
      dvdy = ( 2.0 * y * Ui * Ri * x ) / ( pow(r2,2.) );

      dv2dxx = ( 6.0 * Ui * Ri * x ) / ( pow(r2,2.) ) - ( 8.0 * x * x * x * Ui * Ri ) / ( pow(r2,3.) );
      dv2dyy = ( 2.0 * Ui * Ri * x ) / ( pow(r2,2.) ) - ( 8.0 * y * y * Ui * Ri * x ) / ( pow(r2,3.) );
      dv2dxy = ( 2.0 * y * Ui * Ri ) / ( pow(r2,2.) ) - ( 8.0 * y * Ui * Ri * x * x ) / ( pow(r2,3.) );

      dv3dxxx = (-48.0 * Ui * Ri * x * x ) / ( pow(r2,3.) ) + ( 6.0 * Ui * Ri ) / ( pow(r2,2.) ) + ( 48.0 * x * x * x * x *  Ui * Ri ) / ( pow(r2,4.) );
      dv3dyyy = (-24.0 * y * Ui * Ri * x ) / ( pow(r2,3.) ) + ( 48.0 * y * y * y * Ui * Ri * x ) / ( pow(r2,4.) );
      dv3dxxy = (-24.0 * y * Ui * Ri * x ) / ( pow(r2,3.) ) + ( 48.0 * y * Ui * Ri * x * x * x ) / ( pow(r2,4.) );
      dv3dxyy = ( 2.0 * Ui * Ri ) / ( pow(r2,2.) ) - ( 8.0 * y * y * Ui * Ri ) / ( pow(r2,3.) ) - ( 8.0 * Ui * Ri * x * x ) / ( pow(r2,3.) ) + ( 48.0 * y * y * Ui * Ri * x * x ) / ( pow(r2,4.) );

      // Pressure derivatives.
      
      dPdx = (  2.0 * P * (pow(gm1*0.5,1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * x ) / ( (pow(gm1,1.)) * (pow(r2,2.)) * (pow(term,1.)) );
      dPdy = (  2.0 * P * (pow(gm1*0.5,1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * y ) / ( (pow(gm1,1.)) * (pow(r2,2.)) * (pow(term,1.)) );

      dP2dxx = (  4.0 * P * (pow(gamma,1.)) * (pow(gm1*0.5,2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x * x ) / ( (pow(gm1,2.)) * (pow(r2,4.)) * (pow(term,2.)) )
	     - (  8.0 * P *                   (pow(gm1*0.5,1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * x * x ) / ( (pow(gm1,1.)) * (pow(r2,3.)) * (pow(term,1.)) )
	     + (  2.0 * P *                   (pow(gm1*0.5,1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) )         / ( (pow(gm1,1.)) * (pow(r2,2.)) * (pow(term,1.)) )
	     - (  4.0 * P *                   (pow(gm1*0.5,2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x * x ) / ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,2.)) );

      dP2dyy = (  4.0 * P * (pow(gamma,1.)) * (pow(gm1*0.5,2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y * y ) / ( (pow(gm1,2.)) * (pow(r2,4.)) * (pow(term,2.)) )
	     - (  8.0 * P *                   (pow(gm1*0.5,1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * y * y ) / ( (pow(gm1,1.)) * (pow(r2,3.)) * (pow(term,1.)) )
	     + (  2.0 * P *                   (pow(gm1*0.5,1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) )         / ( (pow(gm1,1.)) * (pow(r2,2.)) * (pow(term,1.)) )
	     - (  4.0 * P *                   (pow(gm1*0.5,2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y * y ) / ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,2.)) );

      dP2dxy = (  4.0 * P * (pow(gamma,1.)) * (pow(gm1*0.5,2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y * x ) / ( (pow(gm1,2.)) * (pow(r2,4.)) * (pow(term,2.)) )
	     - (  8.0 * P *                   (pow(gm1*0.5,1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * x * y ) / ( (pow(gm1,1.)) * (pow(r2,3.)) * (pow(term,1.)) )
	     - (  4.0 * P *                   (pow(gm1*0.5,2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x * y ) / ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,2.)) );

      dP3dxxx = (  8.0 * P * (pow(gamma,2.)) * (pow(gm1*0.5,3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * x * x * x ) / ( (pow(gm1,3.)) * (pow(r2,6.)) * (pow(term,3.)) )
	      - ( 48.0 * P * (pow(gamma,1.)) * (pow(gm1*0.5,2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x * x * x ) / ( (pow(gm1,2.)) * (pow(r2,5.)) * (pow(term,2.)) )
	      + ( 12.0 * P * (pow(gamma,1.)) * (pow(gm1*0.5,2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x )         / ( (pow(gm1,2.)) * (pow(r2,4.)) * (pow(term,2.)) )
	      - ( 24.0 * P * (pow(gamma,1.)) * (pow(gm1*0.5,3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * x * x * x ) / ( (pow(gm1,2.)) * (pow(r2,6.)) * (pow(term,3.)) )              
	      + ( 48.0 * P *                   (pow(gm1*0.5,1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * x * x * x ) / ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,1.)) )
	      - ( 24.0 * P *                   (pow(gm1*0.5,1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * x )         / ( (pow(gm1,1.)) * (pow(r2,3.)) * (pow(term,1.)) )
	      + ( 48.0 * P *                   (pow(gm1*0.5,2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x * x * x ) / ( (pow(gm1,1.)) * (pow(r2,5.)) * (pow(term,2.)) )
	      - ( 12.0 * P *                   (pow(gm1*0.5,2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x )         / ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,2.)) )
	      + ( 16.0 * P *                   (pow(gm1*0.5,3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * x * x * x ) / ( (pow(gm1,1.)) * (pow(r2,6.)) * (pow(term,3.)) );

      dP3dyyy = (  8.0 * P * (pow(gamma,2.)) * (pow(gm1*0.5,3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * y * y * y ) / ( (pow(gm1,3.)) * (pow(r2,6.)) * (pow(term,3.)) )
	      - ( 48.0 * P * (pow(gamma,1.)) * (pow(gm1*0.5,2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y * y * y ) / ( (pow(gm1,2.)) * (pow(r2,5.)) * (pow(term,2.)) )
	      + ( 12.0 * P * (pow(gamma,1.)) * (pow(gm1*0.5,2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y )         / ( (pow(gm1,2.)) * (pow(r2,4.)) * (pow(term,2.)) )
	      - ( 24.0 * P * (pow(gamma,1.)) * (pow(gm1*0.5,3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * y * y * y ) / ( (pow(gm1,2.)) * (pow(r2,6.)) * (pow(term,3.)) )              
	      + ( 48.0 * P *                   (pow(gm1*0.5,1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * y * y * y ) / ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,1.)) )
	      - ( 24.0 * P *                   (pow(gm1*0.5,1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * y )         / ( (pow(gm1,1.)) * (pow(r2,3.)) * (pow(term,1.)) )
	      + ( 48.0 * P *                   (pow(gm1*0.5,2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y * y * y ) / ( (pow(gm1,1.)) * (pow(r2,5.)) * (pow(term,2.)) )
	      - ( 12.0 * P *                   (pow(gm1*0.5,2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y )         / ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,2.)) )
	      + ( 16.0 * P *                   (pow(gm1*0.5,3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * y * y * y ) / ( (pow(gm1,1.)) * (pow(r2,6.)) * (pow(term,3.)) );

      dP3dxxy = (  8.0 * P * (pow(gamma,2.)) * (pow(gm1*0.5,3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * y * x * x ) / ( (pow(gm1,3.)) * (pow(r2,6.)) * (pow(term,3.)) )
	      - ( 48.0 * P * (pow(gamma,1.)) * (pow(gm1*0.5,2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x * x * y ) / ( (pow(gm1,2.)) * (pow(r2,5.)) * (pow(term,2.)) )
	      - ( 24.0 * P * (pow(gamma,1.)) * (pow(gm1*0.5,3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * x * x * y ) / ( (pow(gm1,2.)) * (pow(r2,6.)) * (pow(term,3.)) )
	      + ( 48.0 * P *                   (pow(gm1*0.5,1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * x * x * y ) / ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,1.)) )              
	      + ( 48.0 * P *                   (pow(gm1*0.5,2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x * x * y ) / ( (pow(gm1,1.)) * (pow(r2,5.)) * (pow(term,2.)) )
	      + (  4.0 * P * (pow(gamma,1.)) * (pow(gm1*0.5,2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y )         / ( (pow(gm1,2.)) * (pow(r2,4.)) * (pow(term,2.)) )
	      - (  8.0 * P *                   (pow(gm1*0.5,1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * y )         / ( (pow(gm1,1.)) * (pow(r2,3.)) * (pow(term,1.)) )
	      - (  4.0 * P *                   (pow(gm1*0.5,2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y )         / ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,2.)) )
	      + ( 16.0 * P *                   (pow(gm1*0.5,3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * x * x * y ) / ( (pow(gm1,1.)) * (pow(r2,6.)) * (pow(term,3.)) );

      dP3dxyy = (  8.0 * P * (pow(gamma,2.)) * (pow(gm1*0.5,3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * y * y * x ) / ( (pow(gm1,3.)) * (pow(r2,6.)) * (pow(term,3.)) )
	      - ( 48.0 * P * (pow(gamma,1.)) * (pow(gm1*0.5,2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y * y * x ) / ( (pow(gm1,2.)) * (pow(r2,5.)) * (pow(term,2.)) )
	      + (  4.0 * P * (pow(gamma,1.)) * (pow(gm1*0.5,2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x )         / ( (pow(gm1,2.)) * (pow(r2,4.)) * (pow(term,2.)) )
	      - ( 24.0 * P * (pow(gamma,1.)) * (pow(gm1*0.5,3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * y * y * x ) / ( (pow(gm1,2.)) * (pow(r2,6.)) * (pow(term,3.)) )              
	      + ( 48.0 * P *                   (pow(gm1*0.5,1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * x * y * y ) / ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,1.)) )
	      + ( 48.0 * P *                   (pow(gm1*0.5,2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x * y * y ) / ( (pow(gm1,1.)) * (pow(r2,5.)) * (pow(term,2.)) )
	      - (  8.0 * P *                   (pow(gm1*0.5,1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * x )         / ( (pow(gm1,1.)) * (pow(r2,3.)) * (pow(term,1.)) )
	      + ( 16.0 * P *                   (pow(gm1*0.5,3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * x * y * y ) / ( (pow(gm1,1.)) * (pow(r2,6.)) * (pow(term,3.)) )
	      - (  4.0 * P *                   (pow(gm1*0.5,2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x )         / ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,2.)) );
      
      
      // set the derivatives.
      derivs[0] = drdx;
      derivs[1] = drdy;
      derivs[2] = dr2dxx;
      derivs[3] = dr2dyy;
      derivs[4] = dr2dxy;
      derivs[5] = dr3dxxx;
      derivs[6] = dr3dyyy;
      derivs[7] = dr3dxxy;
      derivs[8] = dr3dxyy;

      derivs[9]  = dudx;
      derivs[10] = dudy;
      derivs[11] = du2dxx;
      derivs[12] = du2dyy;
      derivs[13] = du2dxy;
      derivs[14] = du3dxxx;
      derivs[15] = du3dyyy;
      derivs[16] = du3dxxy;
      derivs[17] = du3dxyy;

      derivs[18] = dvdx;
      derivs[19] = dvdy;
      derivs[20] = dv2dxx;
      derivs[21] = dv2dyy;
      derivs[22] = dv2dxy;
      derivs[23] = dv3dxxx;
      derivs[24] = dv3dyyy;
      derivs[25] = dv3dxxy;
      derivs[26] = dv3dxyy;

      derivs[27] = dPdx;
      derivs[28] = dPdy;
      derivs[29] = dP2dxx;
      derivs[30] = dP2dyy;
      derivs[31] = dP2dxy;
      derivs[32] = dP3dxxx;
      derivs[33] = dP3dyyy;
      derivs[34] = dP3dxxy;
      derivs[35] = dP3dxyy;
    }
  else
    {
      // rho derivatives.
      drdx = ( 2.0 * rho * gm1 * 0.5 * Mi*Mi * Ri2 * x ) / ( gm1 * r2*r2 * term );

      drdy = ( 2.0 * rho * gm1 * 0.5 * Mi*Mi * Ri2 * y ) / ( gm1 * r2*r2 * term );

      dr2dxx =   (  4.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x*x )    / ( (pow(gm1,2.)) * (pow(r2,4.)) * (pow(term,2.)) )
	       - (  8.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * x*x )    / ( (pow(gm1,1.)) * (pow(r2,3.)) * (pow(term,1.)) )
               + (  2.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) )          / ( (pow(gm1,1.)) * (pow(r2,2.)) * (pow(term,1.)) )
               - (  4.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x*x )    / ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,2.)) );

      dr2dyy =   (  4.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y*y )    / ( (pow(gm1,2.)) * (pow(r2,4.)) * (pow(term,2.)) )
	       - (  8.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * y*y )    / ( (pow(gm1,1.)) * (pow(r2,3.)) * (pow(term,1.)) )
               + (  2.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) )          / ( (pow(gm1,1.)) * (pow(r2,2.)) * (pow(term,1.)) )
               - (  4.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y*y )    / ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,2.)) );

      dr2dxy =   (  4.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y*x )    / ( (pow(gm1,2.)) * (pow(r2,4.)) * (pow(term,2.)) )
	       - (  8.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * x*y )    / ( (pow(gm1,1.)) * (pow(r2,3.)) * (pow(term,1.)) )
	       - (  4.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x*y )    / ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,2.)) );

      dr3dxxx =   (  8.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * x*x*x ) / ( (pow(gm1,3.)) * (pow(r2,6.)) * (pow(term,3.)) )
                - ( 48.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x*x*x ) / ( (pow(gm1,2.)) * (pow(r2,5.)) * (pow(term,2.)) )
	        + ( 12.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x )     / ( (pow(gm1,2.)) * (pow(r2,4.)) * (pow(term,2.)) )
                - ( 24.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * x*x*x ) / ( (pow(gm1,2.)) * (pow(r2,6.)) * (pow(term,3.)) )
	        + ( 48.0 * rho *       gm1*0.5       * (pow(Mi,2.)) * (pow(Ri,2.)) * x*x*x ) / (      gm1      * (pow(r2,4.)) * term )
	        - ( 24.0 * rho *       gm1*0.5       * (pow(Mi,2.)) * (pow(Ri,2.)) * x )     / (      gm1      * (pow(r2,3.)) * term )
	        + ( 48.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x*x*x ) / (      gm1      * (pow(r2,5.)) * (pow(term,2.)) )
	        - ( 12.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x )     / (      gm1      * (pow(r2,4.)) * (pow(term,2.)) )
        	+ ( 16.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * x*x*x ) / (      gm1      * (pow(r2,6.)) * (pow(term,3.)) );

      dr3dyyy =   (  8.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * y*y*y ) / ( (pow(gm1,3.)) * (pow(r2,6.)) * (pow(term,3.)) )
                - ( 48.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y*y*y ) / ( (pow(gm1,2.)) * (pow(r2,5.)) * (pow(term,2.)) )
	        + ( 12.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y )     / ( (pow(gm1,2.)) * (pow(r2,4.)) * (pow(term,2.)) )
                - ( 24.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * y*y*y ) / ( (pow(gm1,2.)) * (pow(r2,6.)) * (pow(term,3.)) )
	        + ( 48.0 * rho *       gm1*0.5       * (pow(Mi,2.)) * (pow(Ri,2.)) * y*y*y ) / (      gm1      * (pow(r2,4.)) * term )
	        - ( 24.0 * rho *       gm1*0.5       * (pow(Mi,2.)) * (pow(Ri,2.)) * y )     / (      gm1      * (pow(r2,3.)) * term )
	        + ( 48.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y*y*y ) / (      gm1      * (pow(r2,5.)) * (pow(term,2.)) )
	        - ( 12.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y )     / (      gm1      * (pow(r2,4.)) * (pow(term,2.)) )
        	+ ( 16.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * y*y*y ) / (      gm1      * (pow(r2,6.)) * (pow(term,3.)) );

      dr3dxxy =   (  8.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * y*x*x ) / ( (pow(gm1,3.)) * (pow(r2,6.)) * (pow(term,3.)) )
                - ( 48.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x*x*y ) / ( (pow(gm1,2.)) * (pow(r2,5.)) * (pow(term,2.)) )
	        - ( 24.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * x*x*y ) / ( (pow(gm1,2.)) * (pow(r2,6.)) * (pow(term,3.)) )
                + ( 48.0 * rho *       gm1*0.5       * (pow(Mi,2.)) * (pow(Ri,2.)) * x*x*y ) / (      gm1      * (pow(r2,4.)) * term )
         	+ ( 48.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x*x*y ) / (      gm1      * (pow(r2,5.)) * (pow(term,2.)) )
	        + (  4.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y )     / ( (pow(gm1,2.)) * (pow(r2,4.)) * (pow(term,2.)) )
	        - (  8.0 * rho *       gm1*0.5       * (pow(Mi,2.)) * (pow(Ri,2.)) * y )     / (      gm1      * (pow(r2,3.)) * term )
	        - (  4.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y )     / (      gm1      * (pow(r2,4.)) * (pow(term,2.)) )
        	+ ( 16.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * x*x*y ) / (      gm1      * (pow(r2,6.)) * (pow(term,3.)) );

      
      dr3dxyy =   (  8.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * y*y*x ) / ( (pow(gm1,3.)) * (pow(r2,6.)) * (pow(term,3.)) )
                - ( 48.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y*y*x ) / ( (pow(gm1,2.)) * (pow(r2,5.)) * (pow(term,2.)) )
	        + (  4.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x )     / ( (pow(gm1,2.)) * (pow(r2,4.)) * (pow(term,2.)) )
                - ( 24.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * y*y*x ) / ( (pow(gm1,2.)) * (pow(r2,6.)) * (pow(term,3.)) )
         	+ ( 48.0 * rho *       gm1*0.5       * (pow(Mi,2.)) * (pow(Ri,2.)) * x*y*y ) / (      gm1      * (pow(r2,4.)) * term )
	        + ( 48.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x*y*y ) / (      gm1      * (pow(r2,5.)) * (pow(term,2.)) )
	        - (  8.0 * rho *       gm1*0.5       * (pow(Mi,2.)) * (pow(Ri,2.)) * x )     / (      gm1      * (pow(r2,3.)) * term )
	        + ( 16.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * x*y*y ) / (      gm1      * (pow(r2,6.)) * (pow(term,3.)) )
        	- (  4.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x )     / (      gm1      * (pow(r2,4.)) * (pow(term,2.)) );

      // rho*u derivatives.
      drudx = ( 2.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,3.)) * x * y * Ui ) / ( (pow(gm1,1.)) * (pow(r2,3.)) * (pow(term,1.)) )
	    - ( 2.0 * rho * y * Ui * Ri * x ) / ( (pow(r2,2.)) );

      drudy = (  2.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,3.)) * y * y * Ui ) / ( (pow(gm1,1.)) * (pow(r2,3.)) * (pow(term,1.)) )
	    + (  1.0 * rho * Ui * Ri ) / ( r2 )
	    - (  2.0 * rho * y * y * Ui * Ri ) / ( (pow(r2,2.)) );

      dru2dxx = (  4.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * x * x * y * Ui ) / ( (pow(gm1,2.)) * (pow(r2,5.)) * (pow(term,2.)) )
	      - ( 16.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,3.)) * x * x * y * Ui ) / ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,1.)) )
	      + (  2.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,3.)) * y * Ui )         / ( (pow(gm1,1.)) * (pow(r2,3.)) * (pow(term,1.)) )
	      - (  4.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * x * x * y * Ui ) / ( (pow(gm1,1.)) * (pow(r2,5.)) * (pow(term,2.)) )
	      + (  8.0 * rho * y * Ui * Ri * x * x ) / ( (pow(r2,3.)) )
	      - (  2.0 * rho * y * Ui * Ri )         / ( (pow(r2,2.)) );

      dru2dyy = (  4.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * y * y * y * Ui ) / ( (pow(gm1,2.)) * (pow(r2,5.)) * (pow(term,2.)) )
	      - ( 16.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,3.)) * y * y * y * Ui ) / ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,1.)) )
	      + (  6.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,3.)) * y * Ui )         / ( (pow(gm1,1.)) * (pow(r2,3.)) * (pow(term,1.)) )
	      - (  4.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * y * y * y * Ui ) / ( (pow(gm1,1.)) * (pow(r2,5.)) * (pow(term,2.)) )
	      - (  6.0 * rho * y * Ui * Ri )         / ( (pow(r2,2.)) )
	      + (  8.0 * rho * y * y * y * Ui * Ri ) / ( (pow(r2,3.)) );

      dru2dxy = (  4.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * y * y * x * Ui ) / ( (pow(gm1,2.)) * (pow(r2,5.)) * (pow(term,2.)) )
	      - ( 16.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,3.)) * x * y * y * Ui ) / ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,1.)) )
	      - (  4.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * x * y * y * Ui ) / ( (pow(gm1,1.)) * (pow(r2,5.)) * (pow(term,2.)) )
	      + (  2.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,3.)) * x * Ui )         / ( (pow(gm1,1.)) * (pow(r2,3.)) * (pow(term,1.)) )
	      - (  2.0 * rho * Ui * Ri * x )         / ( (pow(r2,2.)) )
	      + (  8.0 * rho * y * y * Ui * Ri * x ) / ( (pow(r2,3.)) );

      dru3dxxx = (   8.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,7.)) * x * x * x * y * Ui ) / ( (pow(gm1,3.)) * (pow(r2,7.)) * (pow(term,3.)) )
	       - (  72.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * x * x * x * y * Ui ) / ( (pow(gm1,2.)) * (pow(r2,6.)) * (pow(term,2.)) )
	       + (  12.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * x * y * Ui )         / ( (pow(gm1,2.)) * (pow(r2,5.)) * (pow(term,2.)) )
	       - (  24.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,7.)) * x * x * x * y * Ui ) / ( (pow(gm1,2.)) * (pow(r2,7.)) * (pow(term,3.)) )
	       + ( 144.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,3.)) * x * x * x * y * Ui ) / ( (pow(gm1,1.)) * (pow(r2,5.)) * (pow(term,1.)) )
	       - (  48.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,3.)) * x * y * Ui )         / ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,1.)) )
	       + (  72.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * x * x * x * y * Ui ) / ( (pow(gm1,1.)) * (pow(r2,6.)) * (pow(term,2.)) )
	       - (  12.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * y * Ui * x )         / ( (pow(gm1,1.)) * (pow(r2,5.)) * (pow(term,2.)) )
	       + (  16.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,7.)) * x * x * x * y * Ui ) / ( (pow(gm1,1.)) * (pow(r2,7.)) * (pow(term,3.)) )
	       - (  48.0 * rho * y * Ui * Ri * x * x * x ) / ( (pow(r2,4.)) )
	       + (  24.0 * rho * y * Ui * Ri * x )         / ( (pow(r2,3.)) );

      dru3dyyy = (   8.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,7.)) * y * y * y * y * Ui ) / ( (pow(gm1,3.)) * (pow(r2,7.)) * (pow(term,3.)) )
	       - (  72.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * y * y * y * y * Ui ) / ( (pow(gm1,2.)) * (pow(r2,6.)) * (pow(term,2.)) )
	       + (  24.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * y * y * Ui )         / ( (pow(gm1,2.)) * (pow(r2,5.)) * (pow(term,2.)) )
	       - (  24.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,7.)) * y * y * y * y * Ui ) / ( (pow(gm1,2.)) * (pow(r2,7.)) * (pow(term,3.)) )
	       + ( 144.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,3.)) * y * y * y * y * Ui ) / ( (pow(gm1,1.)) * (pow(r2,5.)) * (pow(term,1.)) )
	       - (  96.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,3.)) * y * y * Ui )         / ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,1.)) )
	       + (  72.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * y * y * y * y * Ui ) / ( (pow(gm1,1.)) * (pow(r2,6.)) * (pow(term,2.)) )
	       - (  24.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * y * y * Ui )         / ( (pow(gm1,1.)) * (pow(r2,5.)) * (pow(term,2.)) )
	       + (   6.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,3.)) * Ui )                 / ( (pow(gm1,1.)) * (pow(r2,3.)) * (pow(term,1.)) )
	       + (  16.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,7.)) * y * y * y * y * Ui ) / ( (pow(gm1,1.)) * (pow(r2,7.)) * (pow(term,3.)) )
	       - (   6.0 * rho * Ui * Ri )                 / ( (pow(r2,2.)) )
	       + (  48.0 * rho * y * y * Ui * Ri )         / ( (pow(r2,3.)) )
	       - (  48.0 * rho * y * y * y * y * Ui * Ri ) / ( (pow(r2,4.)) );

      dru3dxxy = (   8.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,7.)) * y * y * x * x * Ui ) / ( (pow(gm1,3.)) * (pow(r2,7.)) * (pow(term,3.)) )
	       - (  72.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * x * x * y * y * Ui ) / ( (pow(gm1,2.)) * (pow(r2,6.)) * (pow(term,2.)) )
	       - (  24.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,7.)) * x * x * y * y * Ui ) / ( (pow(gm1,2.)) * (pow(r2,7.)) * (pow(term,3.)) )
	       + (   4.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * x * x * Ui )         / ( (pow(gm1,2.)) * (pow(r2,5.)) * (pow(term,2.)) )
	       + ( 144.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,3.)) * x * x * y * y * Ui ) / ( (pow(gm1,1.)) * (pow(r2,5.)) * (pow(term,1.)) )
	       + (  72.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * x * x * y * y * Ui ) / ( (pow(gm1,1.)) * (pow(r2,6.)) * (pow(term,2.)) )
	       - (  16.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,3.)) * x * x * Ui )         / ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,1.)) )
	       + (   4.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * y * y * Ui )         / ( (pow(gm1,2.)) * (pow(r2,5.)) * (pow(term,2.)) )
	       - (  16.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,3.)) * y * y * Ui )         / ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,1.)) )
	       - (   4.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * y * y * Ui )         / ( (pow(gm1,1.)) * (pow(r2,5.)) * (pow(term,2.)) )
	       + (   2.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,3.)) * Ui )                 / ( (pow(gm1,1.)) * (pow(r2,3.)) * (pow(term,1.)) )
	       + (  16.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,7.)) * x * x * y * y * Ui ) / ( (pow(gm1,1.)) * (pow(r2,7.)) * (pow(term,3.)) )
	       - (   4.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * x * x * Ui )         / ( (pow(gm1,1.)) * (pow(r2,5.)) * (pow(term,2.)) )
	       + (   8.0 * rho * Ui * Ri * x * x  )        / ( (pow(r2,3.)) )
	       - (  48.0 * rho * y * y * Ui * Ri * x * x ) / ( (pow(r2,4.)) )
	       - (   2.0 * rho * Ui * Ri )                 / ( (pow(r2,2.)) )
	       + (   8.0 * rho * y * y * Ui * Ri )         / ( (pow(r2,3.)) );

      dru3dxyy = (   8.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,7.)) * y * y * y * x * Ui ) / ( (pow(gm1,3.)) * (pow(r2,7.)) * (pow(term,3.)) )
	       - (  72.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * y * y * y * x * Ui ) / ( (pow(gm1,2.)) * (pow(r2,6.)) * (pow(term,2.)) )
	       + (  12.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * x * y * Ui )         / ( (pow(gm1,2.)) * (pow(r2,5.)) * (pow(term,2.)) )
	       - (  24.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,7.)) * y * y * y * x * Ui ) / ( (pow(gm1,2.)) * (pow(r2,7.)) * (pow(term,3.)) )
	       + ( 144.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,3.)) * x * y * y * y * Ui ) / ( (pow(gm1,1.)) * (pow(r2,5.)) * (pow(term,1.)) )
	       + (  72.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * x * y * y * y * Ui ) / ( (pow(gm1,1.)) * (pow(r2,6.)) * (pow(term,2.)) )
	       - (  48.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,3.)) * x * y * Ui )         / ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,1.)) )
	       + (  16.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,7.)) * x * y * y * y * Ui ) / ( (pow(gm1,1.)) * (pow(r2,7.)) * (pow(term,3.)) )
	       - (  12.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * y * Ui * x )         / ( (pow(gm1,1.)) * (pow(r2,5.)) * (pow(term,2.)) )
	       + (  24.0 * rho * y * Ui * Ri * x  )        / ( (pow(r2,3.)) )
	       - (  48.0 * rho * y * y * y * Ui * Ri * x ) / ( (pow(r2,4.)) );

      // rho*v derivatives.
      drvdx = ( -2.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,3.)) * x * x * Ui ) / ( (pow(gm1,1.)) * (pow(r2,3.)) * (pow(term,1.)) )
	    - (  1.0 * rho * Ui * Ri ) / ( r2 )
	    + (  2.0 * rho * x * x * Ui * Ri ) / ( (pow(r2,2.)) );

      drvdy = ( -2.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,3.)) * x * y * Ui ) / ( (pow(gm1,1.)) * (pow(r2,3.)) * (pow(term,1.)) )
	    + ( 2.0 * rho * y * Ui * Ri * x ) / ( (pow(r2,2.)) );

      drv2dxx = ( -4.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * x * x * x * Ui ) / ( (pow(gm1,2.)) * (pow(r2,5.)) * (pow(term,2.)) )
	      + ( 16.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,3.)) * x * x * x * Ui ) / ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,1.)) )
	      - (  6.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,3.)) * x * Ui )         / ( (pow(gm1,1.)) * (pow(r2,3.)) * (pow(term,1.)) )
	      + (  4.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * x * x * x * Ui ) / ( (pow(gm1,1.)) * (pow(r2,5.)) * (pow(term,2.)) )
	      + (  6.0 * rho * Ui * Ri * x )         / ( (pow(r2,2.)) )
	      - (  8.0 * rho * x * x * x * Ui * Ri ) / ( (pow(r2,3.)) );

      drv2dyy = ( -4.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * y * y * x * Ui ) / ( (pow(gm1,2.)) * (pow(r2,5.)) * (pow(term,2.)) )
	      + ( 16.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,3.)) * x * y * y * Ui ) / ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,1.)) )
	      + (  4.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * x * y * y * Ui ) / ( (pow(gm1,1.)) * (pow(r2,5.)) * (pow(term,2.)) )
	      - (  2.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,3.)) * x * Ui )         / ( (pow(gm1,1.)) * (pow(r2,3.)) * (pow(term,1.)) )
	      + (  2.0 * rho * Ui * Ri * x )         / ( (pow(r2,2.)) )
	      - (  8.0 * rho * y * y * Ui * Ri * x ) / ( (pow(r2,3.)) );

      drv2dxy = ( -4.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * x * x * y * Ui ) / ( (pow(gm1,2.)) * (pow(r2,5.)) * (pow(term,2.)) )
	      + ( 16.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,3.)) * x * x * y * Ui ) / ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,1.)) )
	      + (  4.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * x * x * y * Ui ) / ( (pow(gm1,1.)) * (pow(r2,5.)) * (pow(term,2.)) )
	      - (  2.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,3.)) * y * Ui )         / ( (pow(gm1,1.)) * (pow(r2,3.)) * (pow(term,1.)) )
	      + (  2.0 * rho * y * Ui * Ri )         / ( (pow(r2,2.)) )
	      - (  8.0 * rho * y * Ui * Ri * x * x ) / ( (pow(r2,3.)) );

      drv3dxxx = (  -8.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,7.)) * x * x * x * x * Ui ) / ( (pow(gm1,3.)) * (pow(r2,7.)) * (pow(term,3.)) )
	       + (  72.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * x * x * x * x * Ui ) / ( (pow(gm1,2.)) * (pow(r2,6.)) * (pow(term,2.)) )
	       - (  24.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * x * x * Ui )         / ( (pow(gm1,2.)) * (pow(r2,5.)) * (pow(term,2.)) )
	       + (  24.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,7.)) * x * x * x * x * Ui ) / ( (pow(gm1,2.)) * (pow(r2,7.)) * (pow(term,3.)) )
	       - ( 144.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,3.)) * x * x * x * x * Ui ) / ( (pow(gm1,1.)) * (pow(r2,5.)) * (pow(term,1.)) )
	       + (  96.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,3.)) * x * x * Ui )         / ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,1.)) )
	       - (  72.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * x * x * x * x * Ui ) / ( (pow(gm1,1.)) * (pow(r2,6.)) * (pow(term,2.)) )
	       - (   6.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,3.)) * Ui )                 / ( (pow(gm1,1.)) * (pow(r2,3.)) * (pow(term,1.)) )
	       + (  24.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * x * x * Ui )         / ( (pow(gm1,1.)) * (pow(r2,5.)) * (pow(term,2.)) )
	       - (  16.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,7.)) * x * x * x * x * Ui ) / ( (pow(gm1,1.)) * (pow(r2,7.)) * (pow(term,3.)) )
	       - (  48.0 * rho * Ui * Ri * x * x )         / ( (pow(r2,3.)) )
	       + (   6.0 * rho * Ui * Ri )                 / ( (pow(r2,2.)) )
	       + (  48.0 * rho * x * x * x * x * Ui * Ri ) / ( (pow(r2,4.)) );

      drv3dyyy = (  -8.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,7.)) * y * y * y * x * Ui ) / ( (pow(gm1,3.)) * (pow(r2,7.)) * (pow(term,3.)) )
	       + (  72.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * y * y * y * x * Ui ) / ( (pow(gm1,2.)) * (pow(r2,6.)) * (pow(term,2.)) )
	       - (  12.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * x * y * Ui )         / ( (pow(gm1,2.)) * (pow(r2,5.)) * (pow(term,2.)) )
	       + (  24.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,7.)) * y * y * y * x * Ui ) / ( (pow(gm1,2.)) * (pow(r2,7.)) * (pow(term,3.)) )
	       - ( 144.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,3.)) * x * y * y * y * Ui ) / ( (pow(gm1,1.)) * (pow(r2,5.)) * (pow(term,1.)) )
	       - (  72.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * x * y * y * y * Ui ) / ( (pow(gm1,1.)) * (pow(r2,6.)) * (pow(term,2.)) )
	       + (  48.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,3.)) * x * y * Ui )         / ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,1.)) )
	       - (  16.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,7.)) * x * y * y * y * Ui ) / ( (pow(gm1,1.)) * (pow(r2,7.)) * (pow(term,3.)) )
	       + (  12.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * y * Ui * x )         / ( (pow(gm1,1.)) * (pow(r2,5.)) * (pow(term,2.)) )
	       - (  24.0 * rho * y * Ui * Ri * x )         / ( (pow(r2,3.)) )
	       + (  48.0 * rho * y * y * y * Ui * Ri * x ) / ( (pow(r2,4.)) );

      drv3dxxy = (  -8.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,7.)) * x * x * x * y * Ui ) / ( (pow(gm1,3.)) * (pow(r2,7.)) * (pow(term,3.)) )
	       + (  72.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * x * x * x * y * Ui ) / ( (pow(gm1,2.)) * (pow(r2,6.)) * (pow(term,2.)) )
	       + (  24.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,7.)) * x * x * x * y * Ui ) / ( (pow(gm1,2.)) * (pow(r2,7.)) * (pow(term,3.)) )
	       - ( 144.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,3.)) * x * x * x * y * Ui ) / ( (pow(gm1,1.)) * (pow(r2,5.)) * (pow(term,1.)) )
	       - (  72.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * x * x * x * y * Ui ) / ( (pow(gm1,1.)) * (pow(r2,6.)) * (pow(term,2.)) )
	       - (  12.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * x * y * Ui )         / ( (pow(gm1,2.)) * (pow(r2,5.)) * (pow(term,2.)) )
	       + (  48.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,3.)) * x * y * Ui )         / ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,1.)) )
	       + (  12.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * y * Ui * x )         / ( (pow(gm1,1.)) * (pow(r2,5.)) * (pow(term,2.)) )
	       - (  16.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,7.)) * x * x * x * y * Ui ) / ( (pow(gm1,1.)) * (pow(r2,7.)) * (pow(term,3.)) )
	       - (  24.0 * rho * y * Ui * Ri * x  )        / ( (pow(r2,3.)) )
	       + (  48.0 * rho * y * Ui * Ri * x * x * x ) / ( (pow(r2,4.)) );

      drv3dxyy = (  -8.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,7.)) * y * y * x * x * Ui ) / ( (pow(gm1,3.)) * (pow(r2,7.)) * (pow(term,3.)) )
	       + (  72.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * x * x * y * y * Ui ) / ( (pow(gm1,2.)) * (pow(r2,6.)) * (pow(term,2.)) )
	       + (  24.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,7.)) * x * x * y * y * Ui ) / ( (pow(gm1,2.)) * (pow(r2,7.)) * (pow(term,3.)) )
	       - (   4.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * x * x * Ui )         / ( (pow(gm1,2.)) * (pow(r2,5.)) * (pow(term,2.)) )
	       - ( 144.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,3.)) * x * x * y * y * Ui ) / ( (pow(gm1,1.)) * (pow(r2,5.)) * (pow(term,1.)) )
	       - (  72.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * x * x * y * y * Ui ) / ( (pow(gm1,1.)) * (pow(r2,6.)) * (pow(term,2.)) )
	       + (  16.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,3.)) * x * x * Ui )         / ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,1.)) )
	       - (  16.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,7.)) * x * x * y * y * Ui ) / ( (pow(gm1,1.)) * (pow(r2,7.)) * (pow(term,3.)) )
	       + (   4.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * x * x * Ui )         / ( (pow(gm1,1.)) * (pow(r2,5.)) * (pow(term,2.)) )
	       - (   4.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * y * y * Ui )         / ( (pow(gm1,2.)) * (pow(r2,5.)) * (pow(term,2.)) )
	       + (  16.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,3.)) * y * y * Ui )         / ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,1.)) )
	       + (   4.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,5.)) * y * y * Ui )         / ( (pow(gm1,1.)) * (pow(r2,5.)) * (pow(term,2.)) )
	       - (   2.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,3.)) * Ui )                 / ( (pow(gm1,1.)) * (pow(r2,3.)) * (pow(term,1.)) )
	       + (   2.0 * rho * Ui * Ri )                 / ( (pow(r2,2.)) )
	       - (   8.0 * rho * y * y * Ui * Ri )         / ( (pow(r2,3.)) )
	       - (   8.0 * rho * Ui * Ri * x * x )         / ( (pow(r2,3.)) )
	       + (  48.0 * rho * y * y * Ui * Ri * x * x ) / ( (pow(r2,4.)) );

      // energy derivatives.

      dEdx = ( 2.0 * P * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * x ) / ( (pow(gm1,2.)) * (pow(r2,2.)) * (pow(term,1.)) )
	   + ( 1.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * x * ( (y * y * Ui * Ui * Ri * Ri)/(pow(r2,2.)) + (x * x * Ui * Ui * Ri * Ri)/(pow(r2,2.)) ) ) / 
	     ( (pow(gm1,1.)) * (pow(r2,2.)) * (pow(term,1.)) )
       	   + 0.5 * rho * ( (-4.0*y*y*Ui*Ui*Ri*Ri*x)/(pow(r2,3.)) + (2.0*x*Ui*Ui*Ri*Ri)/(pow(r2,2.)) - (4.0*x*x*x*Ui*Ui*Ri*Ri)/(pow(r2,3.)) );

      dEdy = ( 2.0 * P * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * y ) / ( (pow(gm1,2.)) * (pow(r2,2.)) * (pow(term,1.)) )
	   + ( 1.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * y * ( (y * y * Ui * Ui * Ri * Ri)/(pow(r2,2.)) + (x * x * Ui * Ui * Ri * Ri)/(pow(r2,2.)) ) ) / 
	     ( (pow(gm1,1.)) * (pow(r2,2.)) * (pow(term,1.)) )
	   + 0.5 * rho * ( (2.0*y*Ui*Ui*Ri*Ri)/(pow(r2,2.)) - (4.0*y*y*y*Ui*Ui*Ri*Ri)/(pow(r2,3.)) - (4.0*x*x*Ui*Ui*Ri*Ri*y)/(pow(r2,3.)) );

      dE2dxx = ( 4.0 * P * gamma * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x * x ) / ( (pow(gm1,3.)) * (pow(r2,4.)) * (pow(term,2.)) )
	     - ( 8.0 * P *         (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * x * x ) / ( (pow(gm1,2.)) * (pow(r2,3.)) * (pow(term,1.)) )
	     + ( 2.0 * P *         (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) )         / ( (pow(gm1,2.)) * (pow(r2,2.)) * (pow(term,1.)) )
	     - ( 4.0 * P *         (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x * x ) / ( (pow(gm1,2.)) * (pow(r2,4.)) * (pow(term,2.)) )
	     + ( 2.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x * x * term2 ) / ( (pow(gm1,2.)) * (pow(r2,4.)) * (pow(term,2.)) )
	     - ( 4.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * x * x * term2 ) / ( (pow(gm1,1.)) * (pow(r2,3.)) * (pow(term,1.)) )
	     + ( 1.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * term2 )         / ( (pow(gm1,1.)) * (pow(r2,2.)) * (pow(term,1.)) )
	     - ( 2.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x * x * term2 ) / ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,2.)) )
	     + ( 2.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * x * (
										      (-4.0*y*y*Ui*Ui*Ri*Ri*x)/(pow(r2,3.)) + (2.0*x*Ui*Ui*Ri*Ri)/(pow(r2,2.)) - (4.0*x*x*x*Ui*Ui*Ri*Ri)/(pow(r2,3.)) ) ) /
 	       ( (pow(gm1,1.)) * (pow(r2,2.)) * (pow(term,1.)) )
	     + 0.5 * rho * (  (24.0*y*y*Ui*Ui*Ri*Ri*x*x)/(pow(r2,4.)) - (4.0*y*y*Ui*Ui*Ri*Ri)/(pow(r2,3.)) + (2.0*Ui*Ui*Ri*Ri)/(pow(r2,2.)) - (20.0*x*x*Ui*Ui*Ri*Ri)/(pow(r2,3.)) + (24.0*x*x*x*x*Ui*Ui*Ri*Ri)/(pow(r2,4.))  );

      dE2dyy = ( 4.0 * P * gamma * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y * y ) / ( (pow(gm1,3.)) * (pow(r2,4.)) * (pow(term,2.)) )
	     - ( 8.0 * P *         (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * y * y ) / ( (pow(gm1,2.)) * (pow(r2,3.)) * (pow(term,1.)) )
	     + ( 2.0 * P *         (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) )         / ( (pow(gm1,2.)) * (pow(r2,2.)) * (pow(term,1.)) )
	     - ( 4.0 * P *         (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y * y ) / ( (pow(gm1,2.)) * (pow(r2,4.)) * (pow(term,2.)) )
	     + ( 2.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y * y * term2 ) / ( (pow(gm1,2.)) * (pow(r2,4.)) * (pow(term,2.)) )
	     - ( 4.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * y * y * term2 ) / ( (pow(gm1,1.)) * (pow(r2,3.)) * (pow(term,1.)) )
	     + ( 1.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * term2 )         / ( (pow(gm1,1.)) * (pow(r2,2.)) * (pow(term,1.)) )
	     - ( 2.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y * y * term2 ) / ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,2.)) )
	     + ( 2.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * y * (
										      (2.0*y*Ui*Ui*Ri*Ri)/(pow(r2,2.)) - (4.0*y*y*y*Ui*Ui*Ri*Ri)/(pow(r2,3.)) - (4.0*x*x*Ui*Ui*Ri*Ri*y)/(pow(r2,3.)) ) ) /
 	       ( (pow(gm1,1.)) * (pow(r2,2.)) * (pow(term,1.)) )
	     + 0.5 * rho * (  (2.0*Ui*Ui*Ri*Ri)/(pow(r2,2.)) - (20.0*y*y*Ui*Ui*Ri*Ri)/(pow(r2,3.)) + (24.0*y*y*y*y*Ui*Ui*Ri*Ri)/(pow(r2,4.)) + (24.0*y*y*Ui*Ui*Ri*Ri*x*x)/(pow(r2,4.)) - (4.0*x*x*Ui*Ui*Ri*Ri)/(pow(r2,3.))  );

      dE2dxy = ( 4.0 * P * gamma * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y * x ) / ( (pow(gm1,3.)) * (pow(r2,4.)) * (pow(term,2.)) )
	     - ( 8.0 * P *         (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * x * y ) / ( (pow(gm1,2.)) * (pow(r2,3.)) * (pow(term,1.)) )
	     - ( 4.0 * P *         (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x * y)  / ( (pow(gm1,2.)) * (pow(r2,4.)) * (pow(term,2.)) )
	     + ( 2.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y * x * term2 ) / ( (pow(gm1,2.)) * (pow(r2,4.)) * (pow(term,2.)) )
	     - ( 4.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * x * term2 * y ) / ( (pow(gm1,1.)) * (pow(r2,3.)) * (pow(term,1.)) )
	     - ( 2.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x * term2 * y ) / ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,2.)) )
	     + ( 1.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * x * (
										      ( 2.0*y*Ui*Ui*Ri*Ri)/(pow(r2,2.)) - (4.0*y*y*y*Ui*Ui*Ri*Ri)/(pow(r2,3.)) - (4.0*x*x*Ui*Ui*Ri*Ri*y)/(pow(r2,3.)) ) ) /
 	       ( (pow(gm1,1.)) * (pow(r2,2.)) * (pow(term,1.)) )
	     + ( 1.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * y * (
										      (-4.0*y*y*Ui*Ui*Ri*Ri*x)/(pow(r2,3.)) + (2.0*x*Ui*Ui*Ri*Ri)/(pow(r2,2.)) - (4.0*x*x*x*Ui*Ui*Ri*Ri)/(pow(r2,3.)) ) ) /
 	       ( (pow(gm1,1.)) * (pow(r2,2.)) * (pow(term,1.)) )
	     + 0.5 * rho * (  (-16.0*y*Ui*Ui*Ri*Ri*x)/(pow(r2,3.)) + (24.0*y*y*y*Ui*Ui*Ri*Ri*x)/(pow(r2,4.)) + (24.0*x*x*x*Ui*Ui*Ri*Ri*y)/(pow(r2,4.)) );

      dE3dxxx = (  12.0 * P * (pow(gamma,1.)) * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x )         / ( (pow(gm1,3.)) * (pow(r2,4.)) * (pow(term,2.)) )
	      + (   8.0 * P * (pow(gamma,2.)) * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * x * x * x ) / ( (pow(gm1,4.)) * (pow(r2,6.)) * (pow(term,3.)) )
	      - (  48.0 * P * (pow(gamma,1.)) * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x * x * x ) / ( (pow(gm1,3.)) * (pow(r2,5.)) * (pow(term,2.)) )
	      - (  24.0 * P * (pow(gamma,1.)) * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * x * x * x ) / ( (pow(gm1,3.)) * (pow(r2,6.)) * (pow(term,3.)) )
              - (  24.0 * P *                   (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * x )         / ( (pow(gm1,2.)) * (pow(r2,3.)) * (pow(term,1.)) )
              + (  48.0 * P *                   (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * x * x * x ) / ( (pow(gm1,2.)) * (pow(r2,4.)) * (pow(term,1.)) )
              + (  48.0 * P *                   (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x * x * x ) / ( (pow(gm1,2.)) * (pow(r2,5.)) * (pow(term,2.)) )
              - (  12.0 * P *                   (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x )         / ( (pow(gm1,2.)) * (pow(r2,4.)) * (pow(term,2.)) )
              + (  16.0 * P *                   (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * x * x * x ) / ( (pow(gm1,2.)) * (pow(r2,6.)) * (pow(term,3.)) )
	      + 0.5 * rho * (  (-192.0*y*y*Ui*Ui*Ri*Ri*x*x*x)/(pow(r2,5.)) + (72.0*y*y*Ui*Ui*Ri*Ri*x)/(pow(r2,4.)) - (48.0*Ui*Ui*Ri*Ri*x)/(pow(r2,3.))
			       + (216.0*x*x*x*Ui*Ui*Ri*Ri)/(pow(r2,4.)) - (192.0*x*x*x*x*x*Ui*Ui*Ri*Ri)/(pow(r2,5.)) )
              + ( 3.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * (
										   (-4.0*y*y*Ui*Ui*Ri*Ri*x)/(pow(r2,3.)) + (2.0*x*Ui*Ui*Ri*Ri)/(pow(r2,2.)) - (4.0*x*x*x*Ui*Ui*Ri*Ri)/(pow(r2,3.)) ) ) /
 	        ( (pow(gm1,1.)) * (pow(r2,2.)) * (pow(term,1.)) )
              + ( 6.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x * term2 ) / ( (pow(gm1,2.)) * (pow(r2,4.)) * (pow(term,2.)) )
	      + ( 6.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x * x * (
											   (-4.0*y*y*Ui*Ui*Ri*Ri*x)/(pow(r2,3.)) + (2.0*x*Ui*Ui*Ri*Ri)/(pow(r2,2.)) - (4.0*x*x*x*Ui*Ui*Ri*Ri)/(pow(r2,3.)) ) ) /
	        ( (pow(gm1,2.)) * (pow(r2,4.)) * (pow(term,2.)) )
              + ( 4.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * x * x * x * term2 ) / ( (pow(gm1,3.)) * (pow(r2,6.)) * (pow(term,3.)) )
	      - (24.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x * x * x * term2 ) / ( (pow(gm1,2.)) * (pow(r2,5.)) * (pow(term,2.)) )
	      - (12.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * x * x * x * term2 ) / ( (pow(gm1,2.)) * (pow(r2,6.)) * (pow(term,3.)) )
	      - (12.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * x * term2 )         / ( (pow(gm1,1.)) * (pow(r2,3.)) * (pow(term,1.)) )
              - (12.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * x * x * (
										           (-4.0*y*y*Ui*Ui*Ri*Ri*x)/(pow(r2,3.)) + (2.0*x*Ui*Ui*Ri*Ri)/(pow(r2,2.)) - (4.0*x*x*x*Ui*Ui*Ri*Ri)/(pow(r2,3.)) ) ) /
 	        ( (pow(gm1,1.)) * (pow(r2,3.)) * (pow(term,1.)) )
              + (24.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * x * x * x * term2 ) / ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,1.)) )
              + (24.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x * x * x * term2 ) / ( (pow(gm1,1.)) * (pow(r2,5.)) * (pow(term,2.)) )
	      - ( 6.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * term2 * x )         / ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,2.)) )
	      - ( 6.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x * x * (
											   (-4.0*y*y*Ui*Ui*Ri*Ri*x)/(pow(r2,3.)) + (2.0*x*Ui*Ui*Ri*Ri)/(pow(r2,2.)) - (4.0*x*x*x*Ui*Ui*Ri*Ri)/(pow(r2,3.)) ) ) /
	        ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,2.)) )
	      + ( 8.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * x * x * x * term2 ) / ( (pow(gm1,1.)) * (pow(r2,6.)) * (pow(term,3.)) )
	      + ( 3.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * x * (
										       (24.0*y*y*Ui*Ui*Ri*Ri*x*x)/(pow(r2,4.)) - (4.0*y*y*Ui*Ui*Ri*Ri)/(pow(r2,3.)) + (2.0*Ui*Ui*Ri*Ri)/(pow(r2,2.)) 
										       - (20.0*x*x*Ui*Ui*Ri*Ri)/(pow(r2,3.)) + (24.0*x*x*x*x*Ui*Ui*Ri*Ri)/(pow(r2,4.)) ) ) /
	        ( (pow(gm1,1.)) * (pow(r2,2.)) * (pow(term,1.)) );


      dE3dyyy = (  48.0 * P *                   (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * y * y * y ) / ( (pow(gm1,2.)) * (pow(r2,4.)) * (pow(term,1.)) )
	      + (  12.0 * P * (pow(gamma,1.)) * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y )         / ( (pow(gm1,3.)) * (pow(r2,4.)) * (pow(term,2.)) )
	      + (   8.0 * P * (pow(gamma,2.)) * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * y * y * y ) / ( (pow(gm1,4.)) * (pow(r2,6.)) * (pow(term,3.)) )
	      - (  48.0 * P * (pow(gamma,1.)) * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y * y * y ) / ( (pow(gm1,3.)) * (pow(r2,5.)) * (pow(term,2.)) )
              - (  24.0 * P * (pow(gamma,1.)) * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * y * y * y ) / ( (pow(gm1,3.)) * (pow(r2,6.)) * (pow(term,3.)) )
              - (  24.0 * P *                   (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * y )         / ( (pow(gm1,2.)) * (pow(r2,3.)) * (pow(term,1.)) )
              + (  48.0 * P *                   (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y * y * y ) / ( (pow(gm1,2.)) * (pow(r2,5.)) * (pow(term,2.)) )
              - (  12.0 * P *                   (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y )         / ( (pow(gm1,2.)) * (pow(r2,4.)) * (pow(term,2.)) )
              + (  16.0 * P *                   (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * y * y * y ) / ( (pow(gm1,2.)) * (pow(r2,6.)) * (pow(term,3.)) )
	      + 0.5 * rho * ( (-48.0*Ui*Ui*Ri*Ri*y)/(pow(r2,3.)) + (216.0*y*y*y*Ui*Ui*Ri*Ri)/(pow(r2,4.)) - (192.0*y*y*y*y*y*Ui*Ui*Ri*Ri)/(pow(r2,5.))
			      + (72.0*y*Ui*Ui*Ri*Ri*x*x)/(pow(r2,4.)) - (192.0*y*y*y*Ui*Ui*Ri*Ri*x*x)/(pow(r2,5.)) )
	      + ( 3.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * (
										   (2.0*y*Ui*Ui*Ri*Ri)/(pow(r2,2.)) - (4.0*y*y*y*Ui*Ui*Ri*Ri)/(pow(r2,3.)) - (4.0*x*x*Ui*Ui*Ri*Ri*y)/(pow(r2,3.)) ) ) /
 	        ( (pow(gm1,1.)) * (pow(r2,2.)) * (pow(term,1.)) )
              + ( 6.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y * term2 ) / ( (pow(gm1,2.)) * (pow(r2,4.)) * (pow(term,2.)) )
	      + ( 6.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y * y * (
											   (2.0*y*Ui*Ui*Ri*Ri)/(pow(r2,2.)) - (4.0*y*y*y*Ui*Ui*Ri*Ri)/(pow(r2,3.)) - (4.0*x*x*Ui*Ui*Ri*Ri*y)/(pow(r2,3.)) ) ) /
	        ( (pow(gm1,2.)) * (pow(r2,4.)) * (pow(term,2.)) )
              + ( 4.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * y * y * y * term2 ) / ( (pow(gm1,3.)) * (pow(r2,6.)) * (pow(term,3.)) )
	      - (24.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y * y * y * term2 ) / ( (pow(gm1,2.)) * (pow(r2,5.)) * (pow(term,2.)) )
	      - (12.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * y * y * y * term2 ) / ( (pow(gm1,2.)) * (pow(r2,6.)) * (pow(term,3.)) )
	      - (12.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * y * term2 )         / ( (pow(gm1,1.)) * (pow(r2,3.)) * (pow(term,1.)) )
              - (12.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * y * y * (
											   (2.0*y*Ui*Ui*Ri*Ri)/(pow(r2,2.)) - (4.0*y*y*y*Ui*Ui*Ri*Ri)/(pow(r2,3.)) - (4.0*x*x*Ui*Ui*Ri*Ri*y)/(pow(r2,3.)) ) ) /
 	        ( (pow(gm1,1.)) * (pow(r2,3.)) * (pow(term,1.)) )
              + (24.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * y * y * y * term2 ) / ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,1.)) )
              + (24.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y * y * y * term2 ) / ( (pow(gm1,1.)) * (pow(r2,5.)) * (pow(term,2.)) )
	      - ( 6.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * term2 * y )         / ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,2.)) )
	      - ( 6.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y * y * (
											   (2.0*y*Ui*Ui*Ri*Ri)/(pow(r2,2.)) - (4.0*y*y*y*Ui*Ui*Ri*Ri)/(pow(r2,3.)) - (4.0*x*x*Ui*Ui*Ri*Ri*y)/(pow(r2,3.)) ) ) /
	        ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,2.)) )
	      + ( 8.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * y * y * y * term2 ) / ( (pow(gm1,1.)) * (pow(r2,6.)) * (pow(term,3.)) )
	      + ( 3.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * y * (
										       (2.0*Ui*Ui*Ri*Ri)/(pow(r2,2.)) - (20.0*y*y*Ui*Ui*Ri*Ri)/(pow(r2,3.)) + (24.0*y*y*y*y*Ui*Ui*Ri*Ri)/(pow(r2,4.))
										       + (24.0*y*y*Ui*Ui*Ri*Ri*x*x)/(pow(r2,4.)) - (4.0*x*x*Ui*Ui*Ri*Ri)/(pow(r2,3.)) ) ) /
	        ( (pow(gm1,1.)) * (pow(r2,2.)) * (pow(term,1.)) );

      dE3dxxy = (   4.0 * P * (pow(gamma,1.)) * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y )         / ( (pow(gm1,3.)) * (pow(r2,4.)) * (pow(term,2.)) )
	      - (   8.0 * P *                   (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * y )         / ( (pow(gm1,2.)) * (pow(r2,3.)) * (pow(term,1.)) )
	      - (   4.0 * P *                   (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y )         / ( (pow(gm1,2.)) * (pow(r2,4.)) * (pow(term,2.)) )
	      + (   8.0 * P * (pow(gamma,2.)) * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * y * x * x ) / ( (pow(gm1,4.)) * (pow(r2,6.)) * (pow(term,3.)) )
              - (  48.0 * P * (pow(gamma,1.)) * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x * x * y ) / ( (pow(gm1,3.)) * (pow(r2,5.)) * (pow(term,2.)) )
              - (  24.0 * P * (pow(gamma,1.)) * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * x * x * y ) / ( (pow(gm1,3.)) * (pow(r2,6.)) * (pow(term,3.)) )
              + (  48.0 * P *                   (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * x * x * y ) / ( (pow(gm1,2.)) * (pow(r2,4.)) * (pow(term,1.)) )
              + (  48.0 * P *                   (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x * x * y ) / ( (pow(gm1,2.)) * (pow(r2,5.)) * (pow(term,2.)) )
              + (  16.0 * P *                   (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * x * x * y ) / ( (pow(gm1,2.)) * (pow(r2,6.)) * (pow(term,3.)) )
	      + 0.5 * rho * ( (168.0*y*Ui*Ui*Ri*Ri*x*x)/(pow(r2,4.)) - (192.0*y*y*y*Ui*Ui*Ri*Ri*x*x)/(pow(r2,5.)) - (16.0*Ui*Ui*Ri*Ri*y)/(pow(r2,3.))
			      + (24.0*y*y*y*Ui*Ui*Ri*Ri)/(pow(r2,4.)) - (192.0*x*x*x*x*Ui*Ui*Ri*Ri*y)/(pow(r2,5.)) )
	      + ( 1.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * (
										   (2.0*y*Ui*Ui*Ri*Ri)/(pow(r2,2.)) - (4.0*y*y*y*Ui*Ui*Ri*Ri)/(pow(r2,3.)) - (4.0*x*x*Ui*Ui*Ri*Ri*y)/(pow(r2,3.)) ) ) /
 	        ( (pow(gm1,1.)) * (pow(r2,2.)) * (pow(term,1.)) )
              + ( 2.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y * term2 ) / ( (pow(gm1,2.)) * (pow(r2,4.)) * (pow(term,2.)) )
	      - ( 4.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * y * term2 ) / ( (pow(gm1,1.)) * (pow(r2,3.)) * (pow(term,1.)) )
              - ( 2.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * term2 * y ) / ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,2.)) )
	      + ( 2.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x * x * (
											   (2.0*y*Ui*Ui*Ri*Ri)/(pow(r2,2.)) - (4.0*y*y*y*Ui*Ui*Ri*Ri)/(pow(r2,3.)) - (4.0*x*x*Ui*Ui*Ri*Ri*y)/(pow(r2,3.)) ) ) /
 	        ( (pow(gm1,2.)) * (pow(r2,4.)) * (pow(term,2.)) )
              + ( 4.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * y * x * x * term2 ) / ( (pow(gm1,3.)) * (pow(r2,6.)) * (pow(term,3.)) )
              - (24.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x * x * term2 * y ) / ( (pow(gm1,2.)) * (pow(r2,5.)) * (pow(term,2.)) )
              - (12.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * x * x * term2 * y ) / ( (pow(gm1,2.)) * (pow(r2,6.)) * (pow(term,3.)) )
              - ( 4.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * x * x * (
											   (2.0*y*Ui*Ui*Ri*Ri)/(pow(r2,2.)) - (4.0*y*y*y*Ui*Ui*Ri*Ri)/(pow(r2,3.)) - (4.0*x*x*Ui*Ui*Ri*Ri*y)/(pow(r2,3.)) ) ) /
 	        ( (pow(gm1,1.)) * (pow(r2,3.)) * (pow(term,1.)) )
	      + (24.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * x * x * term2 * y ) / ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,1.)) )
              + (24.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x * x * term2 * y ) / ( (pow(gm1,1.)) * (pow(r2,5.)) * (pow(term,2.)) )
              - ( 2.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x * x * (
											   (2.0*y*Ui*Ui*Ri*Ri)/(pow(r2,2.)) - (4.0*y*y*y*Ui*Ui*Ri*Ri)/(pow(r2,3.)) - (4.0*x*x*Ui*Ui*Ri*Ri*y)/(pow(r2,3.)) ) ) /
 	        ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,2.)) )
              + ( 8.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * x * x * term2 * y ) / ( (pow(gm1,1.)) * (pow(r2,6.)) * (pow(term,3.)) )
              + ( 2.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * x * (
										       (-16.0*y*Ui*Ui*Ri*Ri*x)/(pow(r2,3.)) + (24.0*y*y*y*Ui*Ui*Ri*Ri*x)/(pow(r2,4.)) + (24.0*x*x*x*Ui*Ui*Ri*Ri*y)/(pow(r2,4.)) ) ) /
 	        ( (pow(gm1,1.)) * (pow(r2,2.)) * (pow(term,1.)) )
	      + ( 4.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y * x * (
											   (-4.0*y*y*Ui*Ui*Ri*Ri*x)/(pow(r2,3.)) + (2.0*x*Ui*Ui*Ri*Ri)/(pow(r2,2.)) - (4.0*x*x*x*Ui*Ui*Ri*Ri)/(pow(r2,3.)) ) ) /
 	        ( (pow(gm1,2.)) * (pow(r2,4.)) * (pow(term,2.)) )
	      - ( 8.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * x * (
										       (-4.0*y*y*Ui*Ui*Ri*Ri*x)/(pow(r2,3.)) + (2.0*x*Ui*Ui*Ri*Ri)/(pow(r2,2.)) - (4.0*x*x*x*Ui*Ui*Ri*Ri)/(pow(r2,3.)) ) * y ) /
 	        ( (pow(gm1,1.)) * (pow(r2,3.)) * (pow(term,1.)) )
	      - ( 4.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x * (
										       (-4.0*y*y*Ui*Ui*Ri*Ri*x)/(pow(r2,3.)) + (2.0*x*Ui*Ui*Ri*Ri)/(pow(r2,2.)) - (4.0*x*x*x*Ui*Ui*Ri*Ri)/(pow(r2,3.)) ) * y ) /
 	        ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,2.)) )
	      + ( 1.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * y * (
										       (24.0*y*y*Ui*Ui*Ri*Ri*x*x)/(pow(r2,4.)) - (4.0*y*y*Ui*Ui*Ri*Ri)/(pow(r2,3.)) + (2.0*Ui*Ui*Ri*Ri)/(pow(r2,2.))
										       - (20.0*x*x*Ui*Ui*Ri*Ri)/(pow(r2,3.)) + (24.0*x*x*x*x*Ui*Ui*Ri*Ri)/(pow(r2,4.)) ) ) /
	        ( (pow(gm1,1.)) * (pow(r2,2.)) * (pow(term,1.)) );

      dE3dxyy = (   4.0 * P * (pow(gamma,1.)) * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x )         / ( (pow(gm1,3.)) * (pow(r2,4.)) * (pow(term,2.)) )
	      - (   8.0 * P *                   (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * x )         / ( (pow(gm1,2.)) * (pow(r2,3.)) * (pow(term,1.)) )
	      - (   4.0 * P *                   (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x )         / ( (pow(gm1,2.)) * (pow(r2,4.)) * (pow(term,2.)) )
	      + (   8.0 * P * (pow(gamma,2.)) * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * y * y * x ) / ( (pow(gm1,4.)) * (pow(r2,6.)) * (pow(term,3.)) )
              - (  48.0 * P * (pow(gamma,1.)) * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y * y * x ) / ( (pow(gm1,3.)) * (pow(r2,5.)) * (pow(term,2.)) )
              - (  24.0 * P * (pow(gamma,1.)) * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * y * y * x ) / ( (pow(gm1,3.)) * (pow(r2,6.)) * (pow(term,3.)) )
              + (  48.0 * P *                   (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * x * y * y ) / ( (pow(gm1,2.)) * (pow(r2,4.)) * (pow(term,1.)) )
              + (  48.0 * P *                   (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x * y * y ) / ( (pow(gm1,2.)) * (pow(r2,5.)) * (pow(term,2.)) )
              + (  16.0 * P *                   (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * x * y * y ) / ( (pow(gm1,2.)) * (pow(r2,6.)) * (pow(term,3.)) )
	      + ( 1.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * (
										   (-4.0*y*y*Ui*Ui*Ri*Ri*x)/(pow(r2,3.)) + (2.0*x*Ui*Ui*Ri*Ri)/(pow(r2,2.)) - (4.0*x*x*x*Ui*Ui*Ri*Ri)/(pow(r2,3.)) ) ) /
 	        ( (pow(gm1,1.)) * (pow(r2,2.)) * (pow(term,1.)) )
              + ( 2.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x * term2 ) / ( (pow(gm1,2.)) * (pow(r2,4.)) * (pow(term,2.)) )
              - ( 4.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * x * term2 ) / ( (pow(gm1,1.)) * (pow(r2,3.)) * (pow(term,1.)) )
              - ( 2.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * term2 * x ) / ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,2.)) )
              + 0.5 * rho * ( (-16.0*Ui*Ui*Ri*Ri*x)/(pow(r2,3.)) + (168.0*y*y*Ui*Ui*Ri*Ri*x)/(pow(r2,4.)) - (192.0*y*y*y*y*Ui*Ui*Ri*Ri*x)/(pow(r2,5.))
			      - (192.0*y*y*Ui*Ui*Ri*Ri*x*x*x)/(pow(r2,5.)) + (24.0*x*x*x*Ui*Ui*Ri*Ri)/(pow(r2,4.)) )
              + ( 4.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y * x * (
											   (2.0*y*Ui*Ui*Ri*Ri)/(pow(r2,2.)) - (4.0*y*y*y*Ui*Ui*Ri*Ri)/(pow(r2,3.)) - (4.0*x*x*Ui*Ui*Ri*Ri*y)/(pow(r2,3.)) ) ) /
 	        ( (pow(gm1,2.)) * (pow(r2,4.)) * (pow(term,2.)) )
              + ( 4.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * y * y * x * term2 ) / ( (pow(gm1,3.)) * (pow(r2,6.)) * (pow(term,3.)) )
	      - (24.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y * y * x * term2 ) / ( (pow(gm1,2.)) * (pow(r2,5.)) * (pow(term,2.)) )
              - (12.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * y * y * x * term2 ) / ( (pow(gm1,2.)) * (pow(r2,6.)) * (pow(term,3.)) )
              - ( 8.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * x * (
										       (2.0*y*Ui*Ui*Ri*Ri)/(pow(r2,2.)) - (4.0*y*y*y*Ui*Ui*Ri*Ri)/(pow(r2,3.)) - (4.0*x*x*Ui*Ui*Ri*Ri*y)/(pow(r2,3.)) ) * y ) /
 	        ( (pow(gm1,1.)) * (pow(r2,3.)) * (pow(term,1.)) )
              + (24.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * x * term2 * y * y ) / ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,1.)) )
              + (24.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x * term2 * y * y ) / ( (pow(gm1,1.)) * (pow(r2,5.)) * (pow(term,2.)) )
              - ( 4.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * x * (
										       (2.0*y*Ui*Ui*Ri*Ri)/(pow(r2,2.)) - (4.0*y*y*y*Ui*Ui*Ri*Ri)/(pow(r2,3.)) - (4.0*x*x*Ui*Ui*Ri*Ri*y)/(pow(r2,3.)) ) * y ) /
 	        ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,2.)) )
              + ( 8.0 * rho * (pow((gm1*0.5),3.)) * (pow(Mi,6.)) * (pow(Ri,6.)) * x * term2 * y * y ) / ( (pow(gm1,1.)) * (pow(r2,6.)) * (pow(term,3.)) )
              + ( 1.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * x * (
										       (2.0*Ui*Ui*Ri*Ri)/(pow(r2,2.)) - (20.0*y*y*Ui*Ui*Ri*Ri)/(pow(r2,3.)) + (24.0*y*y*y*y*Ui*Ui*Ri*Ri)/(pow(r2,4.))
										       + (24.0*y*y*Ui*Ui*Ri*Ri*x*x)/(pow(r2,4.)) - (4.0*x*x*Ui*Ui*Ri*Ri)/(pow(r2,3.)) ) ) /
	        ( (pow(gm1,1.)) * (pow(r2,2.)) * (pow(term,1.)) )
              + ( 2.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * y * (
										       (-16.0*y*Ui*Ui*Ri*Ri*x)/(pow(r2,3.)) + (24.0*y*y*y*Ui*Ui*Ri*Ri*x)/(pow(r2,4.)) + (24.0*x*x*x*Ui*Ui*Ri*Ri*y)/(pow(r2,4.)) )  ) /
 	        ( (pow(gm1,1.)) * (pow(r2,2.)) * (pow(term,1.)) )
              + ( 2.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y * y * (
											   (-4.0*y*y*Ui*Ui*Ri*Ri*x)/(pow(r2,3.)) + (2.0*x*Ui*Ui*Ri*Ri)/(pow(r2,2.)) - (4.0*x*x*x*Ui*Ui*Ri*Ri)/(pow(r2,3.)) ) ) /
 	        ( (pow(gm1,2.)) * (pow(r2,4.)) * (pow(term,2.)) )
              - ( 4.0 * rho * (pow((gm1*0.5),1.)) * (pow(Mi,2.)) * (pow(Ri,2.)) * y * y * (
										       (-4.0*y*y*Ui*Ui*Ri*Ri*x)/(pow(r2,3.)) + (2.0*x*Ui*Ui*Ri*Ri)/(pow(r2,2.)) - (4.0*x*x*x*Ui*Ui*Ri*Ri)/(pow(r2,3.)) ) ) /
 	        ( (pow(gm1,1.)) * (pow(r2,3.)) * (pow(term,1.)) )
              - ( 2.0 * rho * (pow((gm1*0.5),2.)) * (pow(Mi,4.)) * (pow(Ri,4.)) * y * y * (
											   (-4.0*y*y*Ui*Ui*Ri*Ri*x)/(pow(r2,3.)) + (2.0*x*Ui*Ui*Ri*Ri)/(pow(r2,2.)) - (4.0*x*x*x*Ui*Ui*Ri*Ri)/(pow(r2,3.)) ) ) /
	        ( (pow(gm1,1.)) * (pow(r2,4.)) * (pow(term,2.)) );


      // set the derivatives.
      derivs[0] = drdx;
      derivs[1] = drdy;
      derivs[2] = dr2dxx;
      derivs[3] = dr2dyy;
      derivs[4] = dr2dxy;
      derivs[5] = dr3dxxx;
      derivs[6] = dr3dyyy;
      derivs[7] = dr3dxxy;
      derivs[8] = dr3dxyy;

      derivs[9]  = drudx;
      derivs[10] = drudy;
      derivs[11] = dru2dxx;
      derivs[12] = dru2dyy;
      derivs[13] = dru2dxy;
      derivs[14] = dru3dxxx;
      derivs[15] = dru3dyyy;
      derivs[16] = dru3dxxy;
      derivs[17] = dru3dxyy;

      derivs[18] = drvdx;
      derivs[19] = drvdy;
      derivs[20] = drv2dxx;
      derivs[21] = drv2dyy;
      derivs[22] = drv2dxy;
      derivs[23] = drv3dxxx;
      derivs[24] = drv3dyyy;
      derivs[25] = drv3dxxy;
      derivs[26] = drv3dxyy;

      derivs[27] = dEdx;
      derivs[28] = dEdy;
      derivs[29] = dE2dxx;
      derivs[30] = dE2dyy;
      derivs[31] = dE2dxy;
      derivs[32] = dE3dxxx;
      derivs[33] = dE3dyyy;
      derivs[34] = dE3dxxy;
      derivs[35] = dE3dxyy;
    }
  

  return;
}

//=============================================================
// 
//  Compute_CV_Averages()
//
//  Computes the CV averaged value of the test function for each
//  CV in the grid. This is for helping to debug the Hessian
//  and to get the order of accuracy of the reconstruction.
//
//  GRID *grid;                       // The grid.
//  double *cv_avg;                   // The averaged value.
//
//=============================================================

void Compute_CV_Averages ( GRID *grid, double *cv_avg )
{
  // This procedure mimics the steps used in getting the Moment terms.

  int i, j, k, b, s, n;                  // Loop counters.
  double xc, yc;                         // Values at the centroid.
  double xmid, ymid;                     // Values at face midpoints.
  double len;                            // Length of an edge piece.
  double nx;                             // Normal vector component.
  double ny;
  double xL,xR,yL,yR;
  double t1,t2,t3;                       // Gaussian roots.
  double w1,w2,w3;                       // Gaussian coefficients.
  double xi, yi;                         // Node coordinates.
  double temp;                           // Temporary storage.
  int iedge;                             // Edge ID.
  int num_vert;                          // Number of vertices for an element.
  int nodeL, nodeR;                      // Nodes comprising an edge.
  int node;                              // Generic node.
  int gelem;                             // Global element number.
  int bc;                                // Boundary condition.
  int ghostnode;
  int ghost_node;
  int seg;

  // Curved boundary stuff.
  double XL[2],XR[2],GP[6];
  double w[3];

  // Volume integral variables.
  int e,v,ind;
  int nodes[MAX_NUM_VERT];               // Element vertices.
  int left_id;                           // ID of left edge.
  int right_id;                          // ID of right edge.
  int ct_flag;
  double XC[NDIM];                       // Vectors.
  double xmidL[NDIM];
  double xmidR[NDIM];
  double v_left[NDIM];
  double v_right[NDIM];
  double v_middle[NDIM];
  double bary_coord[7][3];               // The barycentric coordinates.
  double wg[7];                          // The associated weights;
  double b_r = (6. - sqrt(15.))/21.;     // Precompute some the coordinates.
  double b_t = (6. + sqrt(15.))/21.;
  double b_s = (9. + 2.*sqrt(15.))/21.;
  double b_u = (9. - 2.*sqrt(15.))/21.;
  double wA = (155. - sqrt(15.))/1200.;
  double wB = (155. + sqrt(15.))/1200.;
  double ct_gp[7][2];                    // Gauss integration points in (r,s) space (xi,eta).
  double ct_nodes[7][2];                 // The x,y locations of the 6 nodes defining the p2 triangle (curved side).
  double N1,N2,N3,N4,N5,N6;              // Shape functions.
  double N1r,N1s,N2r,N2s,N3r,N3s;        // Shape function derivatives.
  double N4r,N4s,N5r,N5s,N6r,N6s;
  double xr,xs,yr,ys,jacobian;           // Transformation jacobian.
  double XGP[NDIM];
  double QTRI[NUM_VAR];
  double XI,ETA;                         // Triangle coordinates.
  double tri_int;                        // Triangle integration.
  double gp[7][3];                       // Gauss points for integration.
  double dA;                             // Differential area (for a triangle).

  w[0] = 8./9.;
  w[1] = 5./9.;
  w[2] = w[1];

  // Pointers to the Gauss points.
  double *x1 = NULL;
  double *x2 = NULL;
  double *x3 = NULL;

  // Allocate space.
  x1 = (double*)malloc(NDIM*sizeof(double));
  if ( x1 == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'x1' (HESSIAN).\n"); exit(1); }

  x2 = (double*)malloc(NDIM*sizeof(double));
  if ( x2 == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'x2' (HESSIAN).\n"); exit(1); }

  x3 = (double*)malloc(NDIM*sizeof(double));
  if ( x3 == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'x3' (HESSIAN).\n"); exit(1); }

  // Make sure it is clean since length contribtutions are accumulated.
  for ( n=1; n <= grid->nn; n++ )
    {
      cv_avg[n] = 0.;
    }

  // The process is to:
  // 1. Loop over all the subedges which each represent a dual edge.
  // 2. Reconstruct the dual edge from the available information.
  // 3. Apply Gaussian quadrature to calculate the test function for each node (left and right).

  if ( 0 ) // divergence theorem, flux dot n intgral through area
    {
      for ( i=1; i <= grid->nsubedges; i++ )
	{
	  // Find the global edge associated with the subedge.
	  iedge = grid->subedges[i*2];
	  
	  // Find the global element.
	  gelem = grid->subedges[i*2+1];
      
	  // Now get the nodes attached to the edge.
	  nodeL = grid->edges[2*iedge+0];
	  nodeR = grid->edges[2*iedge+1];

	  // Get the node coordinates.
	  xi = grid->x[nodeL*NDIM+0];
	  yi = grid->x[nodeL*NDIM+1];
      
	  // Get the element centroid.
	  xc = grid->el_cent[gelem*2];
	  yc = grid->el_cent[gelem*2+1];

	  // Get the length.
	  len = grid->xn_subedges[i*3+2];

	  // Retrieve the x component of the dual edge normal vector.
	  nx = grid->xn_subedges[i*3+0];
	  ny = grid->xn_subedges[i*3+1];

	  xmid = grid->xm_subedges[i*NDIM+0];
	  ymid = grid->xm_subedges[i*NDIM+1];

	  // For Gaussian quadrature I need to map the points in t-space
	  // onto my dual edge.

	  // I'm also pretty sure that the to be consistent with the right hand
	  // rule the dual edge is defined to be from the edge midpoint to the
	  // element centroid so that the subedge normal vector points away from
	  // the left node.

	  // To get the x,y points from t space use - x = (xL+xR)/2 + (xR-xL)/2*t
	  //                                          y = (yL+yR)/2 + (yR-yL)/2*t
	  //                                          xL,yL = edge midpoint
	  //                                          xR,yR = element centroid

	  t1 = sqrt(15.0)/5.0;
	  w1 = 5./9.;

	  t2 = 0.;
	  w2 = 8./9.;

	  t3 = sqrt(15.0)/(-5.0);
	  w3 = w1;

	  x1[0] = (xmid+xc)*0.5 + (xc-xmid)*0.5*t1;
	  x1[1] = (ymid+yc)*0.5 + (yc-ymid)*0.5*t1;
      
	  x2[0] = (xmid+xc)*0.5 + (xc-xmid)*0.5*t2;
	  x2[1] = (ymid+yc)*0.5 + (yc-ymid)*0.5*t2;
      
	  x3[0] = (xmid+xc)*0.5 + (xc-xmid)*0.5*t3;
	  x3[1] = (ymid+yc)*0.5 + (yc-ymid)*0.5*t3;

	  // Now I can apply the function to the Gauss points.
	  temp = (w1*Test_Function_Divergent(x1)) +
	         (w2*Test_Function_Divergent(x2)) +
	         (w3*Test_Function_Divergent(x3));

	  temp = temp * nx * len/2.;
	  
	  // Accumulate to the node.
	  cv_avg[nodeL] += temp;
	  cv_avg[nodeR] -= temp;
	}

      /*
      printf("CV Averages interior edges only:\n");
      for ( i=1; i <= grid->nn; i++ )
	{
	  printf(" %d = %.15e\n",i,cv_avg[i]);
	}
      printf("\n\n");
      */

      // Now I need to loop over and close off the boundaries.
      
      for ( i=1; i <= grid->nbedges; i++ )
	{
	  // Get the node.
	  node = grid->bedges[i*5+0];
	  ghostnode = grid->bedges[i*5+1] + grid->nn;

	  // Get the boundary.
	  b = grid->bedges[i*5+3];

	  bc = grid->bbc[b];

	  // Get the segment.
	  s = grid->bedges[i*5+4];

	  // Now get the nodes attached to the edge.
	  nodeL = grid->bs[b][s][0];
	  nodeR = grid->bs[b][s][1];

	  // Get the node coordinates.
	  xi = grid->x[node*NDIM+0];
	  yi = grid->x[node*NDIM+1];
	  
	  if ( bc <= 10 )
	    {
	      // Get the edge midpoint.
	      xmid = 0.5*( grid->x[nodeL*NDIM+0] + grid->x[nodeR*NDIM+0] );
	      ymid = 0.5*( grid->x[nodeL*NDIM+1] + grid->x[nodeR*NDIM+1] );
      
	      // Get the length.
	      len = grid->xn_bedges[i*3+2];

	      // Retrieve the x component of the dual edge normal vector.
	      nx = grid->xn_bedges[i*3+0];
	      
	      // For Gaussian quadrature I need to map the points in t-space
	      // onto my dual edge.

	      // I'm also pretty sure that the to be consistent with the right hand
	      // rule the boundary edge is defined from node0 to node1 on the segment.
	      // The normal should thus point out of the mesh.
	  
	      // To get the x,y points from t space use - x = (xL+xR)/2 + (xR-xL)/2*t
	      //                                          y = (yL+yR)/2 + (yR-yL)/2*t
	      // The left and right nodes depend on which node the node is. Is it on the
	      // beginning or end of the edge.

	      if ( node == nodeL )
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
	      w1 = 5./9.;

	      t2 = 0.;
	      w2 = 8./9.;

	      t3 = sqrt(15.0)/(-5.0);
	      w3 = w1;

	      x1[0] = (xL+xR)*0.5 + (xR-xL)*0.5*t1;
	      x1[1] = (yL+yR)*0.5 + (yR-yL)*0.5*t1;
      
	      x2[0] = (xL+xR)*0.5 + (xR-xL)*0.5*t2;
	      x2[1] = (yL+yR)*0.5 + (yR-yL)*0.5*t2;
	  
	      x3[0] = (xL+xR)*0.5 + (xR-xL)*0.5*t3;
	      x3[1] = (yL+yR)*0.5 + (yR-yL)*0.5*t3; 
	  
	      // Now I can apply the function to the Gauss points.
	      temp = (w1*Test_Function_Divergent(x1)) +
		     (w2*Test_Function_Divergent(x2)) +
		     (w3*Test_Function_Divergent(x3));
	  
	      temp = temp * nx * len/2.;
	      
	      // Accumulate to the node.
	      cv_avg[node] += temp;

	      //printf("bedge %d is contributing < %.15e > to node %d\n",i,temp,node);
	    }
	  else
	    {
	      if ( node == nodeL )
		{
		  XL[0] = grid->x[node*NDIM+0];
		  XL[1] = grid->x[node*NDIM+1];
		  XR[0] = grid->x[ghostnode*NDIM+0];
		  XR[1] = grid->x[ghostnode*NDIM+1];
		}
	      else
		{
		  XR[0] = grid->x[node*NDIM+0];
		  XR[1] = grid->x[node*NDIM+1];
		  XL[0] = grid->x[ghostnode*NDIM+0];
		  XL[1] = grid->x[ghostnode*NDIM+1];
		}
	      
	      // Get the Gauss points.
	      curved_boundary_gauss_points ( grid, bc, XL, XR, GP );
	      
	      curved_boundary_arclength ( grid, bc, XL, XR, &len );

	      temp = 0.;
	      
	      for ( j=0; j < 3; j++ )
		{
		  // Get the normal vector.
		  curved_boundary_normal_vector ( grid, bc, &(GP[j*NDIM]), &nx, &ny );
		  
		  // Moved length inside this function. The integral of the test function
		  // acts through a fractional component of the total length.
		  temp += ( w[j] * Test_Function_Divergent( &(GP[j*NDIM]) ) * nx * len);
		}

	      temp = temp * 0.5;
	      
	      cv_avg[node] += temp;
	      //printf("bedge %d is contributing < %.15e > to node %d\n",i,temp,node);
	    }
	}
      
      /*
      printf("CV Averages all edges:\n");
      for ( i=1; i <= grid->nn; i++ )
	{
	  printf(" %d = %.15e\n",i,cv_avg[i]);
	}
      printf("\n\n");
      */
    }
  else // Do a volume integral. Useful if the function is nasty.
    {
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

      // Start the process by looping over all the elements and looping through the internal median dual pieces.
      for ( e=Tri; e <= Quad; e++ )                 // Element types loop.
	{
	  // Find the number of nodes and edges for this element type.
	  num_vert = NumberNodesForElement(e);
	  
	  //printf("  Processing element type %d\n",e);
	  
	  for ( j=1; j <= grid->num_elem[e]; j++ )              // Element loop.
	    {
	      //printf("    Processing element number %d\n",j);
	      
	      // Retrieve the element node ids. This is for convenience.
	      for ( k=0; k < num_vert; k++ )              // Element node extraction loop.
		{
		  nodes[k] = grid->c2n[e][j*num_vert+k];
		}                                         // End element node extraction loop.
	      
	      // Find the element centroid.
	      XC[0] = 0.; XC[1] = 0.;
	      for ( k=0; k < num_vert; k++ )
		{
		  for ( i=0; i < NDIM; i++ )
		    {
		      XC[i] += grid->x[NDIM*nodes[k]+i];
		    }
		}
	      XC[0] /= ((double)num_vert);
	      XC[1] /= ((double)num_vert);

	      // Loop through each node in the element.
	      for ( k=0; k < num_vert; k++ )              // Element node loop.
		{
		  //printf("      Processing node %d\n",k);
		  
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
		      v_middle[i] = XC[i]    - grid->x[NDIM*node+i];
		    }

		  // Process the first triangle.
		  //printf("        Processing triangle 1.\n");

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
		      ct_nodes[1][0] = XC[0];
		      ct_nodes[1][1] = XC[1];

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
			  printf("CRITICAL ERROR: Triangle 1 of element %d of type %d is too concave (compute cv averages: hessian.C)!\n",j,e);
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

			  temp = Test_Function(gp[0]);

			  // Now we accumulate to the triangle.
			  tri_int += ( wg[i] * temp * jacobian * 0.5 );
			}

		      // Add to the running total for the node.
		      cv_avg[node] += tri_int;
      
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
			      gp[v][ind] = (bary_coord[v][0] * grid->x[node*NDIM+ind] +
					    bary_coord[v][1] * xmidL[ind] +
					    bary_coord[v][2] * XC[ind]);
			    }
			  gp[v][2] = 0.;   // Always 0 since 2D vectors.
			}
		      
		      // Now we loop over the Gauss points and get the function value
		      // and perform the numerical integration (the area integral of the triangle).
		      
		      QTRI[0] = 0.;

		      for ( i=0; i < 7; i++ )
			{
			  // Get the value at the point.
			  XGP[0] = gp[i][0];
			  XGP[1] = gp[i][1];
			  
			  temp = Test_Function ( XGP );
			  
			  QTRI[0] += ( wg[i] * temp );
			}
		      
		      QTRI[0] *= dA;
		      
		      // Accumulate the integral to the node.
		      cv_avg[node] += QTRI[0];
		    }

		  // Process the second triangle.
		  
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
		      ct_nodes[1][0] = XC[0];
		      ct_nodes[1][1] = XC[1];
		      
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
			  printf("CRITICAL ERROR: Triangle 2 of element %d of type %d is too concave (compute cv averages: hessian.C)!\n",j,e);
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

			  temp = Test_Function(gp[0]);
			  
			  // Now we accumulate to the triangle.
			  tri_int += ( wg[i] * temp * jacobian * 0.5 );
			}

		      // Add to the running total for the node.
		      cv_avg[node] += tri_int;
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
			      gp[v][ind] = (bary_coord[v][0] * grid->x[node*NDIM+ind] +
					    bary_coord[v][1] * xmidR[ind] +
					    bary_coord[v][2] * XC[ind]);
			    }
			  gp[v][2] = 0.;   // Always 0 since 2D vectors.
			}
		      
		      // Now we loop over the Gauss points and get the function value
		      // and perform the numerical integration (the area integral of the triangle).
		      QTRI[0] = 0.;
		      
		      for ( i=0; i < 7; i++ )
			{
			  // Get the value at the point.
			  XGP[0] = gp[i][0];
			  XGP[1] = gp[i][1];
			  
			  temp = Test_Function ( XGP );
			  
			  QTRI[0] += ( wg[i] * temp );
			}
		      
		      QTRI[0] *= dA;
		      
		      // Accumulate the integral to the node.
		      cv_avg[node] += QTRI[0];
		    }
		}                                       // End element node loop.
	  
	    }                                           // End element loop.
	  
	}                                               // End element type loop.
    }

  if ( 1 )
    {
      temp = 0.;

      for ( i=1; i <= grid->nn; i++ )
	{
	  temp += cv_avg[i];
	}

      printf("TOTAL SUM OF CV NUMERICAL INTEGRALS IS : %.15E\n",temp);
    }

  // Finished the accumulation. Now divide by the control volume area.
  for ( i=1; i<= grid->nn; i++ )
    {
      cv_avg[i] /= grid->cv_area[i];
    }

  // Debug code to print to screen.
  if (0)
    {
      printf("CV AVERAGE DEBUG:\n");
      for ( n=1; n <= grid->nn; n++ )
	{
	  printf("%d :\n",n);
	  printf("  cv_avg  = %f\n",(float)cv_avg[n]);
	}
    }

  // Write out the cv averages to a file.
  FILE *fp;

  fp = fopen("cv_averages.dat","w");
  
  for ( i=1; i <= grid->nn; i++ )
    {
      fprintf(fp,"%.15e\n",cv_avg[i]);
    }
  fclose(fp);

  freenull(x1);
  freenull(x2);
  freenull(x3);

  return;
}

void Annulus_Test_Function ( double *x, double *Q )
{
  double f;                         // Function value.
  double pi = M_PI;
  double r, r2;
  double rho,u,v,P,E;
  double rhoi,Mi,U,Ui,Ri;
  
  rhoi = 1.0;
  Mi = 2.0;
  Ui = 2.0;
  Ri = 2.0;
  
  r2 = (x[0]*x[0]) + (x[1]*x[1]);
  r = sqrt(r2);
      
  rho = pow( ( 1.0 + ((1.4 - 1.0)/2.0)*Mi*Mi*( 1.0 - (Ri*Ri)/r2)) , (1.0/(1.4-1.0)) );
  U   = (Ui*Ri)/r;
  u   = (x[1]*U)/r;
  v   = (-x[0]*U)/r;
  P   = ( pow( rho , 1.4) )/1.4;
  
  E   = P / (1.4 - 1.) + 0.5*rho*(u*u + v*v);

  // conserved: rho,rho*u,rho*v,E
  // primitive: rho,u,v,P

  if ( RECON_PRIM )
    {
      Q[0] = rho;
      Q[1] = u;
      Q[2] = v;
      Q[3] = P;
    }
  else
    {
      Q[0] = rho;
      Q[1] = rho*u;
      Q[2] = rho*v;
      Q[3] = E;
    }

      return;
}

void Test_Annulus_Derivatives ( GRID *grid , PARAMS p )
{
  // Debug variables.
  int i,j,k;                             // Loop counters.
  int v;                                 // Vector loop.
  int ind;                               // Component loop.
  int dmax;                              // Max number of derivatives.
  int nodeL, nodeR;                      // Edge nodes.
  int e;                                 // Element type/loop counter.
  int var;                               // More loop counters.
  int node;                              // Neighboring node.
  double dx[NDIM];                       // dx vector.
  double temp;                           // Intermediate calculation.
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
  double tri_contrib_1, tri_contrib_2;
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
  double *annulus_pv =  NULL;
  double QL[NUM_VAR], QR[NUM_VAR];
  
  // Pointers.
  int nn = grid->nn;
  double *x = grid->x;
  double *hess = grid->hess;
  double *point_value = grid->point_value;

  derivs = (char**)malloc(NUM_MOM*sizeof(char*));
  if ( derivs == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'derivs'.\n"); exit(1); }

  qprime = (double*)malloc((nn+1)*NUM_VAR*NUM_MOM*sizeof(double));
  if ( qprime == NULL ) { printf("MEMORY ERROR: Could not allocate 'qprime'.\n"); exit(0); }

  annulus_pv = (double*)malloc((nn+1)*NUM_VAR*sizeof(double));
  if ( annulus_pv == NULL ) { printf("MEMORY ERROR: Could not allocate 'annulus_pv'.\n"); exit(0); }

  ext_cv_avg = (double*)malloc((nn+1)*NUM_VAR*sizeof(double));
  if ( ext_cv_avg == NULL ) { printf("MEMORY ERROR: Could not allocate 'ext_cv_avg'.\n"); exit(0); }
  
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


  ext_cv_avg = (double*)malloc((nn+1)*NUM_VAR*sizeof(double));
  if ( ext_cv_avg == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'ext_cv_avg'.\n"); exit(1); }

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

  for ( i=1; i <= grid->nn; i++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  ext_cv_avg[i*NUM_VAR+j] = 0.;
	}
    }

  // Initialize the derivatives and the analytical point values.
  for ( i=1; i <= grid->nn; i++ )
    {
      Get_Annulus_Derivatives ( &(grid->x[i*NDIM]) , &(qprime[i*NUM_MOM*NUM_VAR]) );
      Annulus_Test_Function ( &(grid->x[i*NDIM]) , &(annulus_pv[i*NUM_VAR]) );
    }

  if ( p.order == 1 || p.order == 2 ) { dmax = 2; }
  else if ( p.order == 3 ) { dmax = 5; }
  else { dmax = 9; }

  for ( var=0; var < NUM_VAR; var++ )
    {
      if ( var == 0 )
	{
	  printf("TESTING DENSITY:\n\n");
	}
      else if ( var == 1 )
	{
	  if ( RECON_PRIM )
	    printf("TESTING X VELOCITY:\n\n");
	  else
	    printf("TESTING X MOMENTUM:\n\n");
	}
      else if ( var == 2 )
	{
	  if ( RECON_PRIM )
	    printf("TESTING Y VELOCITY:\n\n");
	  else
	    printf("TESTING Y MOMENTUM:\n\n");
	}
      else
	{
	  if ( RECON_PRIM )
	    printf("TESTING PRESSURE:\n\n");
	  else
	    printf("TESTING TOTAL ENERGY:\n\n");
	}

      // test the point values.

      norm = 0.;
      term = 0.;

      for ( i=1; i <= grid->nn; i++ )
	{
	  term = fabs( point_value[i*NUM_VAR+var] - annulus_pv[i*NUM_VAR+var] );
	  norm += term;
	}
      printf("L1 Norm of point value difference = %.15E\n",norm);
      
      norm = 0.;
      for ( i=1; i <= grid->nn; i++ )
	{
	  term = fabs( point_value[i*NUM_VAR+var] - annulus_pv[i*NUM_VAR+var] );
	  norm += (term*term);
	}
      printf("L2 Norm of point value difference = %.15E\n",norm);

      norm = fabs( point_value[1*NUM_VAR+var] - annulus_pv[1*NUM_VAR+var] );
      imax = 1;

      for ( i=2; i <= grid->nn; i++ )
	{
	  term = fabs( point_value[i*NUM_VAR+var] - annulus_pv[i*NUM_VAR+var] );
	  if ( term > norm )
	    {
	      norm = term;
	      imax = i;
	    }
	}
      printf("Infinity norm of the point values difference = %.15E\n",norm);
      printf("This occured at node %d: %f   %f\n",imax,(float)x[imax*NDIM+0],(float)x[imax*NDIM+1]);
      

      norm = 0.;
      term = 0.;
      
      // Find the l1 norm of the differences.
      for ( i=0; i < dmax; i++ )
	{
	  printf("L1 norm of %s derivative difference = ",derivs[i]);
	  norm = 0.;
	  for( j=1; j <= nn; j++ )
	    {
	      term = fabs( qprime[j*NUM_VAR*NUM_MOM+var*NUM_MOM+i] - hess[j*NUM_VAR*NUM_MOM+var*NUM_MOM+i] );
	      norm += term;
	    }
	  printf("%.15E\n",norm);
	}

      // Find the l2 norm of the differences.
      for ( i=0; i < dmax; i++ )
	{
	  printf("L2 norm of %s derivative difference = ",derivs[i]);
	  norm = 0.;
	  for( j=1; j <= nn; j++ )
	    {
	      term = fabs( qprime[j*NUM_VAR*NUM_MOM+var*NUM_MOM+i] - hess[j*NUM_VAR*NUM_MOM+var*NUM_MOM+i] );
	      norm += (term*term);
	    }
	  printf("%.15E\n",sqrt(norm));
	}

      // Find the infinity norm of the differences.
      for( i=0; i < dmax; i++ )
	{
	  printf("Infinty norm of %s derivative difference = ",derivs[i]);
	  norm = fabs( qprime[1*NUM_VAR*NUM_MOM+var*NUM_MOM+i] - hess[1*NUM_VAR*NUM_MOM+var*NUM_MOM+i] );
	  imax = 1;
	  for( j=1; j <= nn; j++ )
	    {
	      term = fabs( qprime[j*NUM_VAR*NUM_MOM+var*NUM_MOM+i] - hess[j*NUM_VAR*NUM_MOM+var*NUM_MOM+i] );
	      if ( term > norm )
		{
		  norm = term;
		  imax = j;
		}
	    }
	  printf("%.15E\n",norm);
	  printf("This occured at node %d: %f   %f\n",imax,(float)x[imax*NDIM+0],(float)x[imax*NDIM+1]);
	}

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

	  Reconstruct_Gauss ( i, nodeL, nodeR, gp[0], QL, QR, grid, p, 1 );

	  recon_val = QL[var];

	  Annulus_Test_Function ( gp[0] , QR );

	  temp = fabs( recon_val - QR[var] );
	  
	  norm_1_s += fabs(temp);
	  norm_2_s += (temp*temp);
	  if ( temp > norm_i_s )
	    {
	      norm_i_s = temp;
	      inf_subedge = i;
	    }
	  
	  // Second gp.
	  dx[0] = x2 - x[nodeL*NDIM+0];
	  dx[1] = y2 - x[nodeL*NDIM+1];
	  gp[0][0] = x2;  gp[0][1] = y2;

	  Reconstruct_Gauss ( i, nodeL, nodeR, gp[0], QL, QR, grid, p, 1 );

	  recon_val = QL[var];

	  Annulus_Test_Function ( gp[0] , QR );

	  temp = fabs( recon_val - QR[var] );

	  norm_1_s += fabs(temp);
	  norm_2_s += (temp*temp);
	  if ( temp > norm_i_s )
	    {
	      norm_i_s = temp;
	      inf_subedge = i;
	    }

	  // Third gp.
	  dx[0] = x3 - x[nodeL*NDIM+0];
	  dx[1] = y3 - x[nodeL*NDIM+1];
	  gp[0][0] = x3;  gp[0][1] = y3;

	  Reconstruct_Gauss ( i, nodeL, nodeR, gp[0], QL, QR, grid, p, 1 );

	  recon_val = QL[var];
	  
	  Annulus_Test_Function ( gp[0] , QR );
	  
	  temp = fabs( recon_val - QR[var] );

	  norm_1_s += fabs(temp);
	  norm_2_s += (temp*temp);
	  if ( temp > norm_i_s )
	    {
	      norm_i_s = temp;
	      inf_subedge = i;
	    }

	  // RIGHT NODE.

	 // First gp.
	  dx[0] = x1 - x[nodeR*NDIM+0];
	  dx[1] = y1 - x[nodeR*NDIM+1];
	  gp[0][0] = x1;  gp[0][1] = y1;

	  Reconstruct_Gauss ( i, nodeL, nodeR, gp[0], QL, QR, grid, p, 1 );

	  recon_val = QR[var];
	  
	  Annulus_Test_Function ( gp[0] , QL );
	  
	  temp = fabs( recon_val - QL[var] );
	  
	  norm_1_s += fabs(temp);
	  norm_2_s += (temp*temp);
	  if ( temp > norm_i_s )
	    {
	      norm_i_s = temp;
	      inf_subedge = i;
	    }

	  // Second gp.
	  dx[0] = x2 - x[nodeR*NDIM+0];
	  dx[1] = y2 - x[nodeR*NDIM+1];
	  gp[0][0] = x2;  gp[0][1] = y2;

	  Reconstruct_Gauss ( i, nodeL, nodeR, gp[0], QL, QR, grid, p, 1 );

	  recon_val = QR[var];
	  
	  Annulus_Test_Function ( gp[0] , QL );
	  
	  temp = fabs( recon_val - QL[var] );

	  norm_1_s += fabs(temp);
	  norm_2_s += (temp*temp);
	  if ( temp > norm_i_s )
	    {
	      norm_i_s = temp;
	      inf_subedge = i;
	    }

	  // Third gp.
	  dx[0] = x3 - x[nodeR*NDIM+0];
	  dx[1] = y3 - x[nodeR*NDIM+1];
	  gp[0][0] = x3;  gp[0][1] = y3;

	  Reconstruct_Gauss ( i, nodeL, nodeR, gp[0], QL, QR, grid, p, 1 );
	  
	  recon_val = QR[var];
	  
	  Annulus_Test_Function ( gp[0] , QL );
	  
	  temp = fabs( recon_val - QL[var] );

	  norm_1_s += fabs(temp);
	  norm_2_s += (temp*temp);
	  if ( temp > norm_i_s )
	    {
	      norm_i_s = temp;
	      inf_subedge = i;
	    }

	}

      // Loop over all the boundary edges and reconstruct to the Gauss points and look at the error created.

      norm_1_b = 0.;
      norm_2_b = 0.;
      norm_i_b = 0.;

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

	      Reconstruct_Gauss_Boundary ( real_node, gp[0], QL, grid, p);

	      recon_val = QL[var];

	      Annulus_Test_Function ( gp[0] , QR );

	      temp = fabs( recon_val - QR[var] );
	  
	      norm_1_b += fabs(temp);
	      norm_2_b += (temp*temp);
	      if ( temp > norm_i_b )
		{
		  norm_i_b = temp;
		  inf_bedge = i;
		}
	      
	      // Second gp.
	      dx[0] = x2 - x[real_node*NDIM+0];
	      dx[1] = y2 - x[real_node*NDIM+1];
	      gp[0][0] = x2;  gp[0][1] = y2;

	      Reconstruct_Gauss_Boundary ( real_node, gp[0], QL, grid, p);

	      recon_val = QL[var];

	      Annulus_Test_Function ( gp[0] , QR );

	      temp = fabs( recon_val - QR[var] );

	      norm_1_b += fabs(temp);
	      norm_2_b += (temp*temp);
	      if ( temp > norm_i_b )
		{
		  norm_i_b = temp;
		  inf_bedge = i;
		}

	      // Third gp.
	      dx[0] = x3 - x[real_node*NDIM+0];
	      dx[1] = y3 - x[real_node*NDIM+1];
	      gp[0][0] = x3;  gp[0][1] = y3;

	      Reconstruct_Gauss_Boundary ( real_node, gp[0], QL, grid, p);

	      recon_val = QL[var];

	      Annulus_Test_Function ( gp[0] , QR );

	      temp = fabs( recon_val - QR[var] );

	      norm_1_b += fabs(temp);
	      norm_2_b += (temp*temp);
	      if ( temp > norm_i_b )
		{
		  norm_i_b = temp;
		  inf_bedge = i;
		}

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

	      Reconstruct_Gauss_Boundary ( real_node, gp[0], QL, grid, p);

	      recon_val = QL[var];

	      Annulus_Test_Function ( gp[0] , QR );

	      temp = fabs( recon_val - QR[var] );
	  
	      norm_1_b += fabs(temp);
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

	      Reconstruct_Gauss_Boundary ( real_node, gp[0], QL, grid, p);

	      recon_val = QL[var];

	      Annulus_Test_Function ( gp[0] , QR );

	      temp = fabs( recon_val - QR[var] );
	  
	      norm_1_b += fabs(temp);
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

	      Reconstruct_Gauss_Boundary ( real_node, gp[0], QL, grid, p);

	      recon_val = QL[var];

	      Annulus_Test_Function ( gp[0] , QR );

	      temp = fabs( recon_val - QR[var] );
	  
	      norm_1_b += fabs(temp);
	      norm_2_b += (temp*temp);
	      if ( temp > norm_i_b )
		{
		  norm_i_b = temp;
		  inf_bedge = i;
		}
	    }

	}

      // Now I can print out all this stuff to the screen.
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

      max_diff = 0.;

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
			  Reconstruct_Gauss_Boundary ( node, gp[0], QL, grid, p);  // one-side reconstruction. It really doesn't care if its on a boundary.

			  ext_val[0] = QL[var];

			  Annulus_Test_Function ( gp[0] , QR );
			  
			  temp = fabs( ext_val[0] - QR[var] );
			  
			  // Compare to the max.
			  max_diff = MAX(max_diff,temp);

			  // Now we accumulate to the triangle.
			  tri_int += ( wg[i] * ext_val[0] * jacobian * 0.5 );
			}
		      tri_contrib_1 = tri_int;
		      // Add to the running total for the node.
		      ext_cv_avg[node*NUM_VAR+var] += tri_int;
		      
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
			  Reconstruct_Gauss_Boundary ( node, gp[i], QL, grid, p);  // one-side reconstruction. It really doesn't care if its on a boundary.

			  ext_val[i] = QL[var];
			  
			  Annulus_Test_Function ( gp[i] , QR );
			  
			  temp = fabs( ext_val[i] - QR[var] );
			  
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
		      ext_cv_avg[node*NUM_VAR+var] += temp;
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
			  Reconstruct_Gauss_Boundary ( node, gp[0], QL, grid, p);  // one-side reconstruction. It really doesn't care if its on a boundary.

			  ext_val[0] = QL[var];

			  Annulus_Test_Function ( gp[0] , QR );
			  
			  temp = fabs( ext_val[0] - QR[var] );
			  
			  // Compare to the max.
			  max_diff = MAX(max_diff,temp);
			  
			  // Now we accumulate to the triangle.
			  tri_int += ( wg[i] * ext_val[0] * jacobian * 0.5 );
			}
		      tri_contrib_2 = tri_int;
		      // Add to the running total for the node.
		      ext_cv_avg[node*NUM_VAR+var] += tri_int;
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
			  Reconstruct_Gauss_Boundary ( node, gp[i], QL, grid, p);  // one-side reconstruction. It really doesn't care if its on a boundary.

			  ext_val[i] = QL[var];

			  Annulus_Test_Function ( gp[i] , QR );
			  
			  temp = fabs( ext_val[i] - QR[var] );
			  
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
		      ext_cv_avg[node*NUM_VAR+var] += temp;
		    }
		}                                       // End element node loop.
	      
	    }                                           // End element loop.

	}                                               // End element type loop.

      // Now get the control volume average of the value.
      for ( i=1; i <= grid->nn; i++ )
	{
	  //ext_cv_avg[i*NUM_VAR+var] = ext_cv_avg[i*NUM_VAR+var] * grid->cv_area[i];
	}

      // Report the largest difference.
      printf("EXTRAPOLATION TEST: The largest error at the Gauss points is: %.15e\n",max_diff);
      
      // Compute the norms of the difference.
      norm = 0.;
      for ( i=1; i <= nn; i++ )
	{
	  norm += fabs( (grid->nQ[i*NUM_VAR+var]*grid->cv_area[i]) - ext_cv_avg[i*NUM_VAR+var] );
	}
      printf("EXTRAPOLATION TEST: L1 norm of difference: %.15e\n",norm);
      
      norm = 0.;
      for ( i=1; i <= nn; i++ )
	{
	  temp = (grid->nQ[i*NUM_VAR+var]*grid->cv_area[i]) - ext_cv_avg[i*NUM_VAR+var];
	  norm += (temp*temp);
	}
      printf("EXTRAPOLATION TEST: L2 norm of difference: %.15e\n",sqrt(norm));
      
      norm = 0.;
      for ( i=1; i <= nn; i++ )
	{
	  temp = fabs( (grid->nQ[i*NUM_VAR+var]*grid->cv_area[i]) - ext_cv_avg[i*NUM_VAR+var] );
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
	  total_integral += ( grid->nQ[i*NUM_VAR+var]*grid->cv_area[i] );
	}

      printf("Total integral of the test function over the domain using divergence and Gaussian Quadrature is %.15E\n",
	     total_integral);
      
      temp = total_integral;
      total_integral = 0.;
      for ( i=1; i <= nn; i++ )
	{
	  total_integral += ( ext_cv_avg[i*NUM_VAR+var] );
	}
      
      printf("Total integral of the test function over the domain using extrapolation and Gaussian Quadrature is %.15E\n",
	     total_integral);
      
      printf("The difference in the integral values is %.15E\n",fabs(total_integral - temp));
      
    } // end variable loop.

  for ( i=0; i < NUM_MOM; i++ )
    freenull(derivs[i]);
  
  freenull(derivs);
 
  freenull(qprime);

  freenull(annulus_pv);
  
  freenull(ext_cv_avg);
 
  return;
}
