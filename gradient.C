//=============================================================
// 
//  gradient.C
//  
//  Functions to calculate the gradient of the conserved variables.
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
#include "gradient.h"
#include "hessian.h"
#include "grid_io.h"


//=============================================================
// 
//  Compute_Gradient_Green()
//
//  Calculates the gradient of the dependent variables using
//  Green-Gauss.
//  
//  GRID *grid                        // The grid.
//
//=============================================================

void Compute_Gradient_Green ( GRID *grid )
{
  int i,j,k;                             // Loop counters.
  int n,b,s;                             // Node,boundary,segment.
  int e;                                 // Element loop counter.
  int var;                               // Variable loop counter.
  int node;                              // Element vertex.
  int num_vert;                          // Number of vertices in an element.
  int gelem;                             // Global element number.
  
  int nodeL, nodeR;                      // Edge nodes.
  double Qav;                            // Average Q value.
  double xc,yc,qc;                       // Values at the centroid.
  double xmid,ymid,qmid;                 // Values at the edge midpoint.
  double nx,ny;                          // Edge normal vector.


  // Gradient debug. Give it a cheese ball function.
  #if DEBUG_GRAD
  for ( i=1; i <= grid->nn; i++ )
    {
      grid->nQ[i*NUM_VAR+0] = 1.*grid->x[i*NDIM] + 2.*grid->x[i*NDIM+1];
      grid->nQ[i*NUM_VAR+1] = 3.*grid->x[i*NDIM] + 4.*grid->x[i*NDIM+1];
      grid->nQ[i*NUM_VAR+2] = 5.*grid->x[i*NDIM] + 6.*grid->x[i*NDIM+1];
      grid->nQ[i*NUM_VAR+3] = 7.*grid->x[i*NDIM] + 8.*grid->x[i*NDIM+1];
    }
  #endif

  // Check for memory.
  if ( grid->grad == NULL )
    {
      grid->grad = (double*)malloc((grid->nn + 1)*NDIM*NUM_VAR*sizeof(double));
      if ( grid->grad == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'grad'.\n"); exit(1); }
    }

  // Make sure that grad is clean.
  for ( n=0; n <= grid->nn; n++ )
    {
      for ( i=0; i < NUM_VAR; i++ )
	{
	  for ( j=0; j < NDIM; j++ )
	    {
	      grid->grad[n*NDIM*NUM_VAR + i*NDIM + j] = 0.;
	    }
	}
    }

  // Loop over the edges in the grid and calculate the contribution.
  if ( 0 )
    {
      for ( i=1; i <= grid->nedges; i++ )
	{
	  nodeL = grid->edges[i*2];
	  nodeR = grid->edges[i*2+1];

	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      Qav = 0.5*(grid->nQ[nodeL*NUM_VAR+j] + grid->nQ[nodeR*NUM_VAR+j]);
	      
	      for ( k=0; k < NDIM; k++ )
		{
		  grid->grad[nodeL*NUM_VAR*NDIM + NDIM*j + k] += (grid->xn_edges[i*3+k] *
								  grid->xn_edges[i*3+2] * Qav);
		  grid->grad[nodeR*NUM_VAR*NDIM + NDIM*j + k] -= (grid->xn_edges[i*3+k] *
								  grid->xn_edges[i*3+2] * Qav);
		}
	    } 
	}
    }
  else
    {
      for ( e=Tri; e <= Quad; e++ )                                                       // Element types loop.
	{
	  // Find the number of nodes and edges for this element type.
	  num_vert = NumberNodesForElement(e);
	  
	  for ( j=1; j <= grid->num_elem[e]; j++ )                                        // Element loop.
	    {
	      gelem = LocalToGlobalID ( grid->num_elem,j,e);
	      
	      xc = grid->el_cent[gelem*NDIM];
	      yc = grid->el_cent[gelem*NDIM+1];
	      
	      // Loop over the nodes and do the computation across the median dual faces.
	      for ( var=0; var < NUM_VAR; var++ )
		{
		  qc = 0.;

		  for ( n=0; n < num_vert; n++ )
		    {
		      node = grid->c2n[e][j*num_vert+n];
		      qc += ( grid->nQ[node*NUM_VAR+var]);
		    }
		  
		  qc /= num_vert;

		  for ( n=0; n < num_vert; n++ )
		    {
		      nodeL = grid->c2n[e][j*num_vert+n];
		      nodeR = grid->c2n[e][j*num_vert+( (n+1)%num_vert) ];

		      xmid = 0.5*(grid->x[nodeL*NDIM+0] + grid->x[nodeR*NDIM+0]);
		      ymid = 0.5*(grid->x[nodeL*NDIM+1] + grid->x[nodeR*NDIM+1]);

		      qmid = 0.5*(grid->nQ[nodeL*NUM_VAR+var] + grid->nQ[nodeR*NUM_VAR+var]);

		      Qav = 0.5*(qc+qmid);

		      nx = yc-ymid;
		      ny = xmid-xc;

		      // Accumulate the gradient.
		      grid->grad[nodeL*NUM_VAR*NDIM + NDIM*var + 0] += (Qav*nx);
		      grid->grad[nodeL*NUM_VAR*NDIM + NDIM*var + 1] += (Qav*ny);

		      grid->grad[nodeR*NUM_VAR*NDIM + NDIM*var + 0] -= (Qav*nx);
		      grid->grad[nodeR*NUM_VAR*NDIM + NDIM*var + 1] -= (Qav*ny);
		    }  // End vertex loop.
		}      // End variable loop.
	    }          // End Element loop.		  
	}              // End Element type loop.
    }  // end else.

  // Loop over the boundaries to close off the control volumes.
  // I'm going to try 5/6 + 1/6 distribution since that's in ux_2d and
  // I believe Dr. Anderson said that was the case. There should be a reason
  // for this in a paper by Barth!

  // Changed back to 3/4 and 1/4. Neither gets it right on the corner boundaries
  // I guess since the discretization is low on the tested grid.
  for ( i=1; i <= grid->nbedges; i++ )
    {
      n = grid->bedges[i*5+0];
      b = grid->bedges[i*5+3];
      s = grid->bedges[i*5+4];

      nodeL = grid->bs[b][s][0];
      nodeR = grid->bs[b][s][1];

      for ( j=0; j < NUM_VAR; j++ )
	{
	  // which node is it?
	  if ( n==nodeL )
	    {
	      Qav = (3./4.)*grid->nQ[nodeL*NUM_VAR+j] + (1./4.)*grid->nQ[nodeR*NUM_VAR+j];
	    }
	  else
	    {
	      Qav = (1./4.)*grid->nQ[nodeL*NUM_VAR+j] + (3./4.)*grid->nQ[nodeR*NUM_VAR+j];
	    }
	  grid->grad[n*NUM_VAR*NDIM + NDIM*j + 0] += grid->xn_bedges[i*3+0] *
		                                     grid->xn_bedges[i*3+2] * Qav;
	  grid->grad[n*NUM_VAR*NDIM + NDIM*j + 1] += grid->xn_bedges[i*3+1] *
	                                             grid->xn_bedges[i*3+2] * Qav;
	}
    }
  
  // Get the cv average value.
  for ( n=1; n <= grid->nn; n++ )
    {
      for ( i=0; i < NUM_VAR; i++ )
	{
	  for ( j=0; j < NDIM; j++ )
	    {
	      grid->grad[n*NDIM*NUM_VAR + i*NDIM + j] /= grid->cv_area[n];
	    }
	}
    }

  #if DEBUG
  printf("GRAD GREEN-GAUSS DEBUG:\n");
  for ( i=1; i <= grid->nn; i++ )
    {
      printf("Node %d:\n",i);
      for ( j=0; j < NUM_VAR; j++ )
	{
	  printf("  Var %d: qx = %f\n",j,(float)grid->grad[i*NDIM*NUM_VAR + j*NDIM + 0]);
	  printf("         qy = %f\n",(float)grid->grad[i*NDIM*NUM_VAR + j*NDIM + 1]);
	}
    }
  #endif

  #if 0
  FILE *fp = NULL;
  fp = fopen("grad_vector.dat","w");
  for ( i=1; i <= grid->nn; i++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  fprintf(fp,"%.15E\n",grid->grad[i*NDIM*NUM_VAR+j*NDIM+0]);
	  fprintf(fp,"%.15E\n",grid->grad[i*NDIM*NUM_VAR+j*NDIM+1]);
	}
    }
  fclose(fp);
  #endif

  return;
}


//=============================================================
// 
//  Compute_Gradient_LeastSquares()
//
//  Calculates the gradient of the dependent variables using
//  the Least Squares method as outlined in Daniel's dissertation.
//  
//  GRID *grid                        // The grid.
//
//=============================================================

void Compute_Gradient_LeastSquares ( GRID *grid )
{
  int i,j,k;                             // Loop counters.
  int n;                                 // Node.
  int nodeL, nodeR;                      // Edge nodes.
  double dQ;                             // Qright - Qleft.

  // Pointers.
  double *wl = NULL;
  double *wr = NULL;

  // Debug stuff.
  int e,v;                               // Element type/loop counter.
  int node;                              // Neighboring node.
  double dx[NDIM];                       // dx vector.
  int nn = grid->nn;
  double *x = grid->x;
  double *cv_avg = NULL;
  int imax,count,ind;
  double *qprime = NULL;
  double term, temp;
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
  double wg[7];                          // The associated weights;
  double b_r = (6. - sqrt(15.))/21.;     // Precompute some the coordinates.
  double b_t = (6. + sqrt(15.))/21.;
  double b_s = (9. + 2.*sqrt(15.))/21.;
  double b_u = (9. - 2.*sqrt(15.))/21.;
  double wA = (155. - sqrt(15.))/1200.;
  double wB = (155. + sqrt(15.))/1200.;

  double gp[7][3];                       // Gauss points for integration.
  double dA;                             // Differential area (for a triangle).
  double ext_val[7];                     // Extrapolation values at Gauss points.
  double max_diff = 0.;                  // Maximum difference between extrapolated value
                                         // and the analytical value.
  double *ext_cv_avg = NULL;             // CV average from the extrapolation.
  double total_integral;

  double fc[NDIM];
  int inf_subedge;
  int iedge, gelem;

  double norm_1_s, norm_2_s, norm_i_s;
  double xmid,ymid;
  double recon_val;

  char buff[100];

  if ( DEBUG_GRAD )
    {
      derivs = (char**)malloc(5*sizeof(char*));
      if ( derivs == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'derivs'.\n"); exit(1); }

      for ( i=0; i < 5; i++ )
	{
	  derivs[i] = NULL;
	  derivs[i] = (char*)malloc(5*sizeof(char));

	  if ( derivs[i] == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'derivs[%d]'.\n",i); exit(1); }
	}

      sprintf((derivs[0]), "Fx");
      sprintf((derivs[1]), "Fy");
      sprintf((derivs[2]), "Fxx");
      sprintf((derivs[3]), "Fyy");
      sprintf((derivs[4]), "Fxy");
    }

  // Gradient debug. Give it a cheese ball function.
  #if DEBUG_GRAD
  for ( i=1; i <= grid->nn; i++ )
    {
      //grid->nQ[i*NUM_VAR+0] = 1.*grid->x[i*NDIM] + 2.*grid->x[i*NDIM+1];

      grid->nQ[i*NUM_VAR+0] = Test_Function(&(grid->x[i*NDIM]));
      grid->nQ[i*NUM_VAR+1] = 3.*grid->x[i*NDIM] + 4.*grid->x[i*NDIM+1];
      grid->nQ[i*NUM_VAR+2] = 5.*grid->x[i*NDIM] + 6.*grid->x[i*NDIM+1];
      grid->nQ[i*NUM_VAR+3] = 7.*grid->x[i*NDIM] + 8.*grid->x[i*NDIM+1];
    }
  #endif

  // Testing stuff
  cv_avg = (double*)malloc((nn+1)*sizeof(double));
  if ( cv_avg == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'cv_avg'.\n"); exit(1); }
  
  ext_cv_avg = (double*)malloc((nn+1)*sizeof(double));
  if ( ext_cv_avg == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'ext_cv_avg'.\n"); exit(1); }

  for ( i=0; i <= grid->nn; i++ )
    ext_cv_avg[i] = 0.;

  for ( i=0; i <= grid->nn; i++ )
    cv_avg[i] = 0.;

  // Get the area averaged values of the test function.
  #if DEBUG_GRAD
  Compute_CV_Averages(grid,cv_avg);
  #endif
  
  // Check for memory.
  if ( grid->grad == NULL )
    {
      grid->grad = (double*)malloc((grid->nn + 1)*NDIM*NUM_VAR*sizeof(double));
      if ( grid->grad == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'grad'.\n"); exit(1); }
    }

  // Get the edge weights.
  if ( grid->edge_weight == NULL )
    {
      Compute_Edge_Weights(grid);
    }

  // Make sure that grad is clean.
  for ( n=1; n <= grid->nn; n++ )
    {
      for ( i=0; i < NUM_VAR; i++ )
	{
	  for ( j=0; j < NDIM; j++ )
	    {
	      grid->grad[n*NDIM*NUM_VAR + i*NDIM + j] = 0.;
	    }
	}
    }

  // Compute the gradient.
  for ( j=0; j < NUM_VAR; j++ )
    {
      for ( i=1; i <= grid->nedges; i++ )
	{
	  // Get the edge nodes.
	  nodeL = grid->edges[2*i+0];
	  nodeR = grid->edges[2*i+1];

	  // Adjust the weight pointers.
	  wl = &(grid->edge_weight[i*2*NDIM]);
	  wr = &(grid->edge_weight[i*2*NDIM+2]);

	  // Get the Q difference.
	  dQ = grid->nQ[nodeR*NUM_VAR+j] - grid->nQ[nodeL*NUM_VAR+j];
	  
	  for ( k=0; k < NDIM; k++ )
	    {
	      grid->grad[nodeL*NUM_VAR*NDIM + NDIM*j + k] += wl[k]*dQ;
	      grid->grad[nodeR*NUM_VAR*NDIM + NDIM*j + k] -= wr[k]*dQ;
	    } 
	} 
    }


  #if DEBUG
  printf("GRAD LEAST SQUARES DEBUG:\n");
  for ( i=1; i <= grid->nn; i++ )
    {
      printf("Node %d:\n",i);
      for ( j=0; j < NUM_VAR; j++ )
	{
	  printf("  Var %d: qx = %f\n",j,(float)grid->grad[i*NDIM*NUM_VAR + j*NDIM + 0]);
	  printf("          qy = %f\n",(float)grid->grad[i*NDIM*NUM_VAR + j*NDIM + 1]);
	}
    }
  #endif

  if ( DEBUG_GRAD )
    {
      // Compare what the Gradient routine returned against
      // what I'm expecting it to be. Print out any large differences.
      // Print out the norm.
      qprime = (double*)malloc((nn+1)*NUM_MOM*sizeof(double));
      if ( qprime == NULL ) { printf("MEMORY ERROR: Could not allocate 'qprime'.\n"); exit(0); }
      
      // Set the analytical values.
      for ( i=1; i <= nn; i++ )
	{
	  Test_Function_Derivatives(&(x[i*NDIM]),&(qprime[i*NUM_MOM]));
	}
      
      for ( i=0; i < 2; i++ )
	{
	  printf("Checking %s derivatives...  ",derivs[i]);
	  count=0;
	  for( j=1; j <= nn; j++ )
	    {
	      if ( fabs( qprime[j*NUM_MOM+i] - grid->grad[j*NUM_VAR*NDIM+0*NDIM+i] ) > 1.E-10 )
		{
		  count++;
		}
	    }
	  printf("Found %d problems!\n",count);
	}
      
      // Find the l1 norm of the differences.
      for ( i=0; i < 2; i++ )
	{
	  printf("L1 norm of %s derivatives = ",derivs[i]);
	  norm = 0.;
	  for( j=1; j <= nn; j++ )
	    {
	      term = fabs( qprime[j*NUM_MOM+i] - grid->grad[j*NUM_VAR*NDIM+0*NDIM+i] );
	      norm += term;
	    }
	  printf("%.15E\n",norm);
	}

      // Find the l2 norm of the differences.
      for ( i=0; i < 2; i++ )
	{
	  printf("L2 norm of %s derivatives = ",derivs[i]);
	  norm = 0.;
	  for( j=1; j <= nn; j++ )
	    {
	      term = fabs( qprime[j*NUM_MOM+i] - grid->grad[j*NUM_VAR*NDIM+0*NDIM+i] );
	      norm += (term*term);
	    }
	  printf("%.15E\n",sqrt(norm));
	}

      // Find the infinity norm of the differences.
      for( i=0; i < 2; i++ )
	{
	  printf("Infinty norm of %s derivatives = ",derivs[i]);
	  norm = fabs( qprime[1*NUM_MOM+i] - grid->grad[1*NUM_VAR*NDIM+0*NDIM+i] );
	  imax = 1;
	  for( j=1; j <= nn; j++ )
	    {
	      norm = fabs( qprime[j*NUM_MOM+i] - grid->grad[j*NUM_VAR*NDIM+0*NDIM+i] );
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
	  for ( i=0; i < 2; i++ )
	    {
	      term = fabs( qprime[j*NUM_MOM+i] - grid->grad[j*NUM_VAR*NDIM+0*NDIM+i] );
	      norm += term;
	    }
	}
      printf("L1 Norm of Derivative Vector = %.15E\n",norm);

      // Find the l2 norm of the derivative vector.
      norm = 0.;
      for( j=1; j <= nn; j++ )
	{
	  for ( i=0; i < 2; i++ )
	    {
	      term = fabs( qprime[j*NUM_MOM+i] - grid->grad[j*NUM_VAR*NDIM+0*NDIM+i] );
	      norm += (term*term);
	    }
	}
      printf("L2 Norm of Derivative Vector = %.15E\n",sqrt(norm));

      // Find the infinity norm of the derivative vector.
      norm = fabs( qprime[1*NUM_MOM+0] - grid->grad[1*NUM_VAR*NDIM+0*NDIM+0] );
      imax = 1;
      k = 0;
      for ( i=1; i < 2; i++ )
	{
	  term = fabs( qprime[1*NUM_MOM+i] - grid->grad[1*NUM_VAR*NDIM+0*NDIM+i] );
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
	  for ( i=0; i < 2; i++ )
	    {
	      term = fabs( qprime[j*NUM_MOM+i] - grid->grad[j*NUM_VAR*NDIM+0*NDIM+i] );
	      if ( term > norm )
		{
		  norm = term;
		  imax = j;
		  k = i;
		}
	    }
	}
      printf("Infinity norm of the derivative vector = %.15E\n",norm);
      printf("This occured at node %d: %f   %f\n",imax,(float)x[imax*NDIM+0],(float)x[imax*NDIM+1]);
      printf("In the %s derivative.\n",derivs[k]);

      printf("\n");
      printf("Reconstruction test...\n");

      // Loop over all the subedges and reconstruct to the subedge center from both directions
      // and look at the error created.
      
      norm_1_s = 0.;
      norm_2_s = 0.;
      norm_i_s = 0.;
      
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

	  fc[0] = (xc[0] + xmid) * 0.5;
	  fc[1] = (xc[1] + ymid) * 0.5;

	  gp[0][0] = fc[0];  gp[0][1] = fc[1];
	  
	  // LEFT SIDE FIRST.

	  dx[0] = fc[0] - x[nodeL*NDIM+0];
	  dx[1] = fc[1] - x[nodeL*NDIM+1];
	  
	  recon_val = (grid->nQ[nodeL*NUM_VAR+0] +
		       grid->grad[nodeL*NUM_VAR*NDIM+0*NDIM+0]*dx[0] +
		       grid->grad[nodeL*NUM_VAR*NDIM+0*NDIM+1]*dx[1] );

	  temp = Test_Function ( gp[0] );
	  
	  temp = fabs( recon_val - temp );
	  
	  norm_1_s += fabs(temp);
	  norm_2_s += (temp*temp);
	  if ( temp > norm_i_s )
	    {
	      norm_i_s = temp;
	      inf_subedge = i;
	    }

	  // RIGHT NODE.

	  dx[0] = fc[0] - x[nodeR*NDIM+0];
	  dx[1] = fc[1] - x[nodeR*NDIM+1];

	  recon_val = (grid->nQ[nodeR*NUM_VAR+0] +
		       grid->grad[nodeR*NUM_VAR*NDIM+0*NDIM+0]*dx[0] +
		       grid->grad[nodeR*NUM_VAR*NDIM+0*NDIM+1]*dx[1] );
	  
	  temp = Test_Function ( gp[0] );
	  
	  temp = fabs( recon_val - temp );
	  
	  norm_1_s += fabs(temp);
	  norm_2_s += (temp*temp);
	  if ( temp > norm_i_s )
	    {
	      norm_i_s = temp;
	      inf_subedge = i;
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
      
      wg[0] = 9./40.;
      wg[1] = wA;
      wg[2] = wA;
      wg[3] = wA;
      wg[4] = wB;
      wg[5] = wB;
      wg[6] = wB;
      
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
		      ext_val[i] = (grid->nQ[node*NUM_VAR+0] +
				    grid->grad[node*NUM_VAR*NDIM+0*NDIM+0]*dx[0] +
				    grid->grad[node*NUM_VAR*NDIM+0*NDIM+1]*dx[1] );

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
		  

		  // Repeat the process for the second triangle.
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
		      ext_val[i] = (grid->nQ[node*NUM_VAR+0] +
				    grid->grad[node*NUM_VAR*NDIM+0*NDIM+0]*dx[0] +
				    grid->grad[node*NUM_VAR*NDIM+0*NDIM+1]*dx[1] );

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
	  grid->nQ[i*NUM_VAR] = fabs(cv_avg[i] - ext_cv_avg[i]);
	}

      // Now write a solution file with this.
      sprintf(buff,"tecplot_diff_cv_ext.dat");
      write_tecplot_solutionC ( buff, grid );
      
    }

  if ( DEBUG_GRAD )
    {
      for ( i=0; i < 5; i++ )
	freenull(derivs[i]);
      
      freenull(derivs);
    }

  freenull(qprime);
  freenull(cv_avg);
  freenull(ext_cv_avg);
  
  return;
}




//=============================================================
// 
//  Compute_Edge_Weights()
//
//  Calculates the edge weights used in the Least Squares
//  gradient calculation.
//  
//  GRID *grid                        // The grid.
//
//=============================================================

void Compute_Edge_Weights ( GRID *grid )
{
  int i,j;                      // Loop counters.
  int nodeL,nodeR;		// Nodes making an edge.
  double dx[NDIM];              // Spatial differences.
  double temp; 			// Dummy variable.
  double *r11 = NULL;           // Terms from the formulation.
  double *r22o = NULL;
  double *r12 = NULL;
  double *s22 = NULL;


  // Allocate the memory.
  r11 = (double*)malloc((grid->nn + 1)*sizeof(double));
  if ( r11 == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'r11'.\n"); exit(1); }

  r22o = (double*)malloc((grid->nn + 1)*sizeof(double));
  if ( r22o == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'r22o'.\n"); exit(1); }

  r12 = (double*)malloc((grid->nn + 1)*sizeof(double));
  if ( r12 == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'r12'.\n"); exit(1); }

  s22 = (double*)malloc((grid->nn + 1)*sizeof(double));
  if ( s22 == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 's22'.\n"); exit(1); }

  if ( grid->edge_weight == NULL )
    {
      grid->edge_weight = (double*)malloc((grid->nedges + 1)*2*NDIM*sizeof(double));
      if ( grid->edge_weight == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'edge_weight'.\n"); exit(1); }
    }

  // Clean out all the allocated memory.
  for ( i=1; i <= grid->nn; i++ )
    {
      r11[i] = 0.;  r22o[i] = 0.;  r12[i] = 0.;  s22[i] = 0.;
    }

  // Loop over the edges to compute the constant coefficients r11 and  r12.
  for ( i=1; i <= grid->nedges; i++ )
    {
      // Get the edge nodes and their spatial distance.
      nodeL = grid->edges[2*i+0];
      nodeR = grid->edges[2*i+1];

      if ( nodeL == 60 || nodeR == 60 )
	{
	  dx[0]=0.;
	}
      
      for ( j=0; j < NDIM; j++ )
	{
	  dx[j] = grid->x[nodeR*NDIM + j] - grid->x[nodeL*NDIM + j];
	}
   
      // Compute r11, which is ||DX||^2
      temp = dx[0]*dx[0];
      r11[nodeL] += temp; 
      r11[nodeR] += temp;

      if ( !(isfinite(temp)) )
	{
	  r11[0] = 0.;
	}
      
      // Compute r12, which is DX . DY 
      temp = dx[0]*dx[1];
      r12[nodeL] += temp;

      if ( !(isfinite(temp)) )
	{
	  r11[0] = 0.;
	}

      temp = (-dx[0])*(-dx[1]);
      r12[nodeR] += temp;

      if ( !(isfinite(temp)) )
	{
	  r11[0] = 0.;
	}

      // Compute r13, which is ||DY||^2
      temp = dx[1]*dx[1];
      s22[nodeL] += temp; 
      s22[nodeR] += temp; 

      if ( !(isfinite(temp)) )
	{
	  r11[0] = 0.;
	}

    }

  // Loop over the nodes to get r22 (r22o).
  for ( i=1; i <= grid->nn; i++ )
    {
      r22o[i] = s22[i] - ( r12[i]*r12[i] )/ r11[i];

      if ( !(isfinite(r22o[i])) )
	{
	  r11[0] = 0.;
	}
    }

  // Debug.
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( !(isfinite(r11[i])) )
	{
	  printf("In Edge weight: r11[%d] is %f.\n",i,(float)r11[i]);
	}
      if ( !(isfinite(r12[i])) )
	{
	  printf("In Edge weight: r12[%d] is %f.\n",i,(float)r12[i]);
	}
      if ( !(isfinite(r22o[i])) )
	{
	  printf("In Edge weight: r22o[%d] is %f.\n",i,(float)r22o[i]);
	}
      if ( !(isfinite(s22[i])) )
	{
	  printf("In Edge weight: s22[%d] is %f.\n",i,(float)s22[i]);
	}
    }

  // Now we can build the edge weights.
  for ( i=1; i <= grid->nedges; i++ )
    {
      nodeL = grid->edges[2*i+0];
      nodeR = grid->edges[2*i+1];
      
      for ( j=0; j < NDIM; j++ )
	{
	  dx[j] = grid->x[nodeR*NDIM + j] - grid->x[nodeL*NDIM + j];
	}

      // Compute Wy and Wx for nodeL.
      // Wy
      grid->edge_weight[i*2*NDIM+1] = (+dx[1] - r12[nodeL]/r11[nodeL]*(+dx[0]))/r22o[nodeL];

      // Wx
      grid->edge_weight[i*2*NDIM+0] = (+dx[0] - r12[nodeL]*grid->edge_weight[i*2*NDIM+1])/r11[nodeL];

      
      // Compute Wy and Wx for nodeR.
      // Wy
      grid->edge_weight[i*2*NDIM+3] = (-dx[1] - r12[nodeR]/r11[nodeR]*(-dx[0]))/r22o[nodeR];

      // Wx
      grid->edge_weight[i*2*NDIM+2] = (-dx[0] - r12[nodeR]*grid->edge_weight[i*2*NDIM+3])/r11[nodeR];
    }

  // Debug
  for ( i=1; i <= grid->nedges; i++ )
    {
      nodeL = grid->edges[2*i+0];
      nodeR = grid->edges[2*i+1];

      if ( !(isfinite(grid->edge_weight[i*2*NDIM+0])) )
	{
	  printf(" In Edge weight: grid->edge_weight[%d*2*NDIM+0] is %f. nodeL = %d  nodeR = %d\n",i,(float)grid->edge_weight[i*2*NDIM+0],nodeL,nodeR);
	}

      if ( !(isfinite(grid->edge_weight[i*2*NDIM+1])) )
	{
	  printf(" In Edge weight: grid->edge_weight[%d*2*NDIM+0] is %f. nodeL = %d  nodeR = %d\n",i,(float)grid->edge_weight[i*2*NDIM+0],nodeL,nodeR);
	}

      if ( !(isfinite(grid->edge_weight[i*2*NDIM+2])) )
	{
	  printf(" In Edge weight: grid->edge_weight[%d*2*NDIM+0] is %f. nodeL = %d  nodeR = %d\n",i,(float)grid->edge_weight[i*2*NDIM+0],nodeL,nodeR);
	}
      
      if ( !(isfinite(grid->edge_weight[i*2*NDIM+3])) )
	{
	  printf(" In Edge weight: grid->edge_weight[%d*2*NDIM+0] is %f. nodeL = %d  nodeR = %d\n",i,(float)grid->edge_weight[i*2*NDIM+0],nodeL,nodeR);
	}
    }

  #if DEBUG
  printf("GRAD LEAST SQUARES EDGE WEIGHTS DEBUG:\n");
  for ( i=1; i <= grid->nedges; i++ )
    {
      printf("Edge %d:\n",i);
      for ( j=0; j < NUM_VAR; j++ )
	{
	  printf("  Edge weight %d : %.15E\n",j+1,grid->edge_weight[i*2*NDIM+j]);
	}
    }
  fflush(stdout);
  #endif

  freenull(r11);
  freenull(r12);
  freenull(s22);
  freenull(r22o);
  
  return;
}
