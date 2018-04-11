//=============================================================
// 
//  limiter.C
//  
//  Functions to compute the reconstruction limiter.
//
//  Written by - Shane Sawyer
//
//=============================================================

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "grid.h"
#include "sys_mem.h"
#include "params.h"
#include "Defined_Variables.h"
#include "residual.h"
#include "flux.h"
#include "mesh_utils.h"
#include "cvbc.h"
#include "jacobian.h"
#include "reconstruct.h"
#include "curved_boundaries.h"
#include "limiter.h"


//=============================================================
// 
//  compute_limiter_barth()
//
//  Computes Barth's limiter for every control volume.
//
//  GRID *grid;                          // The grid.
//  PARAMS p;                            // The parameters.
//
//=============================================================

void compute_limiter_barth ( GRID *grid, PARAMS p )
{
  int i,j;                                      // Loop counters.
  int isubedge;
  int nodeL,nodeR;
  double xc,yc,xmid,ymid;
  double temp;
  int gelem;
  int iedge;
  int gauss_flag = 0;
  double *dQmin = NULL;
  double *dQmax = NULL;
  double qleft[4];                             // Q on the left of the face.
  double qright[4];                            // Q on the right of the face.

  // Variables for higher order flux calculation.
  double x1,x2,x3;                       // Integration points.
  double y1,y2,y3;
  double t1,t2,t3;                       // Gaussian roots.
  double gp[2];                          // Gauss point.
  int DO_GAUSS_QUAD = 0;

  double *phi = NULL;

  if ( grid->phi == NULL )
    {
      grid->phi = (double*)malloc( (grid->nn+1)*2*4*sizeof(double) );
      if ( grid->phi == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'phi'.\n"); exit(1); }
    }

  // The 8 in the memory allocation is for phi and s for all variables.

  phi = grid->phi;

  dQmin = (double*)malloc( (grid->nn+1)*NUM_VAR*sizeof(double));
  dQmax = (double*)malloc( (grid->nn+1)*NUM_VAR*sizeof(double));

  if ( dQmin == NULL )
    { printf("MEMORY ERROR: COULD NOT ALLOCATE 'dQmin'.\n"); exit(1); }

  if ( dQmax == NULL )
    { printf("MEMORY ERROR: COULD NOT ALLOCATE 'dQmax'.\n"); exit(1); }

  // Clean out the memory.
  for ( i=1; i <= grid->nn; i++ )
    {
      for ( j=0; j < 8; j++ )
	phi[i*8+j] = 0.0;

      // Initialize these to 'high' values.
      for ( j=0; j < NUM_VAR; j++ )
	{
	  dQmin[i*NUM_VAR+j] = 0.0;
	  dQmax[i*NUM_VAR+j] = 0.0;
	}
    }

  // Initialize the maximum and minimum CV averages to be the current CV value.
  for ( i=1; i <= grid->nn; i++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  dQmin[i*NUM_VAR+j] = grid->nQ[i*NUM_VAR+j];
	  dQmax[i*NUM_VAR+j] = grid->nQ[i*NUM_VAR+j];
	}
    }

  // Initialize the limiter to 1.0
  for ( i=1; i <= grid->nn; i++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  phi[i*8+j*2+0] = 1.0;
	}
    }

  // The general procedure is to loop over all the nodes in the mesh, then for each variable:
  // 1. find the largest negative Q-Q_i and positive Q-Q_i
  // 2. Get the unconstrained extrapolated value at the quadrature point(s).
  // 3. Compute the maximum value of Phi_ij for the quadrature point(s) j.
  // 4. Set phi(i) = min(Phi_ij).

  // Loop over the edges to find the minimal and maximal values of Qbar in the nodes nearest neighborhood.
  for ( i=1; i <= grid->nedges; i++ )
    {
      nodeL = grid->edges[i*2+0];
      nodeR = grid->edges[i*2+1];

      // Find the minimal/maximal values for each variable.
      for ( j=0; j < NUM_VAR; j++ )
	{
	  // Compare the right node the left node's values and vice versa.
	  dQmin[nodeL*NUM_VAR+j] = MIN(   dQmin[nodeL*NUM_VAR+j] ,
					  grid->nQ[nodeR*NUM_VAR+j] );
	  dQmax[nodeL*NUM_VAR+j] = MAX(   dQmax[nodeL*NUM_VAR+j] ,
					  grid->nQ[nodeR*NUM_VAR+j] );
	  
	  dQmin[nodeR*NUM_VAR+j] = MIN(   dQmin[nodeR*NUM_VAR+j] ,
					  grid->nQ[nodeL*NUM_VAR+j] );
	  dQmax[nodeR*NUM_VAR+j] = MAX(   dQmax[nodeR*NUM_VAR+j] ,
					  grid->nQ[nodeL*NUM_VAR+j] );
	  
	}
    }

  // Now we do an unconstrained reconstruction to every quadrature point and apply the Phi function of
  // Barth and Jespersen.

  // Figure out if we are doing Gaussian Quadrature on the subedges.
  if ( p.order < 3 )
    {
      DO_GAUSS_QUAD = 0;
    }
  else if ( (p.order == 3) && (grid->citer > (p.foits+p.soits)) )
    {
      DO_GAUSS_QUAD = 1;
    }
  else if ( (p.order == 4) && (grid->citer > (p.foits+p.soits)) )
    {
      DO_GAUSS_QUAD = 1;
    }
  
  if ( FULL_GAUSS_FLUX )
    DO_GAUSS_QUAD = 1;
  
  if ( EDGE_BASED )
    DO_GAUSS_QUAD = 0;
  
  for ( isubedge=1; isubedge <= grid->nsubedges; isubedge++ )
    {
      // Extract the real edge this corresponds to.
      iedge = grid->subedges[isubedge*2];

      // Get the left and right nodes.
      nodeL = grid->edges[iedge*2 + 0];
      nodeR = grid->edges[iedge*2 + 1];

      xmid = grid->xm_subedges[isubedge*NDIM+0];
      ymid = grid->xm_subedges[isubedge*NDIM+1];

      // Find the global element.
      gelem = grid->subedges[isubedge*2+1];

      // Get the element centroid.
      xc = grid->el_cent[gelem*2];
      yc = grid->el_cent[gelem*2+1];

      if ( CYLINDER && CYLINDER_HACK )
	{
	  if ( (grid->node_order[nodeL] == 1) &&
	       (grid->node_order[nodeR] == 1) )
	    {
	      if ( DO_GAUSS_QUAD == 1 )
		gauss_flag = 1;
	      else
		gauss_flag = 0;

	      DO_GAUSS_QUAD = 0;
	    }
	}

      if ( DO_GAUSS_QUAD )
	{
	  // Calculate the Gauss points.
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
	  
	  // For each gauss point, we will reconstruct the solution.
	  
	  gp[0] = grid->xg_subedges[isubedge*6+0];
	  gp[1] = grid->xg_subedges[isubedge*6+1];
	  Reconstruct_Gauss_UL ( isubedge, nodeL, nodeR, gp, qleft, qright, grid, p, grid->citer );

	  // Do the left node first.
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      if ( ( qleft[j] - grid->nQ[nodeL*NUM_VAR+j] ) > 0. )
		{
		  temp = MIN( 1.0 , ( dQmax[nodeL*NUM_VAR+j] - grid->nQ[nodeL*NUM_VAR+j] ) / (qleft[j] - grid->nQ[nodeL*NUM_VAR+j]) );
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , temp );
		}
	      else if ( ( qleft[j] - grid->nQ[nodeL*NUM_VAR+j] ) < 0. )
		{
		  temp = MIN( 1.0 , ( dQmin[nodeL*NUM_VAR+j] - grid->nQ[nodeL*NUM_VAR+j] ) / (qleft[j] - grid->nQ[nodeL*NUM_VAR+j]) );
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , temp );
		}
	      else
		{
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , 1.0 );
		}
	    }

	  // Do the right node.
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      if ( ( qright[j] - grid->nQ[nodeR*NUM_VAR+j] ) > 0. )
		{
		  temp = MIN( 1.0 , ( dQmax[nodeR*NUM_VAR+j] - grid->nQ[nodeR*NUM_VAR+j] ) / (qright[j] - grid->nQ[nodeR*NUM_VAR+j]) );
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , temp );
		}
	      else if ( ( qright[j] - grid->nQ[nodeR*NUM_VAR+j] ) < 0. )
		{
		  temp = MIN( 1.0 , ( dQmin[nodeR*NUM_VAR+j] - grid->nQ[nodeR*NUM_VAR+j] ) / (qright[j] - grid->nQ[nodeR*NUM_VAR+j]) );
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , temp );
		}
	      else
		{
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , 1.0 );
		}
	    }


	  gp[0] = grid->xg_subedges[isubedge*6+2];
	  gp[1] = grid->xg_subedges[isubedge*6+3];
	  Reconstruct_Gauss_UL ( isubedge, nodeL, nodeR, gp, qleft, qright, grid, p, grid->citer );

	  // Do the left node first.
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      if ( ( qleft[j] - grid->nQ[nodeL*NUM_VAR+j] ) > 0. )
		{
		  temp = MIN( 1.0 , ( dQmax[nodeL*NUM_VAR+j] - grid->nQ[nodeL*NUM_VAR+j] ) / (qleft[j] - grid->nQ[nodeL*NUM_VAR+j]) );
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , temp );
		}
	      else if ( ( qleft[j] - grid->nQ[nodeL*NUM_VAR+j] ) < 0. )
		{
		  temp = MIN( 1.0 , ( dQmin[nodeL*NUM_VAR+j] - grid->nQ[nodeL*NUM_VAR+j] ) / (qleft[j] - grid->nQ[nodeL*NUM_VAR+j]) );
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , temp );
		}
	      else
		{
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , 1.0 );
		}
	    }

	  // Do the right node.
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      if ( ( qright[j] - grid->nQ[nodeR*NUM_VAR+j] ) > 0. )
		{
		  temp = MIN( 1.0 , ( dQmax[nodeR*NUM_VAR+j] - grid->nQ[nodeR*NUM_VAR+j] ) / (qright[j] - grid->nQ[nodeR*NUM_VAR+j]) );
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , temp );
		}
	      else if ( ( qright[j] - grid->nQ[nodeR*NUM_VAR+j] ) < 0. )
		{
		  temp = MIN( 1.0 , ( dQmin[nodeR*NUM_VAR+j] - grid->nQ[nodeR*NUM_VAR+j] ) / (qright[j] - grid->nQ[nodeR*NUM_VAR+j]) );
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , temp );
		}
	      else
		{
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , 1.0 );
		}
	    }


	  gp[0] = grid->xg_subedges[isubedge*6+4];
	  gp[1] = grid->xg_subedges[isubedge*6+5];
	  Reconstruct_Gauss_UL ( isubedge, nodeL, nodeR, gp, qleft, qright, grid, p, grid->citer );

	  // Do the left node first.
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      if ( ( qleft[j] - grid->nQ[nodeL*NUM_VAR+j] ) > 0. )
		{
		  temp = MIN( 1.0 , ( dQmax[nodeL*NUM_VAR+j] - grid->nQ[nodeL*NUM_VAR+j] ) / (qleft[j] - grid->nQ[nodeL*NUM_VAR+j]) );
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , temp );
		}
	      else if ( ( qleft[j] - grid->nQ[nodeL*NUM_VAR+j] ) < 0. )
		{
		  temp = MIN( 1.0 , ( dQmin[nodeL*NUM_VAR+j] - grid->nQ[nodeL*NUM_VAR+j] ) / (qleft[j] - grid->nQ[nodeL*NUM_VAR+j]) );
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , temp );
		}
	      else
		{
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , 1.0 );
		}
	    }

	  // Do the right node.
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      if ( ( qright[j] - grid->nQ[nodeR*NUM_VAR+j] ) > 0. )
		{
		  temp = MIN( 1.0 , ( dQmax[nodeR*NUM_VAR+j] - grid->nQ[nodeR*NUM_VAR+j] ) / (qright[j] - grid->nQ[nodeR*NUM_VAR+j]) );
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , temp );
		}
	      else if ( ( qright[j] - grid->nQ[nodeR*NUM_VAR+j] ) < 0. )
		{
		  temp = MIN( 1.0 , ( dQmin[nodeR*NUM_VAR+j] - grid->nQ[nodeR*NUM_VAR+j] ) / (qright[j] - grid->nQ[nodeR*NUM_VAR+j]) );
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , temp );
		}
	      else
		{
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , 1.0 );
		}
	    }
	  
	}
      else
	{
	  // Get the Q variables and if high order,  reconstruct the solution.
	  Reconstruct_UL( isubedge, nodeL, nodeR, qleft, qright, grid, p, grid->citer);

	  // Do the left node first.
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      if ( ( qleft[j] - grid->nQ[nodeL*NUM_VAR+j] ) > 0. )
		{
		  temp = MIN( 1.0 , ( dQmax[nodeL*NUM_VAR+j] - grid->nQ[nodeL*NUM_VAR+j] ) / (qleft[j] - grid->nQ[nodeL*NUM_VAR+j]) );
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , temp );
		}
	      else if ( ( qleft[j] - grid->nQ[nodeL*NUM_VAR+j] ) < 0. )
		{
		  temp = MIN( 1.0 , ( dQmin[nodeL*NUM_VAR+j] - grid->nQ[nodeL*NUM_VAR+j] ) / (qleft[j] - grid->nQ[nodeL*NUM_VAR+j]) );
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , temp );
		}
	      else
		{
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , 1.0 );
		}
	    }

	  // Do the right node.
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      if ( ( qright[j] - grid->nQ[nodeR*NUM_VAR+j] ) > 0. )
		{
		  temp = MIN( 1.0 , ( dQmax[nodeR*NUM_VAR+j] - grid->nQ[nodeR*NUM_VAR+j] ) / (qright[j] - grid->nQ[nodeR*NUM_VAR+j]) );
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , temp );
		}
	      else if ( ( qright[j] - grid->nQ[nodeR*NUM_VAR+j] ) < 0. )
		{
		  temp = MIN( 1.0 , ( dQmin[nodeR*NUM_VAR+j] - grid->nQ[nodeR*NUM_VAR+j] ) / (qright[j] - grid->nQ[nodeR*NUM_VAR+j]) );
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , temp );
		}
	      else
		{
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , 1.0 );
		}
	    }

	}

      if ( gauss_flag )
	DO_GAUSS_QUAD = 1;

    }

  // Set s to the same value as phi.
  for ( i=1; i <= grid->nn; i++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  phi[i*8 + j*2 + 1] = phi[i*8 + j*2 + 0];
	}
    }

  // Go ahead and put the minimum phi in the first spot.
  for ( i=1; i <= grid->nn; i++ )
    {
      temp = phi[i*8 + 0*2 + 0];

      for ( j=1; j < NUM_VAR; j++ )
	{
	  temp = MIN ( temp, phi[i*8 + j*2 + 0] );
	}

      phi[i*8 + 0*2 + 0] = temp;
      phi[i*8 + 0*2 + 1] = temp;

      // For high-order, I'll apply the same phi to the Hessian.
    }

  if ( LIMITER_NEIGHBOR_AVERAGE )  // As implemented, I am only looking at the computed minimum phi.
    {
      double *phi_avg = NULL;                      // For neighbor averaging.
      double eps = 1.0E-2;
      double sum;
      int index1, index2, num_nn;
      
      phi_avg = (double*)malloc( (grid->nn+1)*sizeof(double) );
      if ( phi_avg == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'phi_avg'.\n"); exit(1); }
      
      for ( i=1; i <= grid->nn; i++ )
	{
	  // skip those close enough to 0 or to 1.
	  if ( phi[i*8+0] < (0.0+eps) )
	    {
	      phi_avg[i] = phi[i*8];
	      continue;
	    }

	  if ( phi[i*8+0] > (1.0-eps) )
	    {
	      phi_avg[i] = phi[i*8];
	      continue;
	    }

	  sum = 0.;
	  num_nn = 0;

	  index1 = grid->nnsn[i];
	  index2 = grid->nnsn[i+1];

	  for ( j=index1; j < index2; j++ )
	    {
	      nodeR = grid->nsn[j];  // grab the neighbor.

	      sum += phi[nodeR*8+0];
	      num_nn++;
	    }

	  sum += phi[i*8+0];
	  num_nn++;

	  sum = sum / ( (double)num_nn );
	  phi_avg[i] = sum;
	}

      for ( i=1; i <= grid->nn; i++ )
	{
	  phi[i*8+0] = phi_avg[i];
	  phi[i*8+1] = phi_avg[i];
	}
      
      freenull(phi_avg);
    }

  // Sanity check. All values of phi must be in [0,1]!

  for ( i=1; i <= grid->nn; i++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  if ( ( phi[i*8 + j*2 + 0] < 0.0 ) ||
	       ( phi[i*8 + j*2 + 0] > 1.0 )  )
	    {
	      printf("FATAL ERROR: Phi for node %d and variable %d does not have the proper value: < %.15E  >\n",i,j,phi[i*8+j*2+0]);
	      fflush(stdout);
	      exit(1);
	    }
	}
    }

  freenull(dQmin);
  freenull(dQmax);

  return;
}


//=============================================================
// 
//  compute_limiter_venkatakrishnan()
//
//  Computes Venkatakrishnan's modified limiter for every control volume.
//
//  GRID *grid;                          // The grid.
//  PARAMS p;                            // The parameters.
//
//=============================================================

void compute_limiter_venkatakrishnan ( GRID *grid, PARAMS p )
{
  int i,j;                                      // Loop counters.
  int isubedge;
  int nodeL,nodeR;
  double xc,yc,xmid,ymid;
  double temp;
  int gelem;
  int iedge;
  int gauss_flag = 0;
  double *dQmin = NULL;
  double *dQmax = NULL;
  double qleft[4];                             // Q on the left of the face.
  double qright[4];                            // Q on the right of the face.

  // Variables for higher order flux calculation.
  double x1,x2,x3;                       // Integration points.
  double y1,y2,y3;
  double t1,t2,t3;                       // Gaussian roots.
  double gp[2];                          // Gauss point.
  int DO_GAUSS_QUAD = 0;

  double *phi = NULL;

  if ( grid->phi == NULL )
    {
      grid->phi = (double*)malloc( (grid->nn+1)*2*4*sizeof(double) );
      if ( grid->phi == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'phi'.\n"); exit(1); }
    }

  // The 8 in the memory allocation is for phi and s for all variables.

  phi = grid->phi;

  dQmin = (double*)malloc( (grid->nn+1)*NUM_VAR*sizeof(double));
  dQmax = (double*)malloc( (grid->nn+1)*NUM_VAR*sizeof(double));

  if ( dQmin == NULL )
    { printf("MEMORY ERROR: COULD NOT ALLOCATE 'dQmin'.\n"); exit(1); }

  if ( dQmax == NULL )
    { printf("MEMORY ERROR: COULD NOT ALLOCATE 'dQmax'.\n"); exit(1); }

  // Clean out the memory.
  for ( i=1; i <= grid->nn; i++ )
    {
      for ( j=0; j < 8; j++ )
	phi[i*8+j] = 0.0;

      // Initialize these to 'high' values.
      for ( j=0; j < NUM_VAR; j++ )
	{
	  dQmin[i*NUM_VAR+j] = 0.0;
	  dQmax[i*NUM_VAR+j] = 0.0;
	}
    }

  // Initialize the maximum and minimum CV averages to be the current CV value.
  for ( i=1; i <= grid->nn; i++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  dQmin[i*NUM_VAR+j] = grid->nQ[i*NUM_VAR+j];
	  dQmax[i*NUM_VAR+j] = grid->nQ[i*NUM_VAR+j];
	}
    }

  // Initialize the limiter to 1.0
  for ( i=1; i <= grid->nn; i++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  phi[i*8+j*2+0] = 1.0;
	}
    }

  // The general procedure is to loop over all the nodes in the mesh, then for each variable:
  // 1. find the largest negative Q-Q_i and positive Q-Q_i
  // 2. Get the unconstrained extrapolated value at the quadratuer point(s).
  // 3. Compute the maximum value of Phi_ij for the quadrature point(s) j.
  // 4. Set phi(i) = min(Phi_ij).

  // Loop over the edges to find the minimal and maximal values of Qbar in the nodes nearest neighborhood.
  for ( i=1; i <= grid->nedges; i++ )
    {
      nodeL = grid->edges[i*2+0];
      nodeR = grid->edges[i*2+1];

      // Find the minimal/maximal values for each variable.
      for ( j=0; j < NUM_VAR; j++ )
	{
	  // Compare the right node the left node's values and vice versa.
	  dQmin[nodeL*NUM_VAR+j] = MIN(   dQmin[nodeL*NUM_VAR+j] ,
					  grid->nQ[nodeR*NUM_VAR+j] );
	  dQmax[nodeL*NUM_VAR+j] = MAX(   dQmax[nodeL*NUM_VAR+j] ,
					  grid->nQ[nodeR*NUM_VAR+j] );
	  
	  dQmin[nodeR*NUM_VAR+j] = MIN(   dQmin[nodeR*NUM_VAR+j] ,
					  grid->nQ[nodeL*NUM_VAR+j] );
	  dQmax[nodeR*NUM_VAR+j] = MAX(   dQmax[nodeR*NUM_VAR+j] ,
					  grid->nQ[nodeL*NUM_VAR+j] );
	  
	}
    }
  
  // Figure out if we are doing Gaussian Quadrature on the subedges.
  if ( p.order < 3 )
    {
      DO_GAUSS_QUAD = 0;
    }
  else if ( (p.order == 3) && (grid->citer > (p.foits+p.soits)) )
    {
      DO_GAUSS_QUAD = 1;
    }
  else if ( (p.order == 4) && (grid->citer > (p.foits+p.soits)) )
    {
      DO_GAUSS_QUAD = 1;
    }
  
  if ( FULL_GAUSS_FLUX )
    DO_GAUSS_QUAD = 1;
  
  if ( EDGE_BASED )
    DO_GAUSS_QUAD = 0;

  // Now we do an unconstrained reconstruction to every quadrature point and apply the Phi function of
  // Venkatakrishnan.

  for ( isubedge=1; isubedge <= grid->nsubedges; isubedge++ )
    {
      // Extract the real edge this corresponds to.
      iedge = grid->subedges[isubedge*2];

      // Get the left and right nodes.
      nodeL = grid->edges[iedge*2 + 0];
      nodeR = grid->edges[iedge*2 + 1];

      xmid = grid->xm_subedges[isubedge*NDIM+0];
      ymid = grid->xm_subedges[isubedge*NDIM+1];

      // Find the global element.
      gelem = grid->subedges[isubedge*2+1];

      // Get the element centroid.
      xc = grid->el_cent[gelem*2];
      yc = grid->el_cent[gelem*2+1];

      if ( CYLINDER && CYLINDER_HACK )
	{
	  if ( (grid->node_order[nodeL] == 1) &&
	       (grid->node_order[nodeR] == 1) )
	    {
	      if ( DO_GAUSS_QUAD == 1 )
		gauss_flag = 1;
	      else
		gauss_flag = 0;
	      
	      DO_GAUSS_QUAD = 0;
	    }
	}
      
      if ( DO_GAUSS_QUAD )
	{
	  // Calculate the Gauss points.
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
	  
	  // For each gauss point, we will reconstruct the solution.
	  
	  gp[0] = grid->xg_subedges[isubedge*6+0];
	  gp[1] = grid->xg_subedges[isubedge*6+1];
	  Reconstruct_Gauss_UL ( isubedge, nodeL, nodeR, gp, qleft, qright, grid, p, grid->citer );

	  // Do the left node first.
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      if ( ( qleft[j] - grid->nQ[nodeL*NUM_VAR+j] ) > 0. )
		{
		  temp = ( dQmax[nodeL*NUM_VAR+j] - grid->nQ[nodeL*NUM_VAR+j] ) / (qleft[j] - grid->nQ[nodeL*NUM_VAR+j]) ;
		  temp = ( temp*temp + 2.*temp )/( temp*temp + temp + 2. );
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , temp );
		}
	      else if ( ( qleft[j] - grid->nQ[nodeL*NUM_VAR+j] ) < 0. )
		{
		  temp = ( dQmin[nodeL*NUM_VAR+j] - grid->nQ[nodeL*NUM_VAR+j] ) / (qleft[j] - grid->nQ[nodeL*NUM_VAR+j]) ;
		  temp = ( temp*temp + 2.*temp )/( temp*temp + temp + 2. );
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , temp );
		}
	      else
		{
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , 1.0 );
		}
	    }

	  // Do the right node.
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      if ( ( qright[j] - grid->nQ[nodeR*NUM_VAR+j] ) > 0. )
		{
		  temp = ( dQmax[nodeR*NUM_VAR+j] - grid->nQ[nodeR*NUM_VAR+j] ) / (qright[j] - grid->nQ[nodeR*NUM_VAR+j]) ;
		  temp = ( temp*temp + 2.*temp )/( temp*temp + temp + 2. );
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , temp );
		}
	      else if ( ( qright[j] - grid->nQ[nodeR*NUM_VAR+j] ) < 0. )
		{
		  temp = ( dQmin[nodeR*NUM_VAR+j] - grid->nQ[nodeR*NUM_VAR+j] ) / (qright[j] - grid->nQ[nodeR*NUM_VAR+j]) ;
		  temp = ( temp*temp + 2.*temp )/( temp*temp + temp + 2. );
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , temp );
		}
	      else
		{
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , 1.0 );
		}
	    }


	  gp[0] = grid->xg_subedges[isubedge*6+2];
	  gp[1] = grid->xg_subedges[isubedge*6+3];
	  Reconstruct_Gauss_UL ( isubedge, nodeL, nodeR, gp, qleft, qright, grid, p, grid->citer );

	  // Do the left node first.
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      if ( ( qleft[j] - grid->nQ[nodeL*NUM_VAR+j] ) > 0. )
		{
		  temp = ( dQmax[nodeL*NUM_VAR+j] - grid->nQ[nodeL*NUM_VAR+j] ) / (qleft[j] - grid->nQ[nodeL*NUM_VAR+j]) ;
		  temp = ( temp*temp + 2.*temp )/( temp*temp + temp + 2. );
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , temp );
		}
	      else if ( ( qleft[j] - grid->nQ[nodeL*NUM_VAR+j] ) < 0. )
		{
		  temp = ( dQmin[nodeL*NUM_VAR+j] - grid->nQ[nodeL*NUM_VAR+j] ) / (qleft[j] - grid->nQ[nodeL*NUM_VAR+j]) ;
		  temp = ( temp*temp + 2.*temp )/( temp*temp + temp + 2. );
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , temp );
		}
	      else
		{
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , 1.0 );
		}
	    }

	  // Do the right node.
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      if ( ( qright[j] - grid->nQ[nodeR*NUM_VAR+j] ) > 0. )
		{
		  temp = ( dQmax[nodeR*NUM_VAR+j] - grid->nQ[nodeR*NUM_VAR+j] ) / (qright[j] - grid->nQ[nodeR*NUM_VAR+j]) ;
		  temp = ( temp*temp + 2.*temp )/( temp*temp + temp + 2. );
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , temp );
		}
	      else if ( ( qright[j] - grid->nQ[nodeR*NUM_VAR+j] ) < 0. )
		{
		  temp = ( dQmin[nodeR*NUM_VAR+j] - grid->nQ[nodeR*NUM_VAR+j] ) / (qright[j] - grid->nQ[nodeR*NUM_VAR+j]) ;
		  temp = ( temp*temp + 2.*temp )/( temp*temp + temp + 2. );
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , temp );
		}
	      else
		{
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , 1.0 );
		}
	    }
	  

	  gp[0] = grid->xg_subedges[isubedge*6+4];
	  gp[1] = grid->xg_subedges[isubedge*6+5];
	  Reconstruct_Gauss_UL ( isubedge, nodeL, nodeR, gp, qleft, qright, grid, p, grid->citer );

	  // Do the left node first.
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      if ( ( qleft[j] - grid->nQ[nodeL*NUM_VAR+j] ) > 0. )
		{
		  temp = ( dQmax[nodeL*NUM_VAR+j] - grid->nQ[nodeL*NUM_VAR+j] ) / (qleft[j] - grid->nQ[nodeL*NUM_VAR+j]) ;
		  temp = ( temp*temp + 2.*temp )/( temp*temp + temp + 2. );
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , temp );
		}
	      else if ( ( qleft[j] - grid->nQ[nodeL*NUM_VAR+j] ) < 0. )
		{
		  temp = ( dQmin[nodeL*NUM_VAR+j] - grid->nQ[nodeL*NUM_VAR+j] ) / (qleft[j] - grid->nQ[nodeL*NUM_VAR+j]) ;
		  temp = ( temp*temp + 2.*temp )/( temp*temp + temp + 2. );
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , temp );
		}
	      else
		{
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , 1.0 );
		}
	    }

	  // Do the right node.
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      if ( ( qright[j] - grid->nQ[nodeR*NUM_VAR+j] ) > 0. )
		{
		  temp = ( dQmax[nodeR*NUM_VAR+j] - grid->nQ[nodeR*NUM_VAR+j] ) / (qright[j] - grid->nQ[nodeR*NUM_VAR+j]) ;
		  temp = ( temp*temp + 2.*temp )/( temp*temp + temp + 2. );
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , temp );
		}
	      else if ( ( qright[j] - grid->nQ[nodeR*NUM_VAR+j] ) < 0. )
		{
		  temp = ( dQmin[nodeR*NUM_VAR+j] - grid->nQ[nodeR*NUM_VAR+j] ) / (qright[j] - grid->nQ[nodeR*NUM_VAR+j]) ;
		  temp = ( temp*temp + 2.*temp )/( temp*temp + temp + 2. );
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , temp );
		}
	      else
		{
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , 1.0 );
		}
	    }
	  
	}
      else
	{
	  // Get the Q variables and if high order, then reconstruct the solution.
	  Reconstruct_UL( isubedge, nodeL, nodeR, qleft, qright, grid, p, grid->citer);

	  // Do the left node first.
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      if ( ( qleft[j] - grid->nQ[nodeL*NUM_VAR+j] ) > 0. )
		{
		  temp = ( dQmax[nodeL*NUM_VAR+j] - grid->nQ[nodeL*NUM_VAR+j] ) / (qleft[j] - grid->nQ[nodeL*NUM_VAR+j]) ;
		  temp = ( temp*temp + 2.*temp )/( temp*temp + temp + 2. );
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , temp );
		}
	      else if ( ( qleft[j] - grid->nQ[nodeL*NUM_VAR+j] ) < 0. )
		{
		  temp = ( dQmin[nodeL*NUM_VAR+j] - grid->nQ[nodeL*NUM_VAR+j] ) / (qleft[j] - grid->nQ[nodeL*NUM_VAR+j]) ;
		  temp = ( temp*temp + 2.*temp )/( temp*temp + temp + 2. );
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , temp );
		}
	      else
		{
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , 1.0 );
		}
	    }

	  // Do the right node.
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      if ( ( qright[j] - grid->nQ[nodeR*NUM_VAR+j] ) > 0. )
		{
		  temp = ( dQmax[nodeR*NUM_VAR+j] - grid->nQ[nodeR*NUM_VAR+j] ) / (qright[j] - grid->nQ[nodeR*NUM_VAR+j]) ;
		  temp = ( temp*temp + 2.*temp )/( temp*temp + temp + 2. );
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , temp );
		}
	      else if ( ( qright[j] - grid->nQ[nodeR*NUM_VAR+j] ) < 0. )
		{
		  temp = ( dQmin[nodeR*NUM_VAR+j] - grid->nQ[nodeR*NUM_VAR+j] ) / (qright[j] - grid->nQ[nodeR*NUM_VAR+j]) ;
		  temp = ( temp*temp + 2.*temp )/( temp*temp + temp + 2. );
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , temp );
		}
	      else
		{
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , 1.0 );
		}
	    }

	}

      if ( gauss_flag )
	DO_GAUSS_QUAD = 1;
      
    }

  // Set s to the same value as phi.
  for ( i=1; i <= grid->nn; i++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  phi[i*8 + j*2 + 1] = phi[i*8 + j*2 + 0];
	}
    }

  // Go ahead and put the minimum phi in the first spot.
  for ( i=1; i <= grid->nn; i++ )
    {
      temp = phi[i*8 + 0*2 + 0];

      for ( j=1; j < NUM_VAR; j++ )
	{
	  temp = MIN ( temp, phi[i*8 + j*2 + 0] );
	}

      phi[i*8 + 0*2 + 0] = temp;
      phi[i*8 + 0*2 + 1] = temp;
    }

  if ( LIMITER_NEIGHBOR_AVERAGE )  // As implemented, I am only looking at the computed minimum phi.
    {
      double *phi_avg = NULL;                      // For neighbor averaging.
      double eps = 1.0E-2;
      double sum;
      int index1, index2, num_nn;
      
      phi_avg = (double*)malloc( (grid->nn+1)*sizeof(double) );
      if ( phi_avg == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'phi_avg'.\n"); exit(1); }
      
      for ( i=1; i <= grid->nn; i++ )
	{
	  // skip those close enough to 0 or to 1.
	  if ( phi[i*8+0] < (0.0+eps) )
	    {
	      phi_avg[i] = phi[i*8];
	      continue;
	    }

	  if ( phi[i*8+0] > (1.0-eps) )
	    {
	      phi_avg[i] = phi[i*8];
	      continue;
	    }

	  sum = 0.;
	  num_nn = 0;

	  index1 = grid->nnsn[i];
	  index2 = grid->nnsn[i+1];

	  for ( j=index1; j < index2; j++ )
	    {
	      nodeR = grid->nsn[j];  // grab the neighbor.

	      sum += phi[nodeR*8+0];
	      num_nn++;
	    }

	  sum += phi[i*8+0];
	  num_nn++;

	  sum = sum / ( (double)num_nn );
	  phi_avg[i] = sum;
	}

      for ( i=1; i <= grid->nn; i++ )
	{
	  phi[i*8+0] = phi_avg[i];
	  phi[i*8+1] = phi_avg[i];
	}
      
      freenull(phi_avg);
    }

  // Sanity check. All values of phi must be in [0,1]!

  for ( i=1; i <= grid->nn; i++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  if ( ( phi[i*8 + j*2 + 0] < 0.0 ) ||
	       ( phi[i*8 + j*2 + 0] > 1.0 )  )
	    {
	      printf("FATAL ERROR: Phi for node %d and variable %d does not have the proper value: < %.15E  >\n",i,j,phi[i*8+j*2+0]);
	      fflush(stdout);
	      exit(1);
	    }
	}
    }

  freenull(dQmin);
  freenull(dQmax);

  return;
}


//=============================================================
// 
//  compute_limiter_venkatakrishnan2()
//
//  Computes Venkatakrishnan's modified limiter for every control volume using
//  the second approach (based Van Albada).
//
//  GRID *grid;                          // The grid.
//  PARAMS p;                            // The parameters.
//
//=============================================================

void compute_limiter_venkatakrishnan2 ( GRID *grid, PARAMS p )
{
  int i,j;                                      // Loop counters.
  int isubedge;
  int nodeL,nodeR;
  double xc,yc,xmid,ymid;
  double temp;
  int gelem;
  int iedge;
  int start_index, end_index;
  int gauss_flag = 0;
  double dist;
  double eps = 1.0E-12;               // Prevents division by zero.
  double K = 1.0;                     // Some choices are 0.3, 1.0, 10.0 . Too large means no limiting.
  double eps2;
  double del2, del1max, del1min;
  double *dQmin = NULL;
  double *dQmax = NULL;
  double *dx = NULL;
  double qleft[4];                             // Q on the left of the face.
  double qright[4];                            // Q on the right of the face.

  // Variables for higher order flux calculation.
  double x1,x2,x3;                       // Integration points.
  double y1,y2,y3;
  double t1,t2,t3;                       // Gaussian roots.
  double gp[2];                          // Gauss point.
  int DO_GAUSS_QUAD = 0;

  double *phi = NULL;

  if ( grid->phi == NULL )
    {
      grid->phi = (double*)malloc( (grid->nn+1)*2*4*sizeof(double) );
      if ( grid->phi == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'phi'.\n"); exit(1); }
    }

  // The 8 in the memory allocation is for phi and s for all variables.

  phi = grid->phi;

  dQmin = (double*)malloc( (grid->nn+1)*NUM_VAR*sizeof(double));
  dQmax = (double*)malloc( (grid->nn+1)*NUM_VAR*sizeof(double));
  dx = (double*)malloc( (grid->nn+1)*sizeof(double));

  if ( dQmin == NULL )
    { printf("MEMORY ERROR: COULD NOT ALLOCATE 'dQmin'.\n"); exit(1); }

  if ( dQmax == NULL )
    { printf("MEMORY ERROR: COULD NOT ALLOCATE 'dQmax'.\n"); exit(1); }

  if ( dx == NULL )
    { printf("MEMORY ERROR: COULD NOT ALLOCATE 'dx'.\n"); exit(1); }

  // Clean out the memory.
  for ( i=1; i <= grid->nn; i++ )
    {
      dx[i] = 0.0;

      for ( j=0; j < 8; j++ )
	phi[i*8+j] = 0.0;

      // Initialize these to 'high' values.
      for ( j=0; j < NUM_VAR; j++ )
	{
	  dQmin[i*NUM_VAR+j] = 0.0;
	  dQmax[i*NUM_VAR+j] = 0.0;
	}
    }

  // Initialize the maximum and minimum CV averages to be the current CV value.
  for ( i=1; i <= grid->nn; i++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  dQmin[i*NUM_VAR+j] = grid->nQ[i*NUM_VAR+j];
	  dQmax[i*NUM_VAR+j] = grid->nQ[i*NUM_VAR+j];
	}
    }

  // Initialize the limiter to 1.0
  for ( i=1; i <= grid->nn; i++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  phi[i*8+j*2+0] = 1.0;
	}
    }

  // Compute dx, the local mesh spacing, as either the diameter of approximately smallest inscribed circle
  // of the control volume (smallest edge distance time 2) or we assume the control volumes are proportional to
  // an equilateral triangle.

  // inscribed circle.
  for ( i=1; i <= grid->nn; i++ )
    {
      start_index = grid->nnsn[i];
      end_index = grid->nnsn[i+1];

      // The first one is the minimum for now.
      iedge = grid->esn[start_index];
      nodeL = grid->edges[iedge*2+0];
      nodeR = grid->edges[iedge*2+1];

      dist = sqrt( (grid->x[nodeL*NDIM+0] - grid->x[nodeR*NDIM+0])*(grid->x[nodeL*NDIM+0] - grid->x[nodeR*NDIM+0]) +
		       (grid->x[nodeL*NDIM+1] - grid->x[nodeR*NDIM+1])*(grid->x[nodeL*NDIM+1] - grid->x[nodeR*NDIM+1]) );

      dx[i] = dist;

      // Loop over the edges.
      for ( j=(start_index+1); j < end_index; j++ )
	{
	  iedge = grid->esn[j];

	  nodeL = grid->edges[iedge*2+0];
	  nodeR = grid->edges[iedge*2+1];

	  // Compute the distance.
	  dist = sqrt( (grid->x[nodeL*NDIM+0] - grid->x[nodeR*NDIM+0])*(grid->x[nodeL*NDIM+0] - grid->x[nodeR*NDIM+0]) +
		       (grid->x[nodeL*NDIM+1] - grid->x[nodeR*NDIM+1])*(grid->x[nodeL*NDIM+1] - grid->x[nodeR*NDIM+1]) );

	  dx[i] = MIN( dx[i] , dist );
	}
    }

  // Now multiply by two to get diameter.
  for ( i=1; i <= grid->nn; i++ )
    {
      dx[i] = dx[i]*2.0;
    }

  // second method
  /*
  for ( i=1; i <= grid->nn; i++ )
    {
      dx[i] = 2. * sqrt( grid->cv_area[i] / ( 3.*sqrt(3.)) );
    }
  */

  // The general procedure is to loop over all the nodes in the mesh, then for each variable:
  // 1. find the largest negative Q-Q_i and positive Q-Q_i
  // 2. Get the unconstrained extrapolated value at the quadratuer point(s).
  // 3. Compute the maximum value of Phi_ij for the quadrature point(s) j.
  // 4. Set phi(i) = min(Phi_ij).

  // Loop over the edges to find the minimal and maximal values of Qbar in the nodes nearest neighborhood.
  for ( i=1; i <= grid->nedges; i++ )
    {
      nodeL = grid->edges[i*2+0];
      nodeR = grid->edges[i*2+1];

      // Find the minimal/maximal values for each variable.
      for ( j=0; j < NUM_VAR; j++ )
	{
	  // Compare the right node the left node's values and vice versa.
	  dQmin[nodeL*NUM_VAR+j] = MIN(   dQmin[nodeL*NUM_VAR+j] ,
					  grid->nQ[nodeR*NUM_VAR+j] );
	  dQmax[nodeL*NUM_VAR+j] = MAX(   dQmax[nodeL*NUM_VAR+j] ,
					  grid->nQ[nodeR*NUM_VAR+j] );
	  
	  dQmin[nodeR*NUM_VAR+j] = MIN(   dQmin[nodeR*NUM_VAR+j] ,
					  grid->nQ[nodeL*NUM_VAR+j] );
	  dQmax[nodeR*NUM_VAR+j] = MAX(   dQmax[nodeR*NUM_VAR+j] ,
					  grid->nQ[nodeL*NUM_VAR+j] );
	  
	}
    }
  
  // Figure out if we are doing Gaussian Quadrature on the subedges.
  if ( p.order < 3 )
    {
      DO_GAUSS_QUAD = 0;
    }
  else if ( (p.order == 3) && (grid->citer > (p.foits+p.soits)) )
    {
      DO_GAUSS_QUAD = 1;
    }
  else if ( (p.order == 4) && (grid->citer > (p.foits+p.soits)) )
    {
      DO_GAUSS_QUAD = 1;
    }
  
  if ( FULL_GAUSS_FLUX )
    DO_GAUSS_QUAD = 1;
  
  if ( EDGE_BASED )
    DO_GAUSS_QUAD = 0;

  // Now we do an unconstrained reconstruction to every quadrature point and apply the Phi function of
  // Venkatakrishnan.

  for ( isubedge=1; isubedge <= grid->nsubedges; isubedge++ )
    {
      // Extract the real edge this corresponds to.
      iedge = grid->subedges[isubedge*2];

      // Get the left and right nodes.
      nodeL = grid->edges[iedge*2 + 0];
      nodeR = grid->edges[iedge*2 + 1];

      xmid = grid->xm_subedges[isubedge*NDIM+0];
      ymid = grid->xm_subedges[isubedge*NDIM+1];

      // Find the global element.
      gelem = grid->subedges[isubedge*2+1];

      // Get the element centroid.
      xc = grid->el_cent[gelem*2];
      yc = grid->el_cent[gelem*2+1];
      
      if ( CYLINDER && CYLINDER_HACK )
	{
	  if ( (grid->node_order[nodeL] == 1) &&
	       (grid->node_order[nodeR] == 1) )
	    {
	      if ( DO_GAUSS_QUAD == 1 )
		gauss_flag = 1;
	      else
		gauss_flag = 0;
	      
	      DO_GAUSS_QUAD = 0;
	    }
	}

      if ( DO_GAUSS_QUAD )
	{
	  // Calculate the Gauss points.
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
	  
	  // For each gauss point, we will reconstruct the solution.
	  
	  gp[0] = grid->xg_subedges[isubedge*6+0];
	  gp[1] = grid->xg_subedges[isubedge*6+1];
	  Reconstruct_Gauss_UL ( isubedge, nodeL, nodeR, gp, qleft, qright, grid, p, grid->citer );

	  // Do the left node first.
	  eps2 = (K*dx[nodeL]) * (K*dx[nodeL]) * (K*dx[nodeL]);
	  
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      del2 = qleft[j] - grid->nQ[nodeL*NUM_VAR+j];
	      del1min = dQmin[nodeL*NUM_VAR+j] - grid->nQ[nodeL*NUM_VAR+j];
	      del1max = dQmax[nodeL*NUM_VAR+j] - grid->nQ[nodeL*NUM_VAR+j];

	      if ( del2 > 0.0 )
		{
		  temp = (1.0/del2) * ( ( (del1max*del1max + eps2)*del2 + 2.0*del2*del2*del1max ) /
					( del1max*del1max + 2.0*del2*del2 + del1max*del2 + eps ) );
		  
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , temp );
		}
	      else
		{
		  temp = (1.0/del2) * ( ( (del1min*del1min + eps2)*del2 + 2.0*del2*del2*del1min ) /
					( del1min*del1min + 2.0*del2*del2 + del1min*del2 + eps ) );
		  
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , temp );
		}
	    }

	  // Do the right node.
	  eps2 = (K*dx[nodeR]) * (K*dx[nodeR]) * (K*dx[nodeR]);
	  
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      del2 = qright[j] - grid->nQ[nodeR*NUM_VAR+j];
	      del1min = dQmin[nodeR*NUM_VAR+j] - grid->nQ[nodeR*NUM_VAR+j];
	      del1max = dQmax[nodeR*NUM_VAR+j] - grid->nQ[nodeR*NUM_VAR+j];

	      if ( del2 > 0.0 )
		{
		  temp = (1.0/del2) * ( ( (del1max*del1max + eps2)*del2 + 2.0*del2*del2*del1max ) /
					( del1max*del1max + 2.0*del2*del2 + del1max*del2 + eps ) );
		  
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , temp );
		}
	      else
		{
		  temp = (1.0/del2) * ( ( (del1min*del1min + eps2)*del2 + 2.0*del2*del2*del1min ) /
					( del1min*del1min + 2.0*del2*del2 + del1min*del2 + eps ) );
		  
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , temp );
		}
	    }


	  gp[0] = grid->xg_subedges[isubedge*6+2];
	  gp[1] = grid->xg_subedges[isubedge*6+3];
	  Reconstruct_Gauss_UL ( isubedge, nodeL, nodeR, gp, qleft, qright, grid, p, grid->citer );

	  // Do the left node first.
	  eps2 = (K*dx[nodeL]) * (K*dx[nodeL]) * (K*dx[nodeL]);
	  
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      del2 = qleft[j] - grid->nQ[nodeL*NUM_VAR+j];
	      del1min = dQmin[nodeL*NUM_VAR+j] - grid->nQ[nodeL*NUM_VAR+j];
	      del1max = dQmax[nodeL*NUM_VAR+j] - grid->nQ[nodeL*NUM_VAR+j];

	      if ( del2 > 0.0 )
		{
		  temp = (1.0/del2) * ( ( (del1max*del1max + eps2)*del2 + 2.0*del2*del2*del1max ) /
					( del1max*del1max + 2.0*del2*del2 + del1max*del2 + eps ) );
		  
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , temp );
		}
	      else
		{
		  temp = (1.0/del2) * ( ( (del1min*del1min + eps2)*del2 + 2.0*del2*del2*del1min ) /
					( del1min*del1min + 2.0*del2*del2 + del1min*del2 + eps ) );
		  
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , temp );
		}
	    }

	  // Do the right node.
	  eps2 = (K*dx[nodeR]) * (K*dx[nodeR]) * (K*dx[nodeR]);
	  
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      del2 = qright[j] - grid->nQ[nodeR*NUM_VAR+j];
	      del1min = dQmin[nodeR*NUM_VAR+j] - grid->nQ[nodeR*NUM_VAR+j];
	      del1max = dQmax[nodeR*NUM_VAR+j] - grid->nQ[nodeR*NUM_VAR+j];

	      if ( del2 > 0.0 )
		{
		  temp = (1.0/del2) * ( ( (del1max*del1max + eps2)*del2 + 2.0*del2*del2*del1max ) /
					( del1max*del1max + 2.0*del2*del2 + del1max*del2 + eps ) );
		  
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , temp );
		}
	      else
		{
		  temp = (1.0/del2) * ( ( (del1min*del1min + eps2)*del2 + 2.0*del2*del2*del1min ) /
					( del1min*del1min + 2.0*del2*del2 + del1min*del2 + eps ) );
		  
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , temp );
		}
	    }
	  

	  gp[0] = grid->xg_subedges[isubedge*6+4];
	  gp[1] = grid->xg_subedges[isubedge*6+5];
	  Reconstruct_Gauss_UL ( isubedge, nodeL, nodeR, gp, qleft, qright, grid, p, grid->citer );

	  // Do the left node first.
	  eps2 = (K*dx[nodeL]) * (K*dx[nodeL]) * (K*dx[nodeL]);
	  
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      del2 = qleft[j] - grid->nQ[nodeL*NUM_VAR+j];
	      del1min = dQmin[nodeL*NUM_VAR+j] - grid->nQ[nodeL*NUM_VAR+j];
	      del1max = dQmax[nodeL*NUM_VAR+j] - grid->nQ[nodeL*NUM_VAR+j];

	      if ( del2 > 0.0 )
		{
		  temp = (1.0/del2) * ( ( (del1max*del1max + eps2)*del2 + 2.0*del2*del2*del1max ) /
					( del1max*del1max + 2.0*del2*del2 + del1max*del2 + eps ) );
		  
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , temp );
		}
	      else
		{
		  temp = (1.0/del2) * ( ( (del1min*del1min + eps2)*del2 + 2.0*del2*del2*del1min ) /
					( del1min*del1min + 2.0*del2*del2 + del1min*del2 + eps ) );
		  
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , temp );
		}
	    }

	  // Do the right node.
	  eps2 = (K*dx[nodeR]) * (K*dx[nodeR]) * (K*dx[nodeR]);
	  
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      del2 = qright[j] - grid->nQ[nodeR*NUM_VAR+j];
	      del1min = dQmin[nodeR*NUM_VAR+j] - grid->nQ[nodeR*NUM_VAR+j];
	      del1max = dQmax[nodeR*NUM_VAR+j] - grid->nQ[nodeR*NUM_VAR+j];

	      if ( del2 > 0.0 )
		{
		  temp = (1.0/del2) * ( ( (del1max*del1max + eps2)*del2 + 2.0*del2*del2*del1max ) /
					( del1max*del1max + 2.0*del2*del2 + del1max*del2 + eps ) );
		  
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , temp );
		}
	      else
		{
		  temp = (1.0/del2) * ( ( (del1min*del1min + eps2)*del2 + 2.0*del2*del2*del1min ) /
					( del1min*del1min + 2.0*del2*del2 + del1min*del2 + eps ) );
		  
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , temp );
		}
	    }
	  
	}
      else
	{
	  // Get the Q variables and if high order, then reconstruct the solution.
	  Reconstruct_UL( isubedge, nodeL, nodeR, qleft, qright, grid, p, grid->citer);

	  // Do the left node first.
	  eps2 = (K*dx[nodeL]) * (K*dx[nodeL]) * (K*dx[nodeL]);
	  
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      del2 = qleft[j] - grid->nQ[nodeL*NUM_VAR+j];
	      del1min = dQmin[nodeL*NUM_VAR+j] - grid->nQ[nodeL*NUM_VAR+j];
	      del1max = dQmax[nodeL*NUM_VAR+j] - grid->nQ[nodeL*NUM_VAR+j];

	      if ( del2 > 0.0 )
		{
		  temp = (1.0/del2) * ( ( (del1max*del1max + eps2)*del2 + 2.0*del2*del2*del1max ) /
					( del1max*del1max + 2.0*del2*del2 + del1max*del2 + eps ) );
		  
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , temp );
		}
	      else
		{
		  temp = (1.0/del2) * ( ( (del1min*del1min + eps2)*del2 + 2.0*del2*del2*del1min ) /
					( del1min*del1min + 2.0*del2*del2 + del1min*del2 + eps ) );
		  
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , temp );
		}
	    }

	  // Do the right node.
	  eps2 = (K*dx[nodeR]) * (K*dx[nodeR]) * (K*dx[nodeR]);
	  
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      del2 = qright[j] - grid->nQ[nodeR*NUM_VAR+j];
	      del1min = dQmin[nodeR*NUM_VAR+j] - grid->nQ[nodeR*NUM_VAR+j];
	      del1max = dQmax[nodeR*NUM_VAR+j] - grid->nQ[nodeR*NUM_VAR+j];

	      if ( del2 > 0.0 )
		{
		  temp = (1.0/del2) * ( ( (del1max*del1max + eps2)*del2 + 2.0*del2*del2*del1max ) /
					( del1max*del1max + 2.0*del2*del2 + del1max*del2 + eps ) );
		  
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , temp );
		}
	      else
		{
		  temp = (1.0/del2) * ( ( (del1min*del1min + eps2)*del2 + 2.0*del2*del2*del1min ) /
					( del1min*del1min + 2.0*del2*del2 + del1min*del2 + eps ) );
		  
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , temp );
		}
	    }

	}
      
      if ( gauss_flag )
	DO_GAUSS_QUAD = 1;
      
    }
  
  // Set s to the same value as phi.
  for ( i=1; i <= grid->nn; i++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  phi[i*8 + j*2 + 1] = phi[i*8 + j*2 + 0];
	}
    }

  // Go ahead and put the minimum phi in the first spot.
  for ( i=1; i <= grid->nn; i++ )
    {
      temp = phi[i*8 + 0*2 + 0];

      for ( j=1; j < NUM_VAR; j++ )
	{
	  temp = MIN ( temp, phi[i*8 + j*2 + 0] );
	}

      phi[i*8 + 0*2 + 0] = temp;
      phi[i*8 + 0*2 + 1] = temp;
    }

  if ( LIMITER_NEIGHBOR_AVERAGE )  // As implemented, I am only looking at the computed minimum phi.
    {
      double *phi_avg = NULL;                      // For neighbor averaging.
      double eps = 1.0E-2;
      double sum;
      int index1, index2, num_nn;
      
      phi_avg = (double*)malloc( (grid->nn+1)*sizeof(double) );
      if ( phi_avg == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'phi_avg'.\n"); exit(1); }
      
      for ( i=1; i <= grid->nn; i++ )
	{
	  // skip those close enough to 0 or to 1.
	  if ( phi[i*8+0] < (0.0+eps) )
	    {
	      phi_avg[i] = phi[i*8];
	      continue;
	    }

	  if ( phi[i*8+0] > (1.0-eps) )
	    {
	      phi_avg[i] = phi[i*8];
	      continue;
	    }

	  sum = 0.;
	  num_nn = 0;

	  index1 = grid->nnsn[i];
	  index2 = grid->nnsn[i+1];

	  for ( j=index1; j < index2; j++ )
	    {
	      nodeR = grid->nsn[j];  // grab the neighbor.

	      sum += phi[nodeR*8+0];
	      num_nn++;
	    }

	  sum += phi[i*8+0];
	  num_nn++;

	  sum = sum / ( (double)num_nn );
	  phi_avg[i] = sum;
	}

      for ( i=1; i <= grid->nn; i++ )
	{
	  phi[i*8+0] = phi_avg[i];
	  phi[i*8+1] = phi_avg[i];
	}
      
      freenull(phi_avg);
    }

  // Sanity check. All values of phi must be in [0,1]!

  for ( i=1; i <= grid->nn; i++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  if ( ( phi[i*8 + j*2 + 0] < 0.0 ) ||
	       ( phi[i*8 + j*2 + 0] > 1.0 )  )
	    {
	      printf("FATAL ERROR: Phi for node %d and variable %d does not have the proper value: < %.15E  >\n",i,j,phi[i*8+j*2+0]);
	      fflush(stdout);
	      exit(1);
	    }
	}
    }

  freenull(dQmin);
  freenull(dQmax);
  freenull(dx);

  return;
}


//=============================================================
// 
//  compute_limiter_gooch1()
//
//  Computes Venkatakrishnan's modified limiter for every control volume using
//  the second approach (based Van Albada) and computes the high order limiter
//  as described by Nejat in his dissertation.
//
//  GRID *grid;                          // The grid.
//  PARAMS p;                            // The parameters.
//
//=============================================================

void compute_limiter_gooch1 ( GRID *grid, PARAMS p )
{
  int i,j;                                      // Loop counters.
  int isubedge;
  int nodeL,nodeR;
  double xc,yc,xmid,ymid;
  double temp;
  int gelem;
  int iedge;
  int start_index, end_index;
  int gauss_flag = 0;
  double dist;
  double eps = 1.0E-12;               // Prevents division by zero.
  double K = 1.0;                     // Some choices are 0.3, 1.0, 10.0 . Too large means no limiting.
  double eps2;
  double del2, del1max, del1min;
  double *dQmin = NULL;
  double *dQmax = NULL;
  double *dx = NULL;
  double qleft[4];                             // Q on the left of the face.
  double qright[4];                            // Q on the right of the face.

  double S, phi0, sigma, phi_hat;

  // Variables for higher order flux calculation.
  double x1,x2,x3;                       // Integration points.
  double y1,y2,y3;
  double t1,t2,t3;                       // Gaussian roots.
  double gp[2];                          // Gauss point.
  int DO_GAUSS_QUAD = 0;

  double *phi = NULL;

  // tunable parameters.
  phi0 = 0.8;
  S = 20.0;

  if ( grid->phi == NULL )
    {
      grid->phi = (double*)malloc( (grid->nn+1)*2*4*sizeof(double) );
      if ( grid->phi == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'phi'.\n"); exit(1); }
    }

  // The 8 in the memory allocation is for phi and s for all variables.

  phi = grid->phi;

  dQmin = (double*)malloc( (grid->nn+1)*NUM_VAR*sizeof(double));
  dQmax = (double*)malloc( (grid->nn+1)*NUM_VAR*sizeof(double));
  dx = (double*)malloc( (grid->nn+1)*sizeof(double));

  if ( dQmin == NULL )
    { printf("MEMORY ERROR: COULD NOT ALLOCATE 'dQmin'.\n"); exit(1); }

  if ( dQmax == NULL )
    { printf("MEMORY ERROR: COULD NOT ALLOCATE 'dQmax'.\n"); exit(1); }

  if ( dx == NULL )
    { printf("MEMORY ERROR: COULD NOT ALLOCATE 'dx'.\n"); exit(1); }

  // Clean out the memory.
  for ( i=1; i <= grid->nn; i++ )
    {
      dx[i] = 0.0;

      for ( j=0; j < 8; j++ )
	phi[i*8+j] = 0.0;

      // Initialize these to 'high' values.
      for ( j=0; j < NUM_VAR; j++ )
	{
	  dQmin[i*NUM_VAR+j] = 0.0;
	  dQmax[i*NUM_VAR+j] = 0.0;
	}
    }

  // Initialize the maximum and minimum CV averages to be the current CV value.
  for ( i=1; i <= grid->nn; i++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  dQmin[i*NUM_VAR+j] = grid->nQ[i*NUM_VAR+j];
	  dQmax[i*NUM_VAR+j] = grid->nQ[i*NUM_VAR+j];
	}
    }

  // Initialize the limiter to 1.0
  for ( i=1; i <= grid->nn; i++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  phi[i*8+j*2+0] = 1.0;
	}
    }

  // Compute dx, the local mesh spacing, as either the diameter of approximately smallest inscribed circle
  // of the control volume (smallest edge distance time 2) or we assume the control volumes are proportional to
  // an equilateral triangle.

  // inscribed circle.
  for ( i=1; i <= grid->nn; i++ )
    {
      start_index = grid->nnsn[i];
      end_index = grid->nnsn[i+1];

      // The first one is the minimum for now.
      iedge = grid->esn[start_index];
      nodeL = grid->edges[iedge*2+0];
      nodeR = grid->edges[iedge*2+1];

      dist = sqrt( (grid->x[nodeL*NDIM+0] - grid->x[nodeR*NDIM+0])*(grid->x[nodeL*NDIM+0] - grid->x[nodeR*NDIM+0]) +
		       (grid->x[nodeL*NDIM+1] - grid->x[nodeR*NDIM+1])*(grid->x[nodeL*NDIM+1] - grid->x[nodeR*NDIM+1]) );

      dx[i] = dist;

      // Loop over the edges.
      for ( j=(start_index+1); j < end_index; j++ )
	{
	  iedge = grid->esn[j];

	  nodeL = grid->edges[iedge*2+0];
	  nodeR = grid->edges[iedge*2+1];

	  // Compute the distance.
	  dist = sqrt( (grid->x[nodeL*NDIM+0] - grid->x[nodeR*NDIM+0])*(grid->x[nodeL*NDIM+0] - grid->x[nodeR*NDIM+0]) +
		       (grid->x[nodeL*NDIM+1] - grid->x[nodeR*NDIM+1])*(grid->x[nodeL*NDIM+1] - grid->x[nodeR*NDIM+1]) );

	  dx[i] = MIN( dx[i] , dist );
	}
    }

  // Now multiply by two to get diameter.
  for ( i=1; i <= grid->nn; i++ )
    {
      dx[i] = dx[i]*2.0;
    }

  // second method
  /*
  for ( i=1; i <= grid->nn; i++ )
    {
      dx[i] = 2. * sqrt( grid->cv_area[i] / ( 3.*sqrt(3.)) );
    }
  */

  // The general procedure is to loop over all the nodes in the mesh, then for each variable:
  // 1. find the largest negative Q-Q_i and positive Q-Q_i
  // 2. Get the unconstrained extrapolated value at the quadratuer point(s).
  // 3. Compute the maximum value of Phi_ij for the quadrature point(s) j.
  // 4. Set phi(i) = min(Phi_ij).

  // Loop over the edges to find the minimal and maximal values of Qbar in the nodes nearest neighborhood.
  for ( i=1; i <= grid->nedges; i++ )
    {
      nodeL = grid->edges[i*2+0];
      nodeR = grid->edges[i*2+1];

      // Find the minimal/maximal values for each variable.
      for ( j=0; j < NUM_VAR; j++ )
	{
	  // Compare the right node the left node's values and vice versa.
	  dQmin[nodeL*NUM_VAR+j] = MIN(   dQmin[nodeL*NUM_VAR+j] ,
					  grid->nQ[nodeR*NUM_VAR+j] );
	  dQmax[nodeL*NUM_VAR+j] = MAX(   dQmax[nodeL*NUM_VAR+j] ,
					  grid->nQ[nodeR*NUM_VAR+j] );
	  
	  dQmin[nodeR*NUM_VAR+j] = MIN(   dQmin[nodeR*NUM_VAR+j] ,
					  grid->nQ[nodeL*NUM_VAR+j] );
	  dQmax[nodeR*NUM_VAR+j] = MAX(   dQmax[nodeR*NUM_VAR+j] ,
					  grid->nQ[nodeL*NUM_VAR+j] );
	  
	}
    }
  
  // Figure out if we are doing Gaussian Quadrature on the subedges.
  if ( p.order < 3 )
    {
      DO_GAUSS_QUAD = 0;
    }
  else if ( (p.order == 3) && (grid->citer > (p.foits+p.soits)) )
    {
      DO_GAUSS_QUAD = 1;
    }
  else if ( (p.order == 4) && (grid->citer > (p.foits+p.soits)) )
    {
      DO_GAUSS_QUAD = 1;
    }
  
  if ( FULL_GAUSS_FLUX )
    DO_GAUSS_QUAD = 1;
  
  if ( EDGE_BASED )
    DO_GAUSS_QUAD = 0;

  // Now we do an unconstrained reconstruction to every quadrature point and apply the Phi function of
  // Venkatakrishnan.

  for ( isubedge=1; isubedge <= grid->nsubedges; isubedge++ )
    {
      // Extract the real edge this corresponds to.
      iedge = grid->subedges[isubedge*2];

      // Get the left and right nodes.
      nodeL = grid->edges[iedge*2 + 0];
      nodeR = grid->edges[iedge*2 + 1];

      xmid = grid->xm_subedges[isubedge*NDIM+0];
      ymid = grid->xm_subedges[isubedge*NDIM+1];

      // Find the global element.
      gelem = grid->subedges[isubedge*2+1];

      // Get the element centroid.
      xc = grid->el_cent[gelem*2];
      yc = grid->el_cent[gelem*2+1];

      if ( CYLINDER && CYLINDER_HACK )
	{
	  if ( (grid->node_order[nodeL] == 1) &&
	       (grid->node_order[nodeR] == 1) )
	    {
	      if ( DO_GAUSS_QUAD == 1 )
		gauss_flag = 1;
	      else
		gauss_flag = 0;

	      DO_GAUSS_QUAD = 0;
	    }
	}

      if ( DO_GAUSS_QUAD )
	{
	  // Calculate the Gauss points.
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
	  
	  // For each gauss point, we will reconstruct the solution.
	  
	  gp[0] = grid->xg_subedges[isubedge*6+0];
	  gp[1] = grid->xg_subedges[isubedge*6+1];
	  Reconstruct_Gauss_UL ( isubedge, nodeL, nodeR, gp, qleft, qright, grid, p, grid->citer );

	  // Do the left node first.
	  eps2 = (K*dx[nodeL]) * (K*dx[nodeL]) * (K*dx[nodeL]);
	  
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      del2 = qleft[j] - grid->nQ[nodeL*NUM_VAR+j];
	      del1min = dQmin[nodeL*NUM_VAR+j] - grid->nQ[nodeL*NUM_VAR+j];
	      del1max = dQmax[nodeL*NUM_VAR+j] - grid->nQ[nodeL*NUM_VAR+j];

	      if ( del2 > 0.0 )
		{
		  temp = (1.0/del2) * ( ( (del1max*del1max + eps2)*del2 + 2.0*del2*del2*del1max ) /
					( del1max*del1max + 2.0*del2*del2 + del1max*del2 + eps ) );
		  
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , temp );
		}
	      else
		{
		  temp = (1.0/del2) * ( ( (del1min*del1min + eps2)*del2 + 2.0*del2*del2*del1min ) /
					( del1min*del1min + 2.0*del2*del2 + del1min*del2 + eps ) );
		  
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , temp );
		}
	    }

	  // Do the right node.
	  eps2 = (K*dx[nodeR]) * (K*dx[nodeR]) * (K*dx[nodeR]);
	  
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      del2 = qright[j] - grid->nQ[nodeR*NUM_VAR+j];
	      del1min = dQmin[nodeR*NUM_VAR+j] - grid->nQ[nodeR*NUM_VAR+j];
	      del1max = dQmax[nodeR*NUM_VAR+j] - grid->nQ[nodeR*NUM_VAR+j];

	      if ( del2 > 0.0 )
		{
		  temp = (1.0/del2) * ( ( (del1max*del1max + eps2)*del2 + 2.0*del2*del2*del1max ) /
					( del1max*del1max + 2.0*del2*del2 + del1max*del2 + eps ) );
		  
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , temp );
		}
	      else
		{
		  temp = (1.0/del2) * ( ( (del1min*del1min + eps2)*del2 + 2.0*del2*del2*del1min ) /
					( del1min*del1min + 2.0*del2*del2 + del1min*del2 + eps ) );
		  
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , temp );
		}
	    }


	  gp[0] = grid->xg_subedges[isubedge*6+2];
	  gp[1] = grid->xg_subedges[isubedge*6+3];
	  Reconstruct_Gauss_UL ( isubedge, nodeL, nodeR, gp, qleft, qright, grid, p, grid->citer );

	  // Do the left node first.
	  eps2 = (K*dx[nodeL]) * (K*dx[nodeL]) * (K*dx[nodeL]);
	  
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      del2 = qleft[j] - grid->nQ[nodeL*NUM_VAR+j];
	      del1min = dQmin[nodeL*NUM_VAR+j] - grid->nQ[nodeL*NUM_VAR+j];
	      del1max = dQmax[nodeL*NUM_VAR+j] - grid->nQ[nodeL*NUM_VAR+j];

	      if ( del2 > 0.0 )
		{
		  temp = (1.0/del2) * ( ( (del1max*del1max + eps2)*del2 + 2.0*del2*del2*del1max ) /
					( del1max*del1max + 2.0*del2*del2 + del1max*del2 + eps ) );
		  
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , temp );
		}
	      else
		{
		  temp = (1.0/del2) * ( ( (del1min*del1min + eps2)*del2 + 2.0*del2*del2*del1min ) /
					( del1min*del1min + 2.0*del2*del2 + del1min*del2 + eps ) );
		  
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , temp );
		}
	    }

	  // Do the right node.
	  eps2 = (K*dx[nodeR]) * (K*dx[nodeR]) * (K*dx[nodeR]);
	  
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      del2 = qright[j] - grid->nQ[nodeR*NUM_VAR+j];
	      del1min = dQmin[nodeR*NUM_VAR+j] - grid->nQ[nodeR*NUM_VAR+j];
	      del1max = dQmax[nodeR*NUM_VAR+j] - grid->nQ[nodeR*NUM_VAR+j];

	      if ( del2 > 0.0 )
		{
		  temp = (1.0/del2) * ( ( (del1max*del1max + eps2)*del2 + 2.0*del2*del2*del1max ) /
					( del1max*del1max + 2.0*del2*del2 + del1max*del2 + eps ) );
		  
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , temp );
		}
	      else
		{
		  temp = (1.0/del2) * ( ( (del1min*del1min + eps2)*del2 + 2.0*del2*del2*del1min ) /
					( del1min*del1min + 2.0*del2*del2 + del1min*del2 + eps ) );
		  
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , temp );
		}
	    }
	  

	  gp[0] = grid->xg_subedges[isubedge*6+4];
	  gp[1] = grid->xg_subedges[isubedge*6+5];
	  Reconstruct_Gauss_UL ( isubedge, nodeL, nodeR, gp, qleft, qright, grid, p, grid->citer );

	  // Do the left node first.
	  eps2 = (K*dx[nodeL]) * (K*dx[nodeL]) * (K*dx[nodeL]);
	  
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      del2 = qleft[j] - grid->nQ[nodeL*NUM_VAR+j];
	      del1min = dQmin[nodeL*NUM_VAR+j] - grid->nQ[nodeL*NUM_VAR+j];
	      del1max = dQmax[nodeL*NUM_VAR+j] - grid->nQ[nodeL*NUM_VAR+j];

	      if ( del2 > 0.0 )
		{
		  temp = (1.0/del2) * ( ( (del1max*del1max + eps2)*del2 + 2.0*del2*del2*del1max ) /
					( del1max*del1max + 2.0*del2*del2 + del1max*del2 + eps ) );
		  
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , temp );
		}
	      else
		{
		  temp = (1.0/del2) * ( ( (del1min*del1min + eps2)*del2 + 2.0*del2*del2*del1min ) /
					( del1min*del1min + 2.0*del2*del2 + del1min*del2 + eps ) );
		  
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , temp );
		}
	    }

	  // Do the right node.
	  eps2 = (K*dx[nodeR]) * (K*dx[nodeR]) * (K*dx[nodeR]);
	  
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      del2 = qright[j] - grid->nQ[nodeR*NUM_VAR+j];
	      del1min = dQmin[nodeR*NUM_VAR+j] - grid->nQ[nodeR*NUM_VAR+j];
	      del1max = dQmax[nodeR*NUM_VAR+j] - grid->nQ[nodeR*NUM_VAR+j];

	      if ( del2 > 0.0 )
		{
		  temp = (1.0/del2) * ( ( (del1max*del1max + eps2)*del2 + 2.0*del2*del2*del1max ) /
					( del1max*del1max + 2.0*del2*del2 + del1max*del2 + eps ) );
		  
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , temp );
		}
	      else
		{
		  temp = (1.0/del2) * ( ( (del1min*del1min + eps2)*del2 + 2.0*del2*del2*del1min ) /
					( del1min*del1min + 2.0*del2*del2 + del1min*del2 + eps ) );
		  
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , temp );
		}
	    }
	  
	}
      else
	{
	  // Get the Q variables and if high order, then reconstruct the solution.
	  Reconstruct_UL( isubedge, nodeL, nodeR, qleft, qright, grid, p, grid->citer);

	  // Do the left node first.
	  eps2 = (K*dx[nodeL]) * (K*dx[nodeL]) * (K*dx[nodeL]);
	  
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      del2 = qleft[j] - grid->nQ[nodeL*NUM_VAR+j];
	      del1min = dQmin[nodeL*NUM_VAR+j] - grid->nQ[nodeL*NUM_VAR+j];
	      del1max = dQmax[nodeL*NUM_VAR+j] - grid->nQ[nodeL*NUM_VAR+j];

	      if ( del2 > 0.0 )
		{
		  temp = (1.0/del2) * ( ( (del1max*del1max + eps2)*del2 + 2.0*del2*del2*del1max ) /
					( del1max*del1max + 2.0*del2*del2 + del1max*del2 + eps ) );
		  
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , temp );
		}
	      else
		{
		  temp = (1.0/del2) * ( ( (del1min*del1min + eps2)*del2 + 2.0*del2*del2*del1min ) /
					( del1min*del1min + 2.0*del2*del2 + del1min*del2 + eps ) );
		  
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , temp );
		}
	    }

	  // Do the right node.
	  eps2 = (K*dx[nodeR]) * (K*dx[nodeR]) * (K*dx[nodeR]);
	  
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      del2 = qright[j] - grid->nQ[nodeR*NUM_VAR+j];
	      del1min = dQmin[nodeR*NUM_VAR+j] - grid->nQ[nodeR*NUM_VAR+j];
	      del1max = dQmax[nodeR*NUM_VAR+j] - grid->nQ[nodeR*NUM_VAR+j];

	      if ( del2 > 0.0 )
		{
		  temp = (1.0/del2) * ( ( (del1max*del1max + eps2)*del2 + 2.0*del2*del2*del1max ) /
					( del1max*del1max + 2.0*del2*del2 + del1max*del2 + eps ) );
		  
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , temp );
		}
	      else
		{
		  temp = (1.0/del2) * ( ( (del1min*del1min + eps2)*del2 + 2.0*del2*del2*del1min ) /
					( del1min*del1min + 2.0*del2*del2 + del1min*del2 + eps ) );
		  
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , temp );
		}
	    }

	}

      if ( gauss_flag )
	DO_GAUSS_QUAD = 1;
      
    }
  
  if ( LIMITER_NEIGHBOR_AVERAGE )  // As implemented here, I am doing this for the values of phi for all variables.
    {
      double *phi_avg = NULL;                      // For neighbor averaging.
      double epsPHI = 1.0E-2;
      double sum;
      int index1, index2, num_nn;
      int ind;
      
      phi_avg = (double*)malloc( (grid->nn+1)*4*sizeof(double) );
      if ( phi_avg == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'phi_avg'.\n"); exit(1); }
      
      for ( i=1; i <= grid->nn; i++ )
	{
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      // skip those close enough to 0 or to 1.
	      if ( phi[i*8+j*2+0] < (0.0+epsPHI) )
		{
		  phi_avg[i*4+j] = phi[i*8+j*2+0];
		  continue;
		}
	      
	      if ( phi[i*8+j*2+0] > (1.0-epsPHI) )
		{
		  phi_avg[i*4+j] = phi[i*8+j*2+0];
		  continue;
		}
	      
	      sum = 0.;
	      num_nn = 0;
	      
	      index1 = grid->nnsn[i];
	      index2 = grid->nnsn[i+1];
	      
	      for ( ind=index1; ind < index2; ind++ )
		{
		  nodeR = grid->nsn[ind];  // grab the neighbor.
		  
		  sum += phi[nodeR*8+j*2+0];
		  num_nn++;
		}

	      sum += phi[i*8+j*2+0];
	      num_nn++;
	      
	      sum = sum / ( (double)num_nn );
	      phi_avg[i*4+j] = sum;
	    }
	}

      for ( i=1; i <= grid->nn; i++ )
	{
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      phi[i*8+j*2+0] = phi_avg[i*4+j];
	    }
	}      

      freenull(phi_avg);
    }

  // We have computed the values of phi, so now compute the values of sigma and phi^hat.
  for ( i=1; i <= grid->nn; i++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  phi_hat = phi[i*8 + j*2 + 0];

	  sigma = ( 1.0 - tanh( (S*(phi0 - phi_hat) ) ) )*0.5;

	  phi[i*8 + j*2 + 0] = (1.0 - sigma)*phi_hat + sigma;
	  phi[i*8 + j*2 + 1] = sigma;
	}
    }

  // Go ahead and put the minimum phi in the first spot.
  for ( i=1; i <= grid->nn; i++ )
    {
      phi_hat = phi[i*8 + 0*2 + 0];
      sigma = phi[i*8 + 0*2 +1];

      for ( j=1; j < NUM_VAR; j++ )
	{
	  if ( phi[i*8 + j*2 + 0] < phi_hat )
	    {
	      phi_hat = phi[i*8 + 0*2 + 0];
	      sigma = phi[i*8 + 0*2 +1];
	    }
	}
      
      phi[i*8 + 0*2 + 0] = phi_hat;
      phi[i*8 + 0*2 + 1] = sigma;
    }

  // Sanity check. All values of phi must be in [0,1]!

  for ( i=1; i <= grid->nn; i++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  if ( ( phi[i*8 + j*2 + 0] < 0.0 ) ||
	       ( phi[i*8 + j*2 + 0] > 1.0 )  )
	    {
	      printf("FATAL ERROR: Phi for node %d and variable %d does not have the proper value: < %.15E  >\n",i,j,phi[i*8+j*2+0]);
	      fflush(stdout);
	      exit(1);
	    }
	}
    }

  freenull(dQmin);
  freenull(dQmax);
  freenull(dx);

  return;
}



//=============================================================
// 
//  compute_limiter_gooch2()
//
//  Computes the modified new limiter as described by Carl
//  Ollivier-Gooch's team in the literature.
//
//  GRID *grid;                          // The grid.
//  PARAMS p;                            // The parameters.
//
//=============================================================

void compute_limiter_gooch2 ( GRID *grid, PARAMS p )
{
  int i,j;                                      // Loop counters.
  int isubedge;
  int nodeL,nodeR;
  int start_index, end_index;
  double dist;
  double eps = 1.0E-12;               // Prevents division by zero.
  double xc,yc,xmid,ymid;
  double temp;
  int gelem;
  int iedge;
  int gauss_flag = 0;
  double *dQmin = NULL;
  double *dQmax = NULL;
  double *dx = NULL;
  double qleft[4];                             // Q on the left of the face.
  double qright[4];                            // Q on the right of the face.

  // Specific variables to this method.
  double Y1,Y2,Yt;                       // Function values for new min.
  double delta_x;                        // A measure of local grid spacing.
  double K;                              // Tunable parameter.
  double A,B,C,D;
  double area;
  double delta_u, delta_u2;
  double bound,sigma;

  // Variables for higher order flux calculation.
  double x1,x2,x3;                       // Integration points.
  double y1,y2,y3;
  double t1,t2,t3;                       // Gaussian roots.
  double gp[2];                          // Gauss point.
  int DO_GAUSS_QUAD = 0;

  double *phi = NULL;


  // Parameters.
  K = 1.;
  Yt = 1.5;

  A = -4. / 27.;
  B = 0.;
  C = 1.;
  D = 0.;

  /*
  Yt = 1.75;
  A = -16.0 / 343.0;
  B = -8.0 / 49.0;
  C = 1.0;
  D = 0.0;
  */

  if ( grid->phi == NULL )
    {
      grid->phi = (double*)malloc( (grid->nn+1)*2*4*sizeof(double) );
      if ( grid->phi == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'phi'.\n"); exit(1); }
    }

  // The 8 in the memory allocation is for phi and s for all variables.

  phi = grid->phi;

  dQmin = (double*)malloc( (grid->nn+1)*NUM_VAR*sizeof(double));
  dQmax = (double*)malloc( (grid->nn+1)*NUM_VAR*sizeof(double));
  dx = (double*)malloc( (grid->nn+1)*sizeof(double));

  if ( dQmin == NULL )
    { printf("MEMORY ERROR: COULD NOT ALLOCATE 'dQmin'.\n"); exit(1); }

  if ( dQmax == NULL )
    { printf("MEMORY ERROR: COULD NOT ALLOCATE 'dQmax'.\n"); exit(1); }
  
  if ( dx == NULL )
    { printf("MEMORY ERROR: COULD NOT ALLOCATE 'dx'.\n"); exit(1); }

 // Clean out the memory.
  for ( i=1; i <= grid->nn; i++ )
    {
      dx[i] = 0.0;

      for ( j=0; j < 8; j++ )
	phi[i*8+j] = 0.0;

      // Initialize these to 'high' values.
      for ( j=0; j < NUM_VAR; j++ )
	{
	  dQmin[i*NUM_VAR+j] = 0.0;
	  dQmax[i*NUM_VAR+j] = 0.0;
	}
    }

  // Initialize the maximum and minimum CV averages to be the current CV value.
  for ( i=1; i <= grid->nn; i++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  dQmin[i*NUM_VAR+j] = grid->nQ[i*NUM_VAR+j];
	  dQmax[i*NUM_VAR+j] = grid->nQ[i*NUM_VAR+j];
	}
    }

  // Initialize the limiter to 1.0
  for ( i=1; i <= grid->nn; i++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  phi[i*8+j*2+0] = 1.0;
	}
    }

  // Compute dx, the local mesh spacing, as either the diameter of approximately smallest inscribed circle
  // of the control volume (smallest edge distance time 2) or we assume the control volumes are proportional to
  // an equilateral triangle.

  // inscribed circle.
  for ( i=1; i <= grid->nn; i++ )
    {
      start_index = grid->nnsn[i];
      end_index = grid->nnsn[i+1];

      // The first one is the minimum for now.
      iedge = grid->esn[start_index];
      nodeL = grid->edges[iedge*2+0];
      nodeR = grid->edges[iedge*2+1];

      dist = sqrt( (grid->x[nodeL*NDIM+0] - grid->x[nodeR*NDIM+0])*(grid->x[nodeL*NDIM+0] - grid->x[nodeR*NDIM+0]) +
		       (grid->x[nodeL*NDIM+1] - grid->x[nodeR*NDIM+1])*(grid->x[nodeL*NDIM+1] - grid->x[nodeR*NDIM+1]) );

      dx[i] = dist;

      // Loop over the edges.
      for ( j=(start_index+1); j < end_index; j++ )
	{
	  iedge = grid->esn[j];

	  nodeL = grid->edges[iedge*2+0];
	  nodeR = grid->edges[iedge*2+1];

	  // Compute the distance.
	  dist = sqrt( (grid->x[nodeL*NDIM+0] - grid->x[nodeR*NDIM+0])*(grid->x[nodeL*NDIM+0] - grid->x[nodeR*NDIM+0]) +
		       (grid->x[nodeL*NDIM+1] - grid->x[nodeR*NDIM+1])*(grid->x[nodeL*NDIM+1] - grid->x[nodeR*NDIM+1]) );

	  dx[i] = MIN( dx[i] , dist );
	}
    }

  // Now multiply by two to get diameter.
  for ( i=1; i <= grid->nn; i++ )
    {
      dx[i] = dx[i]*2.0;
    }

  // second method
  /*
  for ( i=1; i <= grid->nn; i++ )
    {
      dx[i] = 2. * sqrt( grid->cv_area[i] / ( 3.*sqrt(3.)) );
    }
  */

  // The general procedure is to loop over all the nodes in the mesh, then for each variable:
  // 1. find the largest negative Q-Q_i and positive Q-Q_i
  // 2. Get the unconstrained extrapolated value at the quadrature point(s).
  // 3. Compute the maximum value of Phi_ij for the quadrature point(s) j.
  // 4. Set phi(i) = min(Phi_ij).

  // Loop over the edges to find the minimal and maximal values of Qbar in the nodes nearest neighborhood.
  for ( i=1; i <= grid->nedges; i++ )
    {
      nodeL = grid->edges[i*2+0];
      nodeR = grid->edges[i*2+1];

      // Find the minimal/maximal values for each variable.
      for ( j=0; j < NUM_VAR; j++ )
	{
	  // Compare the right node the left node's values and vice versa.
	  dQmin[nodeL*NUM_VAR+j] = MIN(   dQmin[nodeL*NUM_VAR+j] ,
					  grid->nQ[nodeR*NUM_VAR+j] );
	  dQmax[nodeL*NUM_VAR+j] = MAX(   dQmax[nodeL*NUM_VAR+j] ,
					  grid->nQ[nodeR*NUM_VAR+j] );
	  
	  dQmin[nodeR*NUM_VAR+j] = MIN(   dQmin[nodeR*NUM_VAR+j] ,
					  grid->nQ[nodeL*NUM_VAR+j] );
	  dQmax[nodeR*NUM_VAR+j] = MAX(   dQmax[nodeR*NUM_VAR+j] ,
					  grid->nQ[nodeL*NUM_VAR+j] );
	  
	}
    }


  // Now we do an unconstrained reconstruction to every quadrature point and apply the Phi function.

  // This will return the intermediate phi that can be used in conjuction with a switch to remove
  // the limiter in smooth regions.

  // for yt=3/2 : P(y) = -4/27 y^3 + y

  // Figure out if we are doing Gaussian Quadrature on the subedges.
  if ( p.order < 3 )
    {
      DO_GAUSS_QUAD = 0;
    }
  else if ( (p.order == 3) && (grid->citer > (p.foits+p.soits)) )
    {
      DO_GAUSS_QUAD = 1;
    }
  else if ( (p.order == 4) && (grid->citer > (p.foits+p.soits)) )
    {
      DO_GAUSS_QUAD = 1;
    }
  
  if ( FULL_GAUSS_FLUX )
    DO_GAUSS_QUAD = 1;
  
  if ( EDGE_BASED )
    DO_GAUSS_QUAD = 0;

  for ( isubedge=1; isubedge <= grid->nsubedges; isubedge++ )
    {
      // Extract the real edge this corresponds to.
      iedge = grid->subedges[isubedge*2];

      // Get the left and right nodes.
      nodeL = grid->edges[iedge*2 + 0];
      nodeR = grid->edges[iedge*2 + 1];

      xmid = grid->xm_subedges[isubedge*NDIM+0];
      ymid = grid->xm_subedges[isubedge*NDIM+1];

      // Find the global element.
      gelem = grid->subedges[isubedge*2+1];

      // Get the element centroid.
      xc = grid->el_cent[gelem*2];
      yc = grid->el_cent[gelem*2+1];

      if ( CYLINDER && CYLINDER_HACK )
	{
	  if ( (grid->node_order[nodeL] == 1) &&
	       (grid->node_order[nodeR] == 1) )
	    {
	      if ( DO_GAUSS_QUAD == 1 )
		gauss_flag = 1;
	      else
		gauss_flag = 0;

	      DO_GAUSS_QUAD = 0;
	    }
	}

      if ( DO_GAUSS_QUAD )
	{
	  // Calculate the Gauss points.
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
	  
	  // For each gauss point, we will reconstruct the solution.
	  
	  gp[0] = grid->xg_subedges[isubedge*6+0];
	  gp[1] = grid->xg_subedges[isubedge*6+1];
	  Reconstruct_Gauss_UL ( isubedge, nodeL, nodeR, gp, qleft, qright, grid, p, grid->citer );

	  // Do the left node first.
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      // set y1 and y2.
	      Y1 = ( dQmax[nodeL*NUM_VAR+j] - grid->nQ[nodeL*NUM_VAR+j] ) / (qleft[j] - grid->nQ[nodeL*NUM_VAR+j]);
	      Y2 = ( dQmin[nodeL*NUM_VAR+j] - grid->nQ[nodeL*NUM_VAR+j] ) / (qleft[j] - grid->nQ[nodeL*NUM_VAR+j]);

	      if ( ( qleft[j] - grid->nQ[nodeL*NUM_VAR+j] ) > 0. )
		{
		  // min (1,y1)
		  if ( Y1 < Yt )
		    temp = A*Y1*Y1*Y1 + B*Y1*Y1 + C*Y1 +D;
		  else
		    temp = 1.;
		  
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , temp );
		}
	      else if ( ( qleft[j] - grid->nQ[nodeL*NUM_VAR+j] ) < 0. )
		{
		  // min (1,y2)
		  if ( Y2 < Yt )
		    temp = A*Y2*Y2*Y2 + B*Y2*Y2 + C*Y2 +D;
		  else
		    temp = 1.;

		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , temp );
		}
	      else
		{
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , 1.0 );
		}
	    }

	  // Do the right node.
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      // set y1 and y2
	      Y1 = ( dQmax[nodeR*NUM_VAR+j] - grid->nQ[nodeR*NUM_VAR+j] ) / (qright[j] - grid->nQ[nodeR*NUM_VAR+j]);
	      Y2 = ( dQmin[nodeR*NUM_VAR+j] - grid->nQ[nodeR*NUM_VAR+j] ) / (qright[j] - grid->nQ[nodeR*NUM_VAR+j]);

	      if ( ( qright[j] - grid->nQ[nodeR*NUM_VAR+j] ) > 0. )
		{
		  // min (1,y1)
		  if ( Y1 < Yt )
		    temp = A*Y1*Y1*Y1 + B*Y1*Y1 + C*Y1 +D;
		  else
		    temp = 1.;
		  
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , temp );
		}
	      else if ( ( qright[j] - grid->nQ[nodeR*NUM_VAR+j] ) < 0. )
		{
		  // min (1,y2)
		  if ( Y2 < Yt )
		    temp = A*Y2*Y2*Y2 + B*Y2*Y2 + C*Y2 +D;
		  else
		    temp = 1.;
		  
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , temp );
		}
	      else
		{
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , 1.0 );
		}
	    }


	  gp[0] = grid->xg_subedges[isubedge*6+2];
	  gp[1] = grid->xg_subedges[isubedge*6+3];
	  Reconstruct_Gauss_UL ( isubedge, nodeL, nodeR, gp, qleft, qright, grid, p, grid->citer );

	  // Do the left node first.
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      // set y1 and y2.
	      Y1 = ( dQmax[nodeL*NUM_VAR+j] - grid->nQ[nodeL*NUM_VAR+j] ) / (qleft[j] - grid->nQ[nodeL*NUM_VAR+j]);
	      Y2 = ( dQmin[nodeL*NUM_VAR+j] - grid->nQ[nodeL*NUM_VAR+j] ) / (qleft[j] - grid->nQ[nodeL*NUM_VAR+j]);

	      if ( ( qleft[j] - grid->nQ[nodeL*NUM_VAR+j] ) > 0. )
		{
		  // min (1,y1)
		  if ( Y1 < Yt )
		    temp = A*Y1*Y1*Y1 + B*Y1*Y1 + C*Y1 +D;
		  else
		    temp = 1.;
		  
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , temp );
		}
	      else if ( ( qleft[j] - grid->nQ[nodeL*NUM_VAR+j] ) < 0. )
		{
		  // min (1,y2)
		  if ( Y2 < Yt )
		    temp = A*Y2*Y2*Y2 + B*Y2*Y2 + C*Y2 +D;
		  else
		    temp = 1.;

		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , temp );
		}
	      else
		{
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , 1.0 );
		}
	    }

	  // Do the right node.
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      // set y1 and y2
	      Y1 = ( dQmax[nodeR*NUM_VAR+j] - grid->nQ[nodeR*NUM_VAR+j] ) / (qright[j] - grid->nQ[nodeR*NUM_VAR+j]);
	      Y2 = ( dQmin[nodeR*NUM_VAR+j] - grid->nQ[nodeR*NUM_VAR+j] ) / (qright[j] - grid->nQ[nodeR*NUM_VAR+j]);

	      if ( ( qright[j] - grid->nQ[nodeR*NUM_VAR+j] ) > 0. )
		{
		  // min (1,y1)
		  if ( Y1 < Yt )
		    temp = A*Y1*Y1*Y1 + B*Y1*Y1 + C*Y1 +D;
		  else
		    temp = 1.;
		  
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , temp );
		}
	      else if ( ( qright[j] - grid->nQ[nodeR*NUM_VAR+j] ) < 0. )
		{
		  // min (1,y2)
		  if ( Y2 < Yt )
		    temp = A*Y2*Y2*Y2 + B*Y2*Y2 + C*Y2 +D;
		  else
		    temp = 1.;
		  
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , temp );
		}
	      else
		{
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , 1.0 );
		}
	    }

	  gp[0] = grid->xg_subedges[isubedge*6+4];
	  gp[1] = grid->xg_subedges[isubedge*6+5];
	  Reconstruct_Gauss_UL ( isubedge, nodeL, nodeR, gp, qleft, qright, grid, p, grid->citer );

	  // Do the left node first.
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      // set y1 and y2.
	      Y1 = ( dQmax[nodeL*NUM_VAR+j] - grid->nQ[nodeL*NUM_VAR+j] ) / (qleft[j] - grid->nQ[nodeL*NUM_VAR+j]);
	      Y2 = ( dQmin[nodeL*NUM_VAR+j] - grid->nQ[nodeL*NUM_VAR+j] ) / (qleft[j] - grid->nQ[nodeL*NUM_VAR+j]);

	      if ( ( qleft[j] - grid->nQ[nodeL*NUM_VAR+j] ) > 0. )
		{
		  // min (1,y1)
		  if ( Y1 < Yt )
		    temp = A*Y1*Y1*Y1 + B*Y1*Y1 + C*Y1 +D;
		  else
		    temp = 1.;
		  
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , temp );
		}
	      else if ( ( qleft[j] - grid->nQ[nodeL*NUM_VAR+j] ) < 0. )
		{
		  // min (1,y2)
		  if ( Y2 < Yt )
		    temp = A*Y2*Y2*Y2 + B*Y2*Y2 + C*Y2 +D;
		  else
		    temp = 1.;

		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , temp );
		}
	      else
		{
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , 1.0 );
		}
	    }

	  // Do the right node.
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      // set y1 and y2
	      Y1 = ( dQmax[nodeR*NUM_VAR+j] - grid->nQ[nodeR*NUM_VAR+j] ) / (qright[j] - grid->nQ[nodeR*NUM_VAR+j]);
	      Y2 = ( dQmin[nodeR*NUM_VAR+j] - grid->nQ[nodeR*NUM_VAR+j] ) / (qright[j] - grid->nQ[nodeR*NUM_VAR+j]);

	      if ( ( qright[j] - grid->nQ[nodeR*NUM_VAR+j] ) > 0. )
		{
		  // min (1,y1)
		  if ( Y1 < Yt )
		    temp = A*Y1*Y1*Y1 + B*Y1*Y1 + C*Y1 +D;
		  else
		    temp = 1.;
		  
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , temp );
		}
	      else if ( ( qright[j] - grid->nQ[nodeR*NUM_VAR+j] ) < 0. )
		{
		  // min (1,y2)
		  if ( Y2 < Yt )
		    temp = A*Y2*Y2*Y2 + B*Y2*Y2 + C*Y2 +D;
		  else
		    temp = 1.;
		  
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , temp );
		}
	      else
		{
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , 1.0 );
		}
	    }
	  
	}
      else
	{
	  // Get the Q variables and if high order, then reconstruct the solution.
	  Reconstruct_UL( isubedge, nodeL, nodeR, qleft, qright, grid, p, grid->citer);

	  // Do the left node first.
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      // set y1 and y2.
	      Y1 = ( dQmax[nodeL*NUM_VAR+j] - grid->nQ[nodeL*NUM_VAR+j] ) / (qleft[j] - grid->nQ[nodeL*NUM_VAR+j]);
	      Y2 = ( dQmin[nodeL*NUM_VAR+j] - grid->nQ[nodeL*NUM_VAR+j] ) / (qleft[j] - grid->nQ[nodeL*NUM_VAR+j]);

	      if ( ( qleft[j] - grid->nQ[nodeL*NUM_VAR+j] ) > 0. )
		{
		  // min (1,y1)
		  if ( Y1 < Yt )
		    temp = A*Y1*Y1*Y1 + B*Y1*Y1 + C*Y1 +D;
		  else
		    temp = 1.;
		  
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , temp );
		}
	      else if ( ( qleft[j] - grid->nQ[nodeR*NUM_VAR+j] ) < 0. )
		{
		  // min (1,y2)
		  if ( Y2 < Yt )
		    temp = A*Y2*Y2*Y2 + B*Y2*Y2 + C*Y2 +D;
		  else
		    temp = 1.;

		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , temp );
		}
	      else
		{
		  phi[nodeL*8 + j*2 + 0] = MIN ( phi[nodeL*8 + j*2 + 0] , 1.0 );
		}
	    }

	  // Do the right node.
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      // set y1 and y2
	      Y1 = ( dQmax[nodeR*NUM_VAR+j] - grid->nQ[nodeR*NUM_VAR+j] ) / (qright[j] - grid->nQ[nodeR*NUM_VAR+j]);
	      Y2 = ( dQmin[nodeR*NUM_VAR+j] - grid->nQ[nodeR*NUM_VAR+j] ) / (qright[j] - grid->nQ[nodeR*NUM_VAR+j]);

	      if ( ( qright[j] - grid->nQ[nodeR*NUM_VAR+j] ) > 0. )
		{
		  // min (1,y1)
		  if ( Y1 < Yt )
		    temp = A*Y1*Y1*Y1 + B*Y1*Y1 + C*Y1 +D;
		  else
		    temp = 1.;
		  
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , temp );
		}
	      else if ( ( qright[j] - grid->nQ[nodeR*NUM_VAR+j] ) < 0. )
		{
		  // min (1,y2)
		  if ( Y2 < Yt )
		    temp = A*Y2*Y2*Y2 + B*Y2*Y2 + C*Y2 +D;
		  else
		    temp = 1.;
		  
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , temp );
		}
	      else
		{
		  phi[nodeR*8 + j*2 + 0] = MIN ( phi[nodeR*8 + j*2 + 0] , 1.0 );
		}
	    }

	}

      if ( gauss_flag )
	DO_GAUSS_QUAD = 1;

    }

  if ( LIMITER_NEIGHBOR_AVERAGE )  // As implemented here, I am doing this for the values of phi for all variables.
    {
      double *phi_avg = NULL;                      // For neighbor averaging.
      double epsPHI = 1.0E-2;
      double sum;
      int index1, index2, num_nn;
      int ind;
      
      phi_avg = (double*)malloc( (grid->nn+1)*4*sizeof(double) );
      if ( phi_avg == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'phi_avg'.\n"); exit(1); }
      
      for ( i=1; i <= grid->nn; i++ )
	{
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      // skip those close enough to 0 or to 1.
	      if ( phi[i*8+j*2+0] < (0.0+epsPHI) )
		{
		  phi_avg[i*4+j] = phi[i*8+j*2+0];
		  continue;
		}
	      
	      if ( phi[i*8+j*2+0] > (1.0-epsPHI) )
		{
		  phi_avg[i*4+j] = phi[i*8+j*2+0];
		  continue;
		}
	      
	      sum = 0.;
	      num_nn = 0;
	      
	      index1 = grid->nnsn[i];
	      index2 = grid->nnsn[i+1];
	      
	      for ( ind=index1; ind < index2; ind++ )
		{
		  nodeR = grid->nsn[ind];  // grab the neighbor.
		  
		  sum += phi[nodeR*8+j*2+0];
		  num_nn++;
		}

	      sum += phi[i*8+j*2+0];
	      num_nn++;
	      
	      sum = sum / ( (double)num_nn );
	      phi_avg[i*4+j] = sum;
	    }
	}

      for ( i=1; i <= grid->nn; i++ )
	{
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      phi[i*8+j*2+0] = phi_avg[i*4+j];
	    }
	}      

      freenull(phi_avg);
    }

  // Now for smooth regions of the flow, we can do something a little different.
  // Based on a certain criteria, we can use the switch to control the influence
  // of the higher derivatives' contribution.

  for ( i=1; i <= grid->nn; i++ )
    {
      bound = ( K*dx[i] ) * ( K*dx[i] ) * ( K*dx[i] );

      for ( j=0; j < NUM_VAR; j++ )
	{
	  delta_u = dQmax[i*NUM_VAR+j] - dQmin[i*NUM_VAR+j];
	  delta_u2 = delta_u * delta_u;

	  if ( delta_u2 < bound )
	    {
	      sigma = 1.;
	    }
	  else if ( delta_u2 > 2.*bound )
	    {
	      sigma = 0.;
	    }
	  else
	    {
	      y1 = ( delta_u2 - bound ) / bound;
	      sigma = 2.*y1*y1*y1 - 3.*y1*y1 + 1.;
	    }

	  // store sigma for the variable.
	  phi[i*8+j*2+1] = sigma;

	  // adjust the value of phi.
	  phi[i*8+j*2+0] = sigma + ( 1.0 - sigma ) * phi[i*8+j*2+0];
	  
	}
    }

  // Go ahead and put the minimum phi in the first spot.
  for ( i=1; i <= grid->nn; i++ )
    {
      y1 = phi[i*8 + 0*2 + 0];
      y2 = phi[i*8 + 0*2 + 1];

      for ( j=1; j < NUM_VAR; j++ )
	{
	  if ( phi[i*8 + j*2 + 0] < y1 )
	    {
	      y1 = phi[i*8 + j*2 + 0];
	      y2 = phi[i*8 + j*2 + 1];
	    }
	}

      // Now we have the minimal phi value in the first spot.
      phi[i*8 + 0*2 + 0] = y1;
      phi[i*8 + 0*2 + 1] = y2;
    }

  // Sanity check. All values of phi must be in [0,1]!

  for ( i=1; i <= grid->nn; i++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  if ( ( phi[i*8 + j*2 + 0] < 0.0 ) ||
	       ( phi[i*8 + j*2 + 0] > 1.0 )  )
	    {
	      printf("FATAL ERROR: Phi for node %d and variable %d does not have the proper value: < %.15E  >\n",i,j,phi[i*8+j*2+0]);
	      fflush(stdout);
	      exit(1);
	    }
	}
    }

  freenull(dQmin);
  freenull(dQmax);
  freenull(dx);

  return;
}
