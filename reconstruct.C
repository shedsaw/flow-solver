//=============================================================
// 
//  reconstruct.C
//  
//  Functions to reconstruct the variables.
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
#include "reconstruct.h"

//=============================================================
// 
//  Reconstruct ()
//
//  Applies the reconstruction operator to the left and right
//  state vectors.
//
//  int isedge;                          // Subedge index.
//  int nodeL;                           // The left node.
//  int nodeR;                           // The right node.
//  double QL[4];                        // The left state vector.
//  double QR[4];                        // The right state vector.
//  GRID *grid;                          // The grid.
//  PARAMS p;                            // The parameters.
//  int ts;                              // The current time step.
//
//=============================================================

void Reconstruct ( int isedge, int nodeL, int nodeR, double QL[4], double QR[4], GRID *grid, PARAMS p, int ts )
{
  int j;                                       // Loop counter.
  int rorder;                                  // Reconstruction order.
  int gelem;                                   // Global element index.
  double dx[2];                                // Edge direction.
  double xc[2];                                // Element centroid.
  double ec[2];                                // Edge centroid (center).
  double fc[2];                                // Face coordinates (where to reconstruct).
  double phiL,sigmaL;
  double phiR,sigmaR;

  int rorderL, rorderR;                        // Recostruction order to apply to the individual nodes along the edge.

  // Pointers.
  double *Q = grid->nQ;
  double *grad = grid->grad;
  double *hess = grid->hess;
  double *x = grid->x;

  // Extract the Q values.
  for ( j=0; j < NUM_VAR; j++ )
    {
      QL[j] = Q[nodeL*NUM_VAR+j];
      QR[j] = Q[nodeR*NUM_VAR+j];
    }

  // set phi and sigma.
  phiL = grid->phi[nodeL*8+0*2+0];
  phiR = grid->phi[nodeR*8+0*2+0];

  sigmaL = grid->phi[nodeL*8+0*2+1];
  sigmaR = grid->phi[nodeR*8+0*2+1];

  rorder = p.order;

  if ( p.order == 0 )
    {
      return;
    }
  if ( p.order == 1 )
    {
      if ( ts <= p.foits )
	{
	  rorder = 0;
	  rorderL = 0;
	  rorderR = 0;
	}
      else
	{
	  rorder = p.order;
	  rorderL = 1;
	  rorderR = 1;
	}
    }
  else if ( p.order == 2 ) // This implicitly means we are doing and edge based solve although I'm not checking for the flag here.
    {
      if ( ts <= p.foits )
	{
	  rorder = 0;
	  rorderL = 0;
	  rorderR = 0;
	}
      else
	{
	  rorder = p.order;
	  rorderL = 2;
	  rorderR = 2;
	}
    }
  else if ( p.order == 3 )
    {
      if ( ts <= p.foits )
	{
	  rorder = 0;
	  rorderL = 0;
	  rorderR = 0;
	}
      else if ( ts > p.foits && ts <= (p.foits+p.soits) )
	{
	  rorder = 1;
	  rorderL = 1;
	  rorderR = 1;
	}
      else
	{
	  rorder = 3;
	  
	  rorderL = grid->node_order[nodeL];
	  rorderR = grid->node_order[nodeR];
	}
    }
  else if ( p.order == 4 )
    {
      if ( ts <= p.foits )
	{
	  rorder = 0;
	  rorderL = 0;
	  rorderR = 0;
	}
      else if ( ts > p.foits && ts <= (p.foits+p.soits) )
	{
	  rorder = 1;
	  rorderL = 1;
	  rorderR = 1;
	}
      else if ( ts > (p.foits+p.soits) && ts <= (p.foits+p.soits+p.toits) )
	{
	  rorder = 3;
	  
	  if ( grid->node_order[nodeL] > 2 )
	    rorderL = 3;
	  else
	    rorderL = 1;

	  if ( grid->node_order[nodeR] > 2 )
	    rorderR = 3;
	  else
	    rorderR = 1;
	}
      else
	{
	  rorder = 4;
	  
	  rorderL = grid->node_order[nodeL];
	  rorderR = grid->node_order[nodeR];
	}
    }
  else if ( p.isrestart )
    {
      rorder = p.order;
      rorderL = grid->node_order[nodeL];
      rorderR = grid->node_order[nodeR];
    }

  // Now based on the reconstruction order we apply the correct terms.
  
  // left node.

  if ( rorderL == 0 )
    {
      ; // Nothing to do here.
    }
  else if ( rorderL == 1 )  // Linear reconstruction (2nd order)
    {
      // Compute the coordinates of the reconstruction point.
      gelem = grid->subedges[isedge*2+1];
      xc[0] = grid->el_cent[gelem*NDIM+0];
      xc[1] = grid->el_cent[gelem*NDIM+1];
      
      ec[0] = grid->xm_subedges[isedge*NDIM+0];
      ec[1] = grid->xm_subedges[isedge*NDIM+1];
      
      fc[0] = (xc[0] + ec[0]) * 0.5;
      fc[1] = (xc[1] + ec[1]) * 0.5;
      
      // Do the reconstruction...
      
      dx[0] = fc[0] - x[nodeL*NDIM+0];
      dx[1] = fc[1] - x[nodeL*NDIM+1];
      
      for ( j=0; j < NUM_VAR; j++ )
	{
	  //phiL = grid->phi[nodeL*8+j*2+0];
	  QL[j] += ( phiL*grad[nodeL*NDIM*NUM_VAR+j*NDIM+0]*dx[0] +
		     phiL*grad[nodeL*NDIM*NUM_VAR+j*NDIM+1]*dx[1]);
	}
    }
  else if ( rorderL == 2 )  // Quadratic reconstruction from the linear data. Pseudo 3rd order.
    {
      fc[0] = grid->xm_subedges[isedge*NDIM+0];
      fc[1] = grid->xm_subedges[isedge*NDIM+1];

      // Do the reconstruction...

      dx[0] = x[nodeR*NDIM+0] - x[nodeL*NDIM+0];
      dx[1] = x[nodeR*NDIM+1] - x[nodeL*NDIM+1];

      for ( j=0; j < NUM_VAR; j++ )
	{
	  //phiL = grid->phi[nodeL*8+j*2+0];
	  QL[j] += ( 0.25 * ( Q[nodeR*NUM_VAR+j] - Q[nodeL*NUM_VAR+j] ) +
		     0.25 * (phiL*grad[nodeL*NDIM*NUM_VAR+j*NDIM+0]*dx[0] +
			     phiL*grad[nodeL*NDIM*NUM_VAR+j*NDIM+1]*dx[1]) );
	}
    }
  else if ( rorderL == 3 )
    {
      // Compute the coordinates of the reconstruction point.
      gelem = grid->subedges[isedge*2+1];
      xc[0] = grid->el_cent[gelem*NDIM+0];
      xc[1] = grid->el_cent[gelem*NDIM+1];
      
      ec[0] = grid->xm_subedges[isedge*NDIM+0];
      ec[1] = grid->xm_subedges[isedge*NDIM+1];

      fc[0] = (xc[0] + ec[0]) * 0.5;
      fc[1] = (xc[1] + ec[1]) * 0.5;

      // Do the reconstruction...

      dx[0] = fc[0] - x[nodeL*NDIM+0];
      dx[1] = fc[1] - x[nodeL*NDIM+1];

      if ( PV_EXT )
	{
	  // Extract the Q values.
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      QL[j] = grid->point_value[nodeL*NUM_VAR+j];
	    }
	}
      
      for ( j=0; j < NUM_VAR; j++ )
	{
	  //phiL = grid->phi[nodeL*8+j*2+0];
	  //sigmaL = grid->phi[nodeL*8+j*2+1];
	  QL[j] += ( phiL*hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+0]*dx[0] +
		     phiL*hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+1]*dx[1] +
		     sigmaL*hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+2]*dx[0]*dx[0]*0.5 +
		     sigmaL*hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+3]*dx[1]*dx[1]*0.5 +
		     sigmaL*hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+4]*dx[0]*dx[1] );     
	}
    }
  else if ( rorderL == 4 )
    {
      // Compute the coordinates of the reconstruction point.
      gelem = grid->subedges[isedge*2+1];
      xc[0] = grid->el_cent[gelem*NDIM+0];
      xc[1] = grid->el_cent[gelem*NDIM+1];
      
      ec[0] = grid->xm_subedges[isedge*NDIM+0];
      ec[1] = grid->xm_subedges[isedge*NDIM+1];

      fc[0] = (xc[0] + ec[0]) * 0.5;
      fc[1] = (xc[1] + ec[1]) * 0.5;

      // Do the reconstruction...

      dx[0] = fc[0] - x[nodeL*NDIM+0];
      dx[1] = fc[1] - x[nodeL*NDIM+1];

      if ( PV_EXT )
	{
	  // Extract the Q values.
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      QL[j] = grid->point_value[nodeL*NUM_VAR+j];
	    }
	}

      for ( j=0; j < NUM_VAR; j++ )
	{
	  //phiL = grid->phi[nodeL*8+j*2+0];
	  //sigmaL = grid->phi[nodeL*8+j*2+1];
	  QL[j] += ( phiL*hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+0]*dx[0] +
		     phiL*hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+1]*dx[1] +
		     sigmaL*hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+2]*dx[0]*dx[0]*0.5 +
		     sigmaL*hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+3]*dx[1]*dx[1]*0.5 +
		     sigmaL*hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+4]*dx[0]*dx[1] +
		     sigmaL*hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+5]*dx[0]*dx[0]*dx[0]*(1./6.) +
		     sigmaL*hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+6]*dx[1]*dx[1]*dx[1]*(1./6.) +
		     sigmaL*hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+7]*dx[0]*dx[0]*dx[1]*0.5 +
		     sigmaL*hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+8]*dx[0]*dx[1]*dx[1]*0.5 );
	}
    }
  else
    {
      printf("FATAL ERROR: In Reconstruct(): rorderL is undefined <%d>. subedge %d , nodeL = %d, nodeR = %d, time step %d.\n",rorderL,
	     isedge, nodeL, nodeR, ts);
      fflush(stdout);
      exit(0);
    }

  // right node.

  if ( rorderR == 0 )
    {
      ; // Nothing to do here.
    }
  else if ( rorderR == 1 )  // Linear reconstruction (2nd order)
    {
      // Compute the coordinates of the reconstruction point.
      gelem = grid->subedges[isedge*2+1];
      xc[0] = grid->el_cent[gelem*NDIM+0];
      xc[1] = grid->el_cent[gelem*NDIM+1];
      
      ec[0] = grid->xm_subedges[isedge*NDIM+0];
      ec[1] = grid->xm_subedges[isedge*NDIM+1];
      
      fc[0] = (xc[0] + ec[0]) * 0.5;
      fc[1] = (xc[1] + ec[1]) * 0.5;
      
      // Do the reconstruction...
      dx[0] = fc[0] - x[nodeR*NDIM+0];
      dx[1] = fc[1] - x[nodeR*NDIM+1];

      for ( j=0; j < NUM_VAR; j++ )
	{
	  //phiR = grid->phi[nodeR*8+j*2+0];
	  QR[j] += ( phiR*grad[nodeR*NDIM*NUM_VAR+j*NDIM+0]*dx[0] +
		     phiR*grad[nodeR*NDIM*NUM_VAR+j*NDIM+1]*dx[1]);
	}
    }
  else if ( rorderR == 2 )
    {
      fc[0] = grid->xm_subedges[isedge*NDIM+0];
      fc[1] = grid->xm_subedges[isedge*NDIM+1];

      dx[0] = x[nodeL*NDIM+0] - x[nodeR*NDIM+0];
      dx[1] = x[nodeL*NDIM+1] - x[nodeR*NDIM+1];

      for ( j=0; j < NUM_VAR; j++ )
	{
	  //phiR = grid->phi[nodeR*8+0*2+0];
	  QR[j] += ( 0.25 * ( Q[nodeL*NUM_VAR+j] - Q[nodeR*NUM_VAR+j] ) + 
		     0.25 * ( phiR*grad[nodeR*NDIM*NUM_VAR+j*NDIM+0]*dx[0] +
			      phiR*grad[nodeR*NDIM*NUM_VAR+j*NDIM+1]*dx[1]) );
	}
    }
  else if ( rorderR == 3 )
    {
      // Compute the coordinates of the reconstruction point.
      gelem = grid->subedges[isedge*2+1];
      xc[0] = grid->el_cent[gelem*NDIM+0];
      xc[1] = grid->el_cent[gelem*NDIM+1];
      
      ec[0] = grid->xm_subedges[isedge*NDIM+0];
      ec[1] = grid->xm_subedges[isedge*NDIM+1];

      fc[0] = (xc[0] + ec[0]) * 0.5;
      fc[1] = (xc[1] + ec[1]) * 0.5;

      // Do the reconstruction...
      if ( PV_EXT )
	{
	  // Extract the Q values.
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      QR[j] = grid->point_value[nodeR*NUM_VAR+j];
	    }
	}

      dx[0] = fc[0] - x[nodeR*NDIM+0];
      dx[1] = fc[1] - x[nodeR*NDIM+1];
      
      for ( j=0; j < NUM_VAR; j++ )
	{
	  //phiR = grid->phi[nodeR*8+j*2+0];
	  //sigmaR = grid->phi[nodeR*8+j*2+1];
	  QR[j] += ( phiR*hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+0]*dx[0] +
		     phiR*hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+1]*dx[1] +
		     sigmaR*hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+2]*dx[0]*dx[0]*0.5 +
		     sigmaR*hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+3]*dx[1]*dx[1]*0.5 +
		     sigmaR*hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+4]*dx[0]*dx[1] );     
	}
    }
  else if ( rorderR == 4 )
    {
      // Compute the coordinates of the reconstruction point.
      gelem = grid->subedges[isedge*2+1];
      xc[0] = grid->el_cent[gelem*NDIM+0];
      xc[1] = grid->el_cent[gelem*NDIM+1];
      
      ec[0] = grid->xm_subedges[isedge*NDIM+0];
      ec[1] = grid->xm_subedges[isedge*NDIM+1];

      fc[0] = (xc[0] + ec[0]) * 0.5;
      fc[1] = (xc[1] + ec[1]) * 0.5;

      // Do the reconstruction...

      if ( PV_EXT )
	{
	  // Extract the Q values.
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      QR[j] = grid->point_value[nodeR*NUM_VAR+j];
	    }
	}

      dx[0] = fc[0] - x[nodeR*NDIM+0];
      dx[1] = fc[1] - x[nodeR*NDIM+1];
      
      for ( j=0; j < NUM_VAR; j++ )
	{
	  //phiR = grid->phi[nodeR*8+j*2+0];
	  //sigmaR = grid->phi[nodeR*8+j*2+1];
	  QR[j] += ( phiR*hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+0]*dx[0] +
		     phiR*hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+1]*dx[1] +
		     sigmaR*hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+2]*dx[0]*dx[0]*0.5 +
		     sigmaR*hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+3]*dx[1]*dx[1]*0.5 +
		     sigmaR*hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+4]*dx[0]*dx[1] +
                     sigmaL*hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+5]*dx[0]*dx[0]*dx[0]*(1./6.) +
		     sigmaL*hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+6]*dx[1]*dx[1]*dx[1]*(1./6.) +
		     sigmaL*hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+7]*dx[0]*dx[0]*dx[1]*0.5 +
		     sigmaL*hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+8]*dx[0]*dx[1]*dx[1]*0.5 );     
	}
    }
  else
    {
      printf("FATAL ERROR: In Reconstruct(): rorderR is undefined <%d>. subedge %d , nodeL = %d, nodeR = %d, time step %d.\n",rorderR,
	     isedge, nodeL, nodeR, ts);
      fflush(stdout);
      exit(0);
    }
      
  return;
}


//=============================================================
// 
//  Reconstruct_UL ()
//
//  Applies the reconstruction operator to the left and right
//  state vectors.
//
//  int isedge;                          // Subedge index.
//  int nodeL;                           // The left node.
//  int nodeR;                           // The right node.
//  double QL[4];                        // The left state vector.
//  double QR[4];                        // The right state vector.
//  GRID *grid;                          // The grid.
//  PARAMS p;                            // The parameters.
//  int ts;                              // The current time step.
//
//=============================================================

void Reconstruct_UL ( int isedge, int nodeL, int nodeR, double QL[4], double QR[4], GRID *grid, PARAMS p, int ts )
{
  int j;                                       // Loop counter.
  int rorder;                                  // Reconstruction order.
  int rorderL, rorderR;                        // Recostruction order to apply to the individual nodes along the edge.
  int gelem;                                   // Global element index.
  double dx[2];                                // Edge direction.
  double xc[2];                                // Element centroid.
  double ec[2];                                // Edge centroid (center).
  double fc[2];                                // Face coordinates (where to reconstruct).

  // Pointers.
  double *Q = grid->nQ;
  double *grad = grid->grad;
  double *hess = grid->hess;
  double *x = grid->x;

  // Extract the Q values.
  for ( j=0; j < NUM_VAR; j++ )
    {
      QL[j] = Q[nodeL*NUM_VAR+j];
      QR[j] = Q[nodeR*NUM_VAR+j];
    }

  rorder = p.order;

  if ( p.order == 0 )
    {
      return;
    }
  if ( p.order == 1 )
    {
      if ( ts <= p.foits )
	{
	  rorder = 0;
	  rorderL = 0;
	  rorderR = 0;
	}
      else
	{
	  rorder = p.order;
	  rorderL = 1;
	  rorderR = 1;
	}
    }
  else if ( p.order == 2 ) // This implicitly means we are doing and edge based solve although I'm not checking for the flag here.
    {
      if ( ts <= p.foits )
	{
	  rorder = 0;
	  rorderL = 0;
	  rorderR = 0;
	}
      else
	{
	  rorder = p.order;
	  rorderL = 2;
	  rorderR = 2;
	}
    }
  else if ( p.order == 3 )
    {
      if ( ts <= p.foits )
	{
	  rorder = 0;
	  rorderL = 0;
	  rorderR = 0;
	}
      else if ( ts > p.foits && ts <= (p.foits+p.soits) )
	{
	  rorder = 1;
	  rorderL = 1;
	  rorderR = 1;
	}
      else
	{
	  rorder = 3;
	  
	  rorderL = grid->node_order[nodeL];
	  rorderR = grid->node_order[nodeR];
	}
    }
  else if ( p.order == 4 )
    {
      if ( ts <= p.foits )
	{
	  rorder = 0;
	  rorderL = 0;
	  rorderR = 0;
	}
      else if ( ts > p.foits && ts <= (p.foits+p.soits) )
	{
	  rorder = 1;
	  rorderL = 1;
	  rorderR = 1;
	}
      else if ( ts > (p.foits+p.soits) && ts <= (p.foits+p.soits+p.toits) )
	{
	  rorder = 3;
	  
	  if ( grid->node_order[nodeL] > 2 )
	    rorderL = 3;
	  else
	    rorderL = 1;

	  if ( grid->node_order[nodeR] > 2 )
	    rorderR = 3;
	  else
	    rorderR = 1;
	}
      else
	{
	  rorder = 4;
	  
	  rorderL = grid->node_order[nodeL];
	  rorderR = grid->node_order[nodeR];
	}
    }
  else if ( p.isrestart )
    {
      rorder = p.order;
      rorderL = grid->node_order[nodeL];
      rorderR = grid->node_order[nodeR];
    }

  // Now based on the reconstruction order we apply the correct terms.
  
  // left node.

  if ( rorderL == 0 )
    {
      ; // Nothing to do here.
    }
  else if ( rorderL == 1 )  // Linear reconstruction (2nd order)
    {
      // Compute the coordinates of the reconstruction point.
      gelem = grid->subedges[isedge*2+1];
      xc[0] = grid->el_cent[gelem*NDIM+0];
      xc[1] = grid->el_cent[gelem*NDIM+1];
      
      ec[0] = grid->xm_subedges[isedge*NDIM+0];
      ec[1] = grid->xm_subedges[isedge*NDIM+1];
      
      fc[0] = (xc[0] + ec[0]) * 0.5;
      fc[1] = (xc[1] + ec[1]) * 0.5;
      
      // Do the reconstruction...
      
      dx[0] = fc[0] - x[nodeL*NDIM+0];
      dx[1] = fc[1] - x[nodeL*NDIM+1];
      
      for ( j=0; j < NUM_VAR; j++ )
	{
	  QL[j] += ( grad[nodeL*NDIM*NUM_VAR+j*NDIM+0]*dx[0] +
		     grad[nodeL*NDIM*NUM_VAR+j*NDIM+1]*dx[1]);
	}
    }
  else if ( rorderL == 2 )  // Quadratic reconstruction from the linear data. Pseudo 3rd order.
    {
      fc[0] = grid->xm_subedges[isedge*NDIM+0];
      fc[1] = grid->xm_subedges[isedge*NDIM+1];

      // Do the reconstruction...

      dx[0] = x[nodeR*NDIM+0] - x[nodeL*NDIM+0];
      dx[1] = x[nodeR*NDIM+1] - x[nodeL*NDIM+1];

      for ( j=0; j < NUM_VAR; j++ )
	{
	  //phiL = grid->phi[nodeL*8+j*2+0];
	  QL[j] += ( 0.25 * ( Q[nodeR*NUM_VAR+j] - Q[nodeL*NUM_VAR+j] ) +
		     0.25 * (grad[nodeL*NDIM*NUM_VAR+j*NDIM+0]*dx[0] +
			     grad[nodeL*NDIM*NUM_VAR+j*NDIM+1]*dx[1]) );
	}
    }
  else if ( rorderL == 3 )
    {
      // Compute the coordinates of the reconstruction point.
      gelem = grid->subedges[isedge*2+1];
      xc[0] = grid->el_cent[gelem*NDIM+0];
      xc[1] = grid->el_cent[gelem*NDIM+1];
      
      ec[0] = grid->xm_subedges[isedge*NDIM+0];
      ec[1] = grid->xm_subedges[isedge*NDIM+1];

      fc[0] = (xc[0] + ec[0]) * 0.5;
      fc[1] = (xc[1] + ec[1]) * 0.5;

      // Do the reconstruction...

      dx[0] = fc[0] - x[nodeL*NDIM+0];
      dx[1] = fc[1] - x[nodeL*NDIM+1];

      if ( PV_EXT )
	{
	  // Extract the Q values.
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      QL[j] = grid->point_value[nodeL*NUM_VAR+j];
	    }
	}
      
      for ( j=0; j < NUM_VAR; j++ )
	{
	  QL[j] += ( hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+0]*dx[0] +
		     hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+1]*dx[1] +
		     hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+2]*dx[0]*dx[0]*0.5 +
		     hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+3]*dx[1]*dx[1]*0.5 +
		     hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+4]*dx[0]*dx[1] );     
	}
    }
  else if ( rorderL == 4 )
    {
      // Compute the coordinates of the reconstruction point.
      gelem = grid->subedges[isedge*2+1];
      xc[0] = grid->el_cent[gelem*NDIM+0];
      xc[1] = grid->el_cent[gelem*NDIM+1];
      
      ec[0] = grid->xm_subedges[isedge*NDIM+0];
      ec[1] = grid->xm_subedges[isedge*NDIM+1];

      fc[0] = (xc[0] + ec[0]) * 0.5;
      fc[1] = (xc[1] + ec[1]) * 0.5;

      // Do the reconstruction...

      dx[0] = fc[0] - x[nodeL*NDIM+0];
      dx[1] = fc[1] - x[nodeL*NDIM+1];

      if ( PV_EXT )
	{
	  // Extract the Q values.
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      QL[j] = grid->point_value[nodeL*NUM_VAR+j];
	    }
	}

      for ( j=0; j < NUM_VAR; j++ )
	{
	  QL[j] += ( hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+0]*dx[0] +
		     hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+1]*dx[1] +
		     hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+2]*dx[0]*dx[0]*0.5 +
		     hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+3]*dx[1]*dx[1]*0.5 +
		     hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+4]*dx[0]*dx[1] +
		     hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+5]*dx[0]*dx[0]*dx[0]*(1./6.) +
		     hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+6]*dx[1]*dx[1]*dx[1]*(1./6.) +
		     hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+7]*dx[0]*dx[0]*dx[1]*0.5 +
		     hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+8]*dx[0]*dx[1]*dx[1]*0.5 );
	}
    }
  else
    {
      printf("FATAL ERROR: In Reconstruct_UL(): rorderL is undefined <%d>. subedge %d , nodeL = %d, nodeR = %d, time step %d.\n",rorderL,
	     isedge, nodeL, nodeR, ts);
      fflush(stdout);
      exit(0);
    }

  // right node.

  if ( rorderR == 0 )
    {
      ; // Nothing to do here.
    }
  else if ( rorderR == 1 )  // Linear reconstruction (2nd order)
    {
      // Compute the coordinates of the reconstruction point.
      gelem = grid->subedges[isedge*2+1];
      xc[0] = grid->el_cent[gelem*NDIM+0];
      xc[1] = grid->el_cent[gelem*NDIM+1];
      
      ec[0] = grid->xm_subedges[isedge*NDIM+0];
      ec[1] = grid->xm_subedges[isedge*NDIM+1];
      
      fc[0] = (xc[0] + ec[0]) * 0.5;
      fc[1] = (xc[1] + ec[1]) * 0.5;
      
      // Do the reconstruction...
      dx[0] = fc[0] - x[nodeR*NDIM+0];
      dx[1] = fc[1] - x[nodeR*NDIM+1];

      for ( j=0; j < NUM_VAR; j++ )
	{
	  QR[j] += ( grad[nodeR*NDIM*NUM_VAR+j*NDIM+0]*dx[0] +
		     grad[nodeR*NDIM*NUM_VAR+j*NDIM+1]*dx[1]);
	}
    }
  else if ( rorderR == 2 )
    {
      fc[0] = grid->xm_subedges[isedge*NDIM+0];
      fc[1] = grid->xm_subedges[isedge*NDIM+1];

      dx[0] = x[nodeL*NDIM+0] - x[nodeR*NDIM+0];
      dx[1] = x[nodeL*NDIM+1] - x[nodeR*NDIM+1];

      for ( j=0; j < NUM_VAR; j++ )
	{
	  QR[j] += ( 0.25 * ( Q[nodeL*NUM_VAR+j] - Q[nodeR*NUM_VAR+j] ) + 
		     0.25 * ( grad[nodeR*NDIM*NUM_VAR+j*NDIM+0]*dx[0] +
			      grad[nodeR*NDIM*NUM_VAR+j*NDIM+1]*dx[1]) );
	}
    }
  else if ( rorderR == 3 )
    {
      // Compute the coordinates of the reconstruction point.
      gelem = grid->subedges[isedge*2+1];
      xc[0] = grid->el_cent[gelem*NDIM+0];
      xc[1] = grid->el_cent[gelem*NDIM+1];
      
      ec[0] = grid->xm_subedges[isedge*NDIM+0];
      ec[1] = grid->xm_subedges[isedge*NDIM+1];

      fc[0] = (xc[0] + ec[0]) * 0.5;
      fc[1] = (xc[1] + ec[1]) * 0.5;

      // Do the reconstruction...
      if ( PV_EXT )
	{
	  // Extract the Q values.
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      QR[j] = grid->point_value[nodeR*NUM_VAR+j];
	    }
	}

      dx[0] = fc[0] - x[nodeR*NDIM+0];
      dx[1] = fc[1] - x[nodeR*NDIM+1];
      
      for ( j=0; j < NUM_VAR; j++ )
	{
	  QR[j] += ( hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+0]*dx[0] +
		     hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+1]*dx[1] +
		     hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+2]*dx[0]*dx[0]*0.5 +
		     hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+3]*dx[1]*dx[1]*0.5 +
		     hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+4]*dx[0]*dx[1] );     
	}
    }
  else if ( rorderR == 4 )
    {
      // Compute the coordinates of the reconstruction point.
      gelem = grid->subedges[isedge*2+1];
      xc[0] = grid->el_cent[gelem*NDIM+0];
      xc[1] = grid->el_cent[gelem*NDIM+1];
      
      ec[0] = grid->xm_subedges[isedge*NDIM+0];
      ec[1] = grid->xm_subedges[isedge*NDIM+1];

      fc[0] = (xc[0] + ec[0]) * 0.5;
      fc[1] = (xc[1] + ec[1]) * 0.5;

      // Do the reconstruction...

      if ( PV_EXT )
	{
	  // Extract the Q values.
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      QR[j] = grid->point_value[nodeR*NUM_VAR+j];
	    }
	}

      dx[0] = fc[0] - x[nodeR*NDIM+0];
      dx[1] = fc[1] - x[nodeR*NDIM+1];
      
      for ( j=0; j < NUM_VAR; j++ )
	{
	  QR[j] += ( hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+0]*dx[0] +
		     hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+1]*dx[1] +
		     hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+2]*dx[0]*dx[0]*0.5 +
		     hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+3]*dx[1]*dx[1]*0.5 +
		     hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+4]*dx[0]*dx[1] +
                     hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+5]*dx[0]*dx[0]*dx[0]*(1./6.) +
		     hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+6]*dx[1]*dx[1]*dx[1]*(1./6.) +
		     hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+7]*dx[0]*dx[0]*dx[1]*0.5 +
		     hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+8]*dx[0]*dx[1]*dx[1]*0.5 );     
	}
    }
  else
    {
      printf("FATAL ERROR: In Reconstruct_UL(): rorderR is undefined <%d>. subedge %d , nodeL = %d, nodeR = %d, time step %d.\n",rorderR,
	     isedge, nodeL, nodeR, ts);
      fflush(stdout);
      exit(0);
    }
  

  return;
}


// WEAKLY SUPPORTS THE PSEUDO 3rd ORDER EXTRAPOLATION
//=============================================================
// 
//  generic_reconstruct ()
//
//  Applies the reconstruction operator to a node given the vector.
//
//  I AM NOT DOING ANY CHECKING SO USE THIS METHOD WITH CAUTION.
//  IT IS INTENDED FOR USE ONLY WITH THE VORTEX CONVECTION PROBLEM
//  IN IC.C .
//
//  int node;                            // The node to reconstruct.
//  double Q[4];                         // The state vector.
//  double dx[2];                        // The distance vector.
//  GRID *grid;                          // The grid.
//  PARAMS p;                            // The parameters.
//
//=============================================================

void generic_reconstruct ( int node, double Q[4], double dx[2], GRID *grid, PARAMS p )
{
  int j;                                       // Loop counter.
  int rorder;                                  // Reconstruction order.
  double phi, sigma;

  // Pointers.
  double *grad = grid->grad;
  double *hess = grid->hess;

  phi = grid->phi[node*8+0*2+0];
  sigma = grid->phi[node*8+0*2+1];

  // Extract the Q values.
  for ( j=0; j < NUM_VAR; j++ )
    {
      Q[j] = grid->nQ[node*NUM_VAR+j];
    }

  //rorder = p.order;
  rorder = grid->node_order[node];

  // Now based on the reconstruction order we apply the correct terms.
  if ( rorder == 0 )
    {
      return;
    }

  else if ( rorder == 1 || rorder == 2 )  // Linear reconstruction (2nd order)
    {
      // Do the reconstruction...

      for ( j=0; j < NUM_VAR; j++ )
	{
	  //phi = grid->phi[node*8+j*2+0];
	  Q[j] += ( phi*grad[node*NDIM*NUM_VAR+j*NDIM+0]*dx[0] +
		    phi*grad[node*NDIM*NUM_VAR+j*NDIM+1]*dx[1]);
	}

      return;
    }
  else if ( rorder == 3 )  // Quadratic reconstruction (3rd order).
    {
      // Do the reconstruction...
      
      // Get point values.
      if ( PV_EXT )
	{
	  for ( j=0; j < NUM_VAR; j++ )
	    Q[j] = grid->point_value[node*NUM_VAR+j];
	}

      for ( j=0; j < NUM_VAR; j++ )
	{
	  //phi = grid->phi[node*8+j*2+0];
	  //sigma = grid->phi[node*8+j*2+1];
	  Q[j] += ( phi*hess[node*NUM_MOM*NUM_VAR+j*NUM_MOM+0]*dx[0] +
		    phi*hess[node*NUM_MOM*NUM_VAR+j*NUM_MOM+1]*dx[1] +
		    sigma*hess[node*NUM_MOM*NUM_VAR+j*NUM_MOM+2]*dx[0]*dx[0]*0.5 +
		    sigma*hess[node*NUM_MOM*NUM_VAR+j*NUM_MOM+3]*dx[1]*dx[1]*0.5 +
		    sigma*hess[node*NUM_MOM*NUM_VAR+j*NUM_MOM+4]*dx[0]*dx[1] );     
	}

      return;
    }
  else if ( rorder == 4 )
    {
      // Get point values.
      if ( PV_EXT )
	{
	  for ( j=0; j < NUM_VAR; j++ )
	    Q[j] = grid->point_value[node*NUM_VAR+j];
	}
      
      for ( j=0; j < NUM_VAR; j++ )
	{
	  //phi = grid->phi[node*8+j*2+0];
	  //sigma = grid->phi[node*8+j*2+1];
	  Q[j] += ( phi*hess[node*NUM_MOM*NUM_VAR+j*NUM_MOM+0]*dx[0] +
		    phi*hess[node*NUM_MOM*NUM_VAR+j*NUM_MOM+1]*dx[1] +
		    sigma*hess[node*NUM_MOM*NUM_VAR+j*NUM_MOM+2]*dx[0]*dx[0]*0.5 +
		    sigma*hess[node*NUM_MOM*NUM_VAR+j*NUM_MOM+3]*dx[1]*dx[1]*0.5 +
		    sigma*hess[node*NUM_MOM*NUM_VAR+j*NUM_MOM+4]*dx[0]*dx[1] +
		    sigma*hess[node*NUM_MOM*NUM_VAR+j*NUM_MOM+5]*dx[0]*dx[0]*dx[0]*(1./6.) +
		    sigma*hess[node*NUM_MOM*NUM_VAR+j*NUM_MOM+6]*dx[1]*dx[1]*dx[1]*(1./6.) +
		    sigma*hess[node*NUM_MOM*NUM_VAR+j*NUM_MOM+7]*dx[0]*dx[0]*dx[1]*0.5 +
		    sigma*hess[node*NUM_MOM*NUM_VAR+j*NUM_MOM+8]*dx[0]*dx[1]*dx[1]*0.5 );
	}
      return;
    }
  else
    {
      printf("FATAL ERROR: generic_reconstruct() has an invalid order <%d>.\n",rorder);
      fflush(stdout);
      exit(0);
    }

  return;
}


// DOES NOT SUPPORT THE PSEUDO 3rd ORDER EXTRAPOLATION
//=============================================================
// 
//  generic_reconstruct_UL ()
//
//  Applies the reconstruction operator to a node given the vector.
//
//  I AM NOT DOING ANY CHECKING SO USE THIS METHOD WITH CAUTION.
//  IT IS INTENDED FOR USE ONLY WITH THE VORTEX CONVECTION PROBLEM
//  IN IC.C .
//
//  int node;                            // The node to reconstruct.
//  double Q[4];                         // The state vector.
//  double dx[2];                        // The distance vector.
//  GRID *grid;                          // The grid.
//  PARAMS p;                            // The parameters.
//
//=============================================================

void generic_reconstruct_UL ( int node, double Q[4], double dx[2], GRID *grid, PARAMS p )
{
  int j;                                       // Loop counter.
  int rorder;                                  // Reconstruction order.

  // Pointers.
  double *grad = grid->grad;
  double *hess = grid->hess;

  // Extract the Q values.
  for ( j=0; j < NUM_VAR; j++ )
    {
      Q[j] = grid->nQ[node*NUM_VAR+j];
    }

  //rorder = p.order;
  rorder = grid->node_order[node];

  // Now based on the reconstruction order we apply the correct terms.
  if ( rorder == 0 )
    {
      return;
    }

  else if ( rorder == 1 )  // Linear reconstruction (2nd order)
    {
      // Do the reconstruction...

      for ( j=0; j < NUM_VAR; j++ )
	{
	  Q[j] += ( grad[node*NDIM*NUM_VAR+j*NDIM+0]*dx[0] +
		    grad[node*NDIM*NUM_VAR+j*NDIM+1]*dx[1]);
	}

      return;
    }

  else  if ( rorder == 3 ) // Quadratic reconstruction (3rd order).
    {
      // Do the reconstruction...
      
      // Get point values.
      if ( PV_EXT )
	{
	  for ( j=0; j < NUM_VAR; j++ )
	    Q[j] = grid->point_value[node*NUM_VAR+j];
	}

      for ( j=0; j < NUM_VAR; j++ )
	{
	  Q[j] += ( hess[node*NUM_MOM*NUM_VAR+j*NUM_MOM+0]*dx[0] +
		    hess[node*NUM_MOM*NUM_VAR+j*NUM_MOM+1]*dx[1] +
		    hess[node*NUM_MOM*NUM_VAR+j*NUM_MOM+2]*dx[0]*dx[0]*0.5 +
		    hess[node*NUM_MOM*NUM_VAR+j*NUM_MOM+3]*dx[1]*dx[1]*0.5 +
		    hess[node*NUM_MOM*NUM_VAR+j*NUM_MOM+4]*dx[0]*dx[1] );     
	}

      return;
    }
  else if ( rorder == 4 )
    {
      // Get point values.
      if ( PV_EXT )
	{
	  for ( j=0; j < NUM_VAR; j++ )
	    Q[j] = grid->point_value[node*NUM_VAR+j];
	}
      
      for ( j=0; j < NUM_VAR; j++ )
	{
	  Q[j] += ( hess[node*NUM_MOM*NUM_VAR+j*NUM_MOM+0]*dx[0] +
		    hess[node*NUM_MOM*NUM_VAR+j*NUM_MOM+1]*dx[1] +
		    hess[node*NUM_MOM*NUM_VAR+j*NUM_MOM+2]*dx[0]*dx[0]*0.5 +
		    hess[node*NUM_MOM*NUM_VAR+j*NUM_MOM+3]*dx[1]*dx[1]*0.5 +
		    hess[node*NUM_MOM*NUM_VAR+j*NUM_MOM+4]*dx[0]*dx[1] +
		    hess[node*NUM_MOM*NUM_VAR+j*NUM_MOM+5]*dx[0]*dx[0]*dx[0]*(1./6.) +
		    hess[node*NUM_MOM*NUM_VAR+j*NUM_MOM+6]*dx[1]*dx[1]*dx[1]*(1./6.) +
		    hess[node*NUM_MOM*NUM_VAR+j*NUM_MOM+7]*dx[0]*dx[0]*dx[1]*0.5 +
		    hess[node*NUM_MOM*NUM_VAR+j*NUM_MOM+8]*dx[0]*dx[1]*dx[1]*0.5 );
	}
      return;
    }
  else
    {
      printf("FATAL ERROR: generic_reconstruct_UL() has an invalid order <%d>.\n",rorder);
      fflush(stdout);
      exit(0);
    }

  return;
}


//=============================================================
// 
//  Reconstruct_Gauss ()
//
//  Applies the reconstruction operator to the left and right
//  state vectors.
//
//  int isedge;                          // Subedge index.
//  int nodeL;                           // The left node.
//  int nodeR;                           // The right node.
//  double GP[2];                        // Coordinates of the Gauss point.
//  double QL[4];                        // The left state vector.
//  double QR[4];                        // The right state vector.
//  GRID *grid;                          // The grid.
//  PARAMS p;                            // The parameters.
//  int ts;                              // The current time step.
//
//=============================================================

void Reconstruct_Gauss ( int isedge, int nodeL, int nodeR, double GP[2], double QL[4], double QR[4], GRID *grid, PARAMS p, int ts )
{
  int j;                                       // Loop counter.
  int rorder;                                  // Reconstruction order.
  int rorderL, rorderR;                        // Recostruction order to apply to the individual nodes along the edge.
  double dx[2];                                // Edge direction.
  double phiL,sigmaL;
  double phiR,sigmaR;

  // Pointers.
  double *Q = grid->nQ;
  double *grad = grid->grad;
  double *hess = grid->hess;
  double *x = grid->x;

  // Extract the Q values.
  for ( j=0; j < NUM_VAR; j++ )
    {
      QL[j] = Q[nodeL*NUM_VAR+j];
      QR[j] = Q[nodeR*NUM_VAR+j];
    }

  // set phi and sigma.
  phiL = grid->phi[nodeL*8+0*2+0];
  phiR = grid->phi[nodeR*8+0*2+0];

  sigmaL = grid->phi[nodeL*8+0*2+1];
  sigmaR = grid->phi[nodeR*8+0*2+1];

  if ( p.order == 2 || EDGE_BASED )
    {
      printf("CRITICAL ERROR: The solver is trying to reconstruct to Gaussian quadrature nodes but the solve is specified as edge based!\n");
      printf("In Reconstruct_Gauss().\n");
      fflush(stdout);
      exit(1);
    }

  if ( p.order == 0 )
    {
      rorder = 0;
      return;
    }

  if ( p.order == 1 )
    {
      if ( grid->citer <= p.foits )
	{
	  rorder = 0;
	  rorderL = 0;
	  rorderR = 0;
	}
      else
	{
	  rorder = 1;
	  rorderL = 1;
	  rorderR = 1;
	}
    }
  else if ( p.order == 3 )
    {
      if ( grid->citer <= p.foits )
	{
	  rorder = 0;
	  rorderL = 0;
	  rorderR = 0;
	}
      else if (  ( grid->citer  >   p.foits ) &&
		 ( grid->citer <= ( p.foits + p.soits ) ) )
	{
	  rorder = 1;
	  rorderL = 1;
	  rorderR = 1;
	}
      else
	{
	  rorder = 3;
	  
	  rorderL = grid->node_order[nodeL];
	  rorderR = grid->node_order[nodeR];
	}
    }
  else if ( p.order == 4 )
    {
      if ( grid->citer <= p.foits )
	{
	  rorder = 0;
	  rorderL = 0;
	  rorderR = 0;
	}
      else if ( ( grid->citer >    p.foits && 
		  grid->citer <= ( p.foits + p.soits ) ) )
	{
	  rorder = 1;
	  rorderL = 1;
	  rorderR = 1;
	}
      else if ( ( grid->citer >  ( p.foits + p.soits ) ) && 
		( grid->citer <= ( p.foits + p.soits + p.toits ) ) )
	{
	  rorder = 3;

	  if ( grid->node_order[nodeL] > 2 )
	    rorderL = 3;
	  else
	    rorderL = 1;

	  if ( grid->node_order[nodeR] > 2 )
	    rorderR = 3;
	  else
	    rorderR = 1;
	}
      else
	{
	  rorder = 4;

	  rorderL = grid->node_order[nodeL];
	  rorderR = grid->node_order[nodeR];
	}
    }

  if ( DEBUG_HESS )
    {
      rorder = 3;
      rorderL = 3;
      rorderR = 3;

      for ( j=0; j < NUM_VAR; j++ ) 
	{ 
	  QL[j] = grid->point_value[nodeL*NUM_VAR+j]; 
	  QR[j] = grid->point_value[nodeR*NUM_VAR+j]; 
	} 

    }

  if ( DEBUG_DERIVS ) 
    { 
      rorder = 4;
      rorderL = 4;
      rorderR = 4;
 
      for ( j=0; j < NUM_VAR; j++ )  
        {  
          QL[j] = grid->point_value[nodeL*NUM_VAR+j];  
          QR[j] = grid->point_value[nodeR*NUM_VAR+j];  
        }  
 
    }

  dx[0] = GP[0] - x[nodeL*NDIM+0];
  dx[1] = GP[1] - x[nodeL*NDIM+1];

  //if ( grid->iter <= p.foits ) // do linear reconstruction.
  //  {  return; }

  if ( rorder == 0 )
    return;

  // do the left node first.
  if ( rorderL == 1 )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  //phiL = grid->phi[nodeL*8+j*2+0];
	  QL[j] += ( phiL*grad[nodeL*NDIM*NUM_VAR+j*NDIM+0]*dx[0] +
		     phiL*grad[nodeL*NDIM*NUM_VAR+j*NDIM+1]*dx[1]);
	}
    }
  else if ( rorderL == 3 )
    {
      if ( PV_EXT )
	{
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      QL[j] = grid->point_value[nodeL*NUM_VAR+j];
	    }
	}

      for ( j=0; j < NUM_VAR; j++ )
	{
	  //phiL = grid->phi[nodeL*8+j*2+0];
	  //sigmaL = grid->phi[nodeL*8+j*2+1];
	  QL[j] += ( phiL*hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+0]*dx[0] +
		     phiL*hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+1]*dx[1] +
		     sigmaL*hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+2]*dx[0]*dx[0]*0.5 +
		     sigmaL*hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+3]*dx[1]*dx[1]*0.5 +
		     sigmaL*hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+4]*dx[0]*dx[1] );     
	}
    }
  else if ( rorderL ==  4 )
    {
      if ( PV_EXT )
	{
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      QL[j] = grid->point_value[nodeL*NUM_VAR+j];
	    }
	}

      for ( j=0; j < NUM_VAR; j++ )
	{
	  //phiL = grid->phi[nodeL*8+j*2+0];
	  //sigmaL = grid->phi[nodeL*8+j*2+1];
	  QL[j] += ( phiL*hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+0]*dx[0] +
		     phiL*hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+1]*dx[1] +
		     sigmaL*hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+2]*dx[0]*dx[0]*0.5 +
		     sigmaL*hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+3]*dx[1]*dx[1]*0.5 +
		     sigmaL*hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+4]*dx[0]*dx[1] +
		     sigmaL*hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+5]*dx[0]*dx[0]*dx[0]*(1./6.) +
		     sigmaL*hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+6]*dx[1]*dx[1]*dx[1]*(1./6.) +
		     sigmaL*hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+7]*dx[0]*dx[0]*dx[1]*0.5 +
		     sigmaL*hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+8]*dx[0]*dx[1]*dx[1]*0.5 );
	}
    }
  else
    {
      printf("FATAL ERROR: In Reconstruct_Gauss(): rorderL is undefined <%d>. subedge %d , nodeL = %d, nodeR = %d, time step %d.\n",rorderL,
	     isedge, nodeL, nodeR, ts);
      fflush(stdout);
      exit(0);
    }

  // do the right node.
  
  dx[0] = GP[0] - x[nodeR*NDIM+0];
  dx[1] = GP[1] - x[nodeR*NDIM+1];
  
  if ( rorderR == 1 )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  //phiR = grid->phi[nodeR*8+j*2+0];
	  QR[j] += ( phiR*grad[nodeR*NDIM*NUM_VAR+j*NDIM+0]*dx[0] +
		     phiR*grad[nodeR*NDIM*NUM_VAR+j*NDIM+1]*dx[1]);
	}
    }
  else if ( rorderR == 3 )
    {
      if ( PV_EXT )
	{
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      QR[j] = grid->point_value[nodeR*NUM_VAR+j];
	    }
	}
      
      for ( j=0; j < NUM_VAR; j++ )
	{
	  //phiR = grid->phi[nodeR*8+j*2+0];
	  //sigmaR = grid->phi[nodeR*8+j*2+1];
	  QR[j] += ( phiR*hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+0]*dx[0] +
		     phiR*hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+1]*dx[1] +
		     sigmaR*hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+2]*dx[0]*dx[0]*0.5 +
		     sigmaR*hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+3]*dx[1]*dx[1]*0.5 +
		     sigmaR*hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+4]*dx[0]*dx[1] );     
	}
    }
  else if ( rorderR == 4 )
    {
      if ( PV_EXT )
	{
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      QR[j] = grid->point_value[nodeR*NUM_VAR+j];
	    }
	}

      for ( j=0; j < NUM_VAR; j++ )
	{
	  //phiR = grid->phi[nodeR*8+j*2+0];
	  //sigmaR = grid->phi[nodeR*8+j*2+1];
	  QR[j] += ( phiR*hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+0]*dx[0] +
		     phiR*hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+1]*dx[1] +
		     sigmaR*hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+2]*dx[0]*dx[0]*0.5 +
		     sigmaR*hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+3]*dx[1]*dx[1]*0.5 +
		     sigmaR*hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+4]*dx[0]*dx[1] +
                     sigmaL*hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+5]*dx[0]*dx[0]*dx[0]*(1./6.) +
		     sigmaL*hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+6]*dx[1]*dx[1]*dx[1]*(1./6.) +
		     sigmaL*hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+7]*dx[0]*dx[0]*dx[1]*0.5 +
		     sigmaL*hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+8]*dx[0]*dx[1]*dx[1]*0.5 );     
	}
    }
  else
    {
      printf("FATAL ERROR: In Reconstruct_Gauss(): rorderL is undefined <%d>. subedge %d , nodeR = %d, nodeR = %d, time step %d.\n",rorderR,
	     isedge, nodeL, nodeR, ts);
      fflush(stdout);
      exit(0);
    }
  
  return;
}

//=============================================================
// 
//  Reconstruct_Gauss_UL ()
//
//  Applies the reconstruction operator to the left and right
//  state vectors.
//
//  int isedge;                          // Subedge index.
//  int nodeL;                           // The left node.
//  int nodeR;                           // The right node.
//  double GP[2];                        // Coordinates of the Gauss point.
//  double QL[4];                        // The left state vector.
//  double QR[4];                        // The right state vector.
//  GRID *grid;                          // The grid.
//  PARAMS p;                            // The parameters.
//  int ts;                              // The current time step.
//
//=============================================================

void Reconstruct_Gauss_UL ( int isedge, int nodeL, int nodeR, double GP[2], double QL[4], double QR[4], GRID *grid, PARAMS p, int ts )
{
  int j;                                       // Loop counter.
  int rorder;                                  // Reconstruction order.
  int rorderL, rorderR;                        // Recostruction order to apply to the individual nodes along the edge.
  double dx[2];                                // Edge direction.

  // Pointers.
  double *Q = grid->nQ;
  double *grad = grid->grad;
  double *hess = grid->hess;
  double *x = grid->x;

  // Extract the Q values.
  for ( j=0; j < NUM_VAR; j++ )
    {
      QL[j] = Q[nodeL*NUM_VAR+j];
      QR[j] = Q[nodeR*NUM_VAR+j];
    }

  if ( p.order == 2 || EDGE_BASED )
    {
      printf("CRITICAL ERROR: The solver is trying to reconstruct to Gaussian quadrature nodes but the solve is specified as edge based!\n");
      printf("In Reconstruct_Gauss().\n");
      fflush(stdout);
      exit(1);
    }

  if ( p.order == 0 )
    {
      rorder = 0;
      return;
    }

  if ( p.order == 1 )
    {
      if ( grid->citer <= p.foits )
	{
	  rorder = 0;
	  rorderL = 0;
	  rorderR = 0;
	}
      else
	{
	  rorder = 1;
	  rorderL = 1;
	  rorderR = 1;
	}
    }
  else if ( p.order == 3 )
    {
      if ( grid->citer <= p.foits )
	{
	  rorder = 0;
	  rorderL = 0;
	  rorderR = 0;
	}
      else if (  ( grid->citer  >   p.foits ) &&
		 ( grid->citer <= ( p.foits + p.soits ) ) )
	{
	  rorder = 1;
	  rorderL = 1;
	  rorderR = 1;
	}
      else
	{
	  rorder = 3;
	  
	  rorderL = grid->node_order[nodeL];
	  rorderR = grid->node_order[nodeR];
	}
    }
  else if ( p.order == 4 )
    {
      if ( grid->citer <= p.foits )
	{
	  rorder = 0;
	  rorderL = 0;
	  rorderR = 0;
	}
      else if ( ( grid->citer >    p.foits && 
		  grid->citer <= ( p.foits + p.soits ) ) )
	{
	  rorder = 1;
	  rorderL = 1;
	  rorderR = 1;
	}
      else if ( ( grid->citer >  ( p.foits + p.soits ) ) && 
		( grid->citer <= ( p.foits + p.soits + p.toits ) ) )
	{
	  rorder = 3;

	  if ( grid->node_order[nodeL] > 2 )
	    rorderL = 3;
	  else
	    rorderL = 1;

	  if ( grid->node_order[nodeR] > 2 )
	    rorderR = 3;
	  else
	    rorderR = 1;
	}
      else
	{
	  rorder = 4;

	  rorderL = grid->node_order[nodeL];
	  rorderR = grid->node_order[nodeR];
	}
    }

  if ( DEBUG_HESS )
    {
      rorder = 3;
      rorderL = 3;
      rorderR = 3;

      for ( j=0; j < NUM_VAR; j++ ) 
	{ 
	  QL[j] = grid->point_value[nodeL*NUM_VAR+j]; 
	  QR[j] = grid->point_value[nodeR*NUM_VAR+j]; 
	} 

    }

  if ( DEBUG_DERIVS ) 
    { 
      rorder = 4;
      rorderL = 4;
      rorderR = 4;
 
      for ( j=0; j < NUM_VAR; j++ )  
        {  
          QL[j] = grid->point_value[nodeL*NUM_VAR+j];  
          QR[j] = grid->point_value[nodeR*NUM_VAR+j];  
        }  
 
    }

  dx[0] = GP[0] - x[nodeL*NDIM+0];
  dx[1] = GP[1] - x[nodeL*NDIM+1];


  if ( rorder == 0 )
    return;

  // do the left node first.
  if ( rorderL == 1 )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  QL[j] += ( grad[nodeL*NDIM*NUM_VAR+j*NDIM+0]*dx[0] +
		     grad[nodeL*NDIM*NUM_VAR+j*NDIM+1]*dx[1]);
	}
    }
  else if ( rorderL == 3 )
    {
      if ( PV_EXT )
	{
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      QL[j] = grid->point_value[nodeL*NUM_VAR+j];
	    }
	}

      for ( j=0; j < NUM_VAR; j++ )
	{
	  QL[j] += ( hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+0]*dx[0] +
		     hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+1]*dx[1] +
		     hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+2]*dx[0]*dx[0]*0.5 +
		     hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+3]*dx[1]*dx[1]*0.5 +
		     hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+4]*dx[0]*dx[1] );     
	}
    }
  else if ( rorderL ==  4 )
    {
      if ( PV_EXT )
	{
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      QL[j] = grid->point_value[nodeL*NUM_VAR+j];
	    }
	}

      for ( j=0; j < NUM_VAR; j++ )
	{
	  QL[j] += ( hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+0]*dx[0] +
		     hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+1]*dx[1] +
		     hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+2]*dx[0]*dx[0]*0.5 +
		     hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+3]*dx[1]*dx[1]*0.5 +
		     hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+4]*dx[0]*dx[1] +
		     hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+5]*dx[0]*dx[0]*dx[0]*(1./6.) +
		     hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+6]*dx[1]*dx[1]*dx[1]*(1./6.) +
		     hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+7]*dx[0]*dx[0]*dx[1]*0.5 +
		     hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+8]*dx[0]*dx[1]*dx[1]*0.5 );
	}
    }
  else
    {
      printf("FATAL ERROR: In Reconstruct_Gauss_UL(): rorderL is undefined <%d>. subedge %d , nodeL = %d, nodeR = %d, time step %d.\n",rorderL,
	     isedge, nodeL, nodeR, ts);
      fflush(stdout);
      exit(0);
    }

  // do the right node.
  
  dx[0] = GP[0] - x[nodeR*NDIM+0];
  dx[1] = GP[1] - x[nodeR*NDIM+1];
  
  if ( rorderR == 1 )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  QR[j] += ( grad[nodeR*NDIM*NUM_VAR+j*NDIM+0]*dx[0] +
		     grad[nodeR*NDIM*NUM_VAR+j*NDIM+1]*dx[1]);
	}
    }
  else if ( rorderR == 3 )
    {
      if ( PV_EXT )
	{
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      QR[j] = grid->point_value[nodeR*NUM_VAR+j];
	    }
	}
      
      for ( j=0; j < NUM_VAR; j++ )
	{
	  QR[j] += ( hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+0]*dx[0] +
		     hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+1]*dx[1] +
		     hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+2]*dx[0]*dx[0]*0.5 +
		     hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+3]*dx[1]*dx[1]*0.5 +
		     hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+4]*dx[0]*dx[1] );     
	}
    }
  else if ( rorderR == 4 )
    {
      if ( PV_EXT )
	{
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      QR[j] = grid->point_value[nodeR*NUM_VAR+j];
	    }
	}

      for ( j=0; j < NUM_VAR; j++ )
	{
	  QR[j] += ( hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+0]*dx[0] +
		     hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+1]*dx[1] +
		     hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+2]*dx[0]*dx[0]*0.5 +
		     hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+3]*dx[1]*dx[1]*0.5 +
		     hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+4]*dx[0]*dx[1] +
                     hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+5]*dx[0]*dx[0]*dx[0]*(1./6.) +
		     hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+6]*dx[1]*dx[1]*dx[1]*(1./6.) +
		     hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+7]*dx[0]*dx[0]*dx[1]*0.5 +
		     hess[nodeR*NUM_MOM*NUM_VAR+j*NUM_MOM+8]*dx[0]*dx[1]*dx[1]*0.5 );     
	}
    }
  else
    {
      printf("FATAL ERROR: In Reconstruct_Gauss_UL(): rorderL is undefined <%d>. subedge %d , nodeR = %d, nodeR = %d, time step %d.\n",rorderR,
	     isedge, nodeL, nodeR, ts);
      fflush(stdout);
      exit(0);
    }
  
  return;
}



//=============================================================
// 
//  Reconstruct_Gauss_Boundary ()
//
//  Applies the reconstruction operator to the left
//  state vectors.
//
//  int nodeL;                           // The left node.
//  double GP[2];                        // Coordinates of the Gauss point.
//  double QL[4];                        // The left state vector.
//  GRID *grid;                          // The grid.
//  PARAMS p;                            // The parameters.
//
//=============================================================

void Reconstruct_Gauss_Boundary ( int nodeL, double GP[2], double QL[4], GRID *grid, PARAMS p )
{
  int j;                                       // Loop counter.
  int rorder;                                  // Reconstruction order.
  double dx[2];                                // Edge direction.
  double phi,sigma;
  double xtemp1,xtemp2,xtemp3;
  double ytemp1,ytemp2,ytemp3;

  // Pointers.
  double *Q = grid->nQ;
  double *grad = grid->grad;
  double *hess = grid->hess;
  double *x = grid->x;

  phi = grid->phi[nodeL*8+0*2+0];
  sigma = grid->phi[nodeL*8+0*2+1];

  // Extract the Q values.
  for ( j=0; j < NUM_VAR; j++ )
    {
      QL[j] = Q[nodeL*NUM_VAR+j];
    }
  
  if ( p.order == 0 )
    {
      rorder = 0;
      return;
    }
  
  if ( p.order == 1 )
    {
      if ( grid->citer <= p.foits )
	rorder = 0;
      else
	rorder = 1;
    }
  else if ( p.order == 3 )
    {
      if ( grid->citer <= p.foits )
	rorder = 0;
      else if (  ( grid->citer  >   p.foits ) &&
		 ( grid->citer <= ( p.foits + p.soits ) ) )
	rorder = 1;
      else
	{
	  rorder = grid->node_order[nodeL];
	  
	  if ( rorder > 2 )
	    {
	      if ( PV_EXT )
		{
		  for ( j=0; j < NUM_VAR; j++ )
		    QL[j] = grid->point_value[nodeL*NUM_VAR+j];
		}
	    }
	}
    }
  else if ( p.order == 4 )
    {
      if ( grid->citer <= p.foits )
	rorder = 0;
      else if ( ( grid->citer >    p.foits && 
		  grid->citer <= ( p.foits + p.soits ) ) )
	rorder = 1;
      else if ( ( grid->citer >  ( p.foits + p.soits ) ) && 
		( grid->citer <= ( p.foits + p.soits + p.toits ) ) )
	{
	  if ( grid->node_order[nodeL] > 2 )
	    rorder = 3;
	  else
	    rorder = 1;
	  
	  if (  rorder > 2 )
	    {
	      if ( PV_EXT )
		{
		  for ( j=0; j < NUM_VAR; j++ )
		    {
		      QL[j] = grid->point_value[nodeL*NUM_VAR+j];
		    }
		}
	    }
	}
      else
	{
	  rorder = grid->node_order[nodeL];

	  if ( rorder > 2 )
	    {
	      if ( PV_EXT )
		{
		  for ( j=0; j < NUM_VAR; j++ )
		    {
		      QL[j] = grid->point_value[nodeL*NUM_VAR+j];
		    }
		}
	    }
	}
    }

  if ( DEBUG_HESS ) 
    { 
      rorder = 3; 
      
      for ( j=0; j < NUM_VAR; j++ )  
        {  
          QL[j] = grid->point_value[nodeL*NUM_VAR+j];  
        }  
      
    } 
  
  if ( DEBUG_DERIVS )  
    {  
      rorder = 4;  
      
      for ( j=0; j < NUM_VAR; j++ )   
        {   
          QL[j] = grid->point_value[nodeL*NUM_VAR+j];   
        }   
      
    }
  
  if ( (p.order == 2 || rorder == 2 ) || EDGE_BASED )
    {
      printf("CRITICAL ERROR: The solver is trying to reconstruct to Gaussian quadrature nodes but the solve is specified as edge based!\n");
      printf("In Reconstruct_Gauss_Boundary().\n");
      fflush(stdout);
      exit(1);
    }
  
  dx[0] = GP[0] - x[nodeL*NDIM+0];
  dx[1] = GP[1] - x[nodeL*NDIM+1];

  if ( rorder == 0 )
    return;
  
  if ( rorder == 1 )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  //phi = grid->phi[nodeL*8+j*2+0];
	  QL[j] += ( phi*grad[nodeL*NDIM*NUM_VAR+j*NDIM+0]*dx[0] +
		     phi*grad[nodeL*NDIM*NUM_VAR+j*NDIM+1]*dx[1]);
	}
      
      return;
    }
  else if ( rorder == 3 )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  //phi = grid->phi[nodeL*8+j*2+0];
	  //sigma = grid->phi[nodeL*8+j*2+1];
	  QL[j] += ( phi*hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+0]*dx[0] +
		     phi*hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+1]*dx[1] +
		     sigma*hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+2]*dx[0]*dx[0]*0.5 +
		     sigma*hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+3]*dx[1]*dx[1]*0.5 +
		     sigma*hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+4]*dx[0]*dx[1] );     
	}
      return;
    }
  else // cubic
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  //phi = grid->phi[nodeL*8+j*2+0];
	  //sigma = grid->phi[nodeL*8+j*2+1];
	  QL[j] += ( phi*hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+0]*dx[0] +
		     phi*hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+1]*dx[1] +
		     sigma*hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+2]*dx[0]*dx[0]*0.5 +
		     sigma*hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+3]*dx[1]*dx[1]*0.5 +
		     sigma*hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+4]*dx[0]*dx[1] +
		     sigma*hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+5]*dx[0]*dx[0]*dx[0]*(1./6.) +
		     sigma*hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+6]*dx[1]*dx[1]*dx[1]*(1./6.) +
		     sigma*hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+7]*dx[0]*dx[0]*dx[1]*0.5 +
		     sigma*hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+8]*dx[0]*dx[1]*dx[1]*0.5 );
	}
      return;
    }
  
  return;
}


//=============================================================
// 
//  Reconstruct_Gauss_Boundary_UL ()
//
//  Applies the reconstruction operator to the left
//  state vectors.
//
//  int nodeL;                           // The left node.
//  double GP[2];                        // Coordinates of the Gauss point.
//  double QL[4];                        // The left state vector.
//  GRID *grid;                          // The grid.
//  PARAMS p;                            // The parameters.
//
//=============================================================

void Reconstruct_Gauss_Boundary_UL ( int nodeL, double GP[2], double QL[4], GRID *grid, PARAMS p )
{
  int j;                                       // Loop counter.
  int rorder;                                  // Reconstruction order.
  double dx[2];                                // Edge direction.

  // Pointers.
  double *Q = grid->nQ;
  double *grad = grid->grad;
  double *hess = grid->hess;
  double *x = grid->x;

  // Extract the Q values.
  for ( j=0; j < NUM_VAR; j++ )
    {
      QL[j] = Q[nodeL*NUM_VAR+j];
    }
  
  if ( p.order == 0 )
    {
      rorder = 0;
      return;
    }
  
  if ( p.order == 1 )
    {
      if ( grid->citer <= p.foits )
	rorder = 0;
      else
	rorder = 1;
    }
  else if ( p.order == 3 )
    {
      if ( grid->citer <= p.foits )
	rorder = 0;
      else if (  ( grid->citer  >   p.foits ) &&
		 ( grid->citer <= ( p.foits + p.soits ) ) )
	rorder = 1;
      else
	{
	  rorder = grid->node_order[nodeL];
	  
	  if ( rorder > 2 )
	    {
	      if ( PV_EXT )
		{
		  for ( j=0; j < NUM_VAR; j++ )
		    QL[j] = grid->point_value[nodeL*NUM_VAR+j];
		}
	    }
	}
    }
  else if ( p.order == 4 )
    {
      if ( grid->citer <= p.foits )
	rorder = 0;
      else if ( ( grid->citer >    p.foits && 
		  grid->citer <= ( p.foits + p.soits ) ) )
	rorder = 1;
      else if ( ( grid->citer >  ( p.foits + p.soits ) ) && 
		( grid->citer <= ( p.foits + p.soits + p.toits ) ) )
	{
	  if ( grid->node_order[nodeL] > 2 )
	    rorder = 3;
	  else
	    rorder = 1;
	  
	  if (  rorder > 2 )
	    {
	      if ( PV_EXT )
		{
		  for ( j=0; j < NUM_VAR; j++ )
		    {
		      QL[j] = grid->point_value[nodeL*NUM_VAR+j];
		    }
		}
	    }
	}
      else
	{
	  rorder = grid->node_order[nodeL];

	  if ( rorder > 2 )
	    {
	      if ( PV_EXT )
		{
		  for ( j=0; j < NUM_VAR; j++ )
		    {
		      QL[j] = grid->point_value[nodeL*NUM_VAR+j];
		    }
		}
	    }
	}
    }

  if ( (p.order == 2 || rorder == 2 ) || EDGE_BASED )
    {
      printf("CRITICAL ERROR: The solver is trying to reconstruct to Gaussian quadrature nodes but the solve is specified as edge based!\n");
      printf("In Reconstruct_Gauss_Boundary_UL().\n");
      fflush(stdout);
      exit(1);
    }

  if ( DEBUG_HESS ) 
    { 
      rorder = 3; 
 
      for ( j=0; j < NUM_VAR; j++ )  
        {  
          QL[j] = grid->point_value[nodeL*NUM_VAR+j];  
        }  
 
    } 
 
  if ( DEBUG_DERIVS )  
    {  
      rorder = 4;  
  
      for ( j=0; j < NUM_VAR; j++ )   
        {   
          QL[j] = grid->point_value[nodeL*NUM_VAR+j];   
        }   
  
    }

  dx[0] = GP[0] - x[nodeL*NDIM+0];
  dx[1] = GP[1] - x[nodeL*NDIM+1];

  if ( rorder == 0 )
    return;

  if ( rorder == 1 )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  QL[j] += ( grad[nodeL*NDIM*NUM_VAR+j*NDIM+0]*dx[0] +
		     grad[nodeL*NDIM*NUM_VAR+j*NDIM+1]*dx[1]);
	}
      
      return;
    }
  else if ( rorder == 3 )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  QL[j] += ( hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+0]*dx[0] +
		     hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+1]*dx[1] +
		     hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+2]*dx[0]*dx[0]*0.5 +
		     hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+3]*dx[1]*dx[1]*0.5 +
		     hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+4]*dx[0]*dx[1] );     
	}
      return;
    }
  else // cubic
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  QL[j] += ( hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+0]*dx[0] +
		     hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+1]*dx[1] +
		     hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+2]*dx[0]*dx[0]*0.5 +
		     hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+3]*dx[1]*dx[1]*0.5 +
		     hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+4]*dx[0]*dx[1] +
		     hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+5]*dx[0]*dx[0]*dx[0]*(1./6.) +
		     hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+6]*dx[1]*dx[1]*dx[1]*(1./6.) +
		     hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+7]*dx[0]*dx[0]*dx[1]*0.5 +
		     hess[nodeL*NUM_MOM*NUM_VAR+j*NUM_MOM+8]*dx[0]*dx[1]*dx[1]*0.5 );
	}
      return;
    }
  
  return;
}
