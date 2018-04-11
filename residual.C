//=============================================================
// 
//  residual.C
//  
//  Functions to build the residual.
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

//=============================================================
// 
//  build_interior_residuals()
//
//  Computes the residual over all interior edges/faces.
//
//  GRID *grid;                          // The grid.
//  SYS_MEM *smem;                       // The memory for the system.
//  PARAMS p;                            // The parameters.
//  int jupdate;                         // Jacobian update flag.
//  int ts;                              // The current time step.
//
//=============================================================

void build_interior_residuals ( GRID *grid, SYS_MEM *smem, PARAMS p, int jupdate, int ts )
{
  int i, j, k;                                 // Loop counters.
  int nodeL, nodeR;                            // Edge nodes.
  int gelem;                                   // Global element number.
  int error;                                   // Error code from Flux.
  int isubedge, iedge;                         // Edge tokens.
  double xmid,ymid;                            // Center of the face.
  double xc,yc;                                // Element centroid.
  double nx, ny;                               // Face normal vector.
  double qleft[4];                             // Q on the left of the face.
  double qright[4];                            // Q on the right of the face.
  double flux[4];                              // Flux across the face.
  double len;                                  // Edge length.

  // Variables for higher order flux calculation.
  double x1,x2,x3;                       // Integration points.
  double y1,y2,y3;
  double t1,t2,t3;                       // Gaussian roots.
  double w1,w2,w3;                       // Gaussian coefficients.
  double gp[2];                          // Gauss point.
  double fluxg[4];                       // Flux at the Gauss point.
  int DO_GAUSS_QUAD = 0;

  // Pointers.
  int nn = grid->nn;
  double gamma = p.gamma;
  double *RHS = NULL;                          // Will point to the residual.
  double *RES = grid->R;

  w1 = (5./9.)*0.5;
  w2 = (8./9.)*0.5;
  w3 = w1;

  // Allocate memory full the residual.
  if ( smem->RHS == NULL )
    {
      smem->RHS = (double*)malloc( (nn+1)*NUM_VAR*sizeof(double) );
      if ( smem->RHS == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'RHS'.\n"); exit(1); }
    }

  RHS = smem->RHS;

  if ( grid->R == NULL )
    {
      grid->R = (double*)malloc( (nn+1)*NUM_VAR*sizeof(double) );
      if ( grid->R == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'R' in residual.C .\n"); exit(1); }
    }

  RES = grid->R;

  // Initialize the right hand side and the residual to zero.
  for ( i=1; i <= nn; i++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  RHS[i*NUM_VAR+j] = 0.;
	  RES[i*NUM_VAR+j] = 0.;
	}
    }

  //FILE *fp1 = NULL;
  //fp1 = fopen("debug_gauss_points_solve_density.dat","w");
  //FILE *fp2 = NULL;
  //fp2 = fopen("debug_gauss_points_solve_xmom.dat","w");
  //FILE *fp3 = NULL;
  //fp3 = fopen("debug_gauss_points_solve_ymom.dat","w");
  //FILE *fp4 = NULL;
  //fp4 = fopen("debug_gauss_points_solve_energy.dat","w");
  //double debug_stuff[4][6];

  // Loop over the subedges abd perform the flux calculation.
  for ( isubedge=1; isubedge <= grid->nsubedges; isubedge++ )
    {
      // Extract the real edge this corresponds to.
      iedge = grid->subedges[isubedge*2];

      // Get the normal vector.
      nx = grid->xn_subedges[isubedge*3 + 0];
      ny = grid->xn_subedges[isubedge*3 + 1];
      len = grid->xn_subedges[isubedge*3 + 2];

      // Get the left and right nodes.
      nodeL = grid->edges[iedge*2 + 0];
      nodeR = grid->edges[iedge*2 + 1];

      /*
      if ( nodeL == 60 && ts > p.foits )
	{
	  gelem = 0;
	}

      if ( nodeR == 60 && ts > p.foits )
	{
	  gelem = 0;
	}
      */

      xmid = grid->xm_subedges[isubedge*NDIM+0];
      ymid = grid->xm_subedges[isubedge*NDIM+1];

      // Find the global element.
      gelem = grid->subedges[isubedge*2+1];

      // Get the element centroid.
      xc = grid->el_cent[gelem*2];
      yc = grid->el_cent[gelem*2+1];

      // Get the correct midpoint. General for all edges including those touching 
      // curved boundaries.

      //if ( p.order == 2 && ( ts >= (p.foits+p.soits) ) )
      //if ( p.order == 2 )

      // p.order == 3 will force the 3rd order scheme to do midpoint. WELL NOT ANYMORE!
      //if ( ( (p.order == 3  &&  (ts > (p.foits+p.soits)) )
      //     || FULL_GAUSS_FLUX ) && !EDGE_BASED )  

      // Figure out if we are doing Gaussian Quadrature on the subedges.
      if ( p.order < 3 )
	{
	  DO_GAUSS_QUAD = 0;
	}
      else if ( (p.order == 3) && (ts > (p.foits+p.soits)) )
	{
	  DO_GAUSS_QUAD = 1;
	}
      else if ( (p.order == 4) && (ts > (p.foits+p.soits)) )
	{
	  DO_GAUSS_QUAD = 1;
	}

      if ( FULL_GAUSS_FLUX )
	DO_GAUSS_QUAD = 1;
      
      if ( EDGE_BASED )
	DO_GAUSS_QUAD = 0;

      // HACK !!
      //DO_GAUSS_QUAD = 0;

      if ( CYLINDER && CYLINDER_HACK )
	{
	  if ( (grid->node_order[nodeL] == 1) &&
	       (grid->node_order[nodeR] == 1) )
	    {
	      DO_GAUSS_QUAD = 0;
	    }
	}

	                                            // In other words, the intent is to do Gaussian quadrature on the subedge
                                                    // only if 3rd/4th has been selected and the iteration is appropriate OR if
      if ( DO_GAUSS_QUAD )                          // Gaussian quadrature has been specified BUT NOT if the solver is doing edge based.
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
	  
	  for ( k=0; k < NUM_VAR; k++ )
	    flux[k] = 0.;

	  // For each gauss point, we will reconstruct the solution and then do Roe flux.
	  
	  gp[0] = grid->xg_subedges[isubedge*6+0];
	  gp[1] = grid->xg_subedges[isubedge*6+1];
	  Reconstruct_Gauss ( isubedge, nodeL, nodeR, gp, qleft, qright, grid, p, ts );

	  //debug_stuff[0][0] = qleft[0];
	  //debug_stuff[0][3] = qright[0];
	  //debug_stuff[1][0] = qleft[1];
	  //debug_stuff[1][3] = qright[1];
	  //debug_stuff[2][0] = qleft[2];
	  //debug_stuff[2][3] = qright[2];
	  //debug_stuff[3][0] = qleft[3];
	  //debug_stuff[3][3] = qright[3];
	  
	  error = Roe_flux_centered ( nx, ny, w1*len, gamma, qleft, qright, fluxg );
	  for ( k=0; k < NUM_VAR; k++ )
	    flux[k] = fluxg[k];

	  gp[0] = grid->xg_subedges[isubedge*6+2];
	  gp[1] = grid->xg_subedges[isubedge*6+3];
	  Reconstruct_Gauss ( isubedge, nodeL, nodeR, gp, qleft, qright, grid, p, ts );

	  //debug_stuff[0][1] = qleft[0];
	  //debug_stuff[0][4] = qright[0];
	  //debug_stuff[1][1] = qleft[1];
	  //debug_stuff[1][4] = qright[1];
	  //debug_stuff[2][1] = qleft[2];
	  //debug_stuff[2][4] = qright[2];
	  //debug_stuff[3][1] = qleft[3];
	  //debug_stuff[3][4] = qright[3];
	  
	  error = Roe_flux_centered ( nx, ny, w2*len, gamma, qleft, qright, fluxg );
	  for ( k=0; k < NUM_VAR; k++ )
	    flux[k] += fluxg[k];

	  gp[0] = grid->xg_subedges[isubedge*6+4];
	  gp[1] = grid->xg_subedges[isubedge*6+5];
	  Reconstruct_Gauss ( isubedge, nodeL, nodeR, gp, qleft, qright, grid, p, ts );

	  //debug_stuff[0][2] = qleft[0];
	  //debug_stuff[0][5] = qright[0];
	  //debug_stuff[1][2] = qleft[1];
	  //debug_stuff[1][5] = qright[1];
	  //debug_stuff[2][2] = qleft[2];
	  //debug_stuff[2][5] = qright[2];
	  //debug_stuff[3][2] = qleft[3];
	  //debug_stuff[3][5] = qright[3];

	  error = Roe_flux_centered ( nx, ny, w3*len, gamma, qleft, qright, fluxg );
	  for ( k=0; k < NUM_VAR; k++ )
	    flux[k] += fluxg[k];

	  //for ( k=0; k < 6; k++ )
	  //  {
	  //    fprintf(fp1,"%.15e\n",debug_stuff[0][k]);
	  //    fprintf(fp2,"%.15e\n",debug_stuff[1][k]);
	  //    fprintf(fp3,"%.15e\n",debug_stuff[2][k]);
	  //    fprintf(fp4,"%.15e\n",debug_stuff[3][k]);
	  //  }

	  // end flux

	  // Now I divide by three to do the average of the fluxes.
	  //for ( k=0; k < NUM_VAR; k++ )
	  //  {
	  //    flux[k] = flux[k]*0.5 ;
	  //  }

	  // Extract the Q values for the Jacobian calculation.
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      qleft[j]  = grid->nQ[nodeL*NUM_VAR+j];
	      qright[j] = grid->nQ[nodeR*NUM_VAR+j];
	    }

	  //Reconstruct( isubedge, nodeL, nodeR, qleft, qright, grid, p, ts);

	  // I also need the unpeturbed flux for the jacobian.
	  error = Roe_flux_centered ( nx, ny, len, gamma, qleft, qright, flux );

	}
      else
	{
	  // Get the Q variables and if high order, then reconstruct the solution.
	  Reconstruct( isubedge, nodeL, nodeR, qleft, qright, grid, p, ts);

	  //error = Roe_flux ( nx, ny, len, gamma, qleft, qright, flux );
	  error = Roe_flux_centered ( nx, ny, len, gamma, qleft, qright, flux );

	  // Extract the Q values for the Jacobian calculation.
	  /*
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      qleft[j]  = grid->nQ[nodeL*NUM_VAR+j];
	      qright[j] = grid->nQ[nodeR*NUM_VAR+j];
	    }

	  error = Roe_flux_centered ( nx, ny, len, gamma, qleft, qright, flux );
	  */
	  
	}

      if ( jupdate )
	{
	  compute_numerical_jacobian_roe( grid, smem, p, nodeL,  nodeR, qleft, qright, flux, nx, ny, len );
	  //compute_approximate_jacobian_roe ( grid, smem, p, nodeL, nodeR, qleft, qright, flux, nx, ny, len );
	}
      
      if ( !error )
	{
	  printf("Diagnositics from Roe_flux() in build_residuals.\n");
	  printf("  nodeL = %d\n",nodeL);
	  printf("  nodeR = %d\n",nodeR);
	}
      
      // Add to the 'left' node and subtract from the 'right' node.
      for ( k=0; k < NUM_VAR; k++ )
	{
	  RES[nodeL*NUM_VAR+k] += flux[k];
	  RES[nodeR*NUM_VAR+k] -= flux[k];
	}
    }

  //fclose(fp1);
  //fclose(fp2);
  //fclose(fp3);
  //fclose(fp4);


  return;
}



//=============================================================
// 
//  build_boundary_residuals()
//
//  Computes the residuals on the boundaries.
//
//  GRID *grid;                          // The grid.
//  SYS_MEM *smem;                       // The memory for the system.
//  PARAMS p;                            // The parameters.
//  int jupdate;                         // Boundary jacobian update flag.

//=============================================================

void build_boundary_residuals ( GRID *grid, SYS_MEM *smem, PARAMS p, int jupdate )
{
  int i,j,k,b;                         // Loop counters.
  int real_node, ghost_node;           // Boundary edge nodes.
  int error;
  int bc;                              // Boundary condition.
  double nx,ny,len;                    // Edge normal vector.
  double qleft[4];                     // Q on the left of the face.
  double qright[4];                    // Q on the right of the face.
  double flux[4];                      // Flux across the face.
  double gamma = p.gamma;              // Gamma
  double balance[4];                   // The sum of conserved variables across the boundaries.
  double dfdqL[4][4];                  // Numerical jacobian, temporary storage.
  double pf,uf,vf,rhof;                // Solid surface values.
  double dpdq1,dpdq2,dpdq3,dpdq4;

  double *Q = grid->nQ;
  double *RES = grid->R;
  double *LHS = smem->LHS;
  int *iau = smem->iau;
  
  // Variables for higher order flux.
  double xmid, ymid;                     // Values at face midpoints.
  double x1,x2,x3;                       // Integration points.
  double y1,y2,y3;
  double xL,xR,yL,yR;
  double t1,t2,t3;                       // Gaussian roots.
  double w1,w2,w3;                       // Gaussian coefficients.
  double xi, yi;                         // Node coordinates.
  int nodeL, nodeR;                      // Nodes comprising an edge.
  int s;                                 // Boundary segment.
  double gp[2];                          // Gauss point.
  double qb[4];                          // Q on the boundary (outside values).
  double fluxg[4];                       // Partial flux. At the gauss points.
  int DO_GAUSS_QUAD = 0;

  // Curved boundary stuff.
  double XL[2],XR[2],GP[6];
  double w[3];
  
  w[0] = (8./9.)*0.5;
  w[1] = (5./9.)*0.5;
  w[2] = w[1];

  w1 = (5./9.)*0.5;
  w2 = (8./9.)*0.5;
  w3 = w1;

  // Set the boundary values for this iteration.
  Set_Boundary_Values ( grid,p );

  // Recall that the residual and right hand side vectors have already been set to zero
  // and interior resiudals have been accumulated.

  balance[0] = 0.;
  balance[1] = 0.;
  balance[2] = 0.;
  balance[3] = 0.;
  
  //FILE *fp1 = NULL;
  //fp1 = fopen("debug_gauss_points_solve_density_boundary.dat","w");
  //FILE *fp2 = NULL;
  //fp2 = fopen("debug_gauss_points_solve_xmom_boundary.dat","w");
  //FILE *fp3 = NULL;
  //fp3 = fopen("debug_gauss_points_solve_ymom_boundary.dat","w");
  //FILE *fp4 = NULL;
  //fp4 = fopen("debug_gauss_points_solve_energy_boundary.dat","w");
  //double debug_stuff[4][3];

  // Loop over the boundary edges, do Roe scheme, and accumulate the results to the real node.
  for ( i=1; i <= grid->nbedges; i++ )
    {
      // Grab the data from the structure.
      real_node = grid->bedges[i*5+0];
      ghost_node = grid->bedges[i*5+1] + grid->nn;

      b = grid->bedges[i*5+3];
      bc = grid->bbc[b];

      // Get the normal vector information.
      nx = grid->xn_bedges[i*3+0];
      ny = grid->xn_bedges[i*3+1];
      len = grid->xn_bedges[i*3+2];

      // Add in here the inviscid flux and linearized Jacobian using the CFD2 approach
      // from Dr. Anderson.

      if ( 0 && ( bc == 1 || bc%2 == 1 ) )  // only do for inviscid boundaries.
	{
	  // Zero out the Jacobian.
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      for ( k=0; k < NUM_VAR; k++ )
		{
		  dfdqL[j][k] = 0.;
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

	  // HACK !!
	  //DO_GAUSS_QUAD = 0;

	  if ( CYLINDER && CYLINDER_HACK )
	    {
	      if ( (grid->node_order[real_node] == 1) &&
		   (grid->node_order[ghost_node] == 1) )
		{
		  DO_GAUSS_QUAD = 0;
		}
	    }
	                                                // In other words, the intent is to do Gaussian quadrature on the subedge
	                                                // only if 3rd/4th has been selected and the iteration is appropriate OR if
	  if ( DO_GAUSS_QUAD )                          // Gaussian quadrature has been specified BUT NOT if the solver is doing edge based.
	    {
	      // Get the segment.
	      s = grid->bedges[i*5+4];
	      
	      // Now get the nodes attached to the edge.
	      nodeL = grid->bs[b][s][0];
	      nodeR = grid->bs[b][s][1];
	  
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
	      
		  for ( k=0; k < NUM_VAR; k++ )
		    flux[k] = 0.;
	      
		  gp[0] = grid->xg_bedges[i*6+0];
		  gp[1] = grid->xg_bedges[i*6+1];
		  
		  Reconstruct_Gauss_Boundary ( real_node, gp, qleft, grid, p);

		  if ( RECON_PRIM )
		    {
		      rhof = qleft[0];
		      uf = qleft[1];
		      vf = qleft[2];
		      pf = qleft[3];
		    }
		  else
		    {
		      rhof = qleft[0];
		      uf = qleft[1]/rhof;
		      vf = qleft[2]/rhof;
		      pf = (gamma - 1.0)*( qleft[3] - 0.5*rhof*(uf*uf + vf*vf) );
		    }
		  
		  flux[0] += 0.;
		  flux[1] += nx*pf*len*w1;
		  flux[2] += ny*pf*len*w1;
		  flux[3] += 0.;

		  gp[0] = grid->xg_bedges[i*6+2];
		  gp[1] = grid->xg_bedges[i*6+3];
		  
		  Reconstruct_Gauss_Boundary ( real_node, gp, qleft, grid, p);
		  
		  if ( RECON_PRIM )
		    {
		      rhof = qleft[0];
		      uf = qleft[1];
		      vf = qleft[2];
		      pf = qleft[3];
		    }
		  else
		    {
		      rhof = qleft[0];
		      uf = qleft[1]/rhof;
		      vf = qleft[2]/rhof;
		      pf = (gamma - 1.0)*( qleft[3] - 0.5*rhof*(uf*uf + vf*vf) );
		    }
		  
		  flux[0] += 0.;
		  flux[1] += nx*pf*len*w2;
		  flux[2] += ny*pf*len*w2;
		  flux[3] += 0.;
		  
		  gp[0] = grid->xg_bedges[i*6+4];
		  gp[1] = grid->xg_bedges[i*6+5];
		  
		  Reconstruct_Gauss_Boundary ( real_node, gp, qleft, grid, p);

		  if ( RECON_PRIM )
		    {
		      rhof = qleft[0];
		      uf = qleft[1];
		      vf = qleft[2];
		      pf = qleft[3];
		    }
		  else
		    {
		      rhof = qleft[0];
		      uf = qleft[1]/rhof;
		      vf = qleft[2]/rhof;
		      pf = (gamma - 1.0)*( qleft[3] - 0.5*rhof*(uf*uf + vf*vf) );
		    }
		  
		  flux[0] += 0.;
		  flux[1] += nx*pf*len*w3;
		  flux[2] += ny*pf*len*w3;
		  flux[3] += 0.;

		  // Accumulate the flux.
		  for ( j=0; j < NUM_VAR; j++ )
		    {
		      RES[real_node*NUM_VAR+j] += flux[j];
		      balance[j] += flux[j];
		    }

		  // Compute a low order approximation of the flux jacobian.
		  for ( j=0; j < NUM_VAR; j++ )
		    {
		      qleft[j]  = grid->nQ[real_node*NUM_VAR+j];
		    }

		  if ( RECON_PRIM )
		    {
		      rhof = qleft[0];
		      uf = qleft[1];
		      vf = qleft[2];
		      pf = qleft[3];
		    }
		  else
		    {
		      rhof = qleft[0];
		      uf = qleft[1]/rhof;
		      vf = qleft[2]/rhof;
		      pf = (gamma - 1.0)*( qleft[3] - 0.5*rhof*(uf*uf + vf*vf) );
		    }

		  // Compute the Jacobian.
		  dpdq1 = (gamma-1.)*0.5 * (uf*uf + vf*vf);
		  dpdq2 = -(gamma-1.) * (uf);
		  dpdq3 = -(gamma-1.) * (vf);
		  dpdq4 = gamma - 1.;
		  
		  dfdqL[1][0] = nx*len*dpdq1;
		  dfdqL[1][1] = nx*len*dpdq2;
		  dfdqL[1][2] = nx*len*dpdq3;
		  dfdqL[1][3] = nx*len*dpdq4;
		  
		  dfdqL[2][0] = ny*len*dpdq1;
		  dfdqL[2][1] = ny*len*dpdq2;
		  dfdqL[2][2] = ny*len*dpdq3;
		  dfdqL[2][3] = ny*len*dpdq4;

		  // Store the Jacobian.
		  for ( j=0; j < NUM_VAR; j++ )
		    {
		      for ( k=0; k < NUM_VAR; k++ )
			{
			  LHS[ (iau[real_node])*NUM_VAR*NUM_VAR + j*NUM_VAR + k] += dfdqL[j][k];
			}
		    }
		  
		}
	      else   // Curved boundary.
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
		  //curved_boundary_gauss_points ( grid, bc, XL, XR, GP );
		  
		  //curved_boundary_arclength ( grid, bc, XL, XR, &len );

		  len = grid->xgn_bedges[i*9+2];
		  
		  for ( k=0; k < NUM_VAR; k++ )
		    flux[k] = 0.;

		  for ( j=0; j < 3; j++ )
		    {
		      // Get the normal vector.
		      //curved_boundary_normal_vector ( grid, bc, &(GP[j*NDIM]), &nx, &ny );
		      
		      GP[j*NDIM+0] = grid->xg_bedges[i*6+j*NDIM+0];
		      GP[j*NDIM+1] = grid->xg_bedges[i*6+j*NDIM+1];

		      nx = grid->xgn_bedges[i*9+j*3+0];
		      ny = grid->xgn_bedges[i*9+j*3+1];
		      
		      Reconstruct_Gauss_Boundary ( real_node, &(GP[j*NDIM]), qleft, grid, p);
		      
		      if ( RECON_PRIM )
			{
			  rhof = qleft[0];
			  uf = qleft[1];
			  vf = qleft[2];
			  pf = qleft[3];
			}
		      else
			{
			  rhof = qleft[0];
			  uf = qleft[1]/rhof;
			  vf = qleft[2]/rhof;
			  pf = (gamma - 1.0)*( qleft[3] - 0.5*rhof*(uf*uf + vf*vf) );
			}
		  
		      flux[0] += 0.;
		      flux[1] += nx*pf*len*w[j];
		      flux[2] += ny*pf*len*w[j];
		      flux[3] += 0.;
		    }
		  
		  // Accumulate the flux.
		  for ( j=0; j < NUM_VAR; j++ )
		    {
		      RES[real_node*NUM_VAR+j] += flux[j];
		      balance[j] += flux[j];
		    }

		  // Compute a low order approximation of the flux jacobian.
		  curved_boundary_normal_vector ( grid, bc, XL, &nx, &ny );
		  
		  for ( j=0; j < NUM_VAR; j++ )
		    {
		      qleft[j]  = grid->nQ[real_node*NUM_VAR+j];
		    }

		  if ( RECON_PRIM )
		    {
		      rhof = qleft[0];
		      uf = qleft[1];
		      vf = qleft[2];
		      pf = qleft[3];
		    }
		  else
		    {
		      rhof = qleft[0];
		      uf = qleft[1]/rhof;
		      vf = qleft[2]/rhof;
		      pf = (gamma - 1.0)*( qleft[3] - 0.5*rhof*(uf*uf + vf*vf) );
		    }

		  // Compute the Jacobian.
		  dpdq1 = (gamma-1.)*0.5 * (uf*uf + vf*vf);
		  dpdq2 = -(gamma-1.) * (uf);
		  dpdq3 = -(gamma-1.) * (vf);
		  dpdq4 = gamma - 1.;
		  
		  dfdqL[1][0] = nx*len*dpdq1;
		  dfdqL[1][1] = nx*len*dpdq2;
		  dfdqL[1][2] = nx*len*dpdq3;
		  dfdqL[1][3] = nx*len*dpdq4;
		  
		  dfdqL[2][0] = ny*len*dpdq1;
		  dfdqL[2][1] = ny*len*dpdq2;
		  dfdqL[2][2] = ny*len*dpdq3;
		  dfdqL[2][3] = ny*len*dpdq4;
		  
		  // Store the Jacobian.
		  for ( j=0; j < NUM_VAR; j++ )
		    {
		      for ( k=0; k < NUM_VAR; k++ )
			{
			  LHS[ (iau[real_node])*NUM_VAR*NUM_VAR + j*NUM_VAR + k] += dfdqL[j][k];
			}
		    }
		}
	  
	    }
	  else // Low order, no extrapolation.
	    {
	      // Extact the Q values.
	      for ( j=0; j < NUM_VAR; j++ )
		{
		  qleft[j]  = grid->nQ[real_node*NUM_VAR+j];
		}

	      if ( RECON_PRIM )
		{
		  rhof = qleft[0];
		  uf = qleft[1];
		  vf = qleft[2];
		  pf = qleft[3];
		}
	      else
		{
		  rhof = qleft[0];
		  uf = qleft[1]/rhof;
		  vf = qleft[2]/rhof;
		  pf = (gamma - 1.0)*( qleft[3] - 0.5*rhof*(uf*uf + vf*vf) );
		}
	      
	      flux[0] = 0.;
	      flux[1] = nx*pf*len;
	      flux[2] = ny*pf*len;
	      flux[3] = 0.;
		  
	      // Accumulate the flux.
	      for ( j=0; j < NUM_VAR; j++ )
		{
		  RES[real_node*NUM_VAR+j] += flux[j];
		  balance[j] += flux[j];
		}
		  
	      // Compute the Jacobian.
	      dpdq1 = (gamma-1.)*0.5 * (uf*uf + vf*vf);
	      dpdq2 = -(gamma-1.) * (uf);
	      dpdq3 = -(gamma-1.) * (vf);
	      dpdq4 = gamma - 1.;
	      
	      dfdqL[1][0] = nx*len*dpdq1;
	      dfdqL[1][1] = nx*len*dpdq2;
	      dfdqL[1][2] = nx*len*dpdq3;
	      dfdqL[1][3] = nx*len*dpdq4;
	      
	      dfdqL[2][0] = ny*len*dpdq1;
	      dfdqL[2][1] = ny*len*dpdq2;
	      dfdqL[2][2] = ny*len*dpdq3;
	      dfdqL[2][3] = ny*len*dpdq4;
	      
	      // Store the Jacobian.
	      for ( j=0; j < NUM_VAR; j++ )
		{
		  for ( k=0; k < NUM_VAR; k++ )
		    {
		      LHS[ (iau[real_node])*NUM_VAR*NUM_VAR + j*NUM_VAR + k] += dfdqL[j][k];
		    }
		}
	    }
		  
	  continue;  // Skip the rest of the loop.
	}                                                // End of the CFD 2 style inviscid boundary condition.

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

      // HACK !!
      //DO_GAUSS_QUAD = 0;

      if ( CYLINDER && CYLINDER_HACK )
	{
	  if ( (grid->node_order[real_node] == 1) &&
	       (grid->node_order[ghost_node] == 1) )
	    {
	      DO_GAUSS_QUAD = 0;
	    }
	}
                                                    // In other words, the intent is to do Gaussian quadrature on the subedge
                                                    // only if 3rd/4th has been selected and the iteration is appropriate OR if
      if ( DO_GAUSS_QUAD )                          // Gaussian quadrature has been specified BUT NOT if the solver is doing edge based.
	{
	  // Get the segment.
	  s = grid->bedges[i*5+4];
	  
	  // Now get the nodes attached to the edge.
	  nodeL = grid->bs[b][s][0];
	  nodeR = grid->bs[b][s][1];

	  // Get the node coordinates.
	  xi = grid->x[real_node*NDIM+0];
	  yi = grid->x[real_node*NDIM+1];

	  // Get the edge midpoint.
	  xmid = 0.5*( grid->x[nodeL*NDIM+0] + grid->x[nodeR*NDIM+0] );
	  ymid = 0.5*( grid->x[nodeL*NDIM+1] + grid->x[nodeR*NDIM+1] );

	  if ( bc <= 10 ) // Linear boundary.
	    {
	      // For Gaussian quadrature I need to map the points in t-space
	      // onto my dual edge.

	      // I'm also pretty sure that the to be consistent with the right hand
	      // rule the boundary edge is defined from node0 to node1 on the segment.
	      // The normal should thus point out of the mesh.

	      // To get the x,y points from t space use - x = (xL+xR)/2 + (xR-xL)/2*t
	      //                                          y = (yL+yR)/2 + (yR-yL)/2*t
	      // The left and right nodes depend on which node the node is. Is it on the
	      // beginning or end of the edge.
	  
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
	      
	      // The procedure will be somewhat identical to the one above. However, I will need to adjust the ghost values
	      // at every Gauss point to ensure proper flux treatment.

	      // For each gauss point, we will reconstruct the solution and then do Roe flux.

	      // HACK - this one is really fun. The runs that feature analytical solutions depend on x,y. Get_Boundary_Value()
	      //        uses the node value for x,y. I have to thus change what is stored in grid->x in order to give it the right
	      //        value.

	      if ( (ANNULUS || MMS) || MMS_EXP )
		{
		  t1 = grid->x[real_node*NDIM+0];
		  t2 = grid->x[real_node*NDIM+1];
		}

	      for ( k=0; k < NUM_VAR; k++ )
		flux[k] = 0.;
	      
	      gp[0] = grid->xg_bedges[i*6+0];
	      gp[1] = grid->xg_bedges[i*6+1];
		  
	      Reconstruct_Gauss_Boundary ( real_node, gp, qleft, grid, p);

	      //debug_stuff[0][0] = qleft[0];
	      //debug_stuff[1][0] = qleft[1];
	      //debug_stuff[2][0] = qleft[2];
	      //debug_stuff[3][0] = qleft[3];

	      // Get the new ghost value.
	      
	      // HACK
	      if ( (ANNULUS || MMS) || MMS_EXP )
		{
		  //grid->x[real_node*NDIM+0] = x1;
		  //grid->x[real_node*NDIM+1] = y1;

		  grid->x[real_node*NDIM+0] = gp[0];
		  grid->x[real_node*NDIM+1] = gp[1];
		}

	      Get_Boundary_Value( grid, p, i, qleft, qb );

	      // HACK
	      if ( (ANNULUS || MMS) || MMS_EXP )
		{
		  grid->x[real_node*NDIM+0] = t1;
		  grid->x[real_node*NDIM+1] = t2;
		}

	      error = Roe_flux_centered ( nx, ny, w1*len, gamma, qleft, qb, fluxg );
	      for ( k=0; k < NUM_VAR; k++ )
		flux[k] = fluxg[k];

	      gp[0] = grid->xg_bedges[i*6+2];
	      gp[1] = grid->xg_bedges[i*6+3];
		  
	      Reconstruct_Gauss_Boundary ( real_node, gp, qleft, grid, p);

	      //debug_stuff[0][1] = qleft[0];
	      //debug_stuff[1][1] = qleft[1];
	      //debug_stuff[2][1] = qleft[2];
	      //debug_stuff[3][1] = qleft[3];

	      // Get the new ghost value.

	      // HACK
	      if ( (ANNULUS || MMS) || MMS_EXP )
		{
		  //grid->x[real_node*NDIM+0] = x2;
		  //grid->x[real_node*NDIM+1] = y2;

		  grid->x[real_node*NDIM+0] = gp[0];
		  grid->x[real_node*NDIM+1] = gp[1];
		}
	      
	      Get_Boundary_Value( grid, p, i, qleft, qb );

	      // HACK
	      if ( (ANNULUS || MMS) || MMS_EXP )
		{
		  grid->x[real_node*NDIM+0] = t1;
		  grid->x[real_node*NDIM+1] = t2;
		}

	      error = Roe_flux_centered ( nx, ny, w2*len, gamma, qleft, qb, fluxg );
	      for ( k=0; k < NUM_VAR; k++ )
		flux[k] += fluxg[k];
	      
	      gp[0] = grid->xg_bedges[i*6+4];
	      gp[1] = grid->xg_bedges[i*6+5];
		  
	      Reconstruct_Gauss_Boundary ( real_node, gp, qleft, grid, p);

	      //debug_stuff[0][2] = qleft[0];
	      //debug_stuff[1][2] = qleft[1];
	      //debug_stuff[2][2] = qleft[2];
	      //debug_stuff[3][2] = qleft[3];
	      
	      // Get the new ghost value.

	      // HACK
	      if ( (ANNULUS || MMS) || MMS_EXP )
		{
		  //grid->x[real_node*NDIM+0] = x3;
		  //grid->x[real_node*NDIM+1] = y3;

		  grid->x[real_node*NDIM+0] = gp[0];
		  grid->x[real_node*NDIM+1] = gp[1];
		}
	      
	      Get_Boundary_Value( grid, p, i, qleft, qb );

	      // HACK
	      if ( (ANNULUS || MMS) || MMS_EXP )
		{
		  grid->x[real_node*NDIM+0] = t1;
		  grid->x[real_node*NDIM+1] = t2;
		}
	      
	      error = Roe_flux_centered ( nx, ny, w3*len, gamma, qleft, qb, fluxg );
	      for ( k=0; k < NUM_VAR; k++ )
		flux[k] += fluxg[k];
	      
	      // end flux
	      
	      // Now I divide by three to do the average of the fluxes.
	      for ( k=0; k < NUM_VAR; k++ )
		{
		  //flux[k] = flux[k] / (3.*2.) ;  // 3 for the average and the 2 comes from the mapping of -1 to 1 to the edge.
		  //flux[k] = flux[k]*0.5 ;
		}

	      //for ( k=0; k < 3; k++ )
	      //  {
	      //  fprintf(fp1,"%.15e\n",debug_stuff[0][k]);
	      //  fprintf(fp2,"%.15e\n",debug_stuff[1][k]);
	      //  fprintf(fp3,"%.15e\n",debug_stuff[2][k]);
	      //  fprintf(fp4,"%.15e\n",debug_stuff[3][k]);
	      //	}
	    }
	  else   // Curved boundary.
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
	      //curved_boundary_gauss_points ( grid, bc, XL, XR, GP );

	      //curved_boundary_arclength ( grid, bc, XL, XR, &len );

	      len = grid->xgn_bedges[i*9+2];

	      for ( k=0; k < NUM_VAR; k++ )
		flux[k] = 0.;

	      // HACK - this one is really fun. The runs that feature analytical solutions depend on x,y. Get_Boundary_Value()
	      //        uses the node value for x,y. I have to thus change what is stored in grid->x in order to give it the right
	      //        value.

	      if ( (ANNULUS || MMS) || MMS_EXP )
		{
		  t1 = grid->x[real_node*NDIM+0];
		  t2 = grid->x[real_node*NDIM+1];
		}

	      for ( j=0; j < 3; j++ )
		{
		  // Get the normal vector.
		  //curved_boundary_normal_vector ( grid, bc, &(GP[j*NDIM]), &nx, &ny );

		  GP[j*NDIM+0] = grid->xg_bedges[i*6+j*NDIM+0];
		  GP[j*NDIM+1] = grid->xg_bedges[i*6+j*NDIM+1];

		  nx = grid->xgn_bedges[i*9+j*3+0];
		  ny = grid->xgn_bedges[i*9+j*3+1];

		  Reconstruct_Gauss_Boundary ( real_node, &(GP[j*NDIM]), qleft, grid, p);

		  // HACK - I'm adjusting the normal vector stored in xn_bedges so that Get_Boundary_Value
		  //        will have the normal vector at the Gaussian point.
		  grid->xn_bedges[i*3+0] = nx;
		  grid->xn_bedges[i*3+1] = ny;

		  // HACK
		  if ( (ANNULUS || MMS) || MMS_EXP )
		    {
		      grid->x[real_node*NDIM+0] = GP[j*NDIM+0];
		      grid->x[real_node*NDIM+1] = GP[j*NDIM+1];
		    }
		  
		  // Get the new ghost value.
		  Get_Boundary_Value( grid, p, i, qleft, qb );

		  // HACK
		  if ( (ANNULUS || MMS) || MMS_EXP )
		    {
		      grid->x[real_node*NDIM+0] = t1;
		      grid->x[real_node*NDIM+1] = t2;
		    }
		  
		  error = Roe_flux_centered ( nx, ny, w[j]*len, gamma, qleft, qb, fluxg );
		  for ( k=0; k < NUM_VAR; k++ )
		    flux[k] += fluxg[k];
		}
	      
	      // Now I divide by three to do the average of the fluxes.
	      for ( k=0; k < NUM_VAR; k++ )
		{
		  //flux[k] = flux[k]*0.5 ;
		}

	      // Last thing we do is set up for the approximate jacobian.
	      // Reset the normal vector for the Jacobian calculation below. The arclength should be correct.
	      XL[0] = grid->x[real_node*NDIM+0];  XL[1] = grid->x[real_node*NDIM+1];

	      curved_boundary_normal_vector ( grid, bc, XL, &nx, &ny );

	      // HACK part 2 - reset the normal vector information in xn_bedges with normal vector at the corresponding
	      //               real node.
	      grid->xn_bedges[i*3+0] = nx;
	      grid->xn_bedges[i*3+1] = ny;
	      
	    }
	  
	  // Extract the Q values for the Jacobian calculation.
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      qleft[j]  = grid->nQ[real_node*NUM_VAR+j];
	      qright[j] = grid->nQ[ghost_node*NUM_VAR+j];
	    }

	  // It also needs the unperturbed flux.
	  error = Roe_flux_centered ( nx, ny, len, gamma, qleft, qright, flux );
	}
      else // Low order, no extrapolation.
	{
	  // Extact the Q values.
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      qleft[j]  = grid->nQ[real_node*NUM_VAR+j];
	      qright[j] = grid->nQ[ghost_node*NUM_VAR+j];
	    }
	  
	  error = Roe_flux_centered ( nx, ny, len, gamma, qleft, qright, flux );
	}

      if ( !error )
	{
	  printf("Diagnositics from Roe_flux() in build_boundary_residuals.\n");
	  printf("  nodeL = %d\n",real_node);
	  printf("  nodeR = %d\n",ghost_node);
	  printf("  boundary edge = %d\n",i);
	}

      if ( jupdate )
	{
	  //compute_boundary_numerical_jacobian_roe( grid, smem, p, i, qleft, flux, nx, ny, len );

	  // If I did curved boundaries above the normal vector may be incorrect, so generate the correct
	  // normal vector at the node. Or have the code above do it!

	  compute_boundary_approximate_jacobian_roe( grid, smem, p, i, qleft, qright, flux, nx, ny, len );
	}
      
      // Add to the real node.
      for ( j=0; j < NUM_VAR; j++ )
	{
	  RES[real_node*NUM_VAR+j] += flux[j];
	  balance[j] += flux[j];
	}
    }

  //fclose(fp1);
  //fclose(fp2);
  //fclose(fp3);
  //fclose(fp4);

  // Print the global conservation.
  printf("  Total flux through boundaries:\n");
  printf("    density = %.15E\n",balance[0]); 
  printf("    x-momentum = %.15E\n",balance[1]);
  printf("    y-momentum = %.15E\n",balance[2]);
  printf("    Total Energy = %.15E\n",balance[3]);
  fflush(stdout);
  
  if ( 0 && grid->citer == 1 )  // use grid->citer to get the current time step.
    {
      FILE *fp = fopen("LHS_bdiag.dat","w");
      for ( i=1; i <= grid->nbedges; i++ )
	{
	  real_node = grid->bedges[i*5+0];
	  j = smem->iau[real_node];
	  fprintf(fp,"Node %d, boundary condition = %d\n",real_node,grid->bbc[grid->bedges[i*5+3]]);
	  for ( k=0; k < NUM_VAR; k++ )
	    {
	      for ( b=0; b < NUM_VAR; b++ )
		{
		  fprintf(fp,"%.15E   ",smem->LHS[j*NUM_VAR*NUM_VAR+k*NUM_VAR+b]);
		}
	      fprintf(fp,"\n");
	    }
	  fprintf(fp,"\n\n");
	}
      fclose(fp);
    }

  return;
}


//===================================================================================================================
// NO CURVED BOUNDARY SUPPORT BELOW THIS POINT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//=============================================================
// 
//  fix_farfield_boundary()
//
//  Sets the farfield boundary so the solution is not updated on them.
//
//  GRID *grid;                          // The grid.
//  SYS_MEM *smem;                       // The memory for the system.
//  PARAMS p;                            // The parameters.
//
//=============================================================

void fix_farfield_boundary ( GRID *grid, SYS_MEM *smem, PARAMS p )
{
  int i,j,k,b;                         // Loop counters.
  int real_node;                       // Boundary edge node.
  int bc;                              // Boundary condition.

  // Pointers.
  double *LHS = smem->LHS;
  int *ia = smem->ia;
  int *iau = smem->iau;
  double *RES = grid->R;

  // Loop over the boundary edges and compute the flux and jacobians.
  for ( i=1; i <= grid->nbedges; i++ )
    {
      // Grab the data from the structure.
      real_node = grid->bedges[i*5+0];

      b = grid->bedges[i*5+3];
      bc = grid->bbc[b];

      if ( bc == 1 || bc%2 == 1 )  { continue; }

      // First zero out the flux.
      for ( j=0; j < NUM_VAR; j++ )
	RES[real_node*NUM_VAR+j] = 0.;

      // Now kill the linear system.
      for ( j=ia[real_node]; j < ia[real_node+1]; j++ )
	{
	  for ( k=0; k < NUM_VAR*NUM_VAR; k++ )
	    LHS[j*NUM_VAR*NUM_VAR+k] = 0.;
	}

      // Put the identity matrix back on the diagonal block.
      for ( j=0; j < NUM_VAR; j++ )
	{
	  LHS[iau[real_node] * NUM_VAR * NUM_VAR + j*NUM_VAR + j] = 1.;
	}
      
    }

  return;
}


//=============================================================
// 
//  do_hard_boundary_flux_jac()
//
//  Compute the flux and Jacobians for farfield and solid surface
//  boundaries in the style of CFD 2 (Anderson).
//
//  GRID *grid;                          // The grid.
//  SYS_MEM *smem;                       // The memory for the system.
//  PARAMS p;                            // The parameters.
//
//=============================================================

void do_hard_boundary_flux_jac ( GRID *grid, SYS_MEM *smem, PARAMS p )
{
  int i,j,k,b;                         // Loop counters.
  int real_node;                       // Boundary edge node.
  int bc;                              // Boundary condition.
  double nx,ny,len;                    // Edge normal vector.
  double qleft[4];                     // Q on the left of the face.
  double flux[4];                      // Flux across the face.
  double gamma = p.gamma;              // Gamma
  double dfdqL[4][4];                  // Numerical jacobian, temporary storage.
  double pf,uf,vf,rhof;                // Solid surface values.
  double dpdq1,dpdq2,dpdq3,dpdq4;

  // Pointers.
  double *Q = grid->nQ;
  double *LHS = smem->LHS;
  int *ia = smem->ia;
  int *iau = smem->iau;
  double *RES = grid->R;

  // Loop over the boundary edges and set the freestream values for farfield boundaries.
  for ( i=1; i <= grid->nbedges; i++ )
    {
      // Zero out the Jacobian.
      for ( j=0; j < NUM_VAR; j++ )
	{
	  for ( k=0; k < NUM_VAR; k++ )
	    {
	      dfdqL[j][k] = 0.;
	    }
	}

      // Grab the data from the structure.
      real_node = grid->bedges[i*5+0];

      b = grid->bedges[i*5+3];
      bc = grid->bbc[b];

      // Get the normal vector information.
      nx = grid->xn_bedges[i*3+0];
      ny = grid->xn_bedges[i*3+1];
      len = grid->xn_bedges[i*3+2];

      // Extact the Q values.
      for ( j=0; j < NUM_VAR; j++ )
	{
	  qleft[j]  = Q[real_node*NUM_VAR+j];
	}

      // Solid surface.
      if ( bc == 1 )
	{
	  rhof = qleft[0];
	  uf = qleft[1]/rhof;
	  vf = qleft[2]/rhof;
	  pf = (gamma - 1.0)*( qleft[3] - 0.5*rhof*(uf*uf + vf*vf) );

	  flux[0] = 0.;
	  flux[1] = nx*pf*len;
	  flux[2] = ny*pf*len;
	  flux[3] = 0.;

	  // Accumulate the flux.
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      RES[real_node*NUM_VAR+j] += flux[j];
	    }
	  
	  // Compute the Jacobian.
	  dpdq1 = (gamma-1.)*0.5 * (uf*uf + vf*vf);
	  dpdq2 = -(gamma-1.) * (uf);
	  dpdq3 = -(gamma-1.) * (vf);
	  dpdq4 = gamma - 1.;

	  dfdqL[1][0] = nx*len*dpdq1;
	  dfdqL[1][1] = nx*len*dpdq2;
	  dfdqL[1][2] = nx*len*dpdq3;
	  dfdqL[1][3] = nx*len*dpdq4;

	  dfdqL[2][0] = ny*len*dpdq1;
	  dfdqL[2][1] = ny*len*dpdq2;
	  dfdqL[2][2] = ny*len*dpdq3;
	  dfdqL[2][3] = ny*len*dpdq4;

	  // Store the Jacobian.
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      for ( k=0; k < NUM_VAR; k++ )
		{
		  LHS[ (iau[real_node])*NUM_VAR*NUM_VAR + j*NUM_VAR + k] += dfdqL[j][k];
		}
	    }
	}
      
      if ( bc == 0 )
	{
	  // First zero out the flux.
	  for ( j=0; j < NUM_VAR; j++ )
	    RES[real_node*NUM_VAR+j] = 0.;

	  // Now kill the linear system.
	  for ( j=ia[real_node]; j < ia[real_node+1]; j++ )
	    {
	      for ( k=0; k < NUM_VAR*NUM_VAR; k++ )
		LHS[j*NUM_VAR*NUM_VAR+k] = 0.;
	    }

	  // Put the identity matrix back on the diagonal block.
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      LHS[iau[real_node] * NUM_VAR * NUM_VAR + j*NUM_VAR + j] = 1.;
	    }
	}
    }

  return;
}


//=============================================================
// 
//  do_soft_boundary_flux_jac()
//
//  Compute the flux and Jacobians for farfield and solid surface
//  boundaries in the style of CFD 2 (Anderson) for the solid surface
//  and Roe flux on the farfield (with Qr = infinity conditions ).
//
//  GRID *grid;                          // The grid.
//  SYS_MEM *smem;                       // The memory for the system.
//  PARAMS p;                            // The parameters.
//
//=============================================================

void do_soft_boundary_flux_jac ( GRID *grid, SYS_MEM *smem, PARAMS p )
{
  int i,j,k,b;                         // Loop counters.
  int real_node, ghost_node;           // Boundary edge nodes.
  int error;
  int bc;                              // Boundary condition.
  double nx,ny,len;                    // Edge normal vector.
  double qleft[4];                     // Q on the left of the face.
  double qright[4];                    // Q on the right of the face.
  double Qp[4];                        // Perturbed Q state on the left.
  double flux[4];                      // Flux across the face.
  double fluxp[4];                     // Perturbed flux.
  double gamma = p.gamma;              // Gamma
  double dfdqL[4][4];                  // Numerical jacobian, temporary storage.
  double uinf,vinf,rhoinf,pinf;        // Values at the freestream.
  double pf,uf,vf,rhof;                // Solid surface values.
  double dpdq1,dpdq2,dpdq3,dpdq4;
  double angle;                        // Angle of attack.
  double eps = 1.0E-08;

  // Pointers.
  double *Q = grid->nQ;
  double *LHS = smem->LHS;
  int *iau = smem->iau;
  double *RES = grid->R;

  // Set free stream (infinity) values.
  angle = p.alpha*M_PI/180.;

  rhoinf = 1.0;
  uinf = p.mach*cos(angle);
  vinf = p.mach*sin(angle);
  pinf = 1. / p.gamma;

  // Use these as the right state.
  qright[0] = rhoinf;
  qright[1] = rhoinf*uinf;
  qright[2] = rhoinf*vinf;
  qright[3] = pinf / (gamma-1.) + 0.5*rhoinf*(uinf*uinf + vinf*vinf);

  // Loop over the boundary edges and set the freestream values for farfield boundaries.
  for ( i=1; i <= grid->nbedges; i++ )
    {
      // Zero out the Jacobian.
      for ( j=0; j < NUM_VAR; j++ )
	{
	  for ( k=0; k < NUM_VAR; k++ )
	    {
	      dfdqL[j][k] = 0.;
	    }
	}

      // Grab the data from the structure.
      real_node = grid->bedges[i*5+0];
      ghost_node = grid->bedges[i*5+1] + grid->nn;

      b = grid->bedges[i*5+3];
      bc = grid->bbc[b];

      // Get the normal vector information.
      nx = grid->xn_bedges[i*3+0];
      ny = grid->xn_bedges[i*3+1];
      len = grid->xn_bedges[i*3+2];

      // Extact the Q values.
      for ( j=0; j < NUM_VAR; j++ )
	{
	  qleft[j]  = Q[real_node*NUM_VAR+j];
	}

      // Solid surface.
      if ( bc == 1 )
	{
	  rhof = qleft[0];
	  uf = qleft[1]/rhof;
	  vf = qleft[2]/rhof;
	  pf = (gamma - 1.0)*( qleft[3] - 0.5*rhof*(uf*uf + vf*vf) );

	  flux[0] = 0.;
	  flux[1] = nx*pf*len;
	  flux[2] = ny*pf*len;
	  flux[3] = 0.;

	  // Accumulate the flux.
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      RES[real_node*NUM_VAR+j] += flux[j];
	    }
	  
	  // Compute the Jacobian.
	  dpdq1 = (gamma-1.)*0.5 * (uf*uf + vf*vf);
	  dpdq2 = -(gamma-1.) * (uf);
	  dpdq3 = -(gamma-1.) * (vf);
	  dpdq4 = gamma - 1.;

	  dfdqL[1][0] = nx*len*dpdq1;
	  dfdqL[1][1] = nx*len*dpdq2;
	  dfdqL[1][2] = nx*len*dpdq3;
	  dfdqL[1][3] = nx*len*dpdq4;

	  dfdqL[2][0] = ny*len*dpdq1;
	  dfdqL[2][1] = ny*len*dpdq2;
	  dfdqL[2][2] = ny*len*dpdq3;
	  dfdqL[2][3] = ny*len*dpdq4;

	  // Store the Jacobian.
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      for ( k=0; k < NUM_VAR; k++ )
		{
		  LHS[ (iau[real_node])*NUM_VAR*NUM_VAR + j*NUM_VAR + k] += dfdqL[j][k];
		}
	    }
	}
      
      if ( bc == 0 )
	{
	  // Get the flux.
	  error = Roe_flux_centered ( nx, ny, len, gamma, qleft, qright, flux );

	  if ( !error )
	    {
	      printf("Diagnositics from Roe_flux() in build_boundary_residuals.\n");
	      printf("  nodeL = %d\n",real_node);
	      printf("  nodeR = %d\n",ghost_node);
	      printf("  boundary edge = %d\n",i);
	    }

	  // Accumulate the flux.
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      RES[real_node*NUM_VAR+j] += flux[j];
	    }

	  // Now lets build the flux jacobian leaving the right state unperturbed.
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      // Set the perturbed state vector.
	      for ( k=0; k < NUM_VAR; k++ )
		{
		  Qp[k] = qleft[k];
		}
	      
	      Qp[j] += eps;

	      // Get the numerical flux based on the perturbed state.
	      error = Roe_flux_centered ( nx, ny, len, gamma, Qp, qright, fluxp );

	      // Store the result in the ith column of dfdqL.
	      for ( k=0; k < NUM_VAR; k++ )
		{
		  dfdqL[k][j] = ( fluxp[k] - flux[k] ) / eps;
		}
	    }

	  // Now I can store the result in the matrix.
	  
	  // Now add dfdqL to the diagonal of the physical boundary node.
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      for ( k=0; k < NUM_VAR; k++ )
		{
		  LHS[ (iau[real_node])*NUM_VAR*NUM_VAR + j*NUM_VAR + k] += dfdqL[j][k];
		}
	    }
	}
    }

  return;
}

