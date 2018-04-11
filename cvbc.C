//=============================================================
// 
//  cvbc.C
//  
//  Functions to apply characteristic value boudary conditions
//  to the grid.
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
#include "cvbc.h"
#include "flux.h"
#include "mesh_utils.h"


//=============================================================
// 
//  Set_Boundary_Values()
//
//  Fills in the Q values for the ghost nodes using CVBCs.
//  The method comes from the summer course notes at the SimCenter
//  which are bases on AIAA 1984-1552 by Whitfield and Janus.
//
//  GRID *grid;                          // The grid.
//  PARAMS p;                            // The parameters.
//
//=============================================================

void Set_Boundary_Values ( GRID *grid, PARAMS p )
{
  int i,b,s;                         // Loop counters and generic entities.
  int real_node;                     // Real node on the edge.
  int ghost_node;                    // Ghost node on the edge.
  int bc;                            // Edge boundary condition.

  double uinf,vinf,rhoinf,pinf;      // Values at the freestream.
  double minf;
  double cb;
  double u0,v0,rho0,p0,c0;           // Values at the reference state.
  double theta,ci;                   // Nodal values.
  double rhoi,ui,vi,pi;
  double rhob,ub,vb,pb;              // Boundary values.
  double ev1,ev3,ev4;                // Eigenvalues.

  double nx,ny;                      // Normal vector of edge.

  double angle;                      // Angle of attack

  double x1,y1;
  double PI = M_PI;
  double gm1, gm1inv;
  double r, U;

  angle = p.alpha*M_PI/180.;

  // Get the freestream values.
  minf = p.mach;
  rhoinf = 1.0;
  uinf = p.mach*cos(angle);
  vinf = p.mach*sin(angle);
  pinf = 1. / p.gamma;

  //if ( p.ic == 2 )
  //  pinf = 1.0;

  // Loop over the boundary edges and set the boundary data based on the boundary condition applied.
  for ( i=1; i <= grid->nbedges; i++ )
    {

      // Grab the data from the structure.
      real_node = grid->bedges[i*5+0];
      ghost_node = grid->bedges[i*5+1];
      b = grid->bedges[i*5+3];
      s = grid->bedges[i*5+4];
      
      // Determine the BC.
      bc = grid->bbc[b];

      // Get the normal vector information.
      nx = grid->xn_bedges[i*3+0];
      ny = grid->xn_bedges[i*3+1];

      // Get the local primitive values.

      if ( RECON_PRIM )
	{
	  rhoi = grid->nQ[real_node*NUM_VAR+0];
	  ui =   grid->nQ[real_node*NUM_VAR+1];
	  vi =   grid->nQ[real_node*NUM_VAR+2];
	  pi =   grid->nQ[real_node*NUM_VAR+3];
	  ci = sqrt( (p.gamma*pi)/rhoi );
	}
      else
	{
	  rhoi = grid->nQ[real_node*NUM_VAR+0];
	  ui = grid->nQ[real_node*NUM_VAR+1] / rhoi;
	  vi = grid->nQ[real_node*NUM_VAR+2] / rhoi;
	  pi = (p.gamma - 1.) * ( grid->nQ[real_node*NUM_VAR+3] - 0.5*rhoi*(ui*ui + vi*vi) );
	  ci = sqrt( (p.gamma*pi)/rhoi );
	}

      // Get the local boundary values from the previous state.

      if ( RECON_PRIM )
	{
	  rhob = grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +0 ];
	  ub =   grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +1 ];
	  vb =   grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +2 ];
	  pb =   grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +3 ];
	  cb = sqrt( (p.gamma*pb)/rhob );
	}
      else
	{
	  rhob = grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +0 ];
	  ub = grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +1 ] / rhob;
	  vb = grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +2 ] / rhob;
	  pb = (p.gamma - 1.) * ( grid->nQ[ (grid->nn + ghost_node )*NUM_VAR +3] - 0.5*rhob*(ub*ub + vb*vb) );
	  cb = sqrt( (p.gamma*pb)/rhob );
	}

      // Get the new reference state.

      if ( ANNULUS )
	{
	  //ci = pow( rhoi , (0.5*(p.gamma-1.0)) );
	  //cb = pow( rhob , (0.5*(p.gamma-1.0)) );
	}

      rho0 = 0.5*(rhoi+rhob);
      u0 = 0.5*(ui+ub);
      v0 = 0.5*(vi+vb);
      p0 = 0.5*(pi+pb);
      c0 = 0.5*(ci+cb);

      // Get the local eigenvalues.
      theta = ui*nx + vi*ny;
      ev1 = theta;
      ev3 = theta + ci;
      ev4 = theta - ci;

      // debug the flux integral.
      if ( MMS || (MMS && ANNULUS) )
	{
	  // Every boundary edge will be treated as farfield.
	  x1 = grid->x[real_node*NDIM+0];
	  y1 = grid->x[real_node*NDIM+1];

	  //rhob = 1.0 + sin(PI*x1)*sin(PI*y1);
	  rhob = 1.0 + 0.25*sin(PI*x1)*sin(PI*y1);
	  ub = 0.25 + 0.25*sin(PI*x1)*cos(2.*PI*y1);
	  vb = 0.25 + 0.25*cos(2.*PI*x1)*sin(PI*y1);
	  pb = 1./p.gamma + 0.05*cos(2.*PI*x1)*cos(2.*PI*y1);


	  if ( RECON_PRIM )
	    {
	      grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +0 ] = rhob;
	      grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +1 ] = ub;
	      grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +2 ] = vb;
	      grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +3 ] = pb;
	    }
	  else
	    {
	      grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +0 ] = rhob;
	      grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +1 ] = rhob*ub;
	      grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +2 ] = rhob*vb;
	      grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +3 ] = pb/(p.gamma-1.0) + 0.5*rhob*(ub*ub + vb*vb);
	    }

	  continue;
	}

      if ( MMS_EXP )
	{
	  // Every boundary edge will be treated as farfield.
	  x1 = grid->x[real_node*NDIM+0];
	  y1 = grid->x[real_node*NDIM+1];
	  
	  rhob = exp(x1-1.0) * exp(y1 - 1.0);
	  ub   = exp(x1-1.0) * exp(y1 - 1.0);
	  vb   = exp(x1-1.0) * exp(y1 - 1.0);
	  pb   = ( (1.4-1.0) * 0.5 ) * exp( 3.0 * ( x1 - 1.0) ) * exp( 3.0 * ( y1 - 1.0 ) );
	  
	  if ( RECON_PRIM )
	    {
	      grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +0 ] = rhob;
	      grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +1 ] = ub;
	      grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +2 ] = vb;
	      grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +3 ] = pb;
	    }
	  else
	    {
	      grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +0 ] = rhob;
	      grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +1 ] = rhob*ub;
	      grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +2 ] = rhob*vb;
	      grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +3 ] = pb/(p.gamma-1.0) + 0.5*rhob*(ub*ub+vb*vb);
	    }

	  continue;
	}

      // CODE FOR THE ANNULUS - supersonic vortex problem from Aftosmis.
      if ( ANNULUS )
	{
	  gm1 = p.gamma - 1.0;
	  gm1inv = 1.0 / gm1;

	  minf = 2.0;
	  rhoinf = 1.0;
	  
	  if ( bc == 0 )
	    {
	      r = grid->x[real_node*NDIM + 1];
	      U = (minf*2.)/r;

	      rhob = pow( (1.0 + gm1*2.0*(1.0 - 4.0/(r*r))) , gm1inv );
	      ub = U;
	      vb = 0.;
	      pb = pow( rhob, p.gamma ) / p.gamma;
	    }
	  else if ( bc == 2 )
	    {
	      r = grid->x[real_node*NDIM + 0];
	      U = (minf*2.)/r;

	      if ( MMS )
		{
		  rhob = pow( (1.0 + gm1*2.0*(1.0 - 4.0/(r*r))) , gm1inv );
		  ub = 0.;
		  vb = -U;
		  pb = pow( rhob, p.gamma ) / p.gamma;
		}
	      else
		{
		  rhob = rhoi;
		  ub = ui;
		  vb = vi;
		  pb = pi;
		}
	    }
	  else if ( (bc == 13 || bc == 14) || bc == 1 ) // inner radius
	    {
	      r = 2.0;
	      U = minf;
	      
	      //if ( MMS || 0 )
	      //if ( grid->citer >= 6000 )
	      if ( MMS )
		{
		  rhob = pow( (1.0 + gm1*2.0*(1.0 - 4.0/(r*r))) , gm1inv );
		  ub =  (grid->x[real_node*NDIM+1]*U) / r;
		  vb = -(grid->x[real_node*NDIM+0]*U) / r;
		  pb = pow( rhob, p.gamma ) / p.gamma;
		}
	      else
		{
		  pb = pi + rho0*c0*theta;
		  rhob = rhoi + (pb-pi)/(c0*c0);
		  ub = ui - (pb-pi)/(rho0*c0)*nx;
		  vb = vi - (pb-pi)/(rho0*c0)*ny;
		}
	    }
	  else
	    {
	      r = 3.0;
	      U = minf;
	     
	      if ( MMS )
		{
		  rhob = pow( (1.0 + gm1*2.0*(1.0 - 4.0/(r*r))) , gm1inv );
		  ub =  (grid->x[real_node*NDIM+1]*U) / r;
		  vb = -(grid->x[real_node*NDIM+0]*U) / r;
		  pb = pow( rhob, p.gamma ) / p.gamma;
		}
	      else
		{
		  pb = pi + rho0*c0*theta;
		  rhob = rhoi + (pb-pi)/(c0*c0);
		  ub = ui - (pb-pi)/(rho0*c0)*nx;
		  vb = vi - (pb-pi)/(rho0*c0)*ny;
		}
	    }

	  if ( RECON_PRIM )
	    {
	      grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +0 ] = rhob;
	      grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +1 ] = ub;
	      grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +2 ] = vb;
	      grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +3 ] = pb;
	    }
	  else
	    {
	      grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +0 ] = rhob;
	      grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +1 ] = rhob*ub;
	      grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +2 ] = rhob*vb;
	      grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +3 ] = pb/(p.gamma-1.0) + 0.5*rhob*(ub*ub + vb*vb);
	    }

	  continue;
	}

      // Now I can set the ghost node values based on the applied BC ( currently Farfield or Inviscid Surface ).
      if ( bc == 0 || (bc%2 == 0))  // Farfield Condition
	{
	  // Test the eigenvalues for inflow/outflow.
	  if ( ev1 <= 0.0 && ev4 <= 0.0 ) // Inflow
	    {
	      if ( ev3 >= 0.0 ) //Subsonic
		{
		  pb = 0.5 * ( pinf + pi - rho0*c0*( (uinf-ui)*nx + (vinf-vi)*ny ) );
		  rhob = rhoinf + (pb - pinf)/(c0*c0);
		  ub = uinf + (pb-pinf)/(rho0*c0) * nx;
		  vb = vinf + (pb-pinf)/(rho0*c0) * ny;
		}
	      else if ( ev3 < 0.0 ) // Supersonic  // Edit: changed from ev4
		{
		  rhob = rhoinf;
		  ub = uinf;
		  vb = vinf;
		  pb = pinf;
		}
	      else // Failure.
		{
		  printf("CRITICAL ERROR: In Set_Boundary_Values() the eigenvalues do not correspond to a valid inflow state.\n");
		  printf("  boundary edge = %d, node = %d, ghost node = %d, boundary = %d, segement = %d\n",i,real_node,ghost_node,b,s);
		  printf("  Q[0] = %.15E, Q[1] = %.15E, Q[2] = %.15E, Q[3] = %.15E\n",grid->nQ[real_node*NUM_VAR+0],grid->nQ[real_node*NUM_VAR+1],
			 grid->nQ[real_node*NUM_VAR+2],grid->nQ[real_node*NUM_VAR+3]);
		  printf("  nx = %.15E , ny = %.15E\n",nx,ny);
		  printf("  theta = %.15E , c = %.15E\n",theta,ci);
		  fflush(stdout);
		  exit(1);
		}
	    } // End inflow.

	  else if ( ev1 >= 0.0 && ev3 >= 0.0 ) // Outflow.
	    {
	      if ( ev4 <= 0.0 )  // Subsonic.
		{
		  pb = pinf;
		  rhob = rhoi + (pb - pi)/(c0*c0);
		  ub = ui - (pb-pi)/(rho0*c0) * nx;
		  vb = vi - (pb-pi)/(rho0*c0) * ny;
		}
	      else if ( ev4 > 0.0 )  // Supersonic.
		{
		  rhob = rhoi;
		  ub = ui;
		  vb = vi;
		  pb = pi;
		}
	      else // Failure.
		{
		  printf("CRITICAL ERROR: In Set_Boundary_Values() the eigenvalues do not correspond to a valid outflow state.\n");
		  printf("  boundary edge = %d, node = %d, ghost node = %d, boundary = %d, segement = %d\n",i,real_node,ghost_node,b,s);
		  printf("  Q[0] = %.15E, Q[1] = %.15E, Q[2] = %.15E, Q[3] = %.15E\n",grid->nQ[real_node*NUM_VAR+0],grid->nQ[real_node*NUM_VAR+1],
			 grid->nQ[real_node*NUM_VAR+2],grid->nQ[real_node*NUM_VAR+3]);
		  printf("  nx = %.15E , ny = %.15E\n",nx,ny);
		  printf("  theta = %.15E , c = %.15E\n",theta,ci);
		  fflush(stdout);
		  exit(1);
		}
	    }
	  
	  else
	    {
	      printf("CRITICAL ERROR: In Set_Boundary_Values() the eigenvalues do not correspond to a valid outflow/inflow state but the boundary is tagged farfield.\n");
	      printf("  boundary edge = %d, node = %d, ghost node = %d, boundary = %d, segement = %d\n",i,real_node,ghost_node,b,s);
	      printf("  Q[0] = %.15E, Q[1] = %.15E, Q[2] = %.15E, Q[3] = %.15E\n",grid->nQ[real_node*NUM_VAR+0],grid->nQ[real_node*NUM_VAR+1],
		     grid->nQ[real_node*NUM_VAR+2],grid->nQ[real_node*NUM_VAR+3]);
	      printf("  nx = %.15E , ny = %.15E\n",nx,ny);
	      printf("  theta = %.15E , c = %.15E\n",theta,ci);
	      fflush(stdout);
	      exit(1);
	    }

	} // End Farfield.

      /*
      else if ( bc == 2 ) // Annulus exit. We know it will be supersonic outflow.
	{
	  rhob = rhoi;
	  ub = ui;
	  vb = vi;
	  pb = pi;
	}

      else if ( bc == 3 ) // Annulus entrance. The conditions of entrance are specified.
	{
	  r = grid->x[real_node*NDIM + 1];
	  U = (2.0*2.0)/r;
	  
	  rhob = pow( (1.0 + (1.4-1.0)*2.0*(1.0 - 4.0/(r*r))) , (1.0/(1.4-1.0)) );
	  ub = U;
	  vb = 0.;
	  pb = pow( rhob, p.gamma ) / p.gamma;
	}
      */

      else if ( bc == 1 || (bc%2==1)) // Inviscid surface
	{
	  // I'm not going to look at the eigenvalues, I'm assuming here that they are zero and are set appropriately.
	  pb = pi + rho0*c0*theta;
	  rhob = rhoi + (pb-pi)/(c0*c0);
	  ub = ui - (pb-pi)/(rho0*c0)*nx;
	  vb = vi - (pb-pi)/(rho0*c0)*ny;
	}

      else
	{
	  printf("CRITICAL ERROR: In Set_Boundary_Values() an invalid boundary condition was set.\n");
	  printf("  edge = %d, boundary = %d, segment = %d, bc = %d\n",i,b,s,bc);
	  fflush(stdout);
	  exit(1);
	}

      // Convert the primitive values back to conserved variables.

      if ( RECON_PRIM )
	{
	  grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +0 ] = rhob;
	  grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +1 ] = ub;
	  grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +2 ] = vb;
	  grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +3 ] = pb;
	}
      else
	{
	  grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +0 ] = rhob;
	  grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +1 ] = rhob*ub;
	  grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +2 ] = rhob*vb;
	  grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +3 ] = pb/(p.gamma-1.0) + 0.5*rhob*(ub*ub + vb*vb);
	}

      // Here I need to check my density to make sure its not set below sero. If it is I will do a
      // zero pressure gradient inspired fix to mirror the momentum vector.

      if ( rhob < 0.0 )
	{
	  printf("Density is negative.\n");
	  fflush(stdout);
	  
	  rhob = rhoi;
	  ub = (1. - 2.*nx*nx)*ui - 2.*nx*ny*vi;   // This is a Householder reflector taken from Demmel on page 119
	  vb = (1. - 2.*ny*ny)*vi - 2.*nx*ny*ui;   // of Applied Numerical Linear Algebra.
	  pb = pi;

	  // Convert back to conserved variables.
	  
	  if ( RECON_PRIM )
	    {
	      grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +0 ] = rhob;
	      grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +1 ] = ub;
	      grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +2 ] = vb;
	      grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +3 ] = pb;
	    }
	  else
	    {
	      grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +0 ] = rhob;
	      grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +1 ] = rhob*ub;
	      grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +2 ] = rhob*vb;
	      grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +3 ] = pb/(p.gamma-1.0) + 0.5*rhob*(ub*ub + vb*vb);
	    }
	}
      
    }  // End boundary edge loop.
  
  if ( pb < 0.0 )
    {
      printf("Pressure is negative.\n");
      fflush(stdout);
    }

  return;
}




//=============================================================
// 
//  Get_Boundary_Value()
//
//  Returns the Q values for a ghost node using CVBCs.
//  The method comes from the summer course notes at the SimCenter
//  which are bases on AIAA 1984-1552 by Whitfield and Janus.
//
//  GRID *grid;                          // The grid.
//  PARAMS p;                            // The parameters.
//  int isubedge;                        // The boundary edge index.
//  double Qi[4];                        // The interior cell state.
//  double Qb[4];                        // The ghost state.
//
//=============================================================

void Get_Boundary_Value ( GRID *grid, PARAMS p, int isubedge, double Qi[4], double Qb[4] )
{
  int i,b,s;                         // Loop counters and generic entities.
  int real_node;                     // Real node on the edge.
  int ghost_node;                    // Ghost node on the edge.
  int bc;                            // Edge boundary condition.

  double uinf,vinf,rhoinf,pinf;      // Values at the freestream.
  double minf;
  double u0,v0,rho0,p0,c0;           // Values at the reference state.
  double theta,ci;                   // Nodal values.
  double rhoi,ui,vi,pi;
  double rhob,ub,vb,pb;              // Boundary values.
  double cb;
  double ev1,ev3,ev4;                // Eigenvalues.

  double nx,ny;                      // Normal vector of edge.

  double angle;                      // Angle of attack

  double x1,y1;
  double PI = M_PI;
  double gm1, gm1inv;
  double r, U;

  angle = p.alpha*M_PI/180.;

  // Get the freestream values.
  minf = p.mach;
  rhoinf = 1.0;
  uinf = p.mach*cos(angle);
  vinf = p.mach*sin(angle);
  pinf = 1. / p.gamma;

  //if ( p.ic == 2 )
  //  pinf = 1.;

  i = isubedge;

  // Grab the data from the structure.
  real_node = grid->bedges[i*5+0];
  ghost_node = grid->bedges[i*5+1];
  b = grid->bedges[i*5+3];
  s = grid->bedges[i*5+4];
  
  // Determine the BC.
  bc = grid->bbc[b];

  // Get the normal vector information.
  nx = grid->xn_bedges[i*3+0];
  ny = grid->xn_bedges[i*3+1];
  
  // Get the local primitive values and the local boundary values..

  if ( RECON_PRIM )
    {
      rhoi = Qi[0];
      ui =   Qi[1];
      vi =   Qi[2];
      pi =   Qi[3];
      ci = sqrt( (p.gamma*pi)/rhoi );

      rhob = grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +0 ];
      ub =   grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +1 ];
      vb =   grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +2 ];
      pb =   grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +3 ];
      cb = sqrt( (p.gamma*pb)/rhob );
    }
  else
    {
      rhoi = Qi[0];
      ui =   Qi[1] / rhoi;
      vi =   Qi[2] / rhoi;
      pi = (p.gamma - 1.) * ( Qi[3] - 0.5*rhoi*(ui*ui + vi*vi) );
      ci = sqrt( (p.gamma*pi)/rhoi );

      rhob = grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +0 ];
      ub =   grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +1 ] / rhob;
      vb =   grid->nQ[ (grid->nn + ghost_node)*NUM_VAR +2 ] / rhob;
      pb = (p.gamma - 1.) * ( grid->nQ[ (grid->nn + ghost_node )*NUM_VAR +3] - 0.5*rhob*(ub*ub + vb*vb) );
      cb = sqrt( (p.gamma*pb)/rhob );
    }
  
  // Get the new reference state.

  if ( ANNULUS )
    {
      //ci = pow( rhoi , (0.5*(p.gamma-1.0)) );
      //cb = pow( rhob , (0.5*(p.gamma-1.0)) );
    }
  
  rho0 = 0.5*(rhoi+rhob);
  u0 = 0.5*(ui+ub);
  v0 = 0.5*(vi+vb);
  p0 = 0.5*(pi+pb);
  c0 = 0.5*(ci+cb);

  // Get the local eigenvalues.
  theta = ui*nx + vi*ny;
  ev1 = theta;
  ev3 = theta + ci;
  ev4 = theta - ci;

  // debug the flux integral.
  if ( MMS || (MMS && ANNULUS) )
    {
      // Every boundary edge will be treated as farfield.
      x1 = grid->x[real_node*NDIM+0];
      y1 = grid->x[real_node*NDIM+1];
      
      //rhob = 1.0 + sin(PI*x1)*sin(PI*y1);
      rhob = 1.0 + 0.25*sin(PI*x1)*sin(PI*y1);
      ub = 0.25 + 0.25*sin(PI*x1)*cos(2.*PI*y1);
      vb = 0.25 + 0.25*cos(2.*PI*x1)*sin(PI*y1);
      pb = 1./p.gamma + 0.05*cos(2.*PI*x1)*cos(2.*PI*y1);
      
      if ( RECON_PRIM )
	{
	  Qb[0] = rhob;
	  Qb[1] = ub;
	  Qb[2] = vb;
	  Qb[3] = pb;
	}
      else
	{
	  Qb[0] = rhob;
	  Qb[1] = rhob*ub;
	  Qb[2] = rhob*vb;
	  Qb[3] = pb/(p.gamma-1.0) + 0.5*rhob*(ub*ub + vb*vb);
	}

      return;
    }

  if ( MMS_EXP )
    {
      // Every boundary edge will be treated as farfield.
      x1 = grid->x[real_node*NDIM+0];
      y1 = grid->x[real_node*NDIM+1];
      
      rhob = exp(x1-1.0) * exp(y1 - 1.0);
      ub   = exp(x1-1.0) * exp(y1 - 1.0);
      vb   = exp(x1-1.0) * exp(y1 - 1.0);
      pb   = ( (1.4-1.0) * 0.5 ) * exp( 3.0 * ( x1 - 1.0) ) * exp( 3.0 * ( y1 - 1.0 ) );
      
      if ( RECON_PRIM )
	{
	  Qb[0] = rhob;
	  Qb[1] = ub;
	  Qb[2] = vb;
	  Qb[3] = pb;
	}
      else
	{
	  Qb[0] = rhob;
	  Qb[1] = rhob*ub;
	  Qb[2] = rhob*vb;
	  Qb[3] = pb/(p.gamma-1.0) + 0.5*rhob*(ub*ub + vb*vb);
	}
      
      return;
    }

  // CODE FOR THE ANNULUS - supersonic vortex problem from Aftosmis.
  if ( ANNULUS )
    {
      gm1 = p.gamma - 1.0;
      gm1inv = 1.0 / gm1;

      minf = 2.0;
      rhoinf = 1.0;
	  
      if ( bc == 0 )
	{
	  r = grid->x[real_node*NDIM + 1];
	  U = (minf*2.)/r;
	  
	  rhob = pow( (1.0 + gm1*2.0*(1.0 - 4.0/(r*r))) , gm1inv );
	  ub = U;
	  vb = 0.;
	  pb = pow( rhob, p.gamma ) / p.gamma;
	}
      else if ( bc == 2 )
	{
	  r = grid->x[real_node*NDIM + 0];
	  U = (minf*2.)/r;
	  
	  if ( MMS )
	    {
	      rhob = pow( (1.0 + gm1*2.0*(1.0 - 4.0/(r*r))) , gm1inv );
	      ub = 0.;
	      vb = -U;
	      pb = pow( rhob, p.gamma ) / p.gamma;
	    }
	  else
	    {
	      rhob = rhoi;
	      ub = ui;
	      vb = vi;
	      pb = pi;
	    }
	}
      else if ( (bc == 13 || bc == 14) || bc == 1 ) // inner radius
	{
	  r = 2.0;
	  U = minf;
	  
	  //if ( MMS || 0 )
	  //if ( grid->citer > 6000 )
	  if ( MMS )
	    {
	      rhob = pow( (1.0 + gm1*2.0*(1.0 - 4.0/(r*r))) , gm1inv );
	      ub =  (grid->x[real_node*NDIM+1]*U) / r;
	      vb = -(grid->x[real_node*NDIM+0]*U) / r;
	      pb = pow( rhob, p.gamma ) / p.gamma;
	    }
	  else
	    {
	      pb = pi + rho0*c0*theta;
	      rhob = rhoi + (pb-pi)/(c0*c0);
	      ub = ui - (pb-pi)/(rho0*c0)*nx;
	      vb = vi - (pb-pi)/(rho0*c0)*ny;
	    }
	}
      else
	{
	  r = 3.0;
	  U = minf;
	  
	  if ( MMS )
	    {
	      rhob = pow( (1.0 + gm1*2.0*(1.0 - 4.0/(r*r))) , gm1inv );
	      ub =  (grid->x[real_node*NDIM+1]*U) / r;
	      vb = -(grid->x[real_node*NDIM+0]*U) / r;
	      pb = pow( rhob, p.gamma ) / p.gamma;
	    }
	  else
	    {
	      pb = pi + rho0*c0*theta;
	      rhob = rhoi + (pb-pi)/(c0*c0);
	      ub = ui - (pb-pi)/(rho0*c0)*nx;
	      vb = vi - (pb-pi)/(rho0*c0)*ny;
	    }
	}

      if ( RECON_PRIM )
	{
	  Qb[0] = rhob;
	  Qb[1] = ub;
	  Qb[2] = vb;
	  Qb[3] = pb;
	}
      else
	{
	  Qb[0] = rhob;
	  Qb[1] = rhob*ub;
	  Qb[2] = rhob*vb;
	  Qb[3] = pb/(p.gamma-1.0) + 0.5*rhob*(ub*ub + vb*vb);
	}
      
      return;
    }
  
  // Now I can set the ghost node values based on the applied BC ( currently Farfield or Inviscid Surface ).
  if ( bc == 0 || (bc%2==0) )  // Farfield Condition
    {
      // Test the eigenvalues for inflow/outflow.
      if ( ev1 <= 0.0 && ev4 <= 0.0 ) // Inflow
	{
	  if ( ev3 >= 0.0 ) //Subsonic
	    {
	      pb = 0.5 * ( pinf + pi - rho0*c0*( (uinf-ui)*nx + (vinf-vi)*ny ) );
	      rhob = rhoinf + (pb - pinf)/(c0*c0);
	      ub = uinf + (pb-pinf)/(rho0*c0) * nx;
	      vb = vinf + (pb-pinf)/(rho0*c0) * ny;
	    }
	  else if ( ev3 < 0.0 ) // Supersonic  // Edit: Changed from ev4
	    {
	      rhob = rhoinf;
	      ub = uinf;
	      vb = vinf;
	      pb = pinf;
	    }
	  else // Failure.
	    {
	      printf("CRITICAL ERROR: In Set_Boundary_Values() the eigenvalues do not correspond to a valid inflow state.\n");
	      printf("  boundary edge = %d, node = %d, ghost node = %d, boundary = %d, segement = %d\n",i,real_node,ghost_node,b,s);
	      printf("  Q[0] = %.15E, Q[1] = %.15E, Q[2] = %.15E, Q[3] = %.15E\n",grid->nQ[real_node*NUM_VAR+0],grid->nQ[real_node*NUM_VAR+1],
		     grid->nQ[real_node*NUM_VAR+2],grid->nQ[real_node*NUM_VAR+3]);
	      printf("  nx = %.15E , ny = %.15E\n",nx,ny);
	      printf("  theta = %.15E , c = %.15E\n",theta,ci);
	      fflush(stdout);
	      exit(1);
	    }
	} // End inflow.
      
      else if ( ev1 >= 0.0 && ev3 >= 0.0 ) // Outflow.
	{
	  if ( ev4 <= 0.0 )  // Subsonic.
	    {
	      pb = pinf;
	      rhob = rhoi + (pb - pi)/(c0*c0);
	      ub = ui - (pb-pi)/(rho0*c0) * nx;
	      vb = vi - (pb-pi)/(rho0*c0) * ny;
	    }
	  else if ( ev4 > 0.0 )  // Supersonic.
	    {
	      rhob = rhoi;
	      ub = ui;
	      vb = vi;
	      pb = pi;
	    }
	  else // Failure.
	    {
	      printf("CRITICAL ERROR: In Get_Boundary_Values() the eigenvalues do not correspond to a valid outflow state.\n");
	      printf("  boundary edge = %d, node = %d, ghost node = %d, boundary = %d, segement = %d\n",i,real_node,ghost_node,b,s);
	      printf("  PERTURBED VALUES:\n");
	      printf("  Q[0] = %.15E, Q[1] = %.15E, Q[2] = %.15E, Q[3] = %.15E\n",Qi[0],Qi[1],Qi[2],Qi[3]);
	      printf("  nx = %.15E , ny = %.15E\n",nx,ny);
	      printf("  theta = %.15E , c = %.15E\n",theta,ci);
	      fflush(stdout);
	      exit(1);
	    }
	}
      
      else
	{
	  printf("CRITICAL ERROR: In Get_Boundary_Values() the eigenvalues do not correspond to a valid outflow/inflow state but the boundary is tagged farfield.\n");
	  printf("  boundary edge = %d, node = %d, ghost node = %d, boundary = %d, segement = %d\n",i,real_node,ghost_node,b,s);
	  printf("  PERTURBED VALUES:\n");
	  printf("  Q[0] = %.15E, Q[1] = %.15E, Q[2] = %.15E, Q[3] = %.15E\n",Qi[0],Qi[1],Qi[2],Qi[3]);
	  printf("  nx = %.15E , ny = %.15E\n",nx,ny);
	  printf("  theta = %.15E , c = %.15E\n",theta,ci);
	  fflush(stdout);
	  exit(1);
	}

    } // End Farfield.

  /*
  else if ( bc == 2 ) // Annulus exit. We know it will be supersonic outflow.
    {
      rhob = rhoi;
      ub = ui;
      vb = vi;
      pb = pi;
    }

  else if ( bc == 3 ) // Annulus entrance. The conditions of entrance are specified.
    {
      r = grid->x[real_node*NDIM + 1];
      U = (2.0*2.0)/r;
      
      rhob = pow( (1.0 + (1.4-1.0)*2.0*(1.0 - 4.0/(r*r))) , (1.0/(1.4-1.0)) );
      ub = U;
      vb = 0.;
      pb = pow( rhob, p.gamma ) / p.gamma;
    }
  */

  else if ( bc == 1 || (bc%2==1)) // Inviscid surface
    {
      // I'm not going to look at the eigenvalues, I'm assuming here that they are zero and are set appropriately.
      pb = pi + rho0*c0*theta;
      rhob = rhoi + (pb-pi)/(c0*c0);
      ub = ui - (pb-pi)/(rho0*c0)*nx;
      vb = vi - (pb-pi)/(rho0*c0)*ny;
    }

  else
    {
      printf("CRITICAL ERROR: In Get_Boundary_Values() an invalid boundary condition was set.\n");
      printf("  edge = %d, boundary = %d, segment = %d, bc = %d\n",i,b,s,bc);
      fflush(stdout);
      exit(1);
    }
  
  // Convert the primitive values back to conserved variables.

  if ( RECON_PRIM )
    {
      Qb[0] = rhob;
      Qb[1] = ub;
      Qb[2] = vb;
      Qb[3] = pb;
    }
  else
    {
      Qb[0] = rhob;
      Qb[1] = rhob*ub;
      Qb[2] = rhob*vb;
      Qb[3] = pb/(p.gamma-1.0) + 0.5*rhob*(ub*ub + vb*vb);
    }

  if ( rhob < 0.0 )
    {
      printf("Density is negative in get.\n");
      fflush(stdout);
      
      rhob = rhoi;
      ub = (1. - 2.*nx*nx)*ui - 2.*nx*ny*vi;   // This is a Householder reflector taken from Demmel on page 119
      vb = (1. - 2.*ny*ny)*vi - 2.*nx*ny*ui;   // of Applied Numerical Linear Algebra.
      pb = pi;
      
      // Convert back to conserved variables.

      if ( RECON_PRIM )
	{
	  Qb[0] = rhob;
	  Qb[1] = ub;
	  Qb[2] = vb;
	  Qb[3] = pb;
	}
      else
	{
	  Qb[0] = rhob;
	  Qb[1] = rhob*ub;
	  Qb[2] = rhob*vb;
	  Qb[3] = pb/(p.gamma-1.0) + 0.5*rhob*(ub*ub + vb*vb);
	}
    }

  if ( pb < 0.0 )
    {
      printf("Pressure is negative in get.\n");
      fflush(stdout);
    }
  
  return;
}
