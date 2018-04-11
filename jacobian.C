//=============================================================
// 
//  jacobian.C
//  
//  Functions to compute the Jacobian.
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
#include "jacobian.h"
#include "flux.h"
#include "mesh_utils.h"
#include "cvbc.h"
#include "linear_algebra.h"

//=============================================================
// 
//  compute_numerical_jacobian_roe()
//
//  Computes the flux jacobian contributions over a subedge and
//  accumulates the result to the matrix.
//
//  GRID *grid;                          // The grid.
//  SYS_MEM *smem;                       // The memory for the system.
//  PARAMS p;                            // The parameters.
//  int nodeL;                           // Left node.
//  int nodeR;                           // Right node.
//  double QL[4];                        // The left state.
//  double QR[4];                        // The right state.
//  double flux[4];                      // The unperturbed flux.
//  double nx;                           // Normal vector components.
//  double ny;
//  double len;
//
//=============================================================

void compute_numerical_jacobian_roe ( GRID *grid, SYS_MEM *smem, PARAMS p, int nodeL, int nodeR,
				      double QL[4], double QR[4], double flux[4],
				      double nx, double ny, double len )
{
  int i, j, k, b;                              // Loop counters.
  int error = 0;
  double Qp[4];                                // The perturbed state vector.
  double fluxp[4];                             // The perturbed flux.
  double dfdqL[4][4], dfdqR[4][4];             // Numerical jacobians, temporary storage.
  double eps = 1.0e-08;                        // The perturbation.

  // Pointers.
  double gamma = p.gamma;
  double *LHS = smem->LHS;
  int *ia = smem->ia;
  int *iau = smem->iau;
  int *ja = smem->ja;

  // Compute the first jacobian, with respect to nodeL.
  for ( i=0; i < NUM_VAR; i++ )
    {
      // Set the perturbed state vector.
      for ( j=0; j < NUM_VAR; j++ )
	{
	  Qp[j] = QL[j];
	}
      
      Qp[i] += eps;

      // Get the numerical flux based on the perturbed state.
      //error = Roe_flux ( nx, ny, len, gamma, Qp, QR, fluxp );
      error = Roe_flux_centered ( nx, ny, len, gamma, Qp, QR, fluxp );

      // Store the result in the ith column of dfdqL.
      for ( j=0; j < NUM_VAR; j++ )
	{
	  dfdqL[j][i] = ( fluxp[j] - flux[j] ) / eps;
	}
    }

  
  // Compute the second jacobian, with respect to nodeR.
  for ( i=0; i < NUM_VAR; i++ )
    {
      // Set the perturbed state vector.
      for ( j=0; j < NUM_VAR; j++ )
	{
	  Qp[j] = QR[j];
	}
      
      Qp[i] += eps;

      // Get the numerical flux based on the perturbed state.
      //error = Roe_flux ( nx, ny, len, gamma, QL, Qp, fluxp );
      error = Roe_flux_centered ( nx, ny, len, gamma, QL, Qp, fluxp );

      // Store the result in the ith column of dfdqL.
      for ( j=0; j < NUM_VAR; j++ )
	{
	  dfdqR[j][i] = ( fluxp[j] - flux[j] ) / eps;
	}
    }

  // NaN check.
  
  for ( i=0; i < NUM_VAR; i++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  if ( !(isfinite(dfdqL[i][j])) )
	    {
	      printf("Problem DETECTED in dfdqL[%d][%d]=%f  between nodeL=%d and nodeR=%d.\n",i,j,(float)dfdqL[i][j],nodeL,nodeR);
	    }
	}
    }
  
  for ( i=0; i < NUM_VAR; i++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  if ( !(isfinite(dfdqR[i][j])) )
	    {
	      printf("Problem DETECTED in dfdqR[%d][%d]=%f  between nodeL=%d and nodeR=%d.\n",i,j,(float)dfdqR[i][j],nodeL,nodeR);
	    }
	}
    }

  // Now I can store the results in the matrix.
  
  // Now add dfdqL to the diagonal of the left node.
  for ( i=0; i < NUM_VAR; i++ )
    {
      for ( k=0; k < NUM_VAR; k++ )
	{
	  LHS[ (iau[nodeL])*NUM_VAR*NUM_VAR + i*NUM_VAR + k] += dfdqL[i][k];
	}
    }
	  
  // Now subtract dfdqR from the diagonal of the right node since the normal points torward it.
  for ( i=0; i < NUM_VAR; i++ )
    {
      for ( k=0; k < NUM_VAR; k++ )
	{
	  LHS[ (iau[nodeR])*NUM_VAR*NUM_VAR + i*NUM_VAR + k] -= dfdqR[i][k];
	}
    }
	  
  // Now find the position of the right node in the section of 'ja' associated with the left node
  // and then add dfdqR to that position in LHS.
  for ( i=ia[nodeL]; i < ia[nodeL+1]; i++ )
    {
      if ( ja[i] == nodeR )
	{
	  for ( k=0; k < NUM_VAR; k++ )
	    {
	      for ( b=0; b < NUM_VAR; b++ )
		{
		  LHS[i*NUM_VAR*NUM_VAR+k*NUM_VAR+b] += dfdqR[k][b];
		}
	    }
	}
    }
	  
  // Do the same for the left node in 'ja' for the right node but now subtract dfdqL
  // from that position in LHS.
  for ( i=ia[nodeR]; i < ia[nodeR+1]; i++ )
    {
      if ( ja[i] == nodeL )
	{
	  for ( k=0; k < NUM_VAR; k++ )
	    {
	      for ( b=0; b < NUM_VAR; b++ )
		{
		  LHS[i*NUM_VAR*NUM_VAR+k*NUM_VAR+b] -= dfdqL[k][b];
		}
	    }
	}
    }
  
  return;
}


//=============================================================
// 
//  compute_boundary_numerical_jacobian_roe()
//
//  Computes the flux jacobian contributions over a boundary edge
//  and accumulate the result to the matrix.
//
//  GRID *grid;                          // The grid.
//  SYS_MEM *smem;                       // The memory for the system.
//  PARAMS p;                            // The parameters.
//  int ibedge;                          // The boundary edge.
//  double QL[4];                        // The left state.
//  double flux[4];                      // The unperturbed flux.
//  double nx;                           // Normal vector components.
//  double ny;
//  double len;
//
//=============================================================

void compute_boundary_numerical_jacobian_roe ( GRID *grid, SYS_MEM *smem, PARAMS p, int ibedge,
					       double QL[4], double flux[4], double nx, double ny, double len )
{
  int i, j, k;                                 // Loop counters.
  int n;                                       // Node counter.
  int error = 0;
  double Qp[4];                                // The perturbed state vector.
  double Qb[4];                                // The boundary state vector.
  double fluxp[4];                             // The perturbed flux.
  double dfdqL[4][4];                          // Numerical jacobian, temporary storage.
  double eps = 1.0e-08;                        // The perturbation.

  // Pointers.
  double gamma = p.gamma;
  double *LHS = smem->LHS;
  int *iau = smem->iau;

  n = grid->bedges[ibedge*5+0];    // Get the real node index.

  // Compute the jacobian.
  for ( i=0; i < NUM_VAR; i++ )
    {
      // Set the perturbed state vector.
      for ( j=0; j < NUM_VAR; j++ )
	{
	  Qp[j] = QL[j];
	}
      
      Qp[i] += eps;

      // Update the boundary values to reflect the perturbation.
      Get_Boundary_Value( grid, p, ibedge, Qp, Qb );
      //Get_Boundary_Value( grid, p, ibedge, QL, Qb ); // Freeze the boundary state.

      // Get the numerical flux based on the perturbed state.
      //error = Roe_flux ( nx, ny, len, gamma, Qp, Qb, fluxp );
      error = Roe_flux_centered ( nx, ny, len, gamma, Qp, Qb, fluxp );

      // Store the result in the ith column of dfdqL.
      for ( j=0; j < NUM_VAR; j++ )
	{
	  dfdqL[j][i] = ( fluxp[j] - flux[j] ) / eps;
	}
    }

  // NaN check.
  
  for ( i=0; i < NUM_VAR; i++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  if ( !(isfinite(dfdqL[i][j])) )
	    {
	      printf("Problem DETECTED in dfdqL[%d][%d]=%f on bedge %d, real node = %d.\n",i,j,(float)dfdqL[i][j],ibedge,n);
	    }
	}
    }

  // Now I can store the result in the matrix.
  
  // Now add dfdqL to the diagonal of the physical boundary node.
  for ( i=0; i < NUM_VAR; i++ )
    {
      for ( k=0; k < NUM_VAR; k++ )
	{
	  LHS[ (iau[n])*NUM_VAR*NUM_VAR + i*NUM_VAR + k] += dfdqL[i][k];
	}
    }
  
  return;
}


//=============================================================
// 
//  compute_approximate_jacobian_roe()
//
//  Computes the flux jacobian contributions over a subedge and
//  accumulates the result to the matrix.
//
//  GRID *grid;                          // The grid.
//  SYS_MEM *smem;                       // The memory for the system.
//  PARAMS p;                            // The parameters.
//  int nodeL;                           // Left node.
//  int nodeR;                           // Right node.
//  double QL[4];                        // The left state.
//  double QR[4];                        // The right state.
//  double flux[4];                      // The unperturbed flux.
//  double nx;                           // Normal vector components.
//  double ny;
//  double len;
//
//=============================================================

void compute_approximate_jacobian_roe ( GRID *grid, SYS_MEM *smem, PARAMS p, int nodeL, int nodeR,
					double QL[4], double QR[4], double flux[4],
					double nx, double ny, double len )
{
  int i, j, k, b;                              // Loop counters.
  double dfdqL[4][4], dfdqR[4][4];             // Analytic Euler jacobians, temporary storage.
  double theta,phi;                            // Efficiency variables.

  // Pointers.
  double gamma = p.gamma;
  double *LHS = smem->LHS;
  int *ia = smem->ia;
  int *iau = smem->iau;
  int *ja = smem->ja;

  // Flux variables.
  double rhobar;                     // Roe averaged density.
  double ubar;                       // Roe aveeraged velocities.
  double vbar;
  double Hbar;                       // Roe averaged enthalpy.
  double cbar;                       // Roe averaged speed of sound.
  double pbar;                       // Roe averaged pressure.
  double thetabar;                   // Roe averaged theta.
  double V2bar;                      // Roe averaged V^2.
  double rhoL,uL,vL,HL,pL;           // Left state.
  double rhoR,uR,vR,HR,pR;           // Right state.
  double ER,EL;
  double big_gamma;                  // Multiplication factor to reduce divisions. See Daniel's notes.
  double hL,hR,hbar;
  double ev1,ev2,ev3,ev4;            // Eigenvalues.
  
  double Lambda[NUM_VAR];            // Eigenvalue vector (diagonal matrix).
  double T[NUM_VAR*NUM_VAR];         // Eigenvector matrices.
  double Tinv[NUM_VAR*NUM_VAR];
  double temp[NUM_VAR*NUM_VAR];      // Temporary matrix.
  double Roemat[NUM_VAR*NUM_VAR];    // The Roe matrix.

  // Efficiency variables.
  double oo_2cbar;                   // One over 2 times cbar.
  double oo_2gm1;                    // One over 2 times (gamma - 1).
  double oo_rhobar;                  // One over rhobar;
  double oo_rcbar;                   // One over rhobar times cbar.
  double oo_cbar2;                   // One over cbar squared.

  // Get the right and left states.
  if ( RECON_PRIM )
    {
      rhoL = QL[0];
      uL = QL[1];
      vL = QL[2];
      pL = QL[3];
      EL = pL/(gamma-1.) + 0.5*rhoL*(uL*uL + vL*vL);
      HL = EL + pL;
      hL = HL / rhoL;

      rhoR = QR[0];
      uR = QR[1];
      vR = QR[2];
      pR = QR[3];
      ER = pR/(gamma-1.) + 0.5*rhoR*(uR*uR + vR*vR);
      HR = ER + pR;
      hR = HR / rhoR;
    }
  else
    {
      rhoL = QL[0];
      uL = QL[1]/QL[0];
      vL = QL[2]/QL[0];
      EL = QL[3];
      pL = (gamma - 1.)*(EL - 0.5*rhoL*(uL*uL + vL*vL) );
      HL = EL + pL;
      hL = HL / rhoL;

      rhoR = QR[0];
      uR = QR[1]/QR[0];
      vR = QR[2]/QR[0];
      ER = QR[3];
      pR = (gamma - 1.)*(ER - 0.5*rhoR*(uR*uR + vR*vR) );
      HR = ER + pR;
      hR = HR / rhoR;
    }

  // Generate the Roe averaged values.
  rhobar = sqrt(rhoL*rhoR);
  big_gamma = rhobar / ( rhoL + rhobar ) ;

  ubar = uL + big_gamma*(uR - uL);
  vbar = vL + big_gamma*(vR - vL);
  Hbar = HL + big_gamma*(HR - HL);
  hbar = hL + big_gamma*(hR - hL);

  thetabar = ubar*nx + vbar*ny;

  V2bar = (ubar*ubar + vbar*vbar);

  //pbar = (gamma-1.)/gamma * ( Hbar - rhobar*V2bar*0.5 );
  //cbar = sqrt( (gamma*pbar)/rhobar );

  //cbar = sqrt( (gamma-1.0)*(Hbar - 0.5*V2bar) );
  cbar = sqrt( (gamma-1.0)*(hbar - 0.5*V2bar) );

  // Set the efficiency values.
  oo_2cbar = 1. / ( 2.*cbar );
  oo_2gm1 = 1. / ( 2.*(gamma - 1.) );
  oo_rhobar = 1. / rhobar;
  oo_rcbar = 1. / ( rhobar * cbar );
  oo_cbar2 = 1. / ( cbar * cbar );

  // Set the eigenvalues.
  ev1 = thetabar;
  ev2 = thetabar;
  ev3 = thetabar + cbar;
  ev4 = thetabar - cbar;

  // Get the negative eigenvalues since this is the Roe scheme I'm going to use.
  Lambda[0] = fabs(ev1);
  Lambda[1] = fabs(ev2);
  Lambda[2] = fabs(ev3);
  Lambda[3] = fabs(ev4);

  // Now I can fill in the eigenvector matrices. I have the derivation of them in my notes with help from Hirsch's book (volume 2).
  // Row major format.
  /*
  T[0] = 1.0;
  T[1] = 0.;
  T[2] = rhobar*oo_2cbar;
  T[3] = rhobar*oo_2cbar;
  
  T[4] = ubar;
  T[5] = rhobar*ny;
  T[6] = rhobar*ubar*oo_2cbar + 0.5*rhobar*nx;
  T[7] = rhobar*ubar*oo_2cbar - 0.5*rhobar*nx;
  
  T[8] =  vbar;
  T[9] =  -rhobar*nx;
  T[10] = rhobar*vbar*oo_2cbar + 0.5*rhobar*ny;
  T[11] = rhobar*vbar*oo_2cbar - 0.5*rhobar*ny;

  T[12] = 0.5*V2bar;
  T[13] = rhobar*( ubar*ny - vbar*nx );
  T[14] = 0.5*rhobar*oo_2cbar*V2bar + rhobar*cbar*oo_2gm1 + 0.5*rhobar*thetabar;
  T[15] = 0.5*rhobar*oo_2cbar*V2bar + rhobar*cbar*oo_2gm1 - 0.5*rhobar*thetabar;

  //============================================================================

  Tinv[0] = 1. - 0.5*oo_cbar2*(gamma-1.)*V2bar;
  Tinv[1] = (gamma-1.)*ubar*oo_cbar2;
  Tinv[2] = (gamma-1.)*vbar*oo_cbar2;
  Tinv[3] = -(gamma-1.)*oo_cbar2;

  Tinv[4] = oo_rhobar*( -ubar*ny + vbar*nx );
  Tinv[5] = ny*oo_rhobar;
  Tinv[6] = -nx*oo_rhobar;
  Tinv[7] = 0.;

  Tinv[8]  = -thetabar*oo_rhobar + 0.5*oo_rcbar*(gamma-1.)*V2bar;
  Tinv[9]  = nx*oo_rhobar - (gamma-1.)*ubar*oo_rcbar;
  Tinv[10] = ny*oo_rhobar - (gamma-1.)*vbar*oo_rcbar;
  Tinv[11] = (gamma-1.)*oo_rcbar;

  Tinv[12] = thetabar*oo_rhobar + 0.5*oo_rcbar*(gamma-1.)*V2bar;
  Tinv[13] = -nx*oo_rhobar - (gamma-1.)*ubar*oo_rcbar;
  Tinv[14] = -ny*oo_rhobar - (gamma-1.)*vbar*oo_rcbar;
  Tinv[15] = (gamma-1.)*oo_rcbar;
  */


  // Here I am using the matrices given in Nejat's dissertation. L is T, R is T inverse.
  T[0] = 1.0;
  T[1] = 0.0;
  T[2] = 1.0 / (cbar*cbar);
  T[3] = 1.0 / (cbar*cbar);
  
  T[4] = ubar;
  T[5] = -ny;
  T[6] = ubar / (cbar*cbar) + nx / cbar;
  T[7] = ubar / (cbar*cbar) - nx / cbar;
  
  T[8] = vbar;
  T[9] = nx;
  T[10] = vbar / (cbar*cbar) + ny / cbar;
  T[11] = vbar / (cbar*cbar) - ny / cbar;

  T[12] = 0.5*V2bar;
  T[13] = ubar*(-ny) + vbar*nx;
  T[14] = 1.0/(gamma-1.0) + thetabar/cbar + (0.5*V2bar)/(cbar*cbar);
  T[15] = 1.0/(gamma-1.0) - thetabar/cbar + (0.5*V2bar)/(cbar*cbar);

  //============================================================================

  Tinv[0] = 1.0 - (gamma-1.0)*(0.5*V2bar)/(cbar*cbar);
  Tinv[1] = (gamma-1.0)*ubar/(cbar*cbar);
  Tinv[2] = (gamma-1.0)*vbar/(cbar*cbar);
  Tinv[3] = -1.0*(gamma-1.0)/(cbar*cbar);

  Tinv[4] = -1.0 * ( ubar*(-ny) + vbar*nx );
  Tinv[5] = -ny;
  Tinv[6] = nx;
  Tinv[7] = 0.;

  Tinv[8]  = 0.5 * ( -cbar * thetabar + (gamma-1.0)*0.5*V2bar );
  Tinv[9]  = 0.5 * ( cbar * nx - (gamma-1.0)*ubar );
  Tinv[10] = 0.5 * ( cbar * ny - (gamma-1.0)*vbar );
  Tinv[11] = 0.5 * ( gamma - 1.0 );

  Tinv[12] = 0.5 * ( cbar * thetabar + (gamma-1.0)*0.5*V2bar );
  Tinv[13] = -0.5 * ( cbar * nx + (gamma-1.0)*ubar );
  Tinv[14] = -0.5 * ( cbar * ny + (gamma-1.0)*vbar );
  Tinv[15] = 0.5 * ( gamma - 1.0 );

  
  // Now multiply Lambda times Tinv.
  MMMult_DtG ( Lambda, Tinv, temp, NUM_VAR );

  // And T times the previous result.
  MMMult_GtG ( T, temp, Roemat, NUM_VAR );  // Roe matrix is now filled in.

  // Now we can build the Euler Flux Jacobians from the appropriate states.

  // Left state.
  
  theta = uL*nx + vL*ny;
  phi = 0.5 * ( gamma - 1. ) * ( uL*uL + vL*vL );
  
  dfdqL[0][0] = 0.;
  dfdqL[0][1] = nx;
  dfdqL[0][2] = ny;
  dfdqL[0][3] = 0.;

  dfdqL[1][0] = phi*nx - uL*theta;
  dfdqL[1][1] = (2.-gamma)*uL*nx + theta;
  dfdqL[1][2] = uL*ny - (gamma-1.)*vL*nx;
  dfdqL[1][3] = (gamma-1.)*nx;
  
  dfdqL[2][0] = phi*ny - vL*theta;
  dfdqL[2][1] = vL*nx - (gamma-1.)*uL*ny;
  dfdqL[2][2] = theta + (2.-gamma)*vL*ny;
  dfdqL[2][3] = (gamma-1.)*ny;
  
  dfdqL[3][0] = (2.*phi - gamma*EL/rhoL ) * theta;
  dfdqL[3][1] = (1.-gamma)*uL*theta + nx*(gamma*EL/rhoL - phi);
  dfdqL[3][2] = (1.-gamma)*vL*theta + ny*(gamma*EL/rhoL - phi);
  dfdqL[3][3] = gamma*theta;

  // Right state.
  
  theta = uR*nx + vR*ny;
  phi = 0.5 * ( gamma - 1. ) * ( uR*uR + vR*vR );
  
  dfdqR[0][0] = 0.;
  dfdqR[0][1] = nx;
  dfdqR[0][2] = ny;
  dfdqR[0][3] = 0.;

  dfdqR[1][0] = phi*nx - uR*theta;
  dfdqR[1][1] = (2.-gamma)*uR*nx + theta;
  dfdqR[1][2] = uR*ny - (gamma-1.)*vR*nx;
  dfdqR[1][3] = (gamma-1.)*nx;
  
  dfdqR[2][0] = phi*ny - vR*theta;
  dfdqR[2][1] = vR*nx - (gamma-1.)*uR*ny;
  dfdqR[2][2] = theta + (2.-gamma)*vR*ny;
  dfdqR[2][3] = (gamma-1.)*ny;
  
  dfdqR[3][0] = (2.*phi - gamma*ER/rhoR ) * theta;
  dfdqR[3][1] = (1.-gamma)*uR*theta + nx*(gamma*ER/rhoR - phi);
  dfdqR[3][2] = (1.-gamma)*vR*theta + ny*(gamma*ER/rhoR - phi);
  dfdqR[3][3] = gamma*theta;

  // Loop through and complete the computation before accumulation.
  for ( i=0; i < NUM_VAR; i++ )
    {
      for ( k=0; k < NUM_VAR; k++ )
	{
	  dfdqL[i][k] = 0.5 * ( dfdqL[i][k] + Roemat[i*NUM_VAR + k] ) * len;
	  dfdqR[i][k] = 0.5 * ( dfdqR[i][k] - Roemat[i*NUM_VAR + k] ) * len;
	}
    }


  // NaN check.
  
  for ( i=0; i < NUM_VAR; i++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  if ( !(isfinite(dfdqL[i][j])) )
	    {
	      printf("Problem DETECTED in dfdqL[%d][%d]=%f  between nodeL=%d and nodeR=%d.\n",i,j,(float)dfdqL[i][j],nodeL,nodeR);
	    }
	}
    }
  
  for ( i=0; i < NUM_VAR; i++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  if ( !(isfinite(dfdqR[i][j])) )
	    {
	      printf("Problem DETECTED in dfdqR[%d][%d]=%f  between nodeL=%d and nodeR=%d.\n",i,j,(float)dfdqR[i][j],nodeL,nodeR);
	    }
	}
    }
  
  // Now I can store the results in the matrix.
  
  // Now add dfdqL to the diagonal of the left node.
  for ( i=0; i < NUM_VAR; i++ )
    {
      for ( k=0; k < NUM_VAR; k++ )
	{
	  LHS[ (iau[nodeL])*NUM_VAR*NUM_VAR + i*NUM_VAR + k] += dfdqL[i][k];
	}
    }
	  
  // Now subtract dfdqR from the diagonal of the right node since the normal points torward it.
  for ( i=0; i < NUM_VAR; i++ )
    {
      for ( k=0; k < NUM_VAR; k++ )
	{
	  LHS[ (iau[nodeR])*NUM_VAR*NUM_VAR + i*NUM_VAR + k] -= dfdqR[i][k];
	}
    }
	  
  // Now find the position of the right node in the section of 'ja' associated with the left node
  // and then add dfdqR to that position in LHS.
  for ( i=ia[nodeL]; i < ia[nodeL+1]; i++ )
    {
      if ( ja[i] == nodeR )
	{
	  for ( k=0; k < NUM_VAR; k++ )
	    {
	      for ( b=0; b < NUM_VAR; b++ )
		{
		  LHS[i*NUM_VAR*NUM_VAR+k*NUM_VAR+b] += dfdqR[k][b];
		}
	    }
	}
    }
	  
  // Do the same for the left node in 'ja' for the right node but now subtract dfdqL
  // from that position in LHS.
  for ( i=ia[nodeR]; i < ia[nodeR+1]; i++ )
    {
      if ( ja[i] == nodeL )
	{
	  for ( k=0; k < NUM_VAR; k++ )
	    {
	      for ( b=0; b < NUM_VAR; b++ )
		{
		  LHS[i*NUM_VAR*NUM_VAR+k*NUM_VAR+b] -= dfdqL[k][b];
		}
	    }
	}
    }
  
  return;
}


//=============================================================
// 
//  compute_boundary_approximate_jacobian_roe()
//
//  Computes the flux jacobian contributions over a boundary edge
//  and accumulate the result to the matrix.
//
//  GRID *grid;                          // The grid.
//  SYS_MEM *smem;                       // The memory for the system.
//  PARAMS p;                            // The parameters.
//  int ibedge;                          // The boundary edge.
//  double QL[4];                        // The left state.
//  double QR[4];                        // The ghost state. (Outside)
//  double flux[4];                      // The unperturbed flux.
//  double nx;                           // Normal vector components.
//  double ny;
//  double len;
//
//=============================================================

void compute_boundary_approximate_jacobian_roe ( GRID *grid, SYS_MEM *smem, PARAMS p, int ibedge, double QL[4],
						 double QR[4], double flux[4], double nx, double ny, double len )
{
  int i, j, k;                                 // Loop counters.
  int n;                                       // Node counter.
  double dfdqL[4][4];                          // Analytic Euler jacobian, temporary storage.
  double theta,phi;                            // Efficiency variables.

  // Pointers.
  double gamma = p.gamma;
  double *LHS = smem->LHS;
  int *iau = smem->iau;

  // Flux variables.
  double rhobar;                     // Roe averaged density.
  double ubar;                       // Roe aveeraged velocities.
  double vbar;
  double Hbar;                       // Roe averaged enthalpy.
  double cbar;                       // Roe averaged speed of sound.
  double pbar;                       // Roe averaged pressure.
  double thetabar;                   // Roe averaged theta.
  double V2bar;                      // Roe averaged V^2.
  double rhoL,uL,vL,HL,pL;           // Left state.
  double rhoR,uR,vR,HR,pR;           // Right state.
  double ER,EL;
  double big_gamma;                  // Multiplication factor to reduce divisions. See Daniel's notes.
  double hL,hR,hbar;
  double ev1,ev2,ev3,ev4;            // Eigenvalues.
  
  double Lambda[NUM_VAR];            // Eigenvalue vector (diagonal matrix).
  double T[NUM_VAR*NUM_VAR];         // Eigenvector matrices.
  double Tinv[NUM_VAR*NUM_VAR];
  double temp[NUM_VAR*NUM_VAR];      // Temporary matrix.
  double Roemat[NUM_VAR*NUM_VAR];    // The Roe matrix.

  // Efficiency variables.
  double oo_2cbar;                   // One over 2 times cbar.
  double oo_2gm1;                    // One over 2 times (gamma - 1).
  double oo_rhobar;                  // One over rhobar;
  double oo_rcbar;                   // One over rhobar times cbar.
  double oo_cbar2;                   // One over cbar squared.

  // Get the right and left states.
  if ( RECON_PRIM )
    {
      rhoL = QL[0];
      uL = QL[1];
      vL = QL[2];
      pL = QL[3];
      EL = pL/(gamma-1.) + 0.5*rhoL*(uL*uL + vL*vL);
      HL = EL + pL;
      hL = HL / rhoL;

      rhoR = QR[0];
      uR = QR[1];
      vR = QR[2];
      pR = QR[3];
      ER = pR/(gamma-1.) + 0.5*rhoR*(uR*uR + vR*vR);
      HR = ER + pR;
      hR = HR / rhoR;
    }
  else
    {
      rhoL = QL[0];
      uL = QL[1]/QL[0];
      vL = QL[2]/QL[0];
      EL = QL[3];
      pL = (gamma - 1.)*(EL - 0.5*rhoL*(uL*uL + vL*vL) );
      HL = EL + pL;
      hL = HL / rhoL;

      rhoR = QR[0];
      uR = QR[1]/QR[0];
      vR = QR[2]/QR[0];
      ER = QR[3];
      pR = (gamma - 1.)*(ER - 0.5*rhoR*(uR*uR + vR*vR) );
      HR = ER + pR;
      hR = HR / rhoR;
    }

  // Generate the Roe averaged values.
  rhobar = sqrt(rhoL*rhoR);
  big_gamma = rhobar / ( rhoL + rhobar ) ;

  ubar = uL + big_gamma*(uR - uL);
  vbar = vL + big_gamma*(vR - vL);
  Hbar = HL + big_gamma*(HR - HL);
  hbar = hL + big_gamma*(hR - hL);

  thetabar = ubar*nx + vbar*ny;

  V2bar = (ubar*ubar + vbar*vbar);

  //pbar = (gamma-1.)/gamma * ( Hbar - rhobar*V2bar*0.5 );
  //cbar = sqrt( (gamma*pbar)/rhobar );

  //cbar = sqrt( (gamma-1.0)*(Hbar - 0.5*V2bar) );
  cbar = sqrt( (gamma-1.0)*(hbar - 0.5*V2bar) );
  

  // Set the efficiency values.
  oo_2cbar = 1. / ( 2.*cbar );
  oo_2gm1 = 1. / ( 2.*(gamma - 1.) );
  oo_rhobar = 1. / rhobar;
  oo_rcbar = 1. / ( rhobar * cbar );
  oo_cbar2 = 1. / ( cbar * cbar );

  // Set the eigenvalues.
  ev1 = thetabar;
  ev2 = thetabar;
  ev3 = thetabar + cbar;
  ev4 = thetabar - cbar;

  // Get the negative eigenvalues since this is the Roe scheme I'm going to use.
  Lambda[0] = fabs(ev1);
  Lambda[1] = fabs(ev2);
  Lambda[2] = fabs(ev3);
  Lambda[3] = fabs(ev4);

  // Now I can fill in the eigenvector matrices. I have the derivation of them in my notes with help from Hirsch's book (volume 2).
  // Row major format.
  /*
  T[0] = 1.0;
  T[1] = 0.;
  T[2] = rhobar*oo_2cbar;
  T[3] = rhobar*oo_2cbar;
  
  T[4] = ubar;
  T[5] = rhobar*ny;
  T[6] = rhobar*ubar*oo_2cbar + 0.5*rhobar*nx;
  T[7] = rhobar*ubar*oo_2cbar - 0.5*rhobar*nx;
  
  T[8] =  vbar;
  T[9] =  -rhobar*nx;
  T[10] = rhobar*vbar*oo_2cbar + 0.5*rhobar*ny;
  T[11] = rhobar*vbar*oo_2cbar - 0.5*rhobar*ny;

  T[12] = 0.5*V2bar;
  T[13] = rhobar*( ubar*ny - vbar*nx );
  T[14] = 0.5*rhobar*oo_2cbar*V2bar + rhobar*cbar*oo_2gm1 + 0.5*rhobar*thetabar;
  T[15] = 0.5*rhobar*oo_2cbar*V2bar + rhobar*cbar*oo_2gm1 - 0.5*rhobar*thetabar;

  //============================================================================

  Tinv[0] = 1. - 0.5*oo_cbar2*(gamma-1.)*V2bar;
  Tinv[1] = (gamma-1.)*ubar*oo_cbar2;
  Tinv[2] = (gamma-1.)*vbar*oo_cbar2;
  Tinv[3] = -(gamma-1.)*oo_cbar2;

  Tinv[4] = oo_rhobar*( -ubar*ny + vbar*nx );
  Tinv[5] = ny*oo_rhobar;
  Tinv[6] = -nx*oo_rhobar;
  Tinv[7] = 0.;

  Tinv[8]  = -thetabar*oo_rhobar + 0.5*oo_rcbar*(gamma-1.)*V2bar;
  Tinv[9]  = nx*oo_rhobar - (gamma-1.)*ubar*oo_rcbar;
  Tinv[10] = ny*oo_rhobar - (gamma-1.)*vbar*oo_rcbar;
  Tinv[11] = (gamma-1.)*oo_rcbar;

  Tinv[12] = thetabar*oo_rhobar + 0.5*oo_rcbar*(gamma-1.)*V2bar;
  Tinv[13] = -nx*oo_rhobar - (gamma-1.)*ubar*oo_rcbar;
  Tinv[14] = -ny*oo_rhobar - (gamma-1.)*vbar*oo_rcbar;
  Tinv[15] = (gamma-1.)*oo_rcbar;
  */


  // Here I am using the matrices given in Nejat's dissertation. L is T, R is T inverse.
  T[0] = 1.0;
  T[1] = 0.0;
  T[2] = 1.0 / (cbar*cbar);
  T[3] = 1.0 / (cbar*cbar);
  
  T[4] = ubar;
  T[5] = -ny;
  T[6] = ubar / (cbar*cbar) + nx / cbar;
  T[7] = ubar / (cbar*cbar) - nx / cbar;
  
  T[8] = vbar;
  T[9] = nx;
  T[10] = vbar / (cbar*cbar) + ny / cbar;
  T[11] = vbar / (cbar*cbar) - ny / cbar;

  T[12] = 0.5*V2bar;
  T[13] = ubar*(-ny) + vbar*nx;
  T[14] = 1.0/(gamma-1.0) + thetabar/cbar + (0.5*V2bar)/(cbar*cbar);
  T[15] = 1.0/(gamma-1.0) - thetabar/cbar + (0.5*V2bar)/(cbar*cbar);

  //============================================================================

  Tinv[0] = 1.0 - (gamma-1.0)*(0.5*V2bar)/(cbar*cbar);
  Tinv[1] = (gamma-1.0)*ubar/(cbar*cbar);
  Tinv[2] = (gamma-1.0)*vbar/(cbar*cbar);
  Tinv[3] = -1.0*(gamma-1.0)/(cbar*cbar);

  Tinv[4] = -1.0 * ( ubar*(-ny) + vbar*nx );
  Tinv[5] = -ny;
  Tinv[6] = nx;
  Tinv[7] = 0.;

  Tinv[8]  = 0.5 * ( -cbar * thetabar + (gamma-1.0)*0.5*V2bar );
  Tinv[9]  = 0.5 * ( cbar * nx - (gamma-1.0)*ubar );
  Tinv[10] = 0.5 * ( cbar * ny - (gamma-1.0)*vbar );
  Tinv[11] = 0.5 * ( gamma - 1.0 );

  Tinv[12] = 0.5 * ( cbar * thetabar + (gamma-1.0)*0.5*V2bar );
  Tinv[13] = -0.5 * ( cbar * nx + (gamma-1.0)*ubar );
  Tinv[14] = -0.5 * ( cbar * ny + (gamma-1.0)*vbar );
  Tinv[15] = 0.5 * ( gamma - 1.0 );


  // Now multiply Lambda times Tinv.
  MMMult_DtG ( Lambda, Tinv, temp, NUM_VAR );

  // And T times the previous result.
  MMMult_GtG ( T, temp, Roemat, NUM_VAR );  // Roe matrix is now filled in.

  // Now we can build the Euler Flux Jacobians from the appropriate states.

  // Left state.
  
  theta = uL*nx + vL*ny;
  phi = 0.5 * ( gamma - 1. ) * ( uL*uL + vL*vL );
  
  dfdqL[0][0] = 0.;
  dfdqL[0][1] = nx;
  dfdqL[0][2] = ny;
  dfdqL[0][3] = 0.;

  dfdqL[1][0] = phi*nx - uL*theta;
  dfdqL[1][1] = (2.-gamma)*uL*nx + theta;
  dfdqL[1][2] = uL*ny - (gamma-1.)*vL*nx;
  dfdqL[1][3] = (gamma-1.)*nx;
  
  dfdqL[2][0] = phi*ny - vL*theta;
  dfdqL[2][1] = vL*nx - (gamma-1.)*uL*ny;
  dfdqL[2][2] = theta + (2.-gamma)*vL*ny;
  dfdqL[2][3] = (gamma-1.)*ny;
  
  dfdqL[3][0] = (2.*phi - gamma*EL/rhoL ) * theta;
  dfdqL[3][1] = (1.-gamma)*uL*theta + nx*(gamma*EL/rhoL - phi);
  dfdqL[3][2] = (1.-gamma)*vL*theta + ny*(gamma*EL/rhoL - phi);
  dfdqL[3][3] = gamma*theta;

  // Loop through and complete the computation before accumulation.
  for ( i=0; i < NUM_VAR; i++ )
    {
      for ( k=0; k < NUM_VAR; k++ )
	{
	  dfdqL[i][k] = 0.5 * ( dfdqL[i][k] + Roemat[i*NUM_VAR + k] ) * len;
	}
    }

  // Now I can store the results in the matrix.

  n = grid->bedges[ibedge*5+0];    // Get the real node index.

  // NaN check.
  
  for ( i=0; i < NUM_VAR; i++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  if ( !(isfinite(dfdqL[i][j])) )
	    {
	      printf("NAN DETECTED in dfdqR[%d][%d] on bedge=%d - node = %d.\n",i,j,ibedge,n);
	      printf("Problem DETECTED in dfdqL[%d][%d]=%f  on bedge=%d - node = %d\n",i,j,(float)dfdqL[i][j],ibedge,n);
	    }
	}
    }
  

  // Now add dfdqL to the diagonal of the physical boundary node.
  for ( i=0; i < NUM_VAR; i++ )
    {
      for ( k=0; k < NUM_VAR; k++ )
	{
	  LHS[ (iau[n])*NUM_VAR*NUM_VAR + i*NUM_VAR + k] += dfdqL[i][k];
	}
    }

  return;
}
