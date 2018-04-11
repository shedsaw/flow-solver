//=============================================================
// 
//  flux.cpp
//  
//  Functions to calculate the flux across a face.
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
#include "flux.h"
#include "linear_algebra.h"


//=============================================================
// 
//  Incompressible_Euler_Flux()
//
//  Calculates the euler flux vector from the state vector.
//
//  double *Q                        // The state vector.
//  double *Eulerflux                // The Euler flux vector.
//  double nx                        // Normal vector
//  double ny
//  double len                       // Length of the normal vector.
//  double gamma
//
//=============================================================

void Incompressible_Euler_Flux ( double *Q, double *Eulerflux, double nx, double ny, double len, double gamma )
{
  double P;                   // Pressure
  double rho,u,v,E;           // variables.
  double theta;               // velocity dot n

  if ( RECON_PRIM )
    {
      rho = Q[0];
      u = Q[1];
      v = Q[2];
      P = Q[3];
      E = P/(gamma-1.) + 0.5*rho*(u*u + v*v);
    }
  else
    {
      rho = Q[0];
      u = Q[1] / Q[0];
      v = Q[2] / Q[0];
      E = Q[3];
      P = (gamma - 1.0) * ( E - 0.5*rho*(u*u + v*v) );
    }

  theta = u*nx + v*ny;

  Eulerflux[0] = rho*theta*len;
  Eulerflux[1] = rho*u*theta*len + P*nx*len;
  Eulerflux[2] = rho*v*theta*len + P*ny*len;
  Eulerflux[3] = theta*(E + P)*len;

  return;
}


//=============================================================
// 
//  Roe_flux()
//
//  Calculates the flux through a portion of control volume face
//  using Roe's scheme.
//
//  Returns an error code.
//
//  double nx                        // Normal in the x direction.
//  double ny                        // Normal in the y direction.
//  double len                       // The edge length.
//  double gamma                     // The ratio of specific heats.
//  double qleft[4]                  // Values from the left.
//  double qright[4]                 // Values from the right.
//  double flux[4]                   // The flux to return.
//
//=============================================================

int Roe_flux ( double nx, double ny, double len, double gamma,
	       double qleft[4], double qright[4], double flux[4])
{
  int i;                             // Loop counters.
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
  double EL,ER;
  double big_gamma;                  // Multiplication factor to reduce divisions. See Daniel's notes.

  double ev1,ev2,ev3,ev4;            // Eigenvalues.
  
  double fluxvector[NUM_VAR];        // Incompressible Euler flux vector.
  double diss[NUM_VAR];              // The dissapation added by the Roe scheme.
  double dQ[NUM_VAR];                // Jump in state variables.
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
      rhoL = qleft[0];
      uL = qleft[1];
      vL = qleft[2];
      pL = qleft[3];
      EL = pL/(gamma-1.) + 0.5*rhoL*(uL*uL + vL*vL);
      HL = EL + pL;

      rhoR = qright[0];
      uR = qright[1];
      vR = qright[2];
      pR = qright[3];
      ER = pR/(gamma-1.) + 0.5*rhoR*(uR*uR + vR*vR);
      HR = ER + pR;
    }
  else
    {
      rhoL = qleft[0];
      uL = qleft[1]/qleft[0];
      vL = qleft[2]/qleft[0];
      EL = qleft[3];
      pL = (gamma - 1.)*(EL - 0.5*rhoL*(uL*uL + vL*vL) );
      HL = EL + pL;

      rhoR = qright[0];
      uR = qright[1]/qright[0];
      vR = qright[2]/qright[0];
      ER = qright[3];
      pR = (gamma - 1.)*(ER - 0.5*rhoR*(uR*uR + vR*vR) );
      HR = ER + pR;
    }

  // Generate the Roe averaged values.
  rhobar = sqrt(rhoL*rhoR);
  big_gamma = rhobar / ( rhoL + rhobar ) ;

  ubar = uL + big_gamma*(uR - uL);
  vbar = vL + big_gamma*(vR - vL);
  Hbar = HL + big_gamma*(HR - HL);

  thetabar = ubar*nx + vbar*ny;

  V2bar = (ubar*ubar + vbar*vbar);

  pbar = (gamma-1.)/gamma * ( Hbar - rhobar*V2bar*0.5 );
  cbar = sqrt( (gamma*pbar)/rhobar );

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
  Lambda[0] = 0.5 * ( ev1 - fabs(ev1) );
  Lambda[1] = 0.5 * ( ev2 - fabs(ev2) );
  Lambda[2] = 0.5 * ( ev3 - fabs(ev3) );
  Lambda[3] = 0.5 * ( ev4 - fabs(ev4) );
  
  // Get the jump across the interface.
  if ( RECON_PRIM )
    {
      dQ[0] = rhoR - rhoL;
      dQ[1] = rhoR*uR - rhoL*uL;
      dQ[2] = rhoR*vR - rhoL*vL;
      dQ[3] = ER - EL;
    }
  else
    {
      for ( i=0; i < NUM_VAR; i++ )
	{
	  dQ[i] = qright[i] - qleft[i];
	}
    }

  // Now I can fill in the eigenvector matrices. I have the derivation of them in my notes with help from Hirsch's book (volume 2).
  // Row major format.
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

  // Now multiply Lambda times Tinv.
  MMMult_DtG ( Lambda, Tinv, temp, NUM_VAR );

  // And T times the previous result.
  MMMult_GtG ( T, temp, Roemat, NUM_VAR );

  // Multiply this result with the jump.
  MVMult ( Roemat, dQ, diss, NUM_VAR );

  // Scale this result by the area/length.
  for ( i=0; i < NUM_VAR; i++ )
    {
      diss[i] = diss[i] * len;
    }

  // Now get the Euler flux.
  Incompressible_Euler_Flux( qleft, fluxvector, nx, ny, len, gamma );

  // Add the result.
  for ( i=0; i < NUM_VAR; i++ )
    flux[i] = fluxvector[i] + diss[i];
  
  return 1;
}


//=============================================================
// 
//  Roe_flux_centered()
//
//  Calculates the flux through a portion of control volume face
//  using Roe's scheme.
//
//  Returns an error code.
//
//  double nx                        // Normal in the x direction.
//  double ny                        // Normal in the y direction.
//  double len                       // The edge length.
//  double gamma                     // The ratio of specific heats.
//  double qleft[4]                  // Values from the left.
//  double qright[4]                 // Values from the right.
//  double flux[4]                   // The flux to return.
//
//=============================================================

int Roe_flux_centered ( double nx, double ny, double len, double gamma,
			double qleft[4], double qright[4], double flux[4])
{
  int i;                             // Loop counters.
  double rhobar;                     // Roe averaged density.
  double ubar;                       // Roe aveeraged velocities.
  double vbar;
  double Hbar;                       // Roe averaged enthalpy.
  double cbar;                       // Roe averaged speed of sound.
  double pbar;                       // Roe averaged pressure.
  double thetabar;                   // Roe averaged theta.
  double V2bar;                      // Roe averaged V^2.
  double rhoL,uL,vL,HL,pL,EL;        // Left state.
  double rhoR,uR,vR,HR,pR,ER;        // Right state.
  double big_gamma;                  // Multiplication factor to reduce divisions. See Daniel's notes.
  double hL,hR,hbar;
  double ev1,ev2,ev3,ev4;            // Eigenvalues.
  
  double fluxvectorL[NUM_VAR];       // Incompressible Euler flux vector.
  double fluxvectorR[NUM_VAR];       // Incompressible Euler flux vector.
  double diss[NUM_VAR];              // The dissapation added by the Roe scheme.
  double dQ[NUM_VAR];                // Jump in state variables.
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
      rhoL = qleft[0];
      uL = qleft[1];
      vL = qleft[2];
      pL = qleft[3];
      EL = pL/(gamma-1.) + 0.5*rhoL*(uL*uL + vL*vL);
      HL = EL + pL;
      hL = HL / rhoL;

      rhoR = qright[0];
      uR = qright[1];
      vR = qright[2];
      pR = qright[3];
      ER = pR/(gamma-1.) + 0.5*rhoR*(uR*uR + vR*vR);
      HR = ER + pR;
      hR = HR / rhoR;
    }
  else
    {
      rhoL = qleft[0];
      uL = qleft[1]/qleft[0];
      vL = qleft[2]/qleft[0];
      EL = qleft[3];
      pL = (gamma - 1.)*(EL - 0.5*rhoL*(uL*uL + vL*vL) );
      HL = EL + pL;
      hL = HL / rhoL;

      rhoR = qright[0];
      uR = qright[1]/qright[0];
      vR = qright[2]/qright[0];
      ER = qright[3];
      pR = (gamma - 1.)*(ER - 0.5*rhoR*(uR*uR + vR*vR) );
      HR = ER + pR;
      hR = HR/ rhoR;
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
  
  // Get the jump across the interface.
  if ( RECON_PRIM )
    {
      dQ[0] = rhoR - rhoL;
      dQ[1] = rhoR*uR - rhoL*uL;
      dQ[2] = rhoR*vR - rhoL*vL;
      dQ[3] = ER - EL;
    }
  else
    {
      for ( i=0; i < NUM_VAR; i++ )
	{
	  dQ[i] = qright[i] - qleft[i];
	}
    }
  
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
  MMMult_GtG ( T, temp, Roemat, NUM_VAR );

  // Multiply this result with the jump.
  MVMult ( Roemat, dQ, diss, NUM_VAR );

  // Scale this result by the area/length.
  for ( i=0; i < NUM_VAR; i++ )
    {
      diss[i] = diss[i] * len;
    }

  // Now get the Euler flux.
  Incompressible_Euler_Flux( qleft, fluxvectorL, nx, ny, len, gamma );
  Incompressible_Euler_Flux( qright, fluxvectorR, nx, ny, len, gamma );

  // Add the result.
  for ( i=0; i < NUM_VAR; i++ )
    flux[i] = 0.5 * ( fluxvectorL[i] + fluxvectorR[i] - diss[i] );
  
  return 1;
}


//=============================================================
// 
//  Van_Leer_flux()
//
//  Calculates the flux through a portion of control volume face
//  using Van Leer's Flux Vector Splitting scheme.
//
//  double nx                        // Normal in the x direction.
//  double ny                        // Normal in the y direction.
//  double len                       // The edge length.
//  double gamma                     // The ratio of specific heats.
//  double qleft[4]                  // Values from the left.
//  double qright[4]                 // Values from the right.
//  double flux[4]                   // The flux to return.
//  double dfdqp[4][4]               // The flux Jacobian from F+.
//  double dfdqm[4][4]               // The flux Jacobian from F-.
//
//=============================================================

int Van_Leer_flux ( double nx, double ny, double len, double gamma,
		     double qleft[4], double qright[4],
		     double flux[4], double dfdqp[][4], double dfdqm[][4] )
{
  double pl;                                  // Pressure on the left.
  double pr;                                  // Pressure on the right.
  double cl;                                  // Speed of sound on left.
  double cr;                                  // Speed of sound on right.
  double ubarl;                               // Ubar on the left.
  double ubarr;                               // Ubar on the right.
  double cond;                                // Test condition to determine if supersonic.
  double fluxp[4];                            // F+ flux.
  double fluxm[4];                            // F- flux.
  double rho;                                 // Density.
  double u;                                   // X-velocity.
  double v;                                   // Y-velocity.
  double e;                                   // Energy.
  double drhodq1;                             // Partials of density wrt to the depedent variables.
  double drhodq2;                             
  double drhodq3;
  double drhodq4;
  double dudq1;                               // Partials of x-velocity wrt to the dependent variables.
  double dudq2;
  double dudq3;
  double dudq4;
  double dvdq1;                               // Partials of y-velocity wrt to the dependent variables.
  double dvdq2;
  double dvdq3;
  double dvdq4;
  double detdq1;                              // Partials of total energy wrt to the dependent variables.
  double detdq2;
  double detdq3;
  double detdq4;
  double dpdq1;                               // Partials of pressure wrt to the dependent variables.
  double dpdq2;
  double dpdq3;
  double dpdq4;
  double dcdq1;                               // Partials of speed of sound wrt to the dependent variables.
  double dcdq2;
  double dcdq3;
  double dcdq4;
  double dubardq1;                            // Partials of Ubar wrt to the dependent variables.
  double dubardq2;
  double dubardq3;
  double dubardq4;
  int flag_p_g = 0;                           //|
  int flag_p_l = 0;                           //|
  int flag_m_g = 0;                           //|--->  Test flags.
  int flag_m_l = 0;                           //|

  pl = (gamma-1.)*qleft[3] - ((gamma-1.)/2.)*(qleft[1]*qleft[1] + qleft[2]*qleft[2])/(qleft[0]);
  pr = (gamma-1.)*qright[3] - ((gamma-1.)/2.)*(qright[1]*qright[1] + qright[2]*qright[2])/(qright[0]);

  // Check some values.
  if ( qleft[0] <= 0. )
    {
      printf("FATAL ERROR: In Van Leer Flux.\n");
      printf("  Negative density was supplied to the left side.\n");
      printf("  rhoL = %.15E\n",qleft[0]);
      return 0;
    }

  if ( qright[0] <= 0. )
    {
      printf("FATAL ERROR: In Van Leer Flux.\n");
      printf("  Negative density was supplied to the right side.\n");
      printf("  rhoR = %.15E\n",qright[0]);
      return 0;
    }

  if ( pl <= 0. )
    {
      printf("FATAL ERROR: In Van Leer Flux.\n");
      printf("  Negative pressure was detected on the left side.\n");
      printf("  PL = %.15E\n",pl);
      return 0;
    }

  if ( pr <= 0. )
    {
      printf("FATAL ERROR: In Van Leer Flux.\n");
      printf("  Negative pressure was detected on the right side.\n");
      printf("  PR = %.15E\n",pr);
      return 0;
    }
  
  cl = sqrt( gamma*pl/qleft[0] );
  cr = sqrt( gamma*pr/qright[0] );

  ubarl = (qleft[1]*nx + qleft[2]*ny)/qleft[0];
  ubarr = (qright[1]*nx + qright[2]*ny)/qright[0];

  cond = ubarl/cl;

  // Grab the primitive variables from the left.
  rho = qleft[0];
  u = qleft[1] / qleft[0];
  v = qleft[2] / qleft[0];
  e = rho*( pl/((gamma-1.)*rho) + 0.5*(u*u + v*v) );

  // Calculate the partials to be used in evaulating the flux Jacobians. We start with Qleft.
  drhodq1 = 1.;
  drhodq2 = 0.;
  drhodq3 = 0.;
  drhodq4 = 0.;
  
  detdq1 = 0.;
  detdq2 = 0.;
  detdq3 = 0.;
  detdq4 = 1.;

  dudq1 = -qleft[1]/(qleft[0]*qleft[0]);
  dudq2 = 1./qleft[0];
  dudq3 = 0.;
  dudq4 = 0.;

  dvdq1 = -qleft[2]/(qleft[0]*qleft[0]);
  dvdq2 = 0.;
  dvdq3 = 1./qleft[0];
  dvdq4 = 0.;

  dpdq1 = (gamma-1.)/2. * (qleft[1]*qleft[1] + qleft[2]*qleft[2])/(qleft[0]*qleft[0]);
  dpdq2 = -(gamma-1.) * (qleft[1]/qleft[0]);
  dpdq3 = -(gamma-1.) * (qleft[2]/qleft[0]);
  dpdq4 = gamma - 1.;

  dcdq1 = ( gamma*qleft[0]*dpdq1 - gamma*pl*drhodq1 ) / ( 2.*cl*qleft[0]*qleft[0] );
  dcdq2 = ( gamma*qleft[0]*dpdq2 - gamma*pl*drhodq2 ) / ( 2.*cl*qleft[0]*qleft[0] );
  dcdq3 = ( gamma*qleft[0]*dpdq3 - gamma*pl*drhodq3 ) / ( 2.*cl*qleft[0]*qleft[0] );
  dcdq4 = ( gamma*qleft[0]*dpdq4 - gamma*pl*drhodq4 ) / ( 2.*cl*qleft[0]*qleft[0] );

  dubardq1 = -qleft[1]/(qleft[0]*qleft[0])*nx - qleft[2]/(qleft[0]*qleft[0])*ny;
  dubardq2 = nx / qleft[0];
  dubardq3 = ny / qleft[0];
  dubardq4 = 0.;
  
  //The logic here ensures that only the proper fluxes are initialized for supersonic flows.
  if ( cond > 1. )
    {
      fluxp[0] = qleft[0]*ubarl*len;
      fluxp[1] = (qleft[1]*ubarl + nx*pl)*len;
      fluxp[2] = (qleft[2]*ubarl + ny*pl)*len;
      fluxp[3] = (qleft[3]+pl)*ubarl*len;

      //Compute the fluxes using analytical derivatives.
      //Do DF/DQ and take this to be DF+/DQ:
      //First row
      dfdqp[0][0] = (ubarl*drhodq1 + rho*dubardq1)*len;
      dfdqp[0][1] = (ubarl*drhodq2 + rho*dubardq2)*len;
      dfdqp[0][2] = (ubarl*drhodq3 + rho*dubardq3)*len;
      dfdqp[0][3] = (ubarl*drhodq4 + rho*dubardq4)*len;
      
      //Second row.
      dfdqp[1][0] = (drhodq1*u*ubarl + dudq1*rho*ubarl + dubardq1*rho*u + dpdq1*nx)*len;
      dfdqp[1][1] = (drhodq2*u*ubarl + dudq2*rho*ubarl + dubardq2*rho*u + dpdq2*nx)*len;
      dfdqp[1][2] = (drhodq3*u*ubarl + dudq3*rho*ubarl + dubardq3*rho*u + dpdq3*nx)*len;
      dfdqp[1][3] = (drhodq4*u*ubarl + dudq4*rho*ubarl + dubardq4*rho*u + dpdq4*nx)*len;
      
      //Third row.
      dfdqp[2][0] = (drhodq1*v*ubarl + dvdq1*rho*ubarl + dubardq1*rho*v + dpdq1*ny)*len;
      dfdqp[2][1] = (drhodq2*v*ubarl + dvdq2*rho*ubarl + dubardq2*rho*v + dpdq2*ny)*len;
      dfdqp[2][2] = (drhodq3*v*ubarl + dvdq3*rho*ubarl + dubardq3*rho*v + dpdq3*ny)*len;
      dfdqp[2][3] = (drhodq4*v*ubarl + dvdq4*rho*ubarl + dubardq4*rho*v + dpdq4*ny)*len;
      
      //Fourth Row.
      dfdqp[3][0] = ( (e + pl)*dubardq1 + (detdq1 + dpdq1)*ubarl )*len;
      dfdqp[3][1] = ( (e + pl)*dubardq2 + (detdq2 + dpdq2)*ubarl )*len;
      dfdqp[3][2] = ( (e + pl)*dubardq3 + (detdq3 + dpdq3)*ubarl )*len;
      dfdqp[3][3] = ( (e + pl)*dubardq4 + (detdq4 + dpdq4)*ubarl )*len;
      
      flag_p_g = 1;
    }
  else if ( cond < -1. )
    {
      fluxp[0] = 0.;
      fluxp[1] = 0.;
      fluxp[2] = 0.;
      fluxp[3] = 0.;
      
      dfdqp[0][0] = dfdqp[0][1] = dfdqp[0][2] = dfdqp[0][3] = 0.;
      dfdqp[1][0] = dfdqp[1][1] = dfdqp[1][2] = dfdqp[1][3] = 0.;
      dfdqp[2][0] = dfdqp[2][1] = dfdqp[2][2] = dfdqp[2][3] = 0.;
      dfdqp[3][0] = dfdqp[3][1] = dfdqp[3][2] = dfdqp[3][3] = 0.;
      
      flag_p_l = 1;
    }
  else
    {
      //Calculate the +/- Van Leer flux vectors according to our notes.
      fluxp[0] = 0.25*qleft[0]*cl*( (ubarl/cl + 1.0)*(ubarl/cl + 1.0) ) * len;
      fluxp[1] = fluxp[0] * ( (nx/gamma)*(-ubarl + 2.0*cl) + u );
      fluxp[2] = fluxp[0] * ( (ny/gamma)*(-ubarl + 2.0*cl) + v );
      fluxp[3] = fluxp[0] * ( ((1.0-gamma)*ubarl*ubarl + 2.0*(gamma-1.0)*ubarl*cl + 2.0*cl*cl)/(gamma*gamma-1.0) 
			      + (qleft[1]*qleft[1]+qleft[2]*qleft[2])/(2.*qleft[0]*qleft[0]) );

      //Compute the fluxes using analytical derivatives.
      //First do DF/DQ+:
      //First row
      dfdqp[0][0] = (0.25*cl*(ubarl/cl + 1.)*(ubarl/cl + 1.)*drhodq1 + 0.25*rho*(ubarl/cl + 1.)*(ubarl/cl + 1.)*dcdq1
		     + 0.5*rho*(ubarl/cl + 1.)*dubardq1 - 0.5*((rho*ubarl)/cl)*(ubarl/cl + 1.)*dcdq1)*len;
      dfdqp[0][1] = (0.25*cl*(ubarl/cl + 1.)*(ubarl/cl + 1.)*drhodq2 + 0.25*rho*(ubarl/cl + 1.)*(ubarl/cl + 1.)*dcdq2
		     + 0.5*rho*(ubarl/cl + 1.)*dubardq2 - 0.5*((rho*ubarl)/cl)*(ubarl/cl + 1.)*dcdq2)*len;
      dfdqp[0][2] = (0.25*cl*(ubarl/cl + 1.)*(ubarl/cl + 1.)*drhodq3 + 0.25*rho*(ubarl/cl + 1.)*(ubarl/cl + 1.)*dcdq3
		     + 0.5*rho*(ubarl/cl + 1.)*dubardq3 - 0.5*((rho*ubarl)/cl)*(ubarl/cl + 1.)*dcdq3)*len;
      dfdqp[0][3] = (0.25*cl*(ubarl/cl + 1.)*(ubarl/cl + 1.)*drhodq4 + 0.25*rho*(ubarl/cl + 1.)*(ubarl/cl + 1.)*dcdq4
		     + 0.5*rho*(ubarl/cl + 1.)*dubardq4 - 0.5*((rho*ubarl)/cl)*(ubarl/cl + 1.)*dcdq4)*len;
      
      //Second row.
      dfdqp[1][0] = dfdqp[0][0]*( (nx/gamma)*(-ubarl + 2.0*cl) + u ) + fluxp[0]*( (nx/gamma)*(-dubardq1 + 2.0*dcdq1) +dudq1 );
      dfdqp[1][1] = dfdqp[0][1]*( (nx/gamma)*(-ubarl + 2.0*cl) + u ) + fluxp[0]*( (nx/gamma)*(-dubardq2 + 2.0*dcdq2) +dudq2 );
      dfdqp[1][2] = dfdqp[0][2]*( (nx/gamma)*(-ubarl + 2.0*cl) + u ) + fluxp[0]*( (nx/gamma)*(-dubardq3 + 2.0*dcdq3) +dudq3 );
      dfdqp[1][3] = dfdqp[0][3]*( (nx/gamma)*(-ubarl + 2.0*cl) + u ) + fluxp[0]*( (nx/gamma)*(-dubardq4 + 2.0*dcdq4) +dudq4 );
      
      //Third row.
      dfdqp[2][0] = dfdqp[0][0]*( (ny/gamma)*(-ubarl + 2.0*cl) + v ) + fluxp[0]*( (ny/gamma)*(-dubardq1 + 2.0*dcdq1) +dvdq1 );
      dfdqp[2][1] = dfdqp[0][1]*( (ny/gamma)*(-ubarl + 2.0*cl) + v ) + fluxp[0]*( (ny/gamma)*(-dubardq2 + 2.0*dcdq2) +dvdq2 );
      dfdqp[2][2] = dfdqp[0][2]*( (ny/gamma)*(-ubarl + 2.0*cl) + v ) + fluxp[0]*( (ny/gamma)*(-dubardq3 + 2.0*dcdq3) +dvdq3 );
      dfdqp[2][3] = dfdqp[0][3]*( (ny/gamma)*(-ubarl + 2.0*cl) + v ) + fluxp[0]*( (ny/gamma)*(-dubardq4 + 2.0*dcdq4) +dvdq4 );
      
      //Fourth Row.
      dfdqp[3][0] = dfdqp[0][0]*( ((1.0-gamma)*ubarl*ubarl + 2.0*(gamma-1.0)*ubarl*cl + 2.0*cl*cl)/(gamma*gamma-1.0) + (u*u + v*v)/2.0 )
	+ fluxp[0]*( ((1.0-gamma)*2.*ubarl*dubardq1 + 2.0*(gamma-1.0)*(dubardq1*cl + dcdq1*ubarl) + 4.0*cl*dcdq1)/(gamma*gamma-1.0) 
		     + (u*dudq1 + v*dvdq1) );
      dfdqp[3][1] = dfdqp[0][1]*( ((1.0-gamma)*ubarl*ubarl + 2.0*(gamma-1.0)*ubarl*cl + 2.0*cl*cl)/(gamma*gamma-1.0) + (u*u + v*v)/2.0 )
	+ fluxp[0]*( ((1.0-gamma)*2.*ubarl*dubardq2 + 2.0*(gamma-1.0)*(dubardq2*cl + dcdq2*ubarl) + 4.0*cl*dcdq2)/(gamma*gamma-1.0) 
		     + (u*dudq2 + v*dvdq2) );
      dfdqp[3][2] = dfdqp[0][2]*( ((1.0-gamma)*ubarl*ubarl + 2.0*(gamma-1.0)*ubarl*cl + 2.0*cl*cl)/(gamma*gamma-1.0) + (u*u + v*v)/2.0 )
	+ fluxp[0]*( ((1.0-gamma)*2.*ubarl*dubardq3 + 2.0*(gamma-1.0)*(dubardq3*cl + dcdq3*ubarl) + 4.0*cl*dcdq3)/(gamma*gamma-1.0) 
		     + (u*dudq3 + v*dvdq3) );
      dfdqp[3][3] = dfdqp[0][3]*( ((1.0-gamma)*ubarl*ubarl + 2.0*(gamma-1.0)*ubarl*cl + 2.0*cl*cl)/(gamma*gamma-1.0) + (u*u + v*v)/2.0 )
	+ fluxp[0]*( ((1.0-gamma)*2.*ubarl*dubardq4 + 2.0*(gamma-1.0)*(dubardq4*cl + dcdq4*ubarl) + 4.0*cl*dcdq4)/(gamma*gamma-1.0) 
		     + (u*dudq4 + v*dvdq4) );
    }

  cond = ubarr/cr;

  // Reset the primitive variables with data from the right.
  rho = qright[0];
  u = qright[1] / qright[0];
  v = qright[2] / qright[0];
  e = rho*( pr/((gamma-1.)*rho) + 0.5*(u*u + v*v) );

  // Calculate the partials to be used in evaulating the flux Jacobians using Qright.
  drhodq1 = 1.;
  drhodq2 = 0.;
  drhodq3 = 0.;
  drhodq4 = 0.;
  
  detdq1 = 0.;
  detdq2 = 0.;
  detdq3 = 0.;
  detdq4 = 1.;

  dudq1 = -qright[1]/(qright[0]*qright[0]);
  dudq2 = 1./qright[0];
  dudq3 = 0.;
  dudq4 = 0.;

  dvdq1 = -qright[2]/(qright[0]*qright[0]);
  dvdq2 = 0.;
  dvdq3 = 1./qright[0];
  dvdq4 = 0.;

  dpdq1 = (gamma-1.)/2. * (qright[1]*qright[1] + qright[2]*qright[2])/(qright[0]*qright[0]);
  dpdq2 = -(gamma-1.) * (qright[1]/qright[0]);
  dpdq3 = -(gamma-1.) * (qright[2]/qright[0]);
  dpdq4 = gamma - 1.;

  dcdq1 = ( gamma*qright[0]*dpdq1 - gamma*pr*drhodq1 ) / ( 2.*cr*qright[0]*qright[0] );
  dcdq2 = ( gamma*qright[0]*dpdq2 - gamma*pr*drhodq2 ) / ( 2.*cr*qright[0]*qright[0] );
  dcdq3 = ( gamma*qright[0]*dpdq3 - gamma*pr*drhodq3 ) / ( 2.*cr*qright[0]*qright[0] );
  dcdq4 = ( gamma*qright[0]*dpdq4 - gamma*pr*drhodq4 ) / ( 2.*cr*qright[0]*qright[0] );

  dubardq1 = -qright[1]/(qright[0]*qright[0])*nx - qright[2]/(qright[0]*qright[0])*ny;
  dubardq2 = nx / qright[0];
  dubardq3 = ny / qright[0];
  dubardq4 = 0.;

  if ( cond < -1. )
    {
      fluxm[0] = qright[0]*ubarr*len;
      fluxm[1] = (qright[1]*ubarr + nx*pr)*len;
      fluxm[2] = (qright[2]*ubarr + ny*pr)*len;
      fluxm[3] = (qright[3]+pr)*ubarr*len;

      //Compute the fluxes using analytical derivatives.
      //Do DF/DQ and take this to be DF-/DQ:
      //First row
      dfdqm[0][0] = (ubarr*drhodq1 + rho*dubardq1)*len;
      dfdqm[0][1] = (ubarr*drhodq2 + rho*dubardq2)*len;
      dfdqm[0][2] = (ubarr*drhodq3 + rho*dubardq3)*len;
      dfdqm[0][3] = (ubarr*drhodq4 + rho*dubardq4)*len;
      
      //Second row.
      dfdqm[1][0] = (drhodq1*u*ubarr + dudq1*rho*ubarr + dubardq1*rho*u + dpdq1*nx)*len;
      dfdqm[1][1] = (drhodq2*u*ubarr + dudq2*rho*ubarr + dubardq2*rho*u + dpdq2*nx)*len;
      dfdqm[1][2] = (drhodq3*u*ubarr + dudq3*rho*ubarr + dubardq3*rho*u + dpdq3*nx)*len;
      dfdqm[1][3] = (drhodq4*u*ubarr + dudq4*rho*ubarr + dubardq4*rho*u + dpdq4*nx)*len;
      
      //Third row.
      dfdqm[2][0] = (drhodq1*v*ubarr + dvdq1*rho*ubarr + dubardq1*rho*v + dpdq1*ny)*len;
      dfdqm[2][1] = (drhodq2*v*ubarr + dvdq2*rho*ubarr + dubardq2*rho*v + dpdq2*ny)*len;
      dfdqm[2][2] = (drhodq3*v*ubarr + dvdq3*rho*ubarr + dubardq3*rho*v + dpdq3*ny)*len;
      dfdqm[2][3] = (drhodq4*v*ubarr + dvdq4*rho*ubarr + dubardq4*rho*v + dpdq4*ny)*len;
      
      //Fourth Row.
      dfdqm[3][0] = ( (e + pr)*dubardq1 + (detdq1 + dpdq1)*ubarr )*len;
      dfdqm[3][1] = ( (e + pr)*dubardq2 + (detdq2 + dpdq2)*ubarr )*len;
      dfdqm[3][2] = ( (e + pr)*dubardq3 + (detdq3 + dpdq3)*ubarr )*len;
      dfdqm[3][3] = ( (e + pr)*dubardq4 + (detdq4 + dpdq4)*ubarr )*len;

      flag_m_g = 1;
    }
  else if ( cond > 1. )
    {
      fluxm[0] = 0.;
      fluxm[1] = 0.;
      fluxm[2] = 0.;
      fluxm[3] = 0.;
      
      dfdqm[0][0] = dfdqm[0][1] = dfdqm[0][2] = dfdqm[0][3] = 0.;
      dfdqm[1][0] = dfdqm[1][1] = dfdqm[1][2] = dfdqm[1][3] = 0.;
      dfdqm[2][0] = dfdqm[2][1] = dfdqm[2][2] = dfdqm[2][3] = 0.;
      dfdqm[3][0] = dfdqm[3][1] = dfdqm[3][2] = dfdqm[3][3] = 0.;

      flag_m_l = 1;
    }
  else
    {
      //Calculate the +/- Van Leer flux vectors according to our notes.
      fluxm[0] = -0.25*qright[0]*cr*( (ubarr/cr - 1.0)*(ubarr/cr - 1.0) ) * len;
      fluxm[1] = fluxm[0] * ( (nx/gamma)*(-ubarr - 2.0*cr) + qright[1]/qright[0] );
      fluxm[2] = fluxm[0] * ( (ny/gamma)*(-ubarr - 2.0*cr) + qright[2]/qright[0] );
      fluxm[3] = fluxm[0] * ( ((1.0-gamma)*ubarr*ubarr - 2.0*(gamma-1.0)*ubarr*cr + 2.0*cr*cr)/(gamma*gamma-1.0) 
			      + (qright[1]*qright[1]+qright[2]*qright[2])/(2.*qright[0]*qright[0]) );

     //Now do DF/DQ-:
      //First row
      dfdqm[0][0] = (-0.25*cr*(ubarr/cr - 1.)*(ubarr/cr - 1.)*drhodq1 - 0.25*rho*(ubarr/cr - 1.)*(ubarr/cr - 1.)*dcdq1
		     - 0.5*rho*(ubarr/cr - 1.)*dubardq1 + 0.5*((rho*ubarr)/cr)*(ubarr/cr - 1.)*dcdq1)*len;
      dfdqm[0][1] = (-0.25*cr*(ubarr/cr - 1.)*(ubarr/cr - 1.)*drhodq2 - 0.25*rho*(ubarr/cr - 1.)*(ubarr/cr - 1.)*dcdq2
		     - 0.5*rho*(ubarr/cr - 1.)*dubardq2 + 0.5*((rho*ubarr)/cr)*(ubarr/cr - 1.)*dcdq2)*len;
      dfdqm[0][2] = (-0.25*cr*(ubarr/cr - 1.)*(ubarr/cr - 1.)*drhodq3 - 0.25*rho*(ubarr/cr - 1.)*(ubarr/cr - 1.)*dcdq3
		     - 0.5*rho*(ubarr/cr - 1.)*dubardq3 + 0.5*((rho*ubarr)/cr)*(ubarr/cr - 1.)*dcdq3)*len;
      dfdqm[0][3] = (-0.25*cr*(ubarr/cr - 1.)*(ubarr/cr - 1.)*drhodq4 - 0.25*rho*(ubarr/cr - 1.)*(ubarr/cr - 1.)*dcdq4
		     - 0.5*rho*(ubarr/cr - 1.)*dubardq4 + 0.5*((rho*ubarr)/cr)*(ubarr/cr - 1.)*dcdq4)*len;
      
      //Second row.
      dfdqm[1][0] = dfdqm[0][0]*( (nx/gamma)*(-ubarr - 2.0*cr) + u ) + fluxm[0]*( (nx/gamma)*(-dubardq1 - 2.0*dcdq1) +dudq1 );
      dfdqm[1][1] = dfdqm[0][1]*( (nx/gamma)*(-ubarr - 2.0*cr) + u ) + fluxm[0]*( (nx/gamma)*(-dubardq2 - 2.0*dcdq2) +dudq2 );
      dfdqm[1][2] = dfdqm[0][2]*( (nx/gamma)*(-ubarr - 2.0*cr) + u ) + fluxm[0]*( (nx/gamma)*(-dubardq3 - 2.0*dcdq3) +dudq3 );
      dfdqm[1][3] = dfdqm[0][3]*( (nx/gamma)*(-ubarr - 2.0*cr) + u ) + fluxm[0]*( (nx/gamma)*(-dubardq4 - 2.0*dcdq4) +dudq4 );
      
      //Third row.
      dfdqm[2][0] = dfdqm[0][0]*( (ny/gamma)*(-ubarr - 2.0*cr) + v ) + fluxm[0]*( (ny/gamma)*(-dubardq1 - 2.0*dcdq1) +dvdq1 );
      dfdqm[2][1] = dfdqm[0][1]*( (ny/gamma)*(-ubarr - 2.0*cr) + v ) + fluxm[0]*( (ny/gamma)*(-dubardq2 - 2.0*dcdq2) +dvdq2 );
      dfdqm[2][2] = dfdqm[0][2]*( (ny/gamma)*(-ubarr - 2.0*cr) + v ) + fluxm[0]*( (ny/gamma)*(-dubardq3 - 2.0*dcdq3) +dvdq3 );
      dfdqm[2][3] = dfdqm[0][3]*( (ny/gamma)*(-ubarr - 2.0*cr) + v ) + fluxm[0]*( (ny/gamma)*(-dubardq4 - 2.0*dcdq4) +dvdq4 );
      
      //Fourth Row.
      dfdqm[3][0] = dfdqm[0][0]*( ((1.0-gamma)*ubarr*ubarr - 2.0*(gamma-1.0)*ubarr*cr + 2.0*cr*cr)/(gamma*gamma-1.0) + (u*u + v*v)/2.0 )
	+ fluxm[0]*( ((1.0-gamma)*2.*ubarr*dubardq1 - 2.0*(gamma-1.0)*(dubardq1*cr + dcdq1*ubarr) + 4.0*cr*dcdq1)/(gamma*gamma-1.0) 
		     + (u*dudq1 + v*dvdq1) );
      dfdqm[3][1] = dfdqm[0][1]*( ((1.0-gamma)*ubarr*ubarr - 2.0*(gamma-1.0)*ubarr*cr + 2.0*cr*cr)/(gamma*gamma-1.0) + (u*u + v*v)/2.0 )
	+ fluxm[0]*( ((1.0-gamma)*2.*ubarr*dubardq2 - 2.0*(gamma-1.0)*(dubardq2*cr + dcdq2*ubarr) + 4.0*cr*dcdq2)/(gamma*gamma-1.0) 
		     + (u*dudq2 + v*dvdq2) );
      dfdqm[3][2] = dfdqm[0][2]*( ((1.0-gamma)*ubarr*ubarr - 2.0*(gamma-1.0)*ubarr*cr + 2.0*cr*cr)/(gamma*gamma-1.0) + (u*u + v*v)/2.0 )
	+ fluxm[0]*( ((1.0-gamma)*2.*ubarr*dubardq3 - 2.0*(gamma-1.0)*(dubardq3*cr + dcdq3*ubarr) + 4.0*cr*dcdq3)/(gamma*gamma-1.0) 
		     + (u*dudq3 + v*dvdq3) );
      dfdqm[3][3] = dfdqm[0][3]*( ((1.0-gamma)*ubarr*ubarr - 2.0*(gamma-1.0)*ubarr*cr + 2.0*cr*cr)/(gamma*gamma-1.0) + (u*u + v*v)/2.0 )
	+ fluxm[0]*( ((1.0-gamma)*2.*ubarr*dubardq4 - 2.0*(gamma-1.0)*(dubardq4*cr + dcdq4*ubarr) + 4.0*cr*dcdq4)/(gamma*gamma-1.0) 
		     + (u*dudq4 + v*dvdq4) );
    }

  flux[0] = fluxp[0] + fluxm[0];
  flux[1] = fluxp[1] + fluxm[1];
  flux[2] = fluxp[2] + fluxm[2];
  flux[3] = fluxp[3] + fluxm[3];

  if ( flag_m_g == 1 && flag_p_g == 1 )
    {
      printf("  Flow is Greater than 1.0 at the face in both directions.\n");
      exit(0);
    }

  if ( flag_p_l == 1 && flag_m_l == 1 )
    {
      printf("  Flow is Less than 1.0 at the face in both directions.\n");
      exit(0);
    }

  return 1;
}
