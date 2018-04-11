//=============================================================
// 
//  curved_boundaries.C
//  
//  Functions needed to handle Gaussian quadrature on curved
//  boundaries.
//
//  Written by - Shane Sawyer
//
//=============================================================

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <gsl/gsl_integration.h>  // HACK - i think this is for the installed library. I'm using a static linked library I compile locally now.
#include "gsl/gsl_integration.h"
#include "gsl/gsl_errno.h"
#include "grid.h"
#include "mesh_utils.h"
#include "Defined_Variables.h"
#include "linear_algebra.h"
#include "params.h"
#include "curved_boundaries.h"


// EVERYTHING IS DEFINED WITH OUTWARD POINTING NORMALS.

// LIST OF RESERVED BC's
/*
     0 - Farfield
     1 - Inviscid
     11 - Unit circle @ origin Inviscid
     12 - Unit cirlce @ origin Farfield
     13 - R=2 circle @ origin Inviscid
     14 - R=2 circle @ origin Farfield
     15 - R=3 circle @ origin Inviscid
     16 - R=3 circle @ origin Farfield
     17 - NACA 0012 upper surface Inviscid
     18 - NACA 0012 upper surface Farfield
     19 - NACA 0012 lower surface Inviscid
     20 - NACA 0012 lower surface Farfield
     21 - R=300 circle @ origin Inviscid
     22 - R=300 circle @ origin Farfield

*/

double naca_function ( double x , int bc )
{
  double y = 0.0;
  y = 0.6 * ( 0.2969 * sqrt(x) - 0.1260*x - 0.3516*x*x + 0.2843*x*x*x - 0.1036*x*x*x*x );

  if ( bc == 19 || bc == 20 )
    y = -1.0*y;

  return y;
}

double naca_arclen_function ( double x, void *params )
{
  double y = 0.0;
  y = sqrt( 1.0 + ( 0.08907/sqrt(x) - 0.07560 - 0.42192*x + 0.51174*x*x - 0.24864*x*x*x ) *
                  ( 0.08907/sqrt(x) - 0.07560 - 0.42192*x + 0.51174*x*x - 0.24864*x*x*x ) );

  return y;
}

double naca_arclen_function2 ( double x )
{
  double y = 0.0;
  y = -1.0 * sqrt( 1.0 + ( 0.08907/sqrt(x) - 0.07560 - 0.42192*x + 0.51174*x*x - 0.24864*x*x*x ) *
		         ( 0.08907/sqrt(x) - 0.07560 - 0.42192*x + 0.51174*x*x - 0.24864*x*x*x ) );

  return y;
}

//=============================================================
// 
//  curved_boundary_midpoint()
//
//  Returns the point that lies equidistant from two specified
//  points (which belong to a segment).
//  
//  GRID *grid;                       // The grid.
//  int bc;                           // Determine which shape to apply.
//  double *xL;
//  double *xR;
//  double *xM;
//
//=============================================================

void curved_boundary_midpoint ( GRID *grid, int bc, double *xL, double *xR, double *xM )
{
  double tL,tR,tM;                      // Parameterized coordinates.

  // needed for naca;
  double xn,xnp1,f_xn,fp_xn;
  double Lt, L;                         // lengths
  int MAX_ITERS = 50;
  int i;

  //gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc( 1000 );
  gsl_integration_workspace *work_ptr = NULL;
  gsl_function my_function;
  double error, result;
  double abserror = 1.0e-8;
  double relerror = 1.0e-8;
  double alpha = 1.0;
  double xa,xb;
  void *params_ptr = &alpha;
  my_function.function = &naca_arclen_function;
  my_function.params = params_ptr;

  gsl_set_error_handler_off();
  
  if ( bc == 11 || bc == 12 )  // Half Unit circle at the origin.
    {
      // Do some checks to make sure I get the right angle.
      if ( xL[0] >= 0. && xL[1] >= 0. )
	{
	  tL = acos(xL[0]);
	}
      else if ( xL[0] >= 0.  && xL[1] < 0. )
	{
	  tL = -acos(xL[0]);
	}
      else if ( xL[0] < 0. && xL[1] >= 0. )
	{
	  tL = acos(xL[0]);
	}
      else // -x , -y
	{
	  tL = -acos(xL[0]);
	}

      if ( xR[0] >= 0. && xR[1] >= 0. )
	{
	  tR = acos(xR[0]);
	}
      else if ( xR[0] >= 0.  && xR[1] < 0. )
	{
	  tR = -acos(xR[0]);
	}
      else if ( xR[0] < 0. && xR[1] >= 0. )
	{
	  tR = acos(xR[0]);
	}
      else // -x , -y
	{
	  tR = -acos(xR[0]);
	}

      // Here comes a major hack. For this boundary, I'm going to assume that no edges cross the x axis from nodeL to nodeR.
      // This is easy to ensure from a grid generation perspective. Just mesh half the circle, then mirror about the plane y==0

      /*
      if ( th1 < 0. )
	{
	  thM = 0.5*(th1 - th2);
	}
      else if ( th2 < 0. )
	{
	  thM = 0.5*( th2 - th1);
	}
      */

      if ( tL < 0. && tR > 0. )
	tR = -1.*tR;
      else if ( tL > 0. && tR < 0. )
	tL = -1.*tL;

      tM = 0.5*(tL+tR);

      xM[0] = cos(tM);
      xM[1] = sin(tM);

      return;
    }

  if ( bc == 13 || bc == 14 )
    {
      tL = atan2( xL[1] , xL[0] );

      tR = atan2( xR[1] , xR[0] );

      if ( tL < 0. && tR > 0. )
	tR = -1.*tR;
      else if ( tL > 0. && tR < 0. )
	tL = -1.*tL;

      tM = 0.5*(tL+tR);

      xM[0] = 2.*cos(tM);
      xM[1] = 2.*sin(tM);
      
      return;
    }

  if ( bc == 15 || bc == 16 )
    {
      tL = atan2( xL[1] , xL[0] );

      tR = atan2( xR[1] , xR[0] );

      if ( tL < 0. && tR > 0. )
	tR = -1.*tR;
      else if ( tL > 0. && tR < 0. )
	tL = -1.*tL;

      tM = 0.5*(tL+tR);

      xM[0] = 3.*cos(tM);
      xM[1] = 3.*sin(tM);
      
      return;
    }

  if ( bc == 17 || bc == 18 || bc == 19 || bc == 20 )
    {
      // Initial guess is the half way point across the linear segment joining the two endpoints.
      xM[0] = 0.5*(xL[0] + xR[0]);
      xM[1] = naca_function( xM[0] , bc );

      return;

      if ( xL[0] < xR[0] )
	{
	  xa = xL[0];
	  xb = xR[0];
	}
      else
	{
	  xa = xR[0];
	  xb = xL[0];
	}

      // get the total length.

      work_ptr = gsl_integration_workspace_alloc( 1000 );

      gsl_integration_qags ( &my_function, xa , xb , abserror, relerror, 1000, work_ptr, &result, &error );

      Lt = result;

      // now we start Newton's method to get the midpoint of the curve on the interval [a,b].

      xn = xM[0];

      //printf("Beginning the process:\n");
      //printf(" xa = %.15E , xb = %.15E , x0 = %.15E\n",xa,xb,xn);
      //printf(" Total length is %.15E\n\n",Lt);
      //fflush(stdout);

      for ( i=1; i <= MAX_ITERS; i++ )
	{
	  // compute the current arclength of the guess.
	  gsl_integration_qags ( &my_function, xa , xn , abserror, relerror, 1000, work_ptr, &result, &error );
	  L = result;

	  // compute the functional.
	  f_xn = 0.5*Lt - L;

	  // compute the derivative value.
	  fp_xn = naca_arclen_function2 ( xn );

	  // compute the next iterate.
	  xnp1 = xn - ( f_xn / fp_xn );

	  //printf("Iteration %d:\n",i);
	  //printf("  Current length is L = %.15E\n",L);
	  //printf("  f_xn = 1/2 Lt - L = %.15E\n",f_xn);
	  //printf("  f prime_xn = %.15E\n",fp_xn);
	  //printf("  Next guess: x_n+1 = %.15E\n",xnp1);
	  //fflush(stdout);

	  if ( fabs( f_xn ) < 1.0E-12 )
	    break;

	  xn = xnp1;
	}

      if ( i > MAX_ITERS )
	{
	  printf("FATAL ERROR: In curved_boundary_midpoint(): Newton's method failed to find the midpoint in the specified iterations.\n");
	  printf("  xL = %f , yL = %f ; xR = %f, yR = %f\n",(float)xL[0],(float)xL[1],(float)xR[0],(float)xR[1]);
	  printf("  The current estimate for the midpoint is %.15E , with an error of %.15E\n",xnp1, f_xn);
	  fflush(stdout);

	  gsl_integration_workspace_free(work_ptr);
	  exit(1);
	}
      
      gsl_integration_workspace_free(work_ptr);

      xM[0] = xnp1;
      xM[1] = naca_function( xM[0] , bc );
      
      //if ( bc == 17 || bc == 18 )
      //xM[1] = naca_function( xM[0] , bc );
      //else
      //xM[1] = -1.0 * naca_function( xM[0] , bc );

      //printf("------------------------------------------------------------------------------------------------\n\n");
      //fflush(stdout);
      
      return;
    }

  if ( bc == 21 || bc == 22 )
    {
      tL = atan2( xL[1] , xL[0] );
      
      tR = atan2( xR[1] , xR[0] );

      if ( tL < 0. && tR > 0. )
	tR = -1.*tR;
      else if ( tL > 0. && tR < 0. )
	tL = -1.*tL;
      
      tM = 0.5*(tL+tR);
      
      xM[0] = 300.*cos(tM);
      xM[1] = 300.*sin(tM);
      
      return;
    }

  return;
}

//=============================================================
// 
//  curved_boundary_gauss_points()
//
//  Given two points on a curved boundary, returns the Gaussian
//  quadrature points in x,y space.
//  
//  GRID *grid;                       // The grid.
//  int bc;                           // Determine which shape to apply.
//  double *xL;
//  double *xR;
//  double *GP;
//
//=============================================================

void curved_boundary_gauss_points ( GRID *grid, int bc, double *xL, double *xR, double *GP )
{
  double t1,t2,t3;                      // Gauss points in parameter space.
  double thL,thR;                       // Parameterized coordinates.
  double th1,th2,th3;
  double x1,x2,x3;
  double y1,y2,y3;

  // Set the gauss points coordinates'.
  t1 = 0.;
  t2 = sqrt(0.6);
  t3 = -t2;

  // needed for naca;
  double xn,xnp1,f_xn,fp_xn;
  double Lt, L;                         // lengths
  double beta[3];
  double t[3];
  int MAX_ITERS = 50;
  int i,j;
  double update;

  //gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc( 1000 );
  gsl_integration_workspace *work_ptr = NULL;
  gsl_function my_function;
  double error, result;
  double abserror = 1.0e-8;
  double relerror = 1.0e-8;
  double alpha = 1.0;
  double xa,xb;
  void *params_ptr = &alpha;
  my_function.function = &naca_arclen_function;
  my_function.params = params_ptr;

  //beta[0] = ( 1.0 - sqrt(3./5.))*0.5;
  //beta[1] = 0.5;
  //beta[2] = ( 1.0 + sqrt(3./5.))*0.5;

  beta[0] = 0.5;
  beta[1] = ( 1.0 + sqrt(3./5.))*0.5;
  beta[2] = ( 1.0 - sqrt(3./5.))*0.5;

  t[0] = 0.0;
  t[1] = sqrt(0.6);
  t[2] = -1.0*sqrt(0.6);
    
  if ( bc == 11 || bc == 12 )  // Unit circle at the origin.
    {
      // Get the parameterized coordinates of the points.

      // Do some checks to make sure I get the right angle.
      if ( xL[0] >= 0. && xL[1] >= 0. )
	{
	  thL = acos(xL[0]);
	}
      else if ( xL[0] >= 0.  && xL[1] < 0. )
	{
	  thL = -acos(xL[0]);
	}
      else if ( xL[0] < 0. && xL[1] >= 0. )
	{
	  thL = acos(xL[0]);
	}
      else // -x , -y
	{
	  thL = -acos(xL[0]);
	}

      if ( xR[0] >= 0. && xR[1] >= 0. )
	{
	  thR = acos(xR[0]);
	}
      else if ( xR[0] >= 0.  && xR[1] < 0. )
	{
	  thR = -acos(xR[0]);
	}
      else if ( xR[0] < 0. && xR[1] >= 0. )
	{
	  thR = acos(xR[0]);
	}
      else // -x , -y
	{
	  thR = -acos(xR[0]);
	}

      if ( thL < 0. && thR > 0. )
	thR = -1.*thR;
      else if ( thL > 0. && thR < 0. )
	thL = -1.*thL;
      
      // Gauss points in parameter spcae part 2.
      th1 = (thL+thR)*0.5 + (thR-thL)*0.5*t1;
      th2 = (thL+thR)*0.5 + (thR-thL)*0.5*t2;
      th3 = (thL+thR)*0.5 + (thR-thL)*0.5*t3;

      // Set the x,y coordinates.
      GP[0] = cos(th1);
      GP[1] = sin(th1);

      GP[2] = cos(th2);
      GP[3] = sin(th2);

      GP[4] = cos(th3);
      GP[5] = sin(th3);

      return;
    }

  if ( bc == 13 || bc == 14 )
    {
      thL = atan2( xL[1] , xL[0] );
      
      thR = atan2( xR[1] , xR[0] );
      
      if ( thL < 0. && thR > 0. )
	thR = -1.*thR;
      else if ( thL > 0. && thR < 0. )
	thL = -1.*thL;
      
      // Gauss points in parameter spcae part 2.
      th1 = (thL+thR)*0.5 + (thR-thL)*0.5*t1;
      th2 = (thL+thR)*0.5 + (thR-thL)*0.5*t2;
      th3 = (thL+thR)*0.5 + (thR-thL)*0.5*t3;
      
      // Set the x,y coordinates.
      GP[0] = 2.*cos(th1);
      GP[1] = 2.*sin(th1);

      GP[2] = 2.*cos(th2);
      GP[3] = 2.*sin(th2);

      GP[4] = 2.*cos(th3);
      GP[5] = 2.*sin(th3);
      
      return;
    }

  if ( bc == 15 || bc == 16 )
    {
      thL = atan2( xL[1] , xL[0] );
      
      thR = atan2( xR[1] , xR[0] );
      
      if ( thL < 0. && thR > 0. )
	thR = -1.*thR;
      else if ( thL > 0. && thR < 0. )
	thL = -1.*thL;
      
      // Gauss points in parameter spcae part 2.
      th1 = (thL+thR)*0.5 + (thR-thL)*0.5*t1;
      th2 = (thL+thR)*0.5 + (thR-thL)*0.5*t2;
      th3 = (thL+thR)*0.5 + (thR-thL)*0.5*t3;
      
      // Set the x,y coordinates.
      GP[0] = 3.*cos(th1);
      GP[1] = 3.*sin(th1);

      GP[2] = 3.*cos(th2);
      GP[3] = 3.*sin(th2);

      GP[4] = 3.*cos(th3);
      GP[5] = 3.*sin(th3);
      
      return;
    }

  if ( bc == 17 || bc == 18 || bc == 19 || bc == 20 )
    {
      /*
      x1 = (xL[0] + xR[0])*0.5 + (xR[0] - xL[0])*0.5*t1;
      x2 = (xL[0] + xR[0])*0.5 + (xR[0] - xL[0])*0.5*t2;
      x3 = (xL[0] + xR[0])*0.5 + (xR[0] - xL[0])*0.5*t3;
      
      y1 = naca_function( x1 , bc );
      y2 = naca_function( x2 , bc );
      y3 = naca_function( x3 , bc );

      // Set the x,y coordinates.
      GP[0] = x1;
      GP[1] = y1;

      GP[2] = x2;
      GP[3] = y2;

      GP[4] = x3;
      GP[5] = y3;
      */

      if ( xL[0] < xR[0] )
	{
	  xa = xL[0];
	  xb = xR[0];
	}
      else
	{
	  xa = xR[0];
	  xb = xL[0];
	}

      //printf("Beginning the process of finding gauss points.:\n");
      //printf(" xa = %.15E , xb = %.15E\n",xa,xb);
      //fflush(stdout);
      
      // loop over the needed number of Gauss points and use Newton's method to get the correct location along the curve.

      for ( i=0; i < 3; i++ )
	{
	  //printf("    Working on gauss node %d.\n",i+1);
	  //fflush(stdout);

	  work_ptr = gsl_integration_workspace_alloc( 1000 );

	  gsl_integration_qags ( &my_function, xa , xb , abserror, relerror, 1000, work_ptr, &result, &error );
	  
	  Lt = result;

	  // get an initial estimate.
	  //xn = beta[i] * ( xa + xb );
	  xn = (xa + xb)*0.5 + (xb - xa)*0.5*t[i];

	  //printf("    Total length is %.15E  and my initial guess is %.15E\n",Lt,xn);
	  //printf("    Desired length is %.15E\n",beta[i]*Lt);
	  //fflush(stdout);

	  for ( j=1; j <= MAX_ITERS; j++ )
	    {
	      // compute the current arclength of the guess.
	      gsl_integration_qags ( &my_function, xa , xn , abserror, relerror, 1000, work_ptr, &result, &error );
	      L = result;
	      
	      // compute the functional.
	      f_xn = (beta[i] * Lt) - L;
	      
	      // compute the derivative value.
	      fp_xn = naca_arclen_function2 ( xn );
	      
	      // compute the next iterate.
	      xnp1 = xn - ( f_xn / fp_xn );

	      update = ( f_xn / fp_xn );

	      //printf("Iteration %d:\n",j);
	      //printf("  Current length is L = %.15E\n",L);
	      //printf("  f_xn = 1/2 Lt - L = %.15E\n",f_xn);
	      //printf("  f prime_xn = %.15E\n",fp_xn);
	      //printf("  update is %.15E\n",update);
	      //printf("  Next guess: x_n+1 = %.15E\n",xnp1);
	      //fflush(stdout);

	      while ( xnp1 < xa || xnp1 > xb )
		{
		  //printf("      Update was too severe. Reducing by half from %.15E  to  %.15E\n",update,0.5*update);
		  //fflush(stdout);

		  update = 0.5*update;
		 
		  xnp1 = xn - update;

		  //printf("     New x_n+1 = %.15E\n",xnp1);
		  //fflush(stdout);

		}
	      
	      if ( fabs( f_xn ) < 1.0E-12 )
		break;

	      xn = xnp1;
	    }

	  if ( j > MAX_ITERS )
	    {
	      printf("FATAL ERROR: In curved_boundary_gauss_points(): Newton's method failed to find the point (%d) in the specified iterations.\n",i+1);
	      printf("  xL = %f , yL = %f ; xR = %f, yR = %f\n",(float)xL[0],(float)xL[1],(float)xR[0],(float)xR[1]);
	      printf("  The current estimate for the point is %.15E , with an error of %.15E\n",xnp1, f_xn);
	      fflush(stdout);
	      
	      gsl_integration_workspace_free(work_ptr);
	      exit(1);
	    }

	  
	  gsl_integration_workspace_free(work_ptr);

	  GP[i*2+0] = xnp1;
	  GP[i*2+1] = naca_function( xnp1 , bc );

	  //if ( bc == 19 || bc == 20 )
	    //GP[i*2+1] = -1.0 * GP[i*2+1];

	    //printf("------------------------------------------------------------------------------------------------\n\n");
	    //fflush(stdout);
	}

      // Now if I flipped the nodes on the interval so that xa < xb, then the Gaussian nodes will need to flipped for the case when t= +- sqrt(3/5).
      //if ( xL[0] > xR[0] )  // flip i=1,2
	//{
      //xa = GP[2];
      //  xb = GP[3];

	//  GP[2] = GP[4];
      // GP[3] = GP[5];

      //  GP[4] = xa;
      //  GP[5] = xb;
      //}

      return;
    }

  if ( bc == 21 || bc == 22 )
    {
      thL = atan2( xL[1] , xL[0] );
      
      thR = atan2( xR[1] , xR[0] );
      
      if ( thL < 0. && thR > 0. )
	thR = -1.*thR;
      else if ( thL > 0. && thR < 0. )
	thL = -1.*thL;
      
      // Gauss points in parameter spcae part 2.
      th1 = (thL+thR)*0.5 + (thR-thL)*0.5*t1;
      th2 = (thL+thR)*0.5 + (thR-thL)*0.5*t2;
      th3 = (thL+thR)*0.5 + (thR-thL)*0.5*t3;
      
      // Set the x,y coordinates.
      GP[0] = 300.*cos(th1);
      GP[1] = 300.*sin(th1);

      GP[2] = 300.*cos(th2);
      GP[3] = 300.*sin(th2);

      GP[4] = 300.*cos(th3);
      GP[5] = 300.*sin(th3);
      
      return;
    }

  return;
}

//=============================================================
// 
//  curved_boundary_normal_vector()
//
//  Given a point on a curved boundary, returns the unit normal
//  vector.
//  
//  GRID *grid;                       // The grid.
//  int bc;                           // Determine which shape to apply.
//  double *x;
//  double *nx;
//  double *ny;
//
//=============================================================

void curved_boundary_normal_vector ( GRID *grid, int bc, double *x, double *nx, double *ny )
{
  double t;                               // x in parameterized coordinate space.
  double y,m,len;
    
  if ( bc == 11 || bc == 12 )  // Unit circle at the origin.
    {
      // Get the parameterized coordinate of the point.
      // Do some checks to make sure I get the right angle.
      if ( x[0] >= 0. && x[1] >= 0. )
	{
	  t = acos(x[0]);
	}
      else if ( x[0] >= 0.  && x[1] < 0. )
	{
	  t = -acos(x[0]);
	}
      else if ( x[0] < 0. && x[1] >= 0. )
	{
	  t = acos(x[0]);
	}
      else // -x , -y
	{
	  t = -acos(x[0]);
	}

      if ( bc == 11 )
	{
	  *nx = -cos(t);
	  *ny = -sin(t);
	}
      else  // If the unit circle is on the farfield, this will get the right normal vector direction.
	{
	  *nx = cos(t);
	  *ny = sin(t);
	}
      
      return;
    }

  if ( bc == 13 || bc == 14 )
    {
      t = atan2( x[1] , x[0] );

      if ( bc == 13 )
	{
	  *nx = -cos(t);
	  *ny = -sin(t);
	}
      else
	{
	  //printf("FATAL ERROR: The lower portion of the annulus (Ri) was specified as farfield!\n");
	  //exit(1);

	  *nx = -cos(t);
	  *ny = -sin(t);
	}

      return;
    }

  if ( bc == 15 || bc == 16 )
    {
      t = atan2( x[1] , x[0] );

      if ( bc == 15 )
	{
	  *nx = cos(t);
	  *ny = sin(t);
	}
      else
	{
	  //printf("FATAL ERROR: The upper portion of the annulus (Ro) was specified as farfield!\n");
	  //exit(1);

	  *nx = cos(t);
	  *ny = sin(t);
	}

      return;
    }

  if ( bc == 17 || bc == 18 )
    {
      if ( fabs( x[0] ) < 1.e-10 ) // at the origin.
	{
	  *nx = 1.0;
	  *ny = 0.0;
	}
      else
	{
	  y = 0.08907/sqrt(x[0]) - 0.07560 - 0.42192*x[0] + 0.51174*x[0]*x[0] - 0.24864*x[0]*x[0]*x[0];  // get the tangent line.
	  
	  m = - 1.0 / y;  // get an orthogonal line
	  
	  len = sqrt( 1.0 + m*m );  // get the length

	  *nx = 1.0/len;  // normalize
	  *ny = m / len;
	  
	  if ( *ny > 0.0 )  // flip if the normal points the wrong way.
	    {
	      *nx = -1.0* (*nx);
	      *ny = -1.0* (*ny);
	    }
	} 

      return;
    }

  if ( bc == 19 || bc == 20 )
    {
      if ( fabs( x[0] ) < 1.e-10 ) // at the origin.
	{
	  *nx = 1.0;
	  *ny = 0.0;
	}
      else
	{
	  y = -0.08907/sqrt(x[0]) + 0.07560 + 0.42192*x[0] - 0.51174*x[0]*x[0] + 0.24864*x[0]*x[0]*x[0];  // get the tangent line.
	  
	  m = - 1.0 / y;  // get an orthogonal line
	  
	  len = sqrt( 1.0 + m*m );  // get the length
	  
	  *nx = 1.0/len;  // normalize
	  *ny = m / len;
	  
	  if ( *ny < 0.0 )  // flip if the normal points the wrong way.
	    {
	      *nx = -1.0* (*nx);
	      *ny = -1.0* (*ny);
	    }
	} 

      return;
    }

  if ( bc == 21 || bc == 22 )
    {
      t = atan2( x[1] , x[0] );
      
      if ( bc == 21 )
	{
	  *nx = cos(t);
	  *ny = sin(t);
	}
      else
	{
	  //printf("FATAL ERROR: The upper portion of the annulus (Ro) was specified as farfield!\n");
	  //exit(1);

	  *nx = cos(t);
	  *ny = sin(t);
	}

      return;
    }

  return;
}


//=============================================================
// 
//  curved_boundary_full_normal_vector()
//
//  Given a point on a curved boundary, returns the normal
//  vector.
//  
//  GRID *grid;                       // The grid.
//  int bc;                           // Determine which shape to apply.
//  double *x;
//  double *nx;
//  double *ny;
//
//=============================================================

void curved_boundary_full_normal_vector ( GRID *grid, int bc, double *x, double *nx, double *ny )
{
  // NOT NORMALIZED.

  double t;                               // x in parameterized coordinate space.
    
  if ( bc == 11 || bc == 12 )  // Unit circle at the origin.
    {
      // Get the parameterized coordinate of the point.
      // Do some checks to make sure I get the right angle.
      if ( x[0] >= 0. && x[1] >= 0. )
	{
	  t = acos(x[0]);
	}
      else if ( x[0] >= 0.  && x[1] < 0. )
	{
	  t = -acos(x[0]);
	}
      else if ( x[0] < 0. && x[1] >= 0. )
	{
	  t = acos(x[0]);
	}
      else // -x , -y
	{
	  t = -acos(x[0]);
	}

      if ( bc == 11 )
	{
	  *nx = -cos(t);
	  *ny = -sin(t);
	}
      else // If the unit circle is on the farfield, this will get the right normal vector direction.
	{
	  *nx = cos(t);
	  *ny = sin(t);
	}
      
      return;
    }

  if ( bc == 13 || bc == 14 )
    {
      t = atan2( x[1] , x[0] );

      if ( bc == 13 )
	{
	  *nx = -cos(t);
	  *ny = -sin(t);
	}
      else
	{
	  //printf("FATAL ERROR: The lower portion of the annulus (Ri) was specified as farfield!\n");
	  //exit(1);

	  *nx = -cos(t);
	  *ny = -sin(t);
	}

      return;
    }

  if ( bc == 15 || bc == 16 )
    {
      t = atan2( x[1] , x[0] );

      if ( bc == 15 )
	{
	  *nx = cos(t);
	  *ny = sin(t);
	}
      else
	{
	  //printf("FATAL ERROR: The upper portion of the annulus (Ro) was specified as farfield!\n");
	  //exit(1);

	  *nx = cos(t);
	  *ny = sin(t);
	}

      return;
    }

  if ( bc == 17 || bc == 18 || bc == 19 || bc == 20 )
    {
      printf("FATAL ERROR: This action is not supported by the NACA 0012 boundary condition: curved_boundary_full_normal_vector().\n");
      fflush(stdout);
      exit(0);
    }

  if ( bc == 21 || bc == 22 )
    {
      t = atan2( x[1] , x[0] );

      if ( bc == 21 )
	{
	  *nx = cos(t);
	  *ny = sin(t);
	}
      else
	{
	  //printf("FATAL ERROR: The upper portion of the annulus (Ro) was specified as farfield!\n");
	  //exit(1);

	  *nx = cos(t);
	  *ny = sin(t);
	}

      return;
    }

  return;
}


//=============================================================
// 
//  curved_boundary_arclength()
//
//  Given two points on a curved boundary, return the distance
//  between them as a function of arclength.
//  
//  GRID *grid;                       // The grid.
//  int bc;                           // Determine which shape to apply.
//  double *xL;
//  double *xR;
//  double *arclength;
//
//=============================================================

void curved_boundary_arclength ( GRID *grid, int bc, double *xL, double *xR, double *arclength )
{
  double tL, tR;                            // x in parameterized coordinate space.

  if ( bc == 11 || bc == 12 )  // Unit circle at the origin.
    {
      // Get the parameterized coordinate of the point.

      // Do some checks to make sure I get the right angle.
      if ( xL[0] >= 0. && xL[1] >= 0. )
	{
	  tL = acos(xL[0]);
	}
      else if ( xL[0] >= 0.  && xL[1] < 0. )
	{
	  tL = -acos(xL[0]);
	}
      else if ( xL[0] < 0. && xL[1] >= 0. )
	{
	  tL = acos(xL[0]);
	}
      else // -x , -y
	{
	  tL = -acos(xL[0]);
	}

      if ( xR[0] >= 0. && xR[1] >= 0. )
	{
	  tR = acos(xR[0]);
	}
      else if ( xR[0] >= 0.  && xR[1] < 0. )
	{
	  tR = -acos(xR[0]);
	}
      else if ( xR[0] < 0. && xR[1] >= 0. )
	{
	  tR = acos(xR[0]);
	}
      else // -x , -y
	{
	  tR = -acos(xR[0]);
	}
      
      if ( tL < 0. && tR > 0. )
	tR = -1.*tR;
      else if ( tL > 0. && tR < 0. )
	tL = -1.*tL;

      *arclength = fabs( tR - tL );
      
      return;
    }

  if ( bc == 13 || bc == 14 )
    {
      tL = atan2( xL[1] , xL[0] );

      tR = atan2( xR[1] , xR[0] );

      /*
      if ( tL < 0. && tR > 0. )
	tR = -1.*tR;
      else if ( tL > 0. && tR < 0. )
	tL = -1.*tL;
      */

      *arclength = fabs( tR - tL ) * 2.;
      
      return;
    }

  if ( bc == 15 || bc == 16 )
    {
      tL = atan2( xL[1] , xL[0] );

      tR = atan2( xR[1] , xR[0] );

      /*
      if ( tL < 0. && tR > 0. )
	tR = -1.*tR;
      else if ( tL > 0. && tR < 0. )
	tL = -1.*tL;
      */

      *arclength = fabs( tR - tL ) * 3.;
      
      return;
    }

  if ( bc == 17 || bc == 18 || bc == 19 || bc == 20 )
    {
      // Variables needed for the naca 0012 bc.
      gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc( 1000 );
      gsl_function my_function;
      double error, result;
      double abserror = 1.0e-8;
      double relerror = 1.0e-8;
      double alpha = 1.0;
      double xa,xb;
      void *params_ptr = &alpha;
      my_function.function = &naca_arclen_function;
      my_function.params = params_ptr;

      if ( xL[0] > xR[0] )
	{
	  xa = xR[0];
	  xb = xL[0];
	}
      else
	{
	  xa = xL[0];
	  xb = xR[0];
	}

      gsl_integration_qags ( &my_function, xa , xb , abserror, relerror, 1000, work_ptr, &result, &error );

      *arclength = result;

      gsl_integration_workspace_free(work_ptr);

      if ( result < 0.0 )
	{
	  printf("FATAL ERROR: The arclength along the NACA 0012 airfoil is negative.\n");
	  fflush(stdout);
	  exit(0);
	}
      return;
    }

  if ( bc == 21 || bc == 22 )
    {
      tL = atan2( xL[1] , xL[0] );
      
      tR = atan2( xR[1] , xR[0] );

      if ( tL < 0. && tR > 0. )
	tR = -1.*tR;
      else if ( tL > 0. && tR < 0. )
	tL = -1.*tL;

      *arclength = fabs( tR - tL ) * 300.;
      
      return;
    }  

  return;
}
