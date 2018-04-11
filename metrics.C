//=============================================================
// 
//  metrics.cpp
//  
//  Function to calculate mesh metrics.
//
//  Written by - Shane Sawyer
//
//=============================================================

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "grid.h"
#include "mesh_utils.h"
#include "metrics.h"
#include "Defined_Variables.h"


//=============================================================
// 
//  metrics2D()
//
//  Compute the area of all triangles and quads to ensure positive
//  areas. Check that all quad corner jacobians are positive.
// 
//  GRID *grid                        // The grid.
//
//=============================================================

void metrics2D( GRID *grid )
{
  int i,t;                                    // Loop counters.
  int k=0, m=0;                               // Negative area element counters.
  int n0, n1, n2, n3;                         // Element nodes.
  double area;                                // Element area.
  double mag;                                 // Vector magnitude.
  double j;                                   // Jacobian.
  double javg;                                // Element average jacobian.
  int flag;                                   // Counter for element negative corner jacobians.
  int jstat[4];                               // Jacobian Statistics for the mesh.
                                              //   [0] -> Negative.
                                              //   [1] -> Negative Skewed.
                                              //   [2] -> Positive Skewed.
                                              //   [3] -> Positive.
  double v1[NDIM];                            // Vectors.
  double v2[NDIM];
  double grid_area;                           // Total area of the mesh.

  // Ensure the array is clean at startup.
  jstat[0] = jstat[1] = jstat[2] = jstat[3] = 0;
  grid_area = 0.;

  printf("Starting computation of 2D metrics.................\n");
  if ( grid->num_elem[Tri] > 0 )
    {
      // Compute the area of all triangles.
      for ( t=1; t <= grid->num_elem[Tri]; t++ )
	{
	  // Grab the nodes and create Points and then Vectors.
	  n0 = grid->c2n[Tri][t*3+0];
	  n1 = grid->c2n[Tri][t*3+1];
	  n2 = grid->c2n[Tri][t*3+2];
	  
	  // Create the vectors.
	  for ( i=0; i < NDIM; i++ )
	    {
	      v1[i] = grid->x[NDIM*n1+i] - grid->x[NDIM*n0+i];
	      v2[i] = grid->x[NDIM*n2+i] - grid->x[NDIM*n0+i];
	    }
      
	  // Calculate the area using the cross product.
	  area = 0.5 * ( v1[0]*v2[1] - v1[1]*v2[0]);

	  grid_area += area;
      
	  // If we have non-positive area, increment the counter.
	  if ( area < 0. ) k++;

	  // If we have zero area, we need to post a warning and quit.
	  if ( fabs(area) < 1.e-15 )
	    {
	      printf("  FATAL ERROR - ZERO area detected on triangle %d. AREA = %.15e\n\n",t,area);
	      exit(1);
	    }
	}
      if (k==0) printf("  Triangle Report -> All areas positive.\n");
      else printf("  Triangle Report -> Detected %d negative areas.\n",k);
    }

  // Compute the area of the quads. Use the average of cross products taken from opposite corners
  // to give reasonable results if the quad is non-planar. This shouldn't be an issue for the 2D
  // code, so I may want to use a simpler method later to take advantage of z=0.
  if ( grid->num_elem[Quad] > 0 )
    {
      for ( t=1; t <= grid->num_elem[Quad]; t++ )
	{
	  // Grab the nodes and create Points and then Vectors.
	  n0 = grid->c2n[Quad][t*4+0];
	  n1 = grid->c2n[Quad][t*4+1];
	  n2 = grid->c2n[Quad][t*4+2];
	  n3 = grid->c2n[Quad][t*4+3];

	  // Create the vectors.                          // Using ENCM 516 Notes from Dr. Karman.
	  for ( i=0; i < NDIM; i++ )
	    {
	      v1[i] = grid->x[NDIM*n2+i] - grid->x[NDIM*n0+i];
	      v2[i] = grid->x[NDIM*n3+i] - grid->x[NDIM*n1+i];
	    }

	  // Calculate the area using the cross product.
	  area = 0.5*( v1[0]*v2[1] - v1[1]*v2[0]);  // Assumes planar element in 2D.

	  grid_area += area;

	  // Count up negative areas.
	  if ( area < 0. ) m++;

	  //If we have a zero area, print the warning and then quit.
	  if ( fabs(area) < 1.e-15 )
	    {
	      printf("  FATAL ERROR - ZERO area detected on quad %d. AREA = %.15e\n\n",t,area);
	      exit(1);
	    }
	}

      if (m==0) printf("  Quad Report -> All areas positive.\n");
      else printf("  Quad Report -> Detected %d negative areas.\n",m);

      // Compute the corner jacobians using notes from Dr. Karman in Grid Generation 1.
      // Keep track statistics.
      
      for ( t=1; t <= grid->num_elem[Quad]; t++ )
	{
	  // Set the element average to zero.
	  javg = 0.;

	  // Set the flag to zero - Assume no negative corners.
	  flag = 0;

	  // Grab the nodes and create Points and then Vectors.
	  n0 = grid->c2n[Quad][t*4+0];
	  n1 = grid->c2n[Quad][t*4+1];
	  n2 = grid->c2n[Quad][t*4+2];
	  n3 = grid->c2n[Quad][t*4+3];

	  // First corner.
	  for ( i=0; i < NDIM; i++ )
	    {
	      v1[i] = grid->x[NDIM*n1+i] - grid->x[NDIM*n0+i];
	      v2[i] = grid->x[NDIM*n3+i] - grid->x[NDIM*n0+i];
	    }

	  mag = sqrt(v1[0]*v1[0]+v1[1]*v1[1]);
	  if (mag > 1.0e-20) { v1[0] = v1[0] / mag; v1[1] = v1[1] / mag; }

	  mag = sqrt(v2[0]*v2[0]+v2[1]*v2[1]);
	  if (mag > 1.0e-20) { v2[0] = v2[0] / mag; v2[1] = v2[1] / mag; }
	  
	  j = ( v1[0]*v2[1] - v1[1]*v2[0] );

	  // Keep track of the average.
	  javg += j;
	  
	  // Increment flag if the corner jacobian is negative.
	  if ( j < 0. )
	    {
	      //printf("  FATAL ERROR: Nonpositive corner Jacobian detected on quad %d.\n\n",t);
	      //exit(1);
	      flag++;
	    }

	  // Second corner.
	  for ( i=0; i < NDIM; i++ )
	    {
	      v1[i] = grid->x[NDIM*n2+i] - grid->x[NDIM*n1+i];
	      v2[i] = grid->x[NDIM*n0+i] - grid->x[NDIM*n1+i];
	    }

	  mag = sqrt(v1[0]*v1[0]+v1[1]*v1[1]);
	  if (mag > 1.0e-20) { v1[0] = v1[0] / mag; v1[1] = v1[1] / mag; }

	  mag = sqrt(v2[0]*v2[0]+v2[1]*v2[1]);
	  if (mag > 1.0e-20) { v2[0] = v2[0] / mag; v2[1] = v2[1] / mag; }
	  
	  j = ( v1[0]*v2[1] - v1[1]*v2[0] );

	  javg += j;
	  
	  if ( j < 0. )
	    {
	      //printf("  FATAL ERROR: Nonpositive corner Jacobian detected on quad %d.\n\n",t);
	      //exit(1);
	      flag++;
	    }

	  // Third corner.
	  for ( i=0; i < NDIM; i++ )
	    {
	      v1[i] = grid->x[NDIM*n3+i] - grid->x[NDIM*n2+i];
	      v2[i] = grid->x[NDIM*n1+i] - grid->x[NDIM*n2+i];
	    }

	  mag = sqrt(v1[0]*v1[0]+v1[1]*v1[1]);
	  if (mag > 1.0e-20) { v1[0] = v1[0] / mag; v1[1] = v1[1] / mag; }

	  mag = sqrt(v2[0]*v2[0]+v2[1]*v2[1]);
	  if (mag > 1.0e-20) { v2[0] = v2[0] / mag; v2[1] = v2[1] / mag; }
	  
	  j = ( v1[0]*v2[1] - v1[1]*v2[0] );

	  javg += j;

	  if ( j < 0. )
	    {
	      //printf("  FATAL ERROR: Nonpositive corner Jacobian detected on quad %d.\n\n",t);
	      //exit(1);
	      flag++;
	    }

	  // Fourth corner.
	  for ( i=0; i < NDIM; i++ )
	    {
	      v1[i] = grid->x[NDIM*n0+i] - grid->x[NDIM*n3+i];
	      v2[i] = grid->x[NDIM*n2+i] - grid->x[NDIM*n3+i];
	    }

	  mag = sqrt(v1[0]*v1[0]+v1[1]*v1[1]);
	  if (mag > 1.0e-20) { v1[0] = v1[0] / mag; v1[1] = v1[1] / mag; }

	  mag = sqrt(v2[0]*v2[0]+v2[1]*v2[1]);
	  if (mag > 1.0e-20) { v2[0] = v2[0] / mag; v2[1] = v2[1] / mag; }
	  
	  j = ( v1[0]*v2[1] - v1[1]*v2[0] );

	  javg += j;

	  if ( j < 0. )
	    {
	      //printf("  FATAL ERROR: Nonpositive corner Jacobian detected on quad %d.\n\n",t);
	      //exit(1);
	      flag++;
	    }

	  // Determine the element Jacobian and add it the Jacobian statistics array.
	  
	  // Determine the element average.
	  javg /= 4.;

	  // If all corners are negative, the element must be.
	  if ( flag == 4 ) jstat[0]++;
	  else if ( javg <= 0. && flag < 4 ) jstat[1]++;   // Average is negative, but has at least one positive corner.
	  else if ( javg > 0. && flag > 0 ) jstat[2]++;    // Average is positive, but has at least one negative corner.
	  else jstat[3]++;     // Last resort, the element is positive.
	}
      printf("  Quad Report -> Jacobian Statistics:.\n");
      printf("    Negative.............%d\n",jstat[0]);
      printf("    Negative Skewed......%d\n",jstat[1]);
      printf("    Positive Skewed......%d\n",jstat[2]);
      printf("    Positive.............%d\n\n",jstat[3]);
    }

  printf("Total area spanned by the grid is %.15E\n",grid_area);

  // If all elements have negative areas, then assume the winding is wrong and proceed to rewind the elements.
  if ( k==grid->num_elem[Tri] && m==grid->num_elem[Quad] )
    {
      printf("WARNING: All elements are negative. Mesh elements will be rewound.\n");
      rewind2D(grid);

      // Now validate that all elements have positive area.
      validate2D(grid);
    }
  else if ( m > 0 && k > 0 )
    {
      printf("FATAL ERROR: Negative elements are present in the mesh. Exiting.\n\n");
      exit(1);
    }

  return;
}


//=============================================================
// 
//  validate2D()
//
//  Much the same as metrics2D except that it is checking that
//  all elements have positive area only. Intended to be called
//  only after rewind2D().
//
//  GRID *grid                               // The grid.
//
//=============================================================

void validate2D ( GRID *grid )
{
  int i,t;                                    // Loop counters.
  int k=0, m=0;                               // Negative area element counters.
  int n0, n1, n2, n3;                         // Element nodes.
  double area;                                // Element area.
  double v1[NDIM];                            // Vectors.
  double v2[NDIM];
  double grid_area = 0.;                      // Total grid area.

  // Compute the area of all triangles.
  for ( t=1; t <= grid->num_elem[Tri]; t++ )
    {
      // Grab the nodes and create Points and then Vectors.
      n0 = grid->c2n[Tri][t*3+0];
      n1 = grid->c2n[Tri][t*3+1];
      n2 = grid->c2n[Tri][t*3+2];
      
      // Create the vectors.
      for ( i=0; i < NDIM; i++ )
	{
	  v1[i] = grid->x[NDIM*n1+i] - grid->x[NDIM*n0+i];
	  v2[i] = grid->x[NDIM*n2+i] - grid->x[NDIM*n0+i];
	}
      
      // Calculate the area using the cross product.
      area = 0.5 * ( v1[0]*v2[1] - v1[1]*v2[0]);

      grid_area += area;
      
      // If we have non-positive area, increment the counter.
      if ( area < 0. ) k++;
    }
  
  for ( t=1; t <= grid->num_elem[Quad]; t++ )
    {
      // Grab the nodes and create Points and then Vectors.
      n0 = grid->c2n[Quad][t*4+0];
      n1 = grid->c2n[Quad][t*4+1];
      n2 = grid->c2n[Quad][t*4+2];
      n3 = grid->c2n[Quad][t*4+3];
      
      // Create the vectors.                          // Using ENCM 516 Notes from Dr. Karman.
      for ( i=0; i < NDIM; i++ )
	{
	  v1[i] = grid->x[NDIM*n2+i] - grid->x[NDIM*n0+i];
	  v2[i] = grid->x[NDIM*n3+i] - grid->x[NDIM*n1+i];
	}
      
      // Calculate the area using the cross product.
      area = 0.5*( v1[0]*v2[1] - v1[1]*v2[0]);  // Assumes planar element in 2D.

      grid_area += area;

      // Count up negative areas.
      if ( area < 0. ) m++;
    }

  printf("validate2D: Total area spanned by the grid is %.15E\n",grid_area);

  if ( k != 0 && m != 0 )
    {
      printf("FATAL ERROR: NEGATIVE element areas persist in the mesh. Exiting.\n");
      exit(1);
    }
  else
    printf("validate2D: Mesh is valid.\n");

  return;
}
