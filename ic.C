//=============================================================
// 
//  ic.C
//  
//  Initializes the flow field.
//
//  Written by - Shane Sawyer
//
//=============================================================

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "grid.h"
#include "Defined_Variables.h"
#include "params.h"
#include "reconstruct.h"
#include "mesh_utils.h"
#include "ic.h"
#include "curved_boundaries.h"




//=============================================================
// 
//  init_pv()
//
//  Initializes the flow field of a mesh using the input parameters.
//  Values are calculated as point values.
//
//  GRID *grid;                       // The grid.
//  PARAMS p;                         // The parameters.
//
//=============================================================
void init_pv ( GRID *grid, PARAMS p )
{
  int n;                                 // Node loop counter.
  int i,gn;                              // Edge and ghost node counter.
  int b;

  // Pointers.
  double *x = NULL;
  double *q = NULL;
   
  FILE *fp;                              // Pointer to the file.
  const int bdim = 132;                  // Default array dimension size.
  char buff[bdim];                       // Default buffer for reading in from a stream.

  // Each initial condition function is given the coordinates, parameters, and the Q to store.

  // It is assumed here that memory has been allocated for the solution.
  if ( grid->Q == NULL || grid->dQ == NULL )
    {
      printf("MEMORY ERROR: Q or dQ has not been allocated before initializing the solution!\n");
      exit(0);
    }

  for ( n=1; n <= grid->nn; n++ )
    {
      // Adjust the pointers.
      x = &(grid->x[n*NDIM]);
      q = &(grid->Q[n*NUM_VAR]);

      // Find the appropriate ic function.
      if ( p.isrestart == 1 )
	{ }  // Do nothing for a restart.
      else if ( p.ic == 0 )
	ic_0 ( p,x,q );
      else if ( p.ic == 1 )
	ic_1 ( p,x,q );
      else if ( p.ic == 2 )
	ic_2 ( p,x,q );
      else if ( p.ic == 3 )
	ic_3 ( p,x,q );
      else
	{
	  printf("FATAL ERROR: An unknown ic function was selected: ic = %d\n",p.ic);
	  exit(0);
	}
    }
  
  // Copy the boundary node values to its ghost node.
  for ( i=1; i <= grid->nbedges; i++ )
    {
      // Grab the data from the structure.
      n =  grid->bedges[i*5+0];
      gn = grid->bedges[i*5+1] + grid->nn;
      
      // Copy the variables over.
      grid->Q[gn*NUM_VAR+0] = grid->Q[n*NUM_VAR+0];
      grid->Q[gn*NUM_VAR+1] = grid->Q[n*NUM_VAR+1];
      grid->Q[gn*NUM_VAR+2] = grid->Q[n*NUM_VAR+2];
      grid->Q[gn*NUM_VAR+3] = grid->Q[n*NUM_VAR+3];
      
    }
  
  return;
}


//=============================================================
// 
//  init_cv()
//
//  Initializes the flow field of a mesh using the input parameters.
//  Values are calculated as control volume averages.
//
//  GRID *grid;                       // The grid.
//  PARAMS p;                         // The parameters.
//
//=============================================================
void init_cv ( GRID *grid, PARAMS p )
{
  int n;                                 // Node loop counter.
  int i,j,k,e;                           // Loop counters.
  int v;                                 // Vector loop.
  int ind;                               // Component loop.
  int nodes[MAX_NUM_VERT];               // Element vertices.
  int node;                              // Node index.
  int ghost_node;
  int b,seg,bc;
  int gn;                                // Ghost node index.
  int num_vert;                          // Number of vertices for an element.
  int left_id;                           // ID of left edge.
  int right_id;                          // ID of right edge.
  int ct_flag;
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

  double XGP[NDIM];
  double QGP[NUM_VAR];
  double QTRI[NUM_VAR];

  double gp[7][3];                       // Gauss points for integration.
  double dA;                             // Differential area (for a triangle).

  double ct_gp[7][2];                    // Gauss integration points in (r,s) space (xi,eta).
  double ct_nodes[7][2];                 // The x,y locations of the 6 nodes defining the p2 triangle (curved side).
  double N1,N2,N3,N4,N5,N6;              // Shape functions.
  double N1r,N1s,N2r,N2s,N3r,N3s;        // Shape function derivatives.
  double N4r,N4s,N5r,N5s,N6r,N6s;
  double xr,xs,yr,ys,jacobian;           // Transformation jacobian.
  double xi,eta;                         // Triangle coordinates.
  

  // We need to assert that certain things have been done by checking
  // if the appropriate memory has been allocated.

  // HOLD OFF ON THIS FOR NOW UNLESS WE ABSOLUTELY NEED THIS!


  // Clean out memory because we are accumulating.
  for ( n=1; n <= grid->nn; n++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	grid->Q[n*NUM_VAR+j] = 0.;
    }

  // Now we calculate the area-averaged value of the initial condition using Gaussian quadrature
 
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
	  xc[0] /= ((double)num_vert);
	  xc[1] /= ((double)num_vert);

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

		  xi = ( (ct_nodes[1][0]*ct_nodes[3][1] - ct_nodes[1][1]*ct_nodes[3][0]) -
			 (ct_nodes[5][0]*ct_nodes[3][1] - ct_nodes[5][1]*ct_nodes[3][0]) +
			 (ct_nodes[5][0]*ct_nodes[1][1] - ct_nodes[5][1]*ct_nodes[1][0])  ) /
		       ( (ct_nodes[1][0]*ct_nodes[3][1] - ct_nodes[1][1]*ct_nodes[3][0]) -
			 (ct_nodes[2][0]*ct_nodes[3][1] - ct_nodes[2][1]*ct_nodes[3][0]) +
			 (ct_nodes[2][0]*ct_nodes[1][1] - ct_nodes[2][1]*ct_nodes[1][0])  );

		  eta = ( (ct_nodes[1][0]*ct_nodes[2][1] - ct_nodes[1][1]*ct_nodes[2][0]) -
			  (ct_nodes[5][0]*ct_nodes[2][1] - ct_nodes[5][1]*ct_nodes[2][0]) +
			  (ct_nodes[5][0]*ct_nodes[1][1] - ct_nodes[5][1]*ct_nodes[1][0])  ) /
		        ( (ct_nodes[1][0]*ct_nodes[2][1] - ct_nodes[1][1]*ct_nodes[2][0]) -
			  (ct_nodes[3][0]*ct_nodes[2][1] - ct_nodes[3][1]*ct_nodes[2][0]) +
			  (ct_nodes[3][0]*ct_nodes[1][1] - ct_nodes[1][0]*ct_nodes[3][1])  );

		  if ( eta < 0.25 || xi < 0.25 )  // The transformation will not be invertible.
		    {
		      printf("CRITICAL ERROR: Triangle 1 of element %d of type %d is too concave (ic.C)!\n",j,e);
		      printf("  xi = %.15e     eta = %.15e\n",xi,eta);
		      fflush(stdout);
		      exit(1);
		    }

		  // Do this for each variable.

		  QTRI[0] = 0.; QTRI[1] = 0.; QTRI[2] = 0.; QTRI[3] = 0.;

		  // Now I can loop over the Gaussian integration nodes and apply the quadratic transformation.
		  for ( i=0; i < 7; i++ )
		    {
		      xi = ct_gp[i][0];
		      eta = ct_gp[i][1];
		      
		      // Compute basis functions.
		      N1 = 1. - 3.*xi - 3.*eta + 4.*xi*eta + 2.*xi*xi + 2.*eta*eta;
		      N2 = xi*(2.*xi - 1.);
		      N3 = eta*(2.*eta - 1.);
		      N4 = 4.*xi*(1. - xi - eta);
		      N5 = 4.*xi*eta;
		      N6 = 4.*eta*(1. - xi - eta);

		      // Compute their derivatives.
		      N1r = -3. + 4.*eta + 4.*xi;
		      N1s = -3. + 4.*xi + 4.*eta;
		      N2r = 4.*xi - 1.;
		      N2s = 0.;
		      N3r = 0.;
		      N3s = 4.*eta - 1.;
		      N4r = 4.*(1. - 2.*xi - eta);
		      N4s = -4.*xi;
		      N5r = 4.*eta;
		      N5s = 4.*xi;
		      N6r = -4.*eta;
		      N6s = 4.*(1. - 2.*eta - xi);

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
		      XGP[0] = N1*ct_nodes[1][0] + N2*ct_nodes[2][0] + N3*ct_nodes[3][0] +
			       N4*ct_nodes[4][0] + N5*ct_nodes[5][0] + N6*ct_nodes[6][0];

		      XGP[1] = N1*ct_nodes[1][1] + N2*ct_nodes[2][1] + N3*ct_nodes[3][1] +
			       N4*ct_nodes[4][1] + N5*ct_nodes[5][1] + N6*ct_nodes[6][1];

		      
		      if ( p.ic == 1 )
			ic_1 ( p, XGP, QGP );
		      else if ( p.ic == 2 )
			ic_2 ( p, XGP, QGP );
		      else if ( p.ic == 3 )
			ic_3 ( p, XGP, QGP );
		      else
			ic_0 ( p, XGP, QGP );
		      
		      for ( v=0; v < NUM_VAR; v++ )
			{
			  QTRI[v] += ( wg[i] * QGP[v] *jacobian * 0.5);
			} 
		      
		    }
		  
		  // Add to the running total for the node.
		  for ( v=0; v < NUM_VAR; v++ )
		    {
		      grid->Q[node*NUM_VAR+v] += QTRI[v];
		    }
		  
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
					bary_coord[v][2] * xc[ind]);
			}
		      gp[v][2] = 0.;   // Always 0 since 2D vectors.
		    }

		  // Now we loop over the Gauss points and get the function value
		  // and perform the numerical integration (the area integral of the triangle).
	      
		  // Do this for each variable.
		  
		  QTRI[0] = 0.; QTRI[1] = 0.; QTRI[2] = 0.; QTRI[3] = 0.;
		  
		  for ( i=0; i < 7; i++ )
		    {
		      // Get the value at the point.
		      XGP[0] = gp[i][0];
		      XGP[1] = gp[i][1];
		      
		      if ( p.ic == 1 )
			ic_1 ( p, XGP, QGP );
		      else if ( p.ic == 2 )
			ic_2 ( p, XGP, QGP );
		      else if ( p.ic == 3 )
			ic_3 ( p, XGP, QGP );
		      else
			ic_0 ( p, XGP, QGP );
		      
		      for ( v=0; v < NUM_VAR; v++ )
			{
			  QTRI[v] += ( wg[i] * QGP[v] );
			} 
		    }
		  
		  for ( v=0; v < NUM_VAR; v++ )
		    {
		      QTRI[v] *= dA;
		    }
		  
		  // Accumulate the integral to the node.
		  for ( v=0; v < NUM_VAR; v++ )
		    {
		      grid->Q[node*NUM_VAR+v] += QTRI[v];
		    }
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

		  xi = ( (ct_nodes[1][0]*ct_nodes[3][1] - ct_nodes[1][1]*ct_nodes[3][0]) -
			 (ct_nodes[5][0]*ct_nodes[3][1] - ct_nodes[5][1]*ct_nodes[3][0]) +
			 (ct_nodes[5][0]*ct_nodes[1][1] - ct_nodes[5][1]*ct_nodes[1][0])  ) /
		       ( (ct_nodes[1][0]*ct_nodes[3][1] - ct_nodes[1][1]*ct_nodes[3][0]) -
			 (ct_nodes[2][0]*ct_nodes[3][1] - ct_nodes[2][1]*ct_nodes[3][0]) +
			 (ct_nodes[2][0]*ct_nodes[1][1] - ct_nodes[2][1]*ct_nodes[1][0])  );

		  eta = ( (ct_nodes[1][0]*ct_nodes[2][1] - ct_nodes[1][1]*ct_nodes[2][0]) -
			  (ct_nodes[5][0]*ct_nodes[2][1] - ct_nodes[5][1]*ct_nodes[2][0]) +
			  (ct_nodes[5][0]*ct_nodes[1][1] - ct_nodes[5][1]*ct_nodes[1][0])  ) /
		        ( (ct_nodes[1][0]*ct_nodes[2][1] - ct_nodes[1][1]*ct_nodes[2][0]) -
			  (ct_nodes[3][0]*ct_nodes[2][1] - ct_nodes[3][1]*ct_nodes[2][0]) +
			  (ct_nodes[3][0]*ct_nodes[1][1] - ct_nodes[1][0]*ct_nodes[3][1])  );

		  if ( eta < 0.25 || xi < 0.25 )  // The transformation will not be invertible.
		    {
		      printf("CRITICAL ERROR: Triangle 2 of element %d of type %d is too concave (ic.C)!\n",j,e);
		      printf("  xi = %.15e     eta = %.15e\n",xi,eta);
		      fflush(stdout);
		      exit(1);
		    }

		  // Do this for each variable.

		  QTRI[0] = 0.; QTRI[1] = 0.; QTRI[2] = 0.; QTRI[3] = 0.;

		  // Now I can loop over the Gaussian integration nodes and apply the quadratic transformation.
		  for ( i=0; i < 7; i++ )
		    {
		      xi = ct_gp[i][0];
		      eta = ct_gp[i][1];
		      
		      // Compute basis functions.
		      N1 = 1. - 3.*xi - 3.*eta + 4.*xi*eta + 2.*xi*xi + 2.*eta*eta;
		      N2 = xi*(2.*xi - 1.);
		      N3 = eta*(2.*eta - 1.);
		      N4 = 4.*xi*(1. - xi - eta);
		      N5 = 4.*xi*eta;
		      N6 = 4.*eta*(1. - xi - eta);

		      // Compute their derivatives.
		      N1r = -3. + 4.*eta + 4.*xi;
		      N1s = -3. + 4.*xi + 4.*eta;
		      N2r = 4.*xi - 1.;
		      N2s = 0.;
		      N3r = 0.;
		      N3s = 4.*eta - 1.;
		      N4r = 4.*(1. - 2.*xi - eta);
		      N4s = -4.*xi;
		      N5r = 4.*eta;
		      N5s = 4.*xi;
		      N6r = -4.*eta;
		      N6s = 4.*(1. - 2.*eta - xi);

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
		      XGP[0] = N1*ct_nodes[1][0] + N2*ct_nodes[2][0] + N3*ct_nodes[3][0] +
			       N4*ct_nodes[4][0] + N5*ct_nodes[5][0] + N6*ct_nodes[6][0];

		      XGP[1] = N1*ct_nodes[1][1] + N2*ct_nodes[2][1] + N3*ct_nodes[3][1] +
			       N4*ct_nodes[4][1] + N5*ct_nodes[5][1] + N6*ct_nodes[6][1];

		      
		      if ( p.ic == 1 )
			ic_1 ( p, XGP, QGP );
		      else if ( p.ic == 2 )
			ic_2 ( p, XGP, QGP );
		      else if ( p.ic == 3 )
			ic_3 ( p, XGP, QGP );
		      else
			ic_0 ( p, XGP, QGP );
		      
		      for ( v=0; v < NUM_VAR; v++ )
			{
			  QTRI[v] += ( wg[i] * QGP[v] *jacobian * 0.5);
			} 
		      
		    }
		  
		  // Add to the running total for the node.
		  for ( v=0; v < NUM_VAR; v++ )
		    {
		      grid->Q[node*NUM_VAR+v] += QTRI[v];
		    }
		  
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
					bary_coord[v][2] * xc[ind]);
			}
		      gp[v][2] = 0.;   // Always 0 since 2D vectors.
		    }
	      
		  // Now we loop over the Gauss points and get the function value
		  // and perform the numerical integration (the area integral of the triangle).
		  
		  // Do this for each variable.
		  
		  QTRI[0] = 0.; QTRI[1] = 0.; QTRI[2] = 0.; QTRI[3] = 0.;
	      
		  for ( i=0; i < 7; i++ )
		    {
		      // Get the value at the point.
		      XGP[0] = gp[i][0];
		      XGP[1] = gp[i][1];
		      
		      if ( p.ic == 1 )
			ic_1 ( p, XGP, QGP );
		      else if ( p.ic == 2 )
			ic_2 ( p, XGP, QGP );
		      else if ( p.ic == 3 )
			ic_3 ( p, XGP, QGP );
		      else
			ic_0 ( p, XGP, QGP );
		      
		      for ( v=0; v < NUM_VAR; v++ )
			{
			  QTRI[v] += ( wg[i] * QGP[v] );
			} 
		    }
	      
		  for ( v=0; v < NUM_VAR; v++ )
		    {
		      QTRI[v] *= dA;
		    }
		  
		  // Accumulate the integral to the node.
		  for ( v=0; v < NUM_VAR; v++ )
		    {
		      grid->Q[node*NUM_VAR+v] += QTRI[v];
		    }
		}
	      
	    }                                       // End element node loop.
	      
	}                                           // End element loop.
      
    }                                               // End element type loop.


  // This is a debug check.
  if ( 0 )
    {
      double temp = 0.;
      
      for ( i=1; i<= grid->nn; i++ )
	{
	  temp += grid->Q[i*NUM_VAR+0];
	}
      
      printf("TEMP = %.15E\n",temp);
    }


  // Finished the accumulation. Now divide by the control volume area.
  for ( i=1; i<= grid->nn; i++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  grid->Q[i*NUM_VAR+j] /= grid->cv_area[i];
	}
    }
  
  // Copy the boundary node values to its ghost node.
  for ( i=1; i <= grid->nbedges; i++ )
    {
      // Grab the data from the structure.
      n =  grid->bedges[i*5+0];
      gn = grid->bedges[i*5+1] + grid->nn;

      // Copy the variables over.
      grid->Q[gn*NUM_VAR+0] = grid->Q[n*NUM_VAR+0];
      grid->Q[gn*NUM_VAR+1] = grid->Q[n*NUM_VAR+1];
      grid->Q[gn*NUM_VAR+2] = grid->Q[n*NUM_VAR+2];
      grid->Q[gn*NUM_VAR+3] = grid->Q[n*NUM_VAR+3];
      
    }

  // Write out the data for each variable.
  if ( MMS && RECON_PRIM )
    {
      FILE *fp;

      fp = fopen("density_init.dat","w");

      for ( i=1; i <= grid->nn; i++ )
	{
	  fprintf(fp,"%.15e\n",grid->Q[i*NUM_VAR+0]);
	}
      fclose(fp);

      fp = fopen("xvel_init.dat","w");

      for ( i=1; i <= grid->nn; i++ )
	{
	  fprintf(fp,"%.15e\n",grid->Q[i*NUM_VAR+1]);
	}
      fclose(fp);

      fp = fopen("yvel_init.dat","w");

      for ( i=1; i <= grid->nn; i++ )
	{
	  fprintf(fp,"%.15e\n",grid->Q[i*NUM_VAR+2]);
	}
      fclose(fp);

      fp = fopen("pressure_init.dat","w");

      for ( i=1; i <= grid->nn; i++ )
	{
	  fprintf(fp,"%.15e\n",grid->Q[i*NUM_VAR+3]);
	}
      fclose(fp);
    }

  //exit(0);

  return;
}



//=============================================================
// 
//  ic_0()
//
//  Initializes the flow field of a mesh using the input parameters.
//  Generic flow.
//
//  PARAMS p;                            // The parameters.
//  double *x;                           // Point location.
//  double *q;                           // Solution at point.
//
//=============================================================
void ic_0 ( PARAMS p, double *x, double *q )
{
  double rho;                                   // Fluid density.
  double c;                                     // Local speed of sound.
  double P;                                     // Fluid pressure.
  double u;                                     // X-component of velocity.
  double v;                                     // Y-component of velocity.
  double ei;                                    // Internal energy.
  double et;                                    // Total energy.
  double angle;                                 // Angle of attack in radians.
  double pi = M_PI;
  
  angle = p.alpha*M_PI/180.;

  // Now initialize the dependent and primitive variables to reasonable values. 
  rho = 1.0;
  c = 1.0;
  P = rho*c*c / p.gamma;
  u = p.mach * cos( angle );
  v = p.mach * sin( angle );
  //ei = P / ( (p.gamma - 1.0) * rho );
  //et = rho*( ei + 0.5*(u*u + v*v) );
  et = P/(p.gamma - 1.0) + 0.5*rho*(u*u + v*v);

  q[0] = rho;
  q[1] = rho*u;
  q[2] = rho*v;
  q[3] = et;

  if ( ANNULUS )     // Give it an at rest initial condition.
    {
      rho = 0.6;
      c = 1.0;
      P = rho*c*c / p.gamma;
      //P = 1.0 / p.gamma;
      u = 0.;
      v = 0.;

      q[0] = rho;
      q[1] = 0.;
      q[2] = 0.;
      q[3] = P/(p.gamma-1.0) + 0.5*rho*(u*u + v*v);
    }

  if ( MMS )
    {
      rho = 1.0 + sin(pi*x[0])*sin(pi*x[1]);
      u = 0.25 + 0.25*sin(pi*x[0])*cos(2.*pi*x[1]);
      v = 0.25 + 0.25*cos(2.*pi*x[0])*sin(pi*x[1]);
      P = 1.0/p.gamma + 0.05*cos(2.*pi*x[0])*cos(2.*pi*x[1]);
      et = P / (p.gamma - 1.) + 0.5*rho*(u*u + v*v);

      // Store the conserved variables.
      q[0] = rho;
      q[1] = rho*u;
      q[2] = rho*v;
      q[3] = et;
    }

  if ( ANNULUS && MMS )
    {
      double rhoi,r,r2,Mi,U,Ui,Ri;
      rhoi = 1.0;
      Mi = 2.0;
      Ui = 2.0;
      Ri = 2.0;

      r2 = (x[0]*x[0]) + (x[1]*x[1]);
      r = sqrt(r2);

      rho = pow( ( 1.0 + ((p.gamma - 1.0)/2.0)*Mi*Mi*( 1.0 - (Ri*Ri)/r2)) , (1.0/(p.gamma-1.0)) );
      U   = (Ui*Ri)/r;
      u   = (x[1]*U)/r;
      v   = (-x[0]*U)/r;
      P   = ( pow( rho , p.gamma) )/p.gamma;

      et  = P / (p.gamma - 1.) + 0.5*rho*(u*u + v*v);

      // Store the conserved variables.
      q[0] = rho;
      q[1] = rho*u;
      q[2] = rho*v;
      q[3] = et;
    }

  if ( MMS_EXP )
    {
      rho = exp(x[0]-1.0) * exp(x[1] - 1.0);
      u   = exp(x[0]-1.0) * exp(x[1] - 1.0);
      v   = exp(x[0]-1.0) * exp(x[1] - 1.0);
      et  = 1.5 * exp( 3.0 * ( x[0] - 1.0) ) * exp( 3.0 * ( x[1] - 1.0 ) );
      P   = ( (p.gamma-1.0) * 0.5 ) * exp( 3.0 * ( x[0] - 1.0) ) * exp( 3.0 * ( x[1] - 1.0 ) );

      // Store the conserved variables.
      q[0] = rho;
      q[1] = rho*u;
      q[2] = rho*v;
      q[3] = et;
    }

  // Debug stuff
  
  //q[0] = x[0]*x[0]*x[1]*x[1];
  
  return;
}

//=============================================================
// 
//  ic_1()
//
//  Initializes the flow field of a mesh using the input parameters.
//  Vortex convection problem.
//
//  Based on the code from UX. Domain is:
//                               0 < x < xmax
//                              -5 < y < 5
//
//  PARAMS p;                            // The parameters.
//  double *x;                           // Point location.
//  double *q;                           // Solution at point.
//
//=============================================================
void ic_1 ( PARAMS p, double *x, double *q )
{
  double gamma = p.gamma;                          // Ratio of coefficients of heat.
  double gm1 = gamma - 1.0;
  double mach = p.mach;                            // Mach number.
  double x0;                                       // x location of vortex at time t = 0
  double y0;                                       // y location of vortex at time t = 0
  double r2;                                       // Radius squared.
  double beta;                                     // Vortex strength.
  double fac1,fac2,fac3;                           // Terms in the computation.
  double rho0;                                     // Core density.
  double u0,v0;                                    // Core velocity.
  double P0,T0,s0;                                 // Core pressure,temperature, and entropy.
  double T,u,v,P,rho;                              // Flow properties.


  // Comments from UX=================================================
  // An isentropic perturbation is being added
  
  // This is based on a modification of Yee's test case
  // and comes to us from AEDC (Bobby Nichols) via Roy at UAB
  
  // This is good for the compressible and variable Mach regimes
  //===================================================================

  x0 = 5.;
  y0 = 0.;

  // Initialize the variables.
  r2 = (x[0] - x0)*(x[0] - x0) + (x[1] - y0)*(x[1] - y0);
  beta = 5.0;
  fac1 = beta/(2.0*M_PI);
  fac2 = -0.5*gm1/gamma*fac1*fac1;
  fac3 = exp((1.0 - r2)/2.0);

  rho0 = 1.0;
  u0 = mach*1.0;          // Flow in only the x direction at reference velocity of 1.
  v0 = 0.;
  P0 = 1.0/gamma;
  T0 = P0/rho0;
  s0 = P0/pow(rho0,gamma);

  T = T0 + fac2*fac3*fac3;
  u = u0 - (fac1*fac3)*(x[1] - y0);
  v = v0 + (fac1*fac3)*(x[0] - x0);
  P = pow(pow(T, gamma)/s0,1.0/gm1);
  rho = P/T;

  q[0] = rho;
  q[1] = rho*u;
  q[2] = rho*v;
  q[3] = P/gm1 + 0.5*rho*(u*u + v*v);

  return;
}


//=============================================================
// 
//  ic_2()
//
//  Initializes the flow field of a mesh using the input parameters.
//  Vortex convection problem of Yee. This is derived from Dr. Li
//  Wang's dissertation and the uses the domain from Tenasi.
//
//  Based on the code from UX. Domain is:
//                               0 < x < xmax
//                              -5 < y < 5
//
//  PARAMS p;                            // The parameters.
//  double *x;                           // Point location.
//  double *q;                           // Solution at point.
//
//=============================================================
void ic_2 ( PARAMS p, double *x, double *q )
{
  double gamma = p.gamma;                          // Ratio of coefficients of heat.
  double gm1 = gamma - 1.0;
  double oogm1;                                    // One over (gamma-1)
  double x0,y0;                                    // Vortex center.
  double r2;                                       // Distance from center.
  double phi,sigma;                                // Parameter to control vortex strength.
  double upert,vpert,Tpert;                        // Isentropic vortex perturbations.
  double rho,u,v,P,T,E;                            // Control volume values of conserved variables.
  double Mref,rhoref,Pref,Tref;
  double rhoD,PD;                                  // Dimensional quantities.
  double FC;                                       // Fixed constant.

  Mref = 340.294065;
  Tref = 288.15;
  Pref = 101325.0;
  rhoref = 1.225;

  FC = Pref / pow(rhoref,gamma);

  x0 = 5.;
  y0 = 0.;

  // Initialize the variables.
  r2 = (x[0] - x0)*(x[0] - x0) + (x[1] - y0)*(x[1] - y0);
  oogm1 = 1./gm1;

  sigma = 4.0;
  phi = 1.0;

  upert = (-sigma)/(2.*M_PI) * ( x[1] - y0 ) * exp( phi*(1.-r2) );
  vpert = (sigma)/(2.*M_PI) * ( x[0] - x0 ) * exp( phi*(1.-r2) );
  Tpert = ( -1. * sigma*sigma * gm1 )/( 16. * phi * gamma * M_PI * M_PI ) * exp( 2. * phi * (1.-r2) );

  // Now we perturb the variables.

  // Use nondim T
  T = 1.0 + Tpert;

  // Compute nondim rho.
  rho = pow( ( T ) , oogm1 );

  rhoD = rho*rhoref;
  
  // Compute nondim velocities.
  u = p.mach + upert;
  v = 0. + vpert;

  // Compute the pressure.
  PD = pow( rhoD , gamma ) * FC;

  P = PD / (gamma * Pref );

  // Nondimensionalize the numbers.
  //rho = rho / rhoref;
  //u = u / Mref;
  //v = v / Mref;
  //P = P / (gamma*Pref);

  E = P/gm1 + 0.5*rho*(u*u + v*v);

  q[0] = rho;
  q[1] = rho*u;
  q[2] = rho*v;
  q[3] = E;

  return;
}



//=============================================================
// 
//  ic_3()
//
//  Initialize the flow based on a plausible physical solution
//  as suggested by Ollivier-Gooch in AIAA Journal vol 47, no 9,
//  page 2111.
//
//  PARAMS p;                            // The parameters.
//  double *x;                           // Point location.
//  double *q;                           // Solution at point.
//
//=============================================================

// STORES PRIMITIVE VARIABLES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

void ic_3 ( PARAMS p, double *x, double *q )
{
  double rho,u,v,P;
  double rho0,u0,v0,P0;

  rho0 = 0.1;
  u0 = 0.1;
  v0 = 0.1;
  P0 = 0.1;

  rho = 1. + rho0 * sin(M_PI*x[0]) * sin(M_PI*x[1]);
  u = u0 * sin(M_PI*x[0]) * cos(2.*M_PI*x[1]);
  v = v0 * cos(2.*M_PI*x[0]) * sin(M_PI*x[1]);
  P = 1/p.gamma + P0*sin(2.*M_PI*x[0]) * sin(2.*M_PI*x[1]);

  q[0] = rho;
  q[1] = u;
  q[2] = v;
  q[3] = P;

  return;
}


//=============================================================
// 
//  initialize_vortex_pv()
//
//  Sets flow variables for the Vortex convection problem at a
//  particular time. Used for comparing analytical and numerical
//  solutions.
//
//  Based on the code from UX. Domain is:
//                               0 < x < xmax
//                              -5 < y < 5
//
//  PARAMS p;                            // The parameters.
//  GRID *grid;                          // The grid structure.
//  double *Q;                           // Pointer to the analytical solution.
//  ******  Assumes Q has been initialized. *****
//=============================================================

void initialize_vortex_pv ( PARAMS p, GRID *grid, double *Q )
{
  int i;                                           // Node index.
  double gamma = p.gamma;                          // Ratio of coefficients of heat.
  double gm1 = gamma - 1.0;
  double mach = p.mach;                            // Mach number.
  double x0;                                       // x location of vortex at time t = 0
  double y0;                                       // y location of vortex at time t = 0
  double r2;                                       // Radius squared.
  double beta;                                     // Vortex strength.
  double fac1,fac2,fac3;                           // Terms in the computation.
  double rho0;                                     // Core density.
  double u0,v0;                                    // Core velocity.
  double P0,T0,s0;                                 // Core pressure,temperature, and entropy.
  double T,u,v,P,rho;                              // Flow properties.
  double *x = NULL;                                // Pointer to grid points.
  double *q = NULL;                                // Pointer to state values.


  // Comments from UX=================================================
  // An isentropic perturbation is being added
  
  // This is based on a modification of Yee's test case
  // and comes to us from AEDC (Bobby Nichols) via Roy at UAB
  
  // This is good for the compressible and variable Mach regimes
  //===================================================================

  x0 = 5. + ( grid->citer * p.mintime * p.mach );   // Put the vortex center at the correct spot.
  y0 = 0.;

  for ( i=1; i <= grid->nn; i++ )
    {
      // Adjust the pointers.
      x = &(grid->x[i*NDIM]);
      q = &(Q[i*NUM_VAR]);

      // Initialize the variables.
      r2 = (x[0] - x0)*(x[0] - x0) + (x[1] - y0)*(x[1] - y0);
      beta = 5.0;
      fac1 = beta/(2.0*M_PI);
      fac2 = -0.5*gm1/gamma*fac1*fac1;
      fac3 = exp((1.0 - r2)/2.0);
      
      rho0 = 1.0;
      u0 = mach*1.0;          // Flow in only the x direction at reference velocity of 1.
      v0 = 0.;
      P0 = 1.0/gamma;
      T0 = P0/rho0;
      s0 = P0/pow(rho0,gamma);
      
      T = T0 + fac2*fac3*fac3;
      u = u0 - (fac1*fac3)*(x[1] - y0);
      v = v0 + (fac1*fac3)*(x[0] - x0);
      P = pow(pow(T, gamma)/s0,1.0/gm1);
      rho = P/T;
      
      q[0] = rho;
      q[1] = rho*u;
      q[2] = rho*v;
      q[3] = P/gm1 + 0.5*rho*(u*u + v*v);
    }
  return;
}

//=============================================================
// 
//  initialize_vortex_cv()
//
//  Sets flow variables for the Vortex convection problem at a
//  particular time. Used for comparing analytical and numerical
//  solutions. This is the cv averaged function.
//
//  THIS IS HARD SET FOR STRAIGHT BOUNDARIES.
//
//  Based on the code from UX. Domain is:
//                               0 < x < xmax
//                              -5 < y < 5
//
//  PARAMS p;                            // The parameters.
//  GRID *grid;                          // The grid structure.
//  double *Q;                           // Pointer to the analytical solution.
//  ******  Assumes Q has been initialized. *****
//=============================================================

void initialize_vortex_cv ( PARAMS p, GRID *grid, double *Q )
{
  int n;                                 // Node loop counter.
  int i,j,k,e;                           // Loop counters.
  int v;                                 // Vector loop.
  int ind;                               // Component loop.
  int nodes[MAX_NUM_VERT];               // Element vertices.
  int node;                              // Node index.
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
  double XGP[NDIM];
  double QGP[NUM_VAR];
  double QTRI[NUM_VAR];
  double gp[7][3];                       // Gauss points for integration.
  double dA;                             // Differential area (for a triangle).
  double gamma = p.gamma;                          // Ratio of coefficients of heat.
  double gm1 = gamma - 1.0;
  double mach = p.mach;                            // Mach number.
  double x0;                                       // x location of vortex at time t = 0
  double y0;                                       // y location of vortex at time t = 0
  double r2;                                       // Radius squared.
  double beta;                                     // Vortex strength.
  double fac1,fac2,fac3;                           // Terms in the computation.
  double rho0;                                     // Core density.
  double u0,v0;                                    // Core velocity.
  double P0,T0,s0;                                 // Core pressure,temperature, and entropy.
  double T,U,V,P,rho;                              // Flow properties.


  // Comments from UX=================================================
  // An isentropic perturbation is being added
  
  // This is based on a modification of Yee's test case
  // and comes to us from AEDC (Bobby Nichols) via Roy at UAB
  
  // This is good for the compressible and variable Mach regimes
  //===================================================================

  x0 = 5. + ( grid->citer * p.mintime * p.mach );   // Put the vortex center at the correct spot.
  y0 = 0.;

  // Clean out memory because we are accumulating.
  for ( n=1; n <= grid->nn; n++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	Q[n*NUM_VAR+j] = 0.;
    } 
  
  // Now we calculate the area-averaged value of the initial condition using Gaussian quadrature
 
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
	  xc[0] /= ((double)num_vert);
	  xc[1] /= ((double)num_vert);

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
		      gp[v][ind] = (bary_coord[v][0] * grid->x[node*NDIM+ind] +
				    bary_coord[v][1] * xmidL[ind] +
				    bary_coord[v][2] * xc[ind]);
		    }
		  gp[v][2] = 0.;   // Always 0 since 2D vectors.
		}

	      // Now we loop over the Gauss points and get the function value
	      // and perform the numerical integration (the area integral of the triangle).
	      
	      // Do this for each variable.
		  
	      QTRI[0] = 0.; QTRI[1] = 0.; QTRI[2] = 0.; QTRI[3] = 0.;
		  
	      for ( i=0; i < 7; i++ )
		{
		  // Get the value at the point.
		  XGP[0] = gp[i][0];
		  XGP[1] = gp[i][1];
		      
		  // Initialize the variables.
		  r2 = (XGP[0] - x0)*(XGP[0] - x0) + (XGP[1] - y0)*(XGP[1] - y0);
		  beta = 5.0;
		  fac1 = beta/(2.0*M_PI);
		  fac2 = -0.5*gm1/gamma*fac1*fac1;
		  fac3 = exp((1.0 - r2)/2.0);
      
		  rho0 = 1.0;
		  u0 = mach*1.0;          // Flow in only the x direction at reference velocity of 1.
		  v0 = 0.;
		  P0 = 1.0/gamma;
		  T0 = P0/rho0;
		  s0 = P0/pow(rho0,gamma);
      
		  T = T0 + fac2*fac3*fac3;
		  U = u0 - (fac1*fac3)*(XGP[1] - y0);
		  V = v0 + (fac1*fac3)*(XGP[0] - x0);
		  P = pow(pow(T, gamma)/s0,1.0/gm1);
		  rho = P/T;
      
		  QGP[0] = rho;
		  QGP[1] = rho*U;
		  QGP[2] = rho*V;
		  QGP[3] = P/gm1 + 0.5*rho*(U*U + V*V);
		  
		  for ( v=0; v < NUM_VAR; v++ )
		    {
		      QTRI[v] += ( wg[i] * QGP[v] );
		    } 
		}
	      
	      for ( v=0; v < NUM_VAR; v++ )
		{
		      QTRI[v] *= dA;
		}
	      
	      // Accumulate the integral to the node.
	      for ( v=0; v < NUM_VAR; v++ )
		{
		  Q[node*NUM_VAR+v] += QTRI[v];
		}
	      
	      // Process the second triangle.
	      
	      // Get the area.
	      dA = 0.5*( v_middle[0]*v_right[1] - v_middle[1]*v_right[0] );
	      
	      // Initialize the Gauss points.
	      for ( v=0; v < 7; v++ )
		{
		  for ( ind=0; ind < NDIM; ind++ )
		    {
		      gp[v][ind] = (bary_coord[v][0] * grid->x[node*NDIM+ind] +
				    bary_coord[v][1] * xmidR[ind] +
				    bary_coord[v][2] * xc[ind]);
		    }
		  gp[v][2] = 0.;   // Always 0 since 2D vectors.
		}
	      
	      // Now we loop over the Gauss points and get the function value
	      // and perform the numerical integration (the area integral of the triangle).
	  
	      // Do this for each variable.
	  
	      QTRI[0] = 0.; QTRI[1] = 0.; QTRI[2] = 0.; QTRI[3] = 0.;
	      
	      for ( i=0; i < 7; i++ )
		{
		  // Get the value at the point.
		  XGP[0] = gp[i][0];
		  XGP[1] = gp[i][1];
		      
		  // Initialize the variables.
		  r2 = (XGP[0] - x0)*(XGP[0] - x0) + (XGP[1] - y0)*(XGP[1] - y0);
		  beta = 5.0;
		  fac1 = beta/(2.0*M_PI);
		  fac2 = -0.5*gm1/gamma*fac1*fac1;
		  fac3 = exp((1.0 - r2)/2.0);
		  
		  rho0 = 1.0;
		  u0 = mach*1.0;          // Flow in only the x direction at reference velocity of 1.
		  v0 = 0.;
		  P0 = 1.0/gamma;
		  T0 = P0/rho0;
		  s0 = P0/pow(rho0,gamma);
      
		  T = T0 + fac2*fac3*fac3;
		  U = u0 - (fac1*fac3)*(XGP[1] - y0);
		  V = v0 + (fac1*fac3)*(XGP[0] - x0);
		  P = pow(pow(T, gamma)/s0,1.0/gm1);
		  rho = P/T;
      
		  QGP[0] = rho;
		  QGP[1] = rho*U;
		  QGP[2] = rho*V;
		  QGP[3] = P/gm1 + 0.5*rho*(U*U + V*V);
		  
		  for ( v=0; v < NUM_VAR; v++ )
		    {
		      QTRI[v] += ( wg[i] * QGP[v] );
		    } 
		}
	      
	      for ( v=0; v < NUM_VAR; v++ )
		{
		  QTRI[v] *= dA;
		}
	      
	      // Accumulate the integral to the node.
	      for ( v=0; v < NUM_VAR; v++ )
		{
		  Q[node*NUM_VAR+v] += QTRI[v];
		}
      
	    }                                       // End element node loop.
  
	}                                           // End element loop.

    }                                               // End element type loop.

  // Finished the accumulation. Now divide by the control volume area.
  for ( i=1; i<= grid->nn; i++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  Q[i*NUM_VAR+j] /= grid->cv_area[i];
	}
    }

  return;
}


//=============================================================
// 
//  initialize_vortex_yee_pv()
//
//  Sets flow variables for the Vortex convection problem at a
//  particular time. Used for comparing analytical and numerical
//  solutions.
//
//  Based on the code from UX. Domain is:
//                               0 < x < xmax
//                              -5 < y < 5
//
//  PARAMS p;                            // The parameters.
//  GRID *grid;                          // The grid structure.
//  double *Q;                           // Pointer to the analytical solution.
//  ******  Assumes Q has been initialized. *****
//=============================================================

void initialize_vortex_yee_pv ( PARAMS p, GRID *grid, double *Q )
{
  int i;                                           // Node index.
  double gamma = p.gamma;                          // Ratio of coefficients of heat.
  double gm1 = gamma - 1.0;
  double oogm1;                                    // One over (gamma-1)
  double x0,y0;                                    // Vortex center.
  double r2;                                       // Distance from center.
  double phi,sigma;                                // Parameter to control vortex strength.
  double upert,vpert,Tpert;                        // Isentropic vortex perturbations.
  double rho,u,v,P,T,E;                            // Control volume values of conserved variables.
  double *x = NULL;                                // Pointer to grid points.
  double *q = NULL;                                // Pointer to state values.
  double Mref,rhoref,Pref,Tref;
  double rhoD,PD;                                  // Dimensional quantities.
  double FC;     

  x0 = 5. + ( grid->citer * p.mintime * p.mach );   // Put the vortex center at the correct spot.
  y0 = 0.;

  oogm1 = 1./gm1;
  sigma = 4.0;
  phi = 1.0;

  Mref = 340.294065;
  Tref = 288.15;
  Pref = 101325.0;
  rhoref = 1.225;
  FC = Pref / pow(rhoref,gamma);

  for ( i=1; i <= grid->nn; i++ )
    {
      // Adjust the pointers.
      x = &(grid->x[i*NDIM]);
      q = &(Q[i*NUM_VAR]);

      // Initialize the variables.
      r2 = (x[0] - x0)*(x[0] - x0) + (x[1] - y0)*(x[1] - y0);

      upert = (-sigma)/(2.*M_PI) * ( x[1] - y0 ) * exp( phi*(1.-r2) );
      vpert = (sigma)/(2.*M_PI) * ( x[0] - x0 ) * exp( phi*(1.-r2) );
      Tpert = ( -1. * sigma*sigma * gm1 )/( 16. * phi * gamma * M_PI * M_PI ) * exp( 2. * phi * (1.-r2) );
      
      // Compute the new values for the control volumes.
      T = 1.0 + Tpert;
      rho = pow( ( T ) , oogm1 );
      rhoD = rho*rhoref;
      u = p.mach + upert;
      v = 0. + vpert;
      PD = pow( rhoD , gamma ) * FC;
      P = PD / (gamma * Pref );

      E = P/gm1 + 0.5*rho*(u*u + v*v);
      
      q[0] = rho;
      q[1] = rho*u;
      q[2] = rho*v;
      q[3] = E;
    }
  return;
}

//=============================================================
// 
//  initialize_vortex_yee_cv()
//
//  Sets flow variables for the Vortex convection problem at a
//  particular time. Used for comparing analytical and numerical
//  solutions. This is the cv averaged function.
//
//  THIS IS HARD SET FOR STRAIGHT BOUNDARIES.
//
//  Based on the code from UX. Domain is:
//                               0 < x < xmax
//                              -5 < y < 5
//
//  PARAMS p;                            // The parameters.
//  GRID *grid;                          // The grid structure.
//  double *Q;                           // Pointer to the analytical solution.
//  ******  Assumes Q has been initialized. *****
//=============================================================

void initialize_vortex_yee_cv ( PARAMS p, GRID *grid, double *Q )
{
  int n;                                 // Node loop counter.
  int i,j,k,e;                           // Loop counters.
  int v;                                 // Vector loop.
  int ind;                               // Component loop.
  int nodes[MAX_NUM_VERT];               // Element vertices.
  int node;                              // Node index.
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
  double XGP[NDIM];
  double QGP[NUM_VAR];
  double QTRI[NUM_VAR];
  double gp[7][3];                       // Gauss points for integration.
  double dA;                             // Differential area (for a triangle).
  double gamma = p.gamma;                          // Ratio of coefficients of heat.
  double gm1 = gamma - 1.0;
  double oogm1;                                    // One over (gamma-1)
  double x0,y0;                                    // Vortex center.
  double r2;                                       // Distance from center.
  double phi,sigma;                                // Parameter to control vortex strength.
  double upert,vpert,Tpert;                        // Isentropic vortex perturbations.
  double rho,U,V,P,T,E;                            // Control volume values of conserved variables.
  double Mref,rhoref,Pref,Tref;
  double rhoD,PD;                                  // Dimensional quantities.
  double FC;     

  x0 = 5. + ( grid->citer * p.mintime * p.mach );   // Put the vortex center at the correct spot.
  y0 = 0.;

  oogm1 = 1./gm1;
  sigma = 4.0;
  phi = 1.0;

  Mref = 340.294065;
  Tref = 288.15;
  Pref = 101325.0;
  rhoref = 1.225;
  FC = Pref / pow(rhoref,gamma);

  // Clean out memory because we are accumulating.
  for ( n=1; n <= grid->nn; n++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	Q[n*NUM_VAR+j] = 0.;
    } 
  
  // Now we calculate the area-averaged value of the initial condition using Gaussian quadrature
 
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
	  xc[0] /= ((double)num_vert);
	  xc[1] /= ((double)num_vert);

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
		      gp[v][ind] = (bary_coord[v][0] * grid->x[node*NDIM+ind] +
				    bary_coord[v][1] * xmidL[ind] +
				    bary_coord[v][2] * xc[ind]);
		    }
		  gp[v][2] = 0.;   // Always 0 since 2D vectors.
		}

	      // Now we loop over the Gauss points and get the function value
	      // and perform the numerical integration (the area integral of the triangle).
	      
	      // Do this for each variable.
		  
	      QTRI[0] = 0.; QTRI[1] = 0.; QTRI[2] = 0.; QTRI[3] = 0.;
		  
	      for ( i=0; i < 7; i++ )
		{
		  // Get the value at the point.
		  XGP[0] = gp[i][0];
		  XGP[1] = gp[i][1];
		      
		  // Initialize the variables.
		  r2 = (XGP[0] - x0)*(XGP[0] - x0) + (XGP[1] - y0)*(XGP[1] - y0);

		  upert = (-sigma)/(2.*M_PI) * ( XGP[1] - y0 ) * exp( phi*(1.-r2) );
		  vpert = (sigma)/(2.*M_PI) * ( XGP[0] - x0 ) * exp( phi*(1.-r2) );
		  Tpert = ( -1. * sigma*sigma * gm1 )/( 16. * phi * gamma * M_PI * M_PI ) * exp( 2. * phi * (1.-r2) );
		  
		  // Compute the new values for the control volumes.
		  T = 1.0 + Tpert;
		  rho = pow( ( T ) , oogm1 );
		  rhoD = rho*rhoref;
		  U = p.mach + upert;
		  V = 0. + vpert;
		  PD = pow( rhoD , gamma ) * FC;
		  P = PD / (gamma * Pref );
		  
		  E = P/gm1 + 0.5*rho*(U*U + V*V);
      
		  QGP[0] = rho;
		  QGP[1] = rho*U;
		  QGP[2] = rho*V;
		  QGP[3] = E;
		  
		  for ( v=0; v < NUM_VAR; v++ )
		    {
		      QTRI[v] += ( wg[i] * QGP[v] );
		    } 
		}
	      
	      for ( v=0; v < NUM_VAR; v++ )
		{
		      QTRI[v] *= dA;
		}
	      
	      // Accumulate the integral to the node.
	      for ( v=0; v < NUM_VAR; v++ )
		{
		  Q[node*NUM_VAR+v] += QTRI[v];
		}
	      
	      // Process the second triangle.
	      
	      // Get the area.
	      dA = 0.5*( v_middle[0]*v_right[1] - v_middle[1]*v_right[0] );
	      
	      // Initialize the Gauss points.
	      for ( v=0; v < 7; v++ )
		{
		  for ( ind=0; ind < NDIM; ind++ )
		    {
		      gp[v][ind] = (bary_coord[v][0] * grid->x[node*NDIM+ind] +
				    bary_coord[v][1] * xmidR[ind] +
				    bary_coord[v][2] * xc[ind]);
		    }
		  gp[v][2] = 0.;   // Always 0 since 2D vectors.
		}
	      
	      // Now we loop over the Gauss points and get the function value
	      // and perform the numerical integration (the area integral of the triangle).
	  
	      // Do this for each variable.
	  
	      QTRI[0] = 0.; QTRI[1] = 0.; QTRI[2] = 0.; QTRI[3] = 0.;
	      
	      for ( i=0; i < 7; i++ )
		{
		  // Get the value at the point.
		  XGP[0] = gp[i][0];
		  XGP[1] = gp[i][1];

		  // Initialize the variables.
		  r2 = (XGP[0] - x0)*(XGP[0] - x0) + (XGP[1] - y0)*(XGP[1] - y0);

		  upert = (-sigma)/(2.*M_PI) * ( XGP[1] - y0 ) * exp( phi*(1.-r2) );
		  vpert = (sigma)/(2.*M_PI) * ( XGP[0] - x0 ) * exp( phi*(1.-r2) );
		  Tpert = ( -1. * sigma*sigma * gm1 )/( 16. * phi * gamma * M_PI * M_PI ) * exp( 2. * phi * (1.-r2) );
		  
		  // Compute the new values for the control volumes.
		  T = 1.0 + Tpert;
		  rho = pow( ( T ) , oogm1 );
		  rhoD = rho*rhoref;
		  U = p.mach + upert;
		  V = 0. + vpert;
		  PD = pow( rhoD , gamma ) * FC;
		  P = PD / (gamma * Pref );
		  
		  E = P/gm1 + 0.5*rho*(U*U + V*V);
      
		  QGP[0] = rho;
		  QGP[1] = rho*U;
		  QGP[2] = rho*V;
		  QGP[3] = E;
		  
		  for ( v=0; v < NUM_VAR; v++ )
		    {
		      QTRI[v] += ( wg[i] * QGP[v] );
		    } 
		}
	      
	      for ( v=0; v < NUM_VAR; v++ )
		{
		  QTRI[v] *= dA;
		}
	      
	      // Accumulate the integral to the node.
	      for ( v=0; v < NUM_VAR; v++ )
		{
		  Q[node*NUM_VAR+v] += QTRI[v];
		}
      
	    }                                       // End element node loop.
  
	}                                           // End element loop.

    }                                               // End element type loop.

  // Finished the accumulation. Now divide by the control volume area.
  for ( i=1; i<= grid->nn; i++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  Q[i*NUM_VAR+j] /= grid->cv_area[i];
	}
    }

  return;
}


//=============================================================
// 
//  get_core_pressure()
//
//  Computes the core pressure at the theoretical center of the vortex.
//
//  Based on the code from UX. Domain is:
//                               0 < x < xmax
//                              -5 < y < 5
//
//  GRID *grid;                          // The grid structure.
//  PARAMS p;                            // The parameters.
//  double *time;                        // Current time. Returned.
//  double *core_pressure;               // Current core pressure. Returned.
//
//=============================================================

void get_core_pressure ( GRID *grid, PARAMS p, double *time, double *core_pressure )
{
  int i,j;                                         // Node index.
  double gamma = p.gamma;                          // Ratio of coefficients of heat.
  double gm1 = gamma - 1.0;
  double mach = p.mach;                            // Mach number.
  double x0;                                       // x location of vortex at time t = 0
  double y0;                                       // y location of vortex at time t = 0
  double r2;                                       // Radius squared.
  double beta;                                     // Vortex strength.
  double fac1,fac2,fac3;                           // Terms in the computation.
  double rho0;                                     // Core density.
  double u0,v0;                                    // Core velocity.
  double P0,T0,s0;                                 // Core pressure,temperature, and entropy.
  double T,u,v,P,rho;                              // Flow properties.
  double dx;                                       // Distance in x from center.
  double Q[NUM_VAR];
  double dist[NDIM];

  if ( grid->citer == 0 )  // Calculate the initial core pressure.
    {
      x0 = 5.;
      y0 = 0.;

      // Initialize the variables.
      r2 = 0.;
      beta = 5.0;
      fac1 = beta/(2.0*M_PI);
      fac2 = -0.5*gm1/gamma*fac1*fac1;
      fac3 = exp((1.0 - r2)/2.0);
      
      rho0 = 1.0;
      u0 = mach*1.0;          // Flow in only the x direction at reference velocity of 1.
      v0 = 0.;
      P0 = 1.0/gamma;
      T0 = P0/rho0;
      s0 = P0/pow(rho0,gamma);
      
      T = T0 + fac2*fac3*fac3;
      u = u0 - (fac1*fac3)*(y0 - y0);
      v = v0 + (fac1*fac3)*(x0 - x0);
      P = pow(pow(T, gamma)/s0,1.0/gm1);
      rho = P/T;
      
      *time = 0.;
      *core_pressure = P;

      return;
    }

  // Where should the core be at?
  x0 = 5. + ( grid->citer * p.mintime * p.mach );
  y0 = 0.;

  dx = 1000.0;  // Set it big for now
  j = 0;        // Fictitious node.

  // Now we have to find the closest grid point.
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( fabs(grid->x[i*NDIM+1]) > 10.0E-12 )  // Not on the y=0 line.
	{
	  continue;
	}

      if ( dx > fabs( grid->x[i*NDIM+0] - x0 ) )  // Found a closer point.
	{
	  dx = fabs( grid->x[i*NDIM+0] - x0 );
	  j = i;
	}
      
    }

  // Did I find a point?
  if ( j == 0 )
    {
      printf("CRITICAL ERROR: get_core_pressure() failed to find a closest point!\n");
      exit(0);
    }

  printf("The core is at %f  %f , and the closest node is %d  @  %f  %f\n",(float)x0,(float)y0,j,(float)grid->x[j*NDIM+0],(float)grid->x[j*NDIM+1]);

  // Now extrapolate to get the pressure at the center.

  // set the distance.
  dist[0] = x0 - grid->x[j*NDIM+0];
  dist[1] = y0 - grid->x[j*NDIM+1];

  // Now node j is the magic node. If RECON_PRIM was set then the derivatives and point values are in primitive form but grid->nQ,Q,pQ are in conserved form.
  // reconstruct functions extract from grid->nQ so let's convert it to primitive.

  if ( RECON_PRIM )
    ConvertConservedToPrimitive ( p.gamma , &(grid->nQ[j*NUM_VAR]) );
  
  generic_reconstruct ( j, Q, dist, grid, p );

  // Now I need to get Pressure out of the conserved variables.
  
  if ( RECON_PRIM )
    {
      P = Q[3];
    }
  else
    {
      rho = Q[0];
      u = Q[1] / rho;
      v = Q[2] / rho;
      P = (gamma-1.)*(Q[3] - 0.5*rho*(u*u + v*v));
    }

  // Undo any conversion if necessary.
  if ( RECON_PRIM )
    ConvertPrimitiveToConserved ( p.gamma , &(grid->nQ[j*NUM_VAR]) );

  *time = ( grid->citer * p.mintime );
  *core_pressure = P;

  return;
}

//=============================================================
// 
//  get_core_pressure_yee()
//
//  Computes the core pressure at the theoretical center of the vortex.
//
//  Based on the code from UX. Domain is:
//                               0 < x < xmax
//                              -5 < y < 5
//
//  GRID *grid;                          // The grid structure.
//  PARAMS p;                            // The parameters.
//  double *time;                        // Current time. Returned.
//  double *core_pressure;               // Current core pressure. Returned.
//
//=============================================================

void get_core_pressure_yee ( GRID *grid, PARAMS p, double *time, double *core_pressure )
{
  int i,j;                                         // Node index.
  double gamma = p.gamma;                          // Ratio of coefficients of heat.
  double gm1 = gamma - 1.0;
  double oogm1;                                    // One over (gamma-1)
  double x0,y0;                                    // Vortex center.
  double r2;                                       // Distance from center.
  double phi,sigma;                                // Parameter to control vortex strength.
  double Tpert;                                    // Isentropic vortex perturbations.
  double rho,u,v,P,T;                              // Control volume values of conserved variables.
  double dx;                                       // Distance in x from center.
  double Mref,rhoref,Pref,Tref;
  double rhoD,PD;                                  // Dimensional quantities.
  double FC;
  double Q[NUM_VAR];
  double dist[NDIM];

  Mref = 340.294065;
  Tref = 288.15;
  Pref = 101325.0;
  rhoref = 1.225;
  FC = Pref / pow(rhoref,gamma);

  if ( grid->citer == 0 )  // Calculate the initial core pressure.
    {
      x0 = 5.;
      y0 = 0.;

      // Initialize the variables.
      r2 = 0.;
      oogm1 = 1./gm1;
      sigma = 4.0;
      phi = 1.0;

      Tpert = ( -1. * sigma*sigma * gm1 )/( 16. * phi * gamma * M_PI * M_PI ) * exp( 2. * phi * (1.-r2) );
      
      T = 1.0 + Tpert;
      rho = pow( ( T ) , oogm1 );
      
      rhoD = rho*rhoref;
      PD = pow( rhoD , gamma ) * FC;
      P = PD / (gamma * Pref );
      
      *time = 0.;
      *core_pressure = P;

      return;
    }

  // Where should the core be at?
  x0 = 5. + ( grid->citer * p.mintime * p.mach );
  y0 = 0.;

  dx = 1000.0;  // Set it big for now
  j = 0;        // Fictitious node.

  // Now we have to find the closest grid point.
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( fabs(grid->x[i*NDIM+1]) > 10.0E-12 )  // Not on the y=0 line.
	{
	  continue;
	}

      if ( dx > fabs( grid->x[i*NDIM+0] - x0 ) )  // Found a closer point.
	{
	  dx = fabs( grid->x[i*NDIM+0] - x0 );
	  j = i;
	}
      
    }

  // Did I find a point?
  if ( j == 0 )
    {
      printf("CRITICAL ERROR: get_core_pressure() failed to find a closest point!\n");
      exit(0);
    }

  printf("The core is at %f  %f , and the closest node is %d  @  %f  %f\n",(float)x0,(float)y0,j,(float)grid->x[j*NDIM+0],(float)grid->x[j*NDIM+1]);

  // Now extrapolate to get the pressure at the center.

  // set the distance.
  dist[0] = x0 - grid->x[j*NDIM+0];
  dist[1] = y0 - grid->x[j*NDIM+1];

  // Now node j is the magic node. If RECON_PRIM was set then the derivatives and point values are in primitive form but grid->nQ,Q,pQ are in conserved form.
  // reconstruct functions extract from grid->nQ so let's convert it to primitive.

  if ( RECON_PRIM )
    ConvertConservedToPrimitive ( p.gamma , &(grid->nQ[j*NUM_VAR]) );
  
  generic_reconstruct ( j, Q, dist, grid, p );

  // Now I need to get Pressure out of the conserved variables.

  if ( RECON_PRIM )
    {
      P = Q[3];
    }
  else
    {
      rho = Q[0];
      u = Q[1] / rho;
      v = Q[2] / rho;
      P = (gamma-1.)*(Q[3] - 0.5*rho*(u*u + v*v));
    }

  // Undo any conversion if necessary.
  if ( RECON_PRIM )
    ConvertPrimitiveToConserved ( p.gamma , &(grid->nQ[j*NUM_VAR]) );

  *time = ( grid->citer * p.mintime );
  *core_pressure = P;

  return;
}


//=============================================================
// 
//  write_vortex_pressure_centerline()
//
//  Writes out the pressure along the centerline of the vortex
//  convection problem.
//
//  Based on the code from UX. Domain is:
//                               0 < x < xmax
//                              -5 < y < 5
//
//  GRID *grid;                          // The grid structure.
//  PARAMS p;                            // The parameters.
//  double *time;                        // Current time. Returned.
//
//=============================================================

void write_vortex_pressure_centerline ( GRID *grid, PARAMS p, double *time )
{
  int i;                                           // Node index.
  double gamma = p.gamma;                          // Ratio of coefficients of heat.
  double gm1 = gamma - 1.0;
  double x0;                                       // x location of vortex at time t = 0
  double r2;                                       // Radius squared.
  double beta;                                     // Vortex strength.
  double fac1,fac2,fac3;                           // Terms in the computation.
  double rho0;                                     // Core density.
  double P0,T0,s0;                                 // Core pressure,temperature, and entropy.
  double T,u,v,P,rho;                              // Flow properties.
  double E;
  double dx;                                       // Distance in x from center.
  FILE *vpc_ex = NULL;
  FILE *vpc_num = NULL;
  int num_points;
  double xmax = 1.0E-10;
  double xmin = 1.0E+10;
  double xi;

  // Open the files for writing.
  vpc_ex = fopen("vortex_pressure_centerline_analytical.dat","w");
  vpc_num = fopen("vortex_pressure_centerline_computed.dat","w");

  // Where should the core be at?
  x0 = 5. + ( grid->citer * p.mintime * p.mach);

  num_points = 0;

  // Write the computed solution first.
  for ( i=1; i < grid->nn; i++ )
    {
      // Skip points not at y == 0.
      if ( fabs(grid->x[i*NDIM+1]) > 10.0E-12 )  // Not on the y=0 line.
	{
	  continue;
	}

      xmin = MIN(xmin,grid->x[i*NDIM+0]);
      xmax = MAX(xmax,grid->x[i*NDIM+0]);

      // Grab the variables.
      
      // This is not needed since grid->Q is in conserved form when this function is called from main.C .
      //if ( RECON_PRIM )
      //{
      //  rho = grid->Q[i*NUM_VAR+0];
      //  u =  grid->Q[i*NUM_VAR+1];
      //  v =  grid->Q[i*NUM_VAR+2];
      //  P =  grid->Q[i*NUM_VAR+3];
      //}

      rho = grid->Q[i*NUM_VAR+0];
      u =  grid->Q[i*NUM_VAR+1] / rho;
      v =  grid->Q[i*NUM_VAR+2] / rho;
      E =  grid->Q[i*NUM_VAR+3];

      P = (gamma-1.)*(E - 0.5*rho*(u*u + v*v));
      
      fprintf(vpc_num,"%.15E     %.15E\n",grid->x[i*NDIM+0],P);

      num_points++;
    }
  
  fclose(vpc_num);
  
  // Now do the analytical solution.
  
  beta = 5.0;
  fac1 = beta/(2.0*M_PI);
  fac2 = -0.5*gm1/gamma*fac1*fac1;
  rho0 = 1.0;
  P0 = 1.0/gamma;
  T0 = P0/rho0;
  s0 = P0/pow(rho0,gamma);
  
  dx = (xmax - xmin) / ( (double)(num_points-1));
  
  for ( i=0; i < num_points; i++ )
    {
      xi = xmin + i*dx;
      
      r2 = (xi-x0)*(xi-x0);  // y-y0 = 0.
      fac3 = exp((1.0 - r2)/2.0);
      
      T = T0 + fac2*fac3*fac3;
      P = pow(pow(T, gamma)/s0,1.0/gm1);
      
      fprintf(vpc_ex,"%.15E     %.15E\n",xi,P);
    }
  
  fclose(vpc_ex);
			 

  return;
}


//=============================================================
// 
//  write_vortex_pressure_centerline_yee()
//
//  Writes out the pressure along the centerline of the vortex
//  convection problem.
//
//  Based on the code from UX. Domain is:
//                               0 < x < xmax
//                              -5 < y < 5
//
//  GRID *grid;                          // The grid structure.
//  PARAMS p;                            // The parameters.
//  double *time;                        // Current time. Returned.
//
//=============================================================

void write_vortex_pressure_centerline_yee ( GRID *grid, PARAMS p, double *time )
{
  int i;                                           // Node index.
  double gamma = p.gamma;                          // Ratio of coefficients of heat.
  double gm1 = gamma - 1.0;
  double x0;                                       // x location of vortex at time t = 0
  double r2;                                       // Radius squared.
  double oogm1;                                    // One over (gamma-1)
  double phi,sigma;                                // Parameter to control vortex strength.
  double upert,vpert,Tpert;                        // Isentropic vortex perturbations.
  double rho,u,v,P,T,E;                            // Control volume values of conserved variables.
  FILE *vpc_ex = NULL;
  FILE *vpc_num = NULL;
  int num_points;
  double dx;
  double xmax = 1.0E-10;
  double xmin = 1.0E+10;
  double xi;
  double Mref,rhoref,Pref,Tref;
  double rhoD,PD;                                  // Dimensional quantities.
  double FC; 

  Mref = 340.294065;
  Tref = 288.15;
  Pref = 101325.0;
  rhoref = 1.225;
  FC = Pref / pow(rhoref,gamma);

  // Open the files for writing.
  vpc_ex = fopen("vortex_pressure_centerline_analytical_yee.dat","w");
  vpc_num = fopen("vortex_pressure_centerline_computed_yee.dat","w");

  // Where should the core be at?
  x0 = 5. + ( grid->citer * p.mintime *p.mach );

  oogm1 = 1./gm1;
  sigma = 4.0;
  phi = 1.0;

  num_points = 0;

  // Write the computed solution first.
  for ( i=1; i < grid->nn; i++ )
    {
      // Skip points not at y == 0.
      if ( fabs(grid->x[i*NDIM+1]) > 10.0E-12 )  // Not on the y=0 line.
	{
	  continue;
	}

      xmin = MIN(xmin,grid->x[i*NDIM+0]);
      xmax = MAX(xmax,grid->x[i*NDIM+0]);

      // Grab the variables.

      // This isn't needed since grid->Q is in conserved form when this function is called from main.C .
      //if ( RECON_PRIM )
      //{
      //  rho = grid->Q[i*NUM_VAR+0];
      //  u =  grid->Q[i*NUM_VAR+1];
      //  v =  grid->Q[i*NUM_VAR+2];
      //  P =  grid->Q[i*NUM_VAR+3];
      //}

      rho = grid->Q[i*NUM_VAR+0];
      u =  grid->Q[i*NUM_VAR+1] / rho;
      v =  grid->Q[i*NUM_VAR+2] / rho;
      E =  grid->Q[i*NUM_VAR+3];
      
      P = (gamma-1.)*(E - 0.5*rho*(u*u + v*v));
      
      fprintf(vpc_num,"%.15E     %.15E\n",grid->x[i*NDIM+0],P);

      num_points++;
    }
  
  fclose(vpc_num);
  
  // Now do the analytical solution.
  
  dx = (xmax - xmin) / ( (double)(num_points-1));
  
  for ( i=0; i < num_points; i++ )
    {
      xi = xmin + i*dx;
      
      r2 = (xi-x0)*(xi-x0);  // y-y0 = 0.

      Tpert = ( -1. * sigma*sigma * gm1 )/( 16. * phi * gamma * M_PI * M_PI ) * exp( 2. * phi * (1.-r2) );
      T = 1.0 + Tpert;

      // Compute the new values for the control volumes.
      rho = pow( ( T ) , oogm1 );
      rhoD = rho*rhoref;

      PD = pow( rhoD , gamma ) * FC;

      P = PD / (gamma * Pref );
      
      fprintf(vpc_ex,"%.15E     %.15E\n",xi,P);
    }
  
  fclose(vpc_ex);
  
  return;
}
