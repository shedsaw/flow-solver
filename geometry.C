//=============================================================
// 
//  geometry.C
//  
//  Functions needed building geometry data from the grid.
//
//  Written by - Shane Sawyer
//
//=============================================================

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "grid.h"
#include "geometry.h"
#include "mesh_utils.h"
#include "Defined_Variables.h"
#include "curved_boundaries.h"


//=============================================================
// 
//  Find_Element_Centroids()
//
//  Calculates the centroids for each element in the grid.
//  
//  GRID *grid                        // The grid.
//
//=============================================================

void Find_Element_Centroids ( GRID *grid )
{
  int i, j, k;                           // Loop counters.
  int numelem;                           // Total number of elements.
  int gelem;                             // Global element number.
  int node[MAX_NUM_VERT];                // Element vertices.
  double xc, yc;                         // Values at the centroid.
  int num_vert;                          // Number of vertices for an element.
  
  numelem = 0;
  for ( i=Tri; i <= Quad; i++ )
    {
      numelem += grid->num_elem[i];
    }

  // Make sure space has been allocated.
  if ( grid->el_cent == NULL )
    {
      grid->el_cent = (double*)malloc((numelem+1)*NDIM*sizeof(double));
      if ( grid->el_cent == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'el_cent'.\n"); exit(1); }
    }

  // Loop over the element types and elements to get the algebraic centroid.
  for ( i=Tri; i <= Quad; i++ )                                                       // Element types loop.
    {
      // Find the number of nodes and edges for this element type.
      num_vert = NumberNodesForElement(i);

      for ( j=1; j <= grid->num_elem[i]; j++ )                                        // Element loop.
	{
	  // Retrieve the element node ids.
	  for ( k=0; k < num_vert; k++ )              // Element node extraction loop.
	    {
	      node[k] = grid->c2n[i][j*num_vert+k];
	    }                                         // End element node extraction loop.

	  xc = 0.;  yc = 0.;
	  for ( k=0; k < num_vert; k++ )              // Vertex sum loop.
	    {
	      xc += grid->x[NDIM*node[k]];
	      yc += grid->x[NDIM*node[k]+1];
	    }                                         // End vertes sum loop.
	  xc /= (num_vert*1.);
	  yc /= (num_vert*1.);

	  // Get the global element id.
	  gelem = LocalToGlobalID ( grid->num_elem,
				    j,
				    i );

	  // Store the centroid.
	  grid->el_cent[gelem*NDIM] = xc;
	  grid->el_cent[gelem*NDIM+1] = yc;
	}                                                                           // End element loop.
    }                                                                              // End element type loop.

  // Debug code to print to screen.
  #if DEBUG
  printf("EL_CENT DEBUG:\n");
  for ( i=Tri; i <= Quad; i++ )
    {
      printf("Element Type %d:\n",i);
      for ( j=1; j <= grid->num_elem[i]; j++ )
	{
	  gelem = LocalToGlobalID ( grid->num_elem,j,i );
	  printf(" Local Elem %d , Global Elem %d, xc = %lf  yc = %lf \n",j,
		 gelem,grid->el_cent[gelem*NDIM],grid->el_cent[gelem*NDIM+1]);
	}
      printf("\n");
    }
  #endif

  return;
}



//=============================================================
// 
//  Build_Edge_Data()
//
//  Calculates the normal vector for the edges. This also fills
//  in the subedge data.
//  
//  GRID *grid                        // The grid.
//
//=============================================================

void Build_Edge_Data ( GRID *grid )
{
  int i,j,k,n;                           // Loop counters.
  int tsubedges;                         // Total number of subedges.
  int iedge;                             // Global edge.
  int isubedge;                          // Subedge being worked on.
  int gelem;                             // Global element.
  int start_idx, end_idx;                // Array bounds in esn.
  int nodeL, nodeR;                      // Edge/element nodes.
  int num_vert;                          // Number of vertices for an element.
  int *edgehits = NULL;                  // Number of elements that border the edge, 1 or 2.
  double xc, yc;                         // Values at the centroid.
  double xmid, ymid;                     // Values at face midpoints.
  double nx,ny;                          // Normal vector components.
  double mag;                            // Magnitude of a vector.

  iedge = 0;
  isubedge = 0;

  edgehits = (int*)malloc((grid->nedges+1)*sizeof(int));
  if ( edgehits == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'edgehits'.\n"); exit(1); }

  // Set the edge hits to 0 initially.
  for ( i=1; i <= grid->nedges; i++ )
    {
      edgehits[i] = 0;
    }

  // Get the total number of subedges.
  tsubedges = grid->num_elem[Tri]*3 + grid->num_elem[Quad]*4;
  
  // Now we can allocate space. Assume that edges has been built!
  
  if ( grid->subedges == NULL )
    {
      grid->subedges = (int*)malloc((tsubedges+1)*2*sizeof(int));
      if ( grid->subedges == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'subedges'.\n"); exit(1); }
    }

  if ( grid->xn_edges == NULL )
    {
      grid->xn_edges = (double*)malloc((grid->nedges+1)*3*sizeof(double));
      if ( grid->xn_edges == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'xn_edges'.\n"); exit(1); }
    }

  if ( grid->xn_subedges == NULL )
    {
      grid->xn_subedges = (double*)malloc((tsubedges+1)*3*sizeof(double));
      if ( grid->xn_subedges == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'xn_subedges'.\n"); exit(1); }
    }

  if ( grid->xm_subedges == NULL )
    {
      grid->xm_subedges = (double*)malloc((tsubedges+1)*NDIM*sizeof(double));
      if ( grid->xm_subedges == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'xm_subedges'.\n"); exit(1); }
    }

  // Clean out memory.
  for ( i=1; i <= grid->nedges; i++ )
    {
      for ( j=0; j < 3; j++ )
	{
	  grid->xn_edges[i*3+j] = 0.;
	}
    }

  // Clean out memory.
  for ( i=1; i <= tsubedges; i++ )
    {
      for ( j=0; j < 3; j++ )
	{
	  grid->xn_subedges[i*3+j] = 0.;
	}
    }

  for ( i=1; i <= tsubedges; i++ )
    {
      for ( j=0; j < NDIM; j++ )
	{
	  grid->xm_subedges[i*NDIM+j] = 0.;
	}
    }

  // I think the best way to get everything I need is to loop over the elements and get my subedges (median dual)
  // first then construct the edges data from that.

  // Loop over the element types and elements to get the median dual pieces.
  for ( i=Tri; i <= Quad; i++ )                                                       // Element types loop.
    {
      // Find the number of nodes and edges for this element type.
      num_vert = NumberNodesForElement(i);

      for ( j=1; j <= grid->num_elem[i]; j++ )                                        // Element loop.
	{
	  // Get the global element id.
	  gelem = LocalToGlobalID ( grid->num_elem,
				    j,
				    i );

	  // Get the centroid.
	  xc = grid->el_cent[gelem*NDIM];
          yc = grid->el_cent[gelem*NDIM+1];

	  // Loop over the element nodes (edges).
	  for ( k=0; k < num_vert; k++ )                                            // Element node loop.
	    {
	      isubedge++;

	      nodeL = k;
	      nodeR = (k+1)%num_vert;

	      // Get the global node numbers.
	      nodeL = grid->c2n[i][j*num_vert + nodeL];
	      nodeR = grid->c2n[i][j*num_vert + nodeR];
	      
	      // The edge is defined as the left node having the lower index and points towards the
	      // right node away from the left node. So lets swap the nodes now if they need it.
	      Swap(&nodeL,&nodeR);

	      // Calculate the edge midpoint.
	      xmid = 0.5*( grid->x[nodeL*NDIM+0] + grid->x[nodeR*NDIM+0] );
	      ymid = 0.5*( grid->x[nodeL*NDIM+1] + grid->x[nodeR*NDIM+1] );

	      // Now we need to find the global edge that connects the two nodes.
	      start_idx = grid->nnsn[nodeL];
	      end_idx = grid->nnsn[nodeL+1];

	      for ( n=start_idx; n < end_idx; n++ )        // Surrounding edge loop.
		{
		  iedge = grid->esn[n];

		  if ( grid->edges[iedge*2 + 1] == nodeR )
		    {  break; }  // Found the right edge.

		}                                          // End surrounding edge loop.

	      nx = yc - ymid;
	      ny = xmid - xc;

	      // If I've swapped I need to adjust the normal vector.
	      if ( nodeR == grid->c2n[i][j*num_vert + k] )
		{  nx = -nx;  ny = -ny; }

	      // Store the data.
	      grid->subedges[isubedge*2] = iedge;
	      grid->subedges[isubedge*2+1] = gelem;

	      grid->xn_subedges[isubedge*3] = nx;
	      grid->xn_subedges[isubedge*3+1] = ny;

	      grid->xm_subedges[isubedge*NDIM+0] = xmid;
	      grid->xm_subedges[isubedge*NDIM+1] = ymid;

	      grid->xn_edges[iedge*3] += nx;
	      grid->xn_edges[iedge*3+1] += ny;

	      edgehits[iedge]++;

	    }                                                                       // End element node loop. 
	}                                                                           // End element loop.
    }                                                                               // End element type loop.

  // Sanity check.
  if ( isubedge != tsubedges )
    {
      printf("FATAL ERROR: In Build_Edges_Data() a topological error was encountered.\n");
      printf("  isubedge = %d , tsubedges = %d\n",isubedge,tsubedges);
      exit(0);
    }

  grid->nsubedges = isubedge;

  // Sanity Check.
  // Every edge should have been hit at least once, but not more than twice.
  for ( i=1; i <= grid->nedges; i++ )
    {
      if ( edgehits[i] < 1 || edgehits[i] > 2 )
	{
	  printf("FATAL ERROR: In Build_Edges_Data() a topological error was encountered.\n");
	  printf("  Listing Bad Edges:\n");
	  for ( j=1; j <= grid->nedges; j++ )
	    {
	      if ( edgehits[i] < 1 || edgehits[i] > 2 )
		{
		  printf("   Edge %d - %d hits!\n",j,edgehits[j]);
		}
	    }
	  exit(0);
	}
    }

  // Now I need to average the data in xn_edges. - Not sure this is correct!
  /*
  for ( i=1; i <= grid->nedges; i++ )
    {
      grid->xn_edges[i*3] /= ( (double)edgehits[i] );
      grid->xn_edges[i*3+1] /= ( (double)edgehits[i] );
    }
  */

  // Now we normalize the the normal vectors and store the length.
  for ( i=1; i <= grid->nedges; i++ )
    {
      mag = sqrt( grid->xn_edges[i*3]*grid->xn_edges[i*3] +
		  grid->xn_edges[i*3+1]*grid->xn_edges[i*3+1]);

      if ( mag > 1.0E-20 )
	{
	  grid->xn_edges[i*3] /= mag;
	  grid->xn_edges[i*3+1] /= mag;
	}
      grid->xn_edges[i*3+2] = mag;
    }

  for ( i=1; i <= grid->nsubedges; i++ )
    {
      mag = sqrt( grid->xn_subedges[i*3]*grid->xn_subedges[i*3] +
		  grid->xn_subedges[i*3+1]*grid->xn_subedges[i*3+1]);

      if ( mag > 1.0E-20 )
	{
	  grid->xn_subedges[i*3] /= mag;
	  grid->xn_subedges[i*3+1] /= mag;
	}
      grid->xn_subedges[i*3+2] = mag;
    }
  

  freenull(edgehits);

  // Debug code to print to screen.
  #if DEBUG
  printf("XN_EDGES DEBUG:\n");
  for ( i=1; i <= grid->nedges; i++ )
    {
      printf("Edge %d:\n",i);
      printf("  nx = %lf  ny = %lf  length = %lf\n",grid->xn_edges[i*3],
	     grid->xn_edges[i*3+1],grid->xn_edges[i*3+2]);
    }
  #endif

  #if DEBUG
  printf("SUBEDGES DEBUG:\n");
  for ( i=1; i <= grid->nsubedges; i++ )
    {
      printf("Subedge %d:\n",i);
      printf("  iedge = %d , gelem = %d\n",grid->subedges[i*2],
	     grid->subedges[i*2+1]);
    }
  #endif

  #if DEBUG
  printf("XN_SUBEDGES DEBUG:\n");
  for ( i=1; i <= grid->nsubedges; i++ )
    {
      printf("Subedge %d:\n",i);
      printf("  nx = %lf  ny = %lf  length = %lf\n",grid->xn_subedges[i*3],
	     grid->xn_subedges[i*3+1],grid->xn_subedges[i*3+2]);
    }
  #endif

  #if DEBUG
  printf("XM_SUBEDGES DEBUG:\n");
  for ( i=1; i <= grid->nsubedges; i++ )
    {
      printf("Subedge %d:\n",i);
      printf("  xmid = %.15E  ymid = %.15E\n",grid->xm_subedges[i*NDIM+0],
	                                      grid->xm_subedges[i*NDIM+1]);
    }
  #endif
  
  return;
}




//=============================================================
// 
//  Build_Boundary_Edge_Data()
//
//  Calculates the normal vector for the boundary edges. This also fills
//  in the bedge data, xn_bedge data, and constructs the ghost nodes for the grid.
//  
//  GRID *grid                        // The grid.
//
//=============================================================

void Build_Boundary_Edge_Data ( GRID *grid )
{
  int i,j,k;                             // Loop counters.
  int b,s;                               // Loop counters.
  int tbedges;                           // Total number of boundary edges.
  int ibedge;                            // Boundary edge being worked on.
  int gelem;                             // Global element.
  int n0,n1;                             // Boundary edge nodes.
  int num_vert;                          // Number of vertices for an element.
  double xmid, ymid;                     // Values at face midpoints.
  double nx,ny;                          // Normal vector components.
  double mag;                            // Magnitude of a vector.
  double *temp = NULL;                   // Temporary pointer.
  int *itemp = NULL;                     // Temporary integer array pointer.
  int current_ghost_node = 1;
  int bc;
  double xL[2],xR[2],xM[2];
  int node, ghostnode;

  // Count the total number of boundary edges.
  tbedges = 0;
  for ( b=1; b <= grid->nb; b++ )
    {
      tbedges += grid->nbs[b];
    }

  tbedges = tbedges*2;
  grid->nbedges = tbedges;

  // This also gives me the number of ghost nodes present in the mesh.
  grid->nn_ghost = tbedges;

  // Adjust space for the ghost nodes in the x,Q,and nQ arrays. Assumes these arrays have already
  // been allocated (as they should have been when reading in the grid).
  
  itemp = (int*)realloc( (void*)grid->node_state , (grid->nn + grid->nn_ghost + 1)*sizeof(int) );

  if ( itemp == NULL )
    {
      printf("MEMORY ERROR: COULD NOT REALLOCATE 'node_state' while adjusting for ghost nodes.\n");
      exit(1);
    }
  else
    {
      grid->node_state = itemp;
    }

  itemp = (int*)realloc( (void*)grid->node_order , (grid->nn + grid->nn_ghost + 1)*sizeof(int) );

  if ( itemp == NULL )
    {
      printf("MEMORY ERROR: COULD NOT REALLOCATE 'node_order' while adjusting for ghost nodes.\n");
      exit(1);
    }
  else
    {
      grid->node_order = itemp;
    }
  

  temp = (double*)realloc( (void*)grid->x , (grid->nn + grid->nn_ghost + 1)*NDIM*sizeof(double) );

  if ( temp == NULL )
    {
      printf("MEMORY ERROR: COULD NOT REALLOCATE 'x' while adjusting for ghost nodes.\n");
      exit(1);
    }
  else
    {
      grid->x = temp;
    }
  
  temp = (double*)realloc( (void*)grid->Q , (grid->nn + grid->nn_ghost + 1)*NUM_VAR*sizeof(double) );

  if ( temp == NULL )
    {
      printf("MEMORY ERROR: COULD NOT REALLOCATE 'Q' while adjusting for ghost nodes.\n");
      exit(1);
    }
  else
    {
      grid->Q = temp;
    }

  temp = (double*)realloc( (void*)grid->nQ , (grid->nn + grid->nn_ghost + 1)*NUM_VAR*sizeof(double) );

  if ( temp == NULL )
    {
      printf("MEMORY ERROR: COULD NOT REALLOCATE 'nQ' while adjusting for ghost nodes.\n");
      exit(1);
    }
  else
    {
      grid->nQ = temp;
    }

  temp = (double*)realloc( (void*)grid->pQ , (grid->nn + grid->nn_ghost + 1)*NUM_VAR*sizeof(double) );

  if ( temp == NULL )
    {
      printf("MEMORY ERROR: COULD NOT REALLOCATE 'pQ' while adjusting for ghost nodes.\n");
      exit(1);
    }
  else
    {
      grid->pQ = temp;
    }

  // Allocate space for the bedges.
  if ( grid->bedges == NULL )
    {
      grid->bedges = (int*)malloc((tbedges+1)*5*sizeof(int));
      if ( grid->bedges == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'bedges'.\n"); exit(1); }
    }

  if ( grid->xn_bedges == NULL )
    {
      grid->xn_bedges = (double*)malloc((tbedges+1)*3*sizeof(double));
      if ( grid->xn_bedges == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'xn_bedges'.\n"); exit(1); }
    }

  // Now lets fill in the node_state array. First assume everything is marked as interior.
  for ( i=1; i <= grid->nn; i++ )
    {
      grid->node_state[i] = INTERIOR;
    }

  // Loop over the boundary segments and create the bedge data.
  ibedge = 1;  current_ghost_node = 1;
  for ( b=1; b <= grid->nb; b++ )
    {                                                       // Boundary loop.
      for ( s=1; s <= grid->nbs[b]; s++ )
	{                                  
	  n0 = grid->bs[b][s][0];
	  n1 = grid->bs[b][s][1];
	  bc = grid->bbc[b];

	  // Adjust the node state for the boundary nodes on this segment.
	  grid->node_state[n0] = BOUNDARY;
	  grid->node_state[n1] = BOUNDARY;

	  if ( bc <= 10 )
	    {
	      xmid = (grid->x[NDIM*n0]+grid->x[NDIM*n1])/2.;
	      ymid = (grid->x[NDIM*n0+1]+grid->x[NDIM*n1+1])/2.;
	  
	      nx = grid->x[NDIM*n1+1] - ymid;
	      ny = xmid - grid->x[NDIM*n1+0];

	      //Store the data.
	      grid->bedges[ibedge*5+0] = n0;
	      grid->bedges[ibedge*5+1] = current_ghost_node;
	      grid->bedges[ibedge*5+3] = b;
	      grid->bedges[ibedge*5+4] = s;

	      grid->xn_bedges[ibedge*3+0] = nx;
	      grid->xn_bedges[ibedge*3+1] = ny;

	      grid->x[NDIM*(grid->nn + current_ghost_node)+0] = ( grid->x[NDIM*n0] + xmid )/2.;
	      grid->x[NDIM*(grid->nn + current_ghost_node)+1] = ( grid->x[NDIM*n0+1] + ymid )/2.;

	      // HACK!
	      if ( 0 )
		{
		  if ( b == 2 )
		    {
		      xL[0] = grid->x[NDIM*n0+0];
		      xL[1] = grid->x[NDIM*n0+1];
		      
		      xR[0] = grid->x[NDIM*n1+0];
		      xR[1] = grid->x[NDIM*n1+1];
		      
		      curved_boundary_midpoint ( grid, 14, xL, xR, xM );
		      
		      grid->x[NDIM*(grid->nn + current_ghost_node)+0] = xM[0];
		      grid->x[NDIM*(grid->nn + current_ghost_node)+1] = xM[1];
		    }

		  if ( b == 4 )
		    {
		      xL[0] = grid->x[NDIM*n0+0];
		      xL[1] = grid->x[NDIM*n0+1];
		      
		      xR[0] = grid->x[NDIM*n1+0];
		      xR[1] = grid->x[NDIM*n1+1];
		      
		      curved_boundary_midpoint ( grid, 16, xL, xR, xM );
		      
		      grid->x[NDIM*(grid->nn + current_ghost_node)+0] = xM[0];
		      grid->x[NDIM*(grid->nn + current_ghost_node)+1] = xM[1];
		    }
		}
	      
	      // Mark the ghost node as such.
	      grid->node_state[grid->nn + current_ghost_node] = GHOST;
	  
	      ibedge++;
	      current_ghost_node++;

	      //Store the data.
	      grid->bedges[ibedge*5+0] = n1;
	      grid->bedges[ibedge*5+1] = current_ghost_node;
	      grid->bedges[ibedge*5+3] = b;
	      grid->bedges[ibedge*5+4] = s;
	      
	      grid->xn_bedges[ibedge*3] = nx;
	      grid->xn_bedges[ibedge*3+1] = ny;
	      
	      grid->x[NDIM*(grid->nn + current_ghost_node)+0] = ( grid->x[NDIM*n1] + xmid )/2.;
	      grid->x[NDIM*(grid->nn + current_ghost_node)+1] = ( grid->x[NDIM*n1+1] + ymid )/2.;

	      if ( 0 )
		{
		  if ( b == 2 || b == 4 )
		    {
		      grid->x[NDIM*(grid->nn + current_ghost_node)+0] = xM[0];
		      grid->x[NDIM*(grid->nn + current_ghost_node)+1] = xM[1];
		    }
		}

	      // Mark the ghost node as such.
	      grid->node_state[grid->nn + current_ghost_node] = GHOST;
	  
	      ibedge++;
	      current_ghost_node++;
	    }
	  else   // Do curved boundaries.
	    {
	      xL[0] = grid->x[NDIM*n0+0];
	      xL[1] = grid->x[NDIM*n0+1];
	      
	      xR[0] = grid->x[NDIM*n1+0];
	      xR[1] = grid->x[NDIM*n1+1];

	      //printf("Preparing to call curved_boundary_midpoint() with left node %d and right node %d.\n",n0,n1);
	      //fflush(stdout);
	      
	      // Get the midpoint.
	      curved_boundary_midpoint ( grid, bc, xL, xR, xM );

	      //Store the data.
	      grid->bedges[ibedge*5+0] = n0;
	      grid->bedges[ibedge*5+1] = current_ghost_node;
	      grid->bedges[ibedge*5+3] = b;
	      grid->bedges[ibedge*5+4] = s;

	      // Lets store the normal at the actual node.
	      //curved_boundary_normal_vector ( grid, bc, xM, &nx, &ny );
	      curved_boundary_normal_vector ( grid, bc, xL, &nx, &ny );

	      grid->xn_bedges[ibedge*3+0] = nx;
	      grid->xn_bedges[ibedge*3+1] = ny;

	      grid->x[NDIM*(grid->nn + current_ghost_node)+0] = xM[0];
	      grid->x[NDIM*(grid->nn + current_ghost_node)+1] = xM[1];

	      // Mark the ghost node as such.
	      grid->node_state[grid->nn + current_ghost_node] = GHOST;
	  
	      ibedge++;
	      current_ghost_node++;

	      //Store the data.
	      grid->bedges[ibedge*5+0] = n1;
	      grid->bedges[ibedge*5+1] = current_ghost_node;
	      grid->bedges[ibedge*5+3] = b;
	      grid->bedges[ibedge*5+4] = s;

	      // Reuse normal from first node.
	      //grid->xn_bedges[ibedge*3] = nx;
	      //grid->xn_bedges[ibedge*3+1] = ny;

	      // Lets store the normal at the actual node.
	      curved_boundary_normal_vector ( grid, bc, xR, &nx, &ny );

	      grid->xn_bedges[ibedge*3] = nx;
	      grid->xn_bedges[ibedge*3+1] = ny;
	      
	      grid->x[NDIM*(grid->nn + current_ghost_node)+0] = xM[0];
	      grid->x[NDIM*(grid->nn + current_ghost_node)+1] = xM[1];

	      // Mark the ghost node as such.
	      grid->node_state[grid->nn + current_ghost_node] = GHOST;
	  
	      ibedge++;
	      current_ghost_node++;
	    }
	}
    }                                                     // End boundary loop.

  // Sanity check.
  if ( tbedges != (ibedge-1) )
    {
      printf("FATAL ERROR: In Build_Boundary_Edge_Data(), total number of edges greated does not match expected number.\n");
      printf("  tbedges = %d , ibedge = %d\n",tbedges,ibedge);
      exit(0);
    }
  if ( tbedges != (current_ghost_node-1) )
    {
      printf("FATAL ERROR: In Build_Boundary_Edge_Data(), total number of ghost nodes greated does not match expected number.\n");
      printf("  tbedges = %d , number of ghost nodes = %d\n",tbedges,current_ghost_node);
      exit(0);
    }

  printf("  Created %d ghost nodes.\n",grid->nn_ghost);

  // I still need to stitch on the correct interior element to the boundary edge in question.
  
  // Loop over all boundary segments.
  for ( b=1; b <= grid->nb; b++ )
    {                                                       // Boundary loop.
      for ( s=1; s <= grid->nbs[b]; s++ )
	{                                  
	  n0 = grid->bs[b][s][0];
	  n1 = grid->bs[b][s][1];
	  
	  // Loop over the element types and elements.
	  for ( i=Tri; i <= Quad; i++ )                                                       // Element types loop.
	    {
	      // Find the number of nodes and edges for this element type.
	      num_vert = NumberNodesForElement(i);
	      
	      for ( j=1; j <= grid->num_elem[i]; j++ )                                        // Element loop.
		{
		  // Match the first node in the boundary segment.
		  for ( k=0; k < num_vert; k++ )
		    {
		      if ( n0 == grid->c2n[i][j*num_vert+k] )
			break;
		    }

		  // k is set to correct local node or == num_vert.
		  if ( k == num_vert )
		    continue;

		  // Now is the next node in the element's list n1?
		  if ( n1 != grid->c2n[i][j*num_vert+ ((k+1)%num_vert) ] )
		    continue;

		  // We have found our element.
		  break;
		}                                                                           // End element loop.

	      // Did we find the element, or do we proceed to the next element type?
	      if ( j <= grid->num_elem[i] )
		break;
	    }                                                                           // End element type loop.

	  // At this point we should the correct element.
	  if ( i > Quad )
	    {
	      printf("FATAL ERROR: In Build_Boundary_Edge_Data(), a segment was not able to find a volume element.\n");
	      exit(1);
	    }
	  
	  // Get the global element id.
	  gelem = LocalToGlobalID ( grid->num_elem,
				    j,
				    i );
	  
	  // Now I just zip through the bedges data structure and find the ones with the correct entries for boundary (b) and segment(s).
	  for ( i=1; i <= grid->nbedges; i++ )
	    {
	      if ( b == grid->bedges[i*5+3] &&
		   s == grid->bedges[i*5+4] )
		{
		  grid->bedges[i*5+2] = gelem;
		}
	    }

	}   // End segment loop.
    }    // End boundary loop.

  // Normalize the vectors and store the length.
  for ( i=1; i <= grid->nbedges; i++ )
    {
      mag = sqrt( grid->xn_bedges[i*3]*grid->xn_bedges[i*3] +
		  grid->xn_bedges[i*3+1]*grid->xn_bedges[i*3+1]);

      bc = grid->bbc[ grid->bedges[i*5+3] ];

      if ( bc > 10 )
	{
	  node = grid->bedges[i*5+0];
	  ghostnode = grid->bedges[i*5+1] + grid->nn;

	  xL[0] = grid->x[node*NDIM+0];
	  xL[1] = grid->x[node*NDIM+1];
	  xR[0] = grid->x[ghostnode*NDIM+0];
	  xR[1] = grid->x[ghostnode*NDIM+1];

	  curved_boundary_arclength ( grid, bc, xL, xR, &mag );

	  // For the curved boundaries, the unit normal is already stored.
	  grid->xn_bedges[i*3+2] = mag;
	}
      else
	{
	  if ( mag > 1.0E-20 )
	    {
	      grid->xn_bedges[i*3] /= mag;
	      grid->xn_bedges[i*3+1] /= mag;
	    }
	  grid->xn_bedges[i*3+2] = mag;
	}
    }

  #if DEBUG
  printf("BEDGES DEBUG:\n");
  for ( i=1; i <= grid->nbedges; i++ )
    {
      printf("Bedge %d:\n",i);
      printf("  node = %d , ghost node = %d, gelem = %d, boundary = %d , segment = %d\n",grid->bedges[i*5],
	     grid->bedges[i*5+1],grid->bedges[i*5+2],grid->bedges[i*5+3],grid->bedges[i*5+4]);
    }
  #endif

  #if DEBUG
  printf("XN_BEDGES DEBUG:\n");
  for ( i=1; i <= grid->nbedges; i++ )
    {
      printf("Bedge %d:\n",i);
      printf("  nx = %lf  ny = %lf  length = %lf\n",grid->xn_bedges[i*3],
	     grid->xn_bedges[i*3+1],grid->xn_bedges[i*3+2]);
    }
  #endif

  return;
}

// Now I need a function to go back in and correct the subedges that connect to the ghost nodes
// lying on curved boundaries.


//=============================================================
// 
//  Fix_Subedge_Data()
//
//  Corrects the definition of subedges that are connected to
//  ghost nodes on curved boundaries.
//  
//  GRID *grid                        // The grid.
//
//=============================================================

void Fix_Subedge_Data ( GRID *grid )
{
  int i,j;                               // Loop counters.
  int b,s;                               // Loop counters.
  int iedge;                             // Global edge.
  int gelem;                             // Global element.
  int igelem;
  double xmid, ymid;                     // Values at face midpoints.
  double nx,ny;                          // Normal vector components.
  int ghostnode;
  int nL,nR;
  int nodeL,nodeR;
  double len;
  int bc;
  double XC[NDIM],XM[NDIM];

  i = 0;
  for ( b=1; b <= grid->nb; b++ )
    {
      bc = grid->bbc[b];
      if ( bc > 10 )
	i = 1;
    }

  if ( i == 0 )
    {
      printf("Fix_Subedge_Data(): No curved boundaries are detected.\n");
      return;
    }

  // Loop over the boundary edges and process curved ones.

  for ( i=1; i <= grid->nbedges; i++ )
    {
      // Get the data.
      ghostnode = grid->bedges[i*5+1];
      gelem = grid->bedges[i*5+2];
      b = grid->bedges[i*5+3];
      s = grid->bedges[i*5+4];

      bc = grid->bbc[b];

      if ( bc <= 10 )
	continue;   // Skip linear edges.

      // Get the edge data.
      nL = grid->bs[b][s][0];
      nR = grid->bs[b][s][1];

      // Now lets find the subedge.
      for ( j=1; j <= grid->nsubedges; j++ )
	{
	  iedge = grid->subedges[j*2+0];
	  igelem = grid->subedges[j*2+1];

	  nodeL = grid->edges[iedge*2+0];
	  nodeR = grid->edges[iedge*2+1];

	  if (  ( nodeL == nL || nodeL == nR ) &&
		( nodeR == nL || nodeR == nR ) )
	    {
	      break;  // Found the correct edge.
	    }
	}

      if ( j > grid->nsubedges )
	{
	  printf(" FATAL ERROR: In Fix_Subedge_Data(), the corresponding global edge could not be found.\n");
	  exit(1);
	}

      if ( gelem != igelem )
	{
	  printf(" FATAL ERROR: In Fix_Subedge_Data(), the boundary edge and subedge do not have the same global element.\n");
	  exit(1);
	}
      

      // Now lets fix the subedge.
      
      // Get the centroid.
      XC[0] = grid->el_cent[gelem*NDIM];
      XC[1] = grid->el_cent[gelem*NDIM+1];

      // Get the midpoint on the curve.
      XM[0] = grid->x[ (grid->nn + ghostnode)*NDIM + 0];
      XM[1] = grid->x[ (grid->nn + ghostnode)*NDIM + 1];

      // Find the new subedge midpoint.
      //xmid = 0.5*( XC[0] + XM[0] );
      //ymid = 0.5*( XC[1] + XM[1] );

      xmid = XM[0];
      ymid = XM[1];

      nx = XC[1] - XM[1];
      ny = XM[0] - XC[0];

      len = sqrt( nx*nx + ny*ny );

      // Normalize the vector.
      nx /= len;
      ny /= len;

      if ( nL > nR )
	{
	  nx = -nx;
	  ny = -ny;
	}

      if ( DEBUG )
	{
	  printf("  FIXING SUBEDGE %d : nx = %.15E  ny = %.15E  len = %.15E\n",j,grid->xn_subedges[j*3+0],grid->xn_subedges[j*3+1],grid->xn_subedges[j*3+2]);
	  printf("             Now its  nx = %.15E  ny = %.15E  len = %.15E\n",nx,ny,len);
	  printf("             Old: xmid = %.15E     ymid = %.15E\n",grid->xm_subedges[j*NDIM+0],grid->xm_subedges[j*NDIM+1]);
	  printf("             New: xmid = %.15E     ymid = %.15E\n",xmid,ymid);
	  fflush(stdout);
	}

      // Store the new info.
      grid->xn_subedges[j*3+0] = nx;
      grid->xn_subedges[j*3+1] = ny;
      grid->xn_subedges[j*3+2] = len;

      grid->xm_subedges[j*NDIM+0] = xmid;
      grid->xm_subedges[j*NDIM+1] = ymid;
    }
  
  return;
}

//=============================================================
// 
//  Generate_Gaussian_Data_For_Edges()
//
//  Fills in the Gaussian node data for the subedges and bedges.
//  
//  GRID *grid                        // The grid.
//
//=============================================================

void Generate_Gaussian_Data_For_Edge ( GRID *grid )
{
  int i,j,k,n;                           // Loop counters.
  int tsubedges;                         // Total number of subedges.
  int tbedges;                           // Total number of bedges.
  int iedge;                             // Global edge.
  int isubedge;                          // Subedge being worked on.
  int ibedge;                            // Bedge being worked on.
  int gelem;                             // Global element.
  int nodeL, nodeR;                      // Edge/element nodes.
  int b,bc,s;                            // Boundary segment data.
  int real_node,ghost_node;              
  double xi,yi;
  double xL,yL,xR,yR;
  double XL[NDIM],XR[NDIM];
  double GP[6];
  double xc, yc;                         // Values at the centroid.
  double xmid, ymid;                     // Values at face midpoints.
  double nx,ny;                          // Normal vector components.
  double len;                            // Normal vector length.
  double x1,x2,x3,y1,y2,y3;              // Gaussian nodes.
  double t1,t2,t3;                       // Parameter nodes.

  tsubedges = grid->nsubedges;
  tbedges = grid->nbedges;

  // Allocate space for the data.
  if ( grid->xg_subedges == NULL )
    {
      grid->xg_subedges = (double*)malloc((tsubedges+1)*3*NDIM*sizeof(double));
      if ( grid->xg_subedges == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'xg_subedges'.\n"); exit(1); }
    }
  
  if ( grid->xg_bedges == NULL )
    {
      grid->xg_bedges = (double*)malloc((tbedges+1)*3*NDIM*sizeof(double));
      if ( grid->xg_bedges == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'xg_bedges'.\n"); exit(1); }
    }
  
  if ( grid->xgn_bedges == NULL )
    {
      grid->xgn_bedges = (double*)malloc((tbedges+1)*3*(NDIM+1)*sizeof(double));
      if ( grid->xgn_bedges == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'xgn_bedges'.\n"); exit(1); }
    }

  // Now we can initialize the Gaussian nodes of the subedges. This assumes that all subedges are linear.

  for ( isubedge=1; isubedge <= grid->nsubedges; isubedge++ )
    {
      // Extract the real edge this corresponds to.
      iedge = grid->subedges[isubedge*2];

      xmid = grid->xm_subedges[isubedge*NDIM+0];
      ymid = grid->xm_subedges[isubedge*NDIM+1];

      // Find the global element.
      gelem = grid->subedges[isubedge*2+1];

      // Get the element centroid.
      xc = grid->el_cent[gelem*2];
      yc = grid->el_cent[gelem*2+1];
  
      t1 = sqrt(15.0)/5.0;
      t2 = 0.;
      t3 = sqrt(15.0)/(-5.0);

      x1 = (xmid+xc)*0.5 + (xc-xmid)*0.5*t1;
      y1 = (ymid+yc)*0.5 + (yc-ymid)*0.5*t1;
      
      x2 = (xmid+xc)*0.5 + (xc-xmid)*0.5*t2;
      y2 = (ymid+yc)*0.5 + (yc-ymid)*0.5*t2;
      
      x3 = (xmid+xc)*0.5 + (xc-xmid)*0.5*t3;
      y3 = (ymid+yc)*0.5 + (yc-ymid)*0.5*t3;

      // Store the data.
      grid->xg_subedges[isubedge*6+0] = x1;
      grid->xg_subedges[isubedge*6+1] = y1;
      grid->xg_subedges[isubedge*6+2] = x2;
      grid->xg_subedges[isubedge*6+3] = y2;
      grid->xg_subedges[isubedge*6+4] = x3;
      grid->xg_subedges[isubedge*6+5] = y3;
    }
  // Important to note here is that I don't store normal components since the subedges are assumed linear,
  // they are the same as the subedge.


  // Now we can initialize the Gaussian nodes of the bedges.
  for ( i=1; i <= grid->nbedges; i++ )
    {
      // Grab the data from the structure.
      real_node = grid->bedges[i*5+0];
      ghost_node = grid->bedges[i*5+1] + grid->nn;

      // Get the normal vector information.
      nx = grid->xn_bedges[i*3+0];
      ny = grid->xn_bedges[i*3+1];
      len = grid->xn_bedges[i*3+2];

      b = grid->bedges[i*5+3];
      bc = grid->bbc[b];

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
  
	  t1 = sqrt(15.0)/5.0;
	  t2 = 0.;
	  t3 = sqrt(15.0)/(-5.0);

	  x1 = (xL+xR)*0.5 + (xR-xL)*0.5*t1;
	  y1 = (yL+yR)*0.5 + (yR-yL)*0.5*t1;
	  
	  x2 = (xL+xR)*0.5 + (xR-xL)*0.5*t2;
	  y2 = (yL+yR)*0.5 + (yR-yL)*0.5*t2;

	  x3 = (xL+xR)*0.5 + (xR-xL)*0.5*t3;
	  y3 = (yL+yR)*0.5 + (yR-yL)*0.5*t3;

	  // Store the data.
	  grid->xg_bedges[i*6+0] = x1;
	  grid->xg_bedges[i*6+1] = y1;
	  grid->xg_bedges[i*6+2] = x2;
	  grid->xg_bedges[i*6+3] = y2;
	  grid->xg_bedges[i*6+4] = x3;
	  grid->xg_bedges[i*6+5] = y3;

	  // store the normal components.
	  grid->xgn_bedges[i*9+0] = nx;
	  grid->xgn_bedges[i*9+1] = ny;
	  grid->xgn_bedges[i*9+2] = len;

	  grid->xgn_bedges[i*9+3] = nx;
	  grid->xgn_bedges[i*9+4] = ny;
	  grid->xgn_bedges[i*9+5] = len;

	  grid->xgn_bedges[i*9+6] = nx;
	  grid->xgn_bedges[i*9+7] = ny;
	  grid->xgn_bedges[i*9+8] = len;
	}
      else // Curved boundary.
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
	  curved_boundary_gauss_points ( grid, bc, XL, XR, GP );
	  
	  curved_boundary_arclength ( grid, bc, XL, XR, &len );

	  for ( j=0; j < 3; j++ )
	    {
	      // Get the normal vector.
	      curved_boundary_normal_vector ( grid, bc, &(GP[j*NDIM]), &nx, &ny );

	      // Store the data.
	      grid->xg_bedges[i*6+j*NDIM]   = GP[j*NDIM];
	      grid->xg_bedges[i*6+j*NDIM+1] = GP[j*NDIM+1];

	      grid->xgn_bedges[i*9+j*3+0] = nx;
	      grid->xgn_bedges[i*9+j*3+1] = ny;
	      grid->xgn_bedges[i*9+j*3+2] = len;
	    }
	}
    }
	  

  return;
}
