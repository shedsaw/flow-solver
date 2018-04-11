//=============================================================
// 
//  mesh_utils.C
//  
//  Functions that are needed to modify/fix the grid or retrieve
//  grid information.
//
//  Written by - Shane Sawyer
//
//=============================================================

#include <stdio.h>
#include <stdlib.h>
#include "grid.h"
#include "maps.h"
#include "mesh_utils.h"
#include "Defined_Variables.h"


//=============================================================
// 
//  rewind2D()
//
//  Rewinds all elements in a 2D mesh. Reorients boundary
//  segments if necessary.
// 
//  GRID *grid                          // The grid.
//
//=============================================================

void rewind2D ( GRID *grid )
{
  int i,j,b,s;                           // Element (segment) loop counter.
  int n1,n2;                             // Element (segment) node id placeholders.
  double nx,ny;                          // Segment normal vector.
  double vx, vy;                         // Segment to centroid vector.
  double xmid, ymid;                     // Mid-segment point.
  double xc, yc;                         // Element centroid.
  int start_idx;                         // Start index in csn.
  int end_idx;                           // End index in csn.
  int elem_id;                           // Element id.
  int type;                              // Element type.
  int num_vert;                          // Number of vertices for an element.
  int flag;                              // Logic flag.
  int num_seg_fixed=0;                   // Number of segments that required swapping.
  double dot;                            // Vector dot product.

  type = 0;

  // Rewind all triangles.
  for ( i=1; i <= grid->num_elem[Tri]; i++ )
    {
      // For triangles, only need to hold onto one node and then switch the other one with it. Node 1 will stay fixed.
      n1 = grid->c2n[Tri][i*3+1];
      grid->c2n[Tri][i*3+1] = grid->c2n[Tri][i*3+2];
      grid->c2n[Tri][i*3+2] = n1;
    }

  // Rewind all quads.
  for ( i=1; i <= grid->num_elem[Quad]; i++ )
    {
      // For quads, nodes 1 and 3 (Fortran indexing) stay fixed. Swap nodes 2 and 4.
      n1 = grid->c2n[Quad][i*4+1];
      grid->c2n[Quad][i*4+1] = grid->c2n[Quad][i*4+3];
      grid->c2n[Quad][i*4+3] = n1;
    }

  // The boundary segments need verification to make sure that the edge normal points out of the
  // mesh. To do this, we need to find which element the edge belongs, so the node to cell map is needed.
  printf("  Verifying boundary segment winding...\n");
  Cells_Surrounding_Node(grid->nn,
			 grid->num_elem,
			 grid->c2n,
			 &(grid->csn),
			 &(grid->ncsn));

  // Now loop over the boundary segments and make sure that the edge normal points out of the mesh.
  for ( b=1; b <= grid->nb; b++ )
    {                                                       // Boundary loop.
      for ( s=1; s <= grid->nbs[b]; s++ )
	{                                                   // Segment loop.

	  // Reset the flag.
	  flag = 0;

	  // Retrieve the edge nodes and form the normal vector.
	  n1 = grid->bs[b][s][0];
	  n2 = grid->bs[b][s][1];
	  
	  // Find the segment midpoint.
	  xmid = (grid->x[NDIM*n1]+grid->x[NDIM*n2])/2.;
	  ymid = (grid->x[NDIM*n1+1]+grid->x[NDIM*n2+1])/2.;

	  // Set up the search indices.
	  start_idx = grid->ncsn[n1];
	  end_idx = grid->ncsn[n1+1];

	  // Loop the csn list and find the correct element.
	  for ( i=start_idx; i < end_idx; i++ )            // Element search loop.
	    {
	      // Grab the element id (global).
	      elem_id = grid->csn[i];

	      // Get the local element id.
	      elem_id = GlobalToLocalID( grid->num_elem , elem_id , &type );

	      // Get the number of vertices.
	      num_vert = NumberNodesForElement(type);

	      // See if n2 is in the node list for this element.
	      for ( j=0; j < num_vert; j++ )              // Element node loop.
		{
		  if ( n2 == grid->c2n[type][elem_id*num_vert+j] )
		    {
		      flag = 1;
		      break;
		    }
		}                                        // End element node loop.

	      // If we found the element, break the search loop.
	      if ( flag == 1 )
		break;
	    }                                            // End element search loop.

	  // At this point flag should be 1, if not that means an error occured in the search or map
	  // construction.
	  if ( flag != 1 )
	    {
	      printf("FATAL ERROR: In boundary %d, segment %d:\n",b,s);
	      printf("  Nodes %d and %d form an edge but do not share element connectivity.\n",n1,n2);
	      exit(1);
	    }

	  // Find the element centroid.
	  xc = 0.;  yc = 0.;
	  for ( i=0; i < num_vert; i++ )
	    {
	      xc += grid->x[NDIM*(grid->c2n[type][elem_id*num_vert+i])];
	      yc += grid->x[NDIM*(grid->c2n[type][elem_id*num_vert+i])+1];
	    }
	  xc /= (num_vert*1.);
	  yc /= (num_vert*1.);

	  // Create the vector from the mid-edge point to the centroid.
	  vx = xc - xmid;
	  vy = yc - ymid;

	  // Create the edge normal vector.
	  nx = grid->x[NDIM*n2+1] - ymid;
	  ny = xmid - grid->x[NDIM*n2];

	  // Take the dot product of the two vectors.
	  dot = vx*nx + vy*ny;

	  // If the dot product is positive, then the angle is less than 90 degrees and points into the mesh, so swap.
	  if ( dot > 0. )
	    {
	      grid->bs[b][s][0] = n2;
	      grid->bs[b][s][1] = n1;
	      num_seg_fixed++;	      
	    }

	}                                                   // End segment loop.
    }                                                    // End boundary loop.


  printf("rewind2D(): %d boundary segments swapped nodes.\n",num_seg_fixed);

  // Free up the memory for csn,ncsn since they are remade later.
  freenull(grid->ncsn);
  freenull(grid->csn);

  return;
}


//=============================================================
// 
//  NumberNodesForElement()
//
//  Returns the number of vertices an element type has.
//
//  int e;                 // Element type.
//  
//=============================================================

int NumberNodesForElement ( int e )
{
  switch(e)
    {
    case(Edge):
      return (2);
      break;
    case(Tri):
      return (3);
      break;
    case(Quad):
      return (4);
      break;
    default:
      printf("FATAL ERROR: 'NumberNodesForElement' received an invalid element type ID : %d\n",e);
      exit(1);
    }
}



//=============================================================
// 
//  NumberEdgesForElement()
//
//  Returns the number of edges an element type has.
//
//  int e;                 // Element type.
//  
//=============================================================

int NumberEdgesForElement ( int e )
{
  switch(e)
    {
    case(Edge):
      return (1);
      break;
    case(Tri):
      return (3);
      break;
    case(Quad):
      return (4);
      break;
    default:
      printf("FATAL ERROR: 'NumberEdgesForElement' received an invalid element type ID : %d\n",e);
      exit(1);
    }
}


//=============================================================
// 
//  NumberNeighborsForElement()
//
//  Returns the number of neighbors an element type has.
//
//  int e;                 // Element type.
//  
//=============================================================

int NumberNeighborsForElement ( int e )
{
  switch(e)
    {
    case(Edge):
      return (1);
      break;
    case(Tri):
      return (3);
      break;
    case(Quad):
      return (4);
      break;
    default:
      printf("FATAL ERROR: 'NumberNeighborsForElement' received an invalid element type ID : %d\n",e);
      exit(1);
    }
}




//=============================================================
// 
//  RetrieveElementNodesSurroundingVertex()
//
//  Given an element and a vertex with element numbering (0,1,..),
//  finds the list of element vertices that connect with this vertex.
//
//  int **c2n;                        // The cell to node table.
//  int elem_id;                      // The global element id.
//  int node_id;                      // The element vertex in question.
//  int *list;                        // The list of connected nodes (in local numbering).
//  int *length;                      // Length of the above list.  
//
//=============================================================

void RetrieveElementNodesSurroundingVertex ( int **c2n,
					     int *num_elem,
					     int elem_id,
					     int node_id,
					     int *list,
					     int *length )
{
  int i;                            // Loop counter.
  int type;                         // Element type.
  int vert_id;                      // The local id of the node in the element vertex list.
  int num_vert;                     // Number of vertices.
  
  type = 0;  vert_id = 0;

  // Get the local element id and more importantly the type.
  elem_id = GlobalToLocalID ( num_elem, elem_id , &type );

  // Determine which vertex node_id is.

  num_vert = NumberNodesForElement(type);    // Get the number of vertices.
  for ( i=0; i < num_vert; i++ )
    {
      if ( node_id == c2n[type][num_vert*elem_id + i] )
	{
	  vert_id = i;                      // Break early when we find the id.
	  break;
	}
    }

  // Use switch statements to pick out the element type first, then switch statements to find the proper vertex.
  switch (type)
    {
    case Edge:
      printf("FATAL ERROR: Edge is not a supported type in 'RetrieveElementNodesSurroundingVertex'.\n");
      exit(1);
      break;
    case Tri:

      // Since the length is the same for any vertex on a triangle set it here.
      *length = 2;

      switch (vert_id)                      // Vertex offset in the last column before bracket.  
	{                                   //   \/
	case 0:
	  list[0] = c2n[type][num_vert*elem_id + 1];
	  list[1] = c2n[type][num_vert*elem_id + 2];
	  break;
	case 1:
	  list[0] = c2n[type][num_vert*elem_id + 0];
	  list[1] = c2n[type][num_vert*elem_id + 2];
	  break;
	case 2:
	  list[0] = c2n[type][num_vert*elem_id + 0];
	  list[1] = c2n[type][num_vert*elem_id + 1];
	  break;
	default:
	  printf("FATAL ERROR: Vertex id is out of range in 'RetrieveElementNodesSurroundingVertex'. ID=%d,type=%d\n",vert_id,type);
	  exit(1);
	}
      break;
    case Quad:

      // Since the length is the same for any vertex on a quad, set it here..
      *length = 2;

      switch (vert_id)
	{
	case 0:
	  list[0] = c2n[type][num_vert*elem_id + 1];
	  list[1] = c2n[type][num_vert*elem_id + 3];
	  break;
	case 1:
	  list[0] = c2n[type][num_vert*elem_id + 0];
	  list[1] = c2n[type][num_vert*elem_id + 2];
	  break;
	case 2:
	  list[0] = c2n[type][num_vert*elem_id + 1];
	  list[1] = c2n[type][num_vert*elem_id + 3];
	  break;
	case 3:
	  list[0] = c2n[type][num_vert*elem_id + 0];
	  list[1] = c2n[type][num_vert*elem_id + 2];
	  break;
	default:
	  printf("FATAL ERROR: Vertex id is out of range in 'RetrieveElementNodesSurroundingVertex'. ID=%d,type=%d\n",vert_id,type);
	  exit(1);
	}
      break;
    default:
      printf("FATAL ERROR: Invalid Element type in 'RetrieveElementNodesSurroundingVertex'. Element type=%d\n",type);
      exit(1);
    }
  return;
  
}



//=============================================================
// 
//  RetrieveElementNodesFormingEdge()
//
//  Given an element and which edge under consideration, returns the 
//  left and right nodes of that edge with the condition that the
//  left node is of lower order than the right node if desired by setting
//  order to 1.
//
//  int **c2n;                        // The cell to node table.
//  int elem_id;                      // The element id.
//  int type;                         // The element type.
//  int edge_id;                      // The element edge in question.
//  int *left;                        // The left node of the edge.
//  int *right;                       // The right node of the edge.
//  int order;                        // Logic flag. 0 = do not sort the nodes by id.
//                                    //             1 = Sort the nodes such that left node has the lower index.
//
//=============================================================

void RetrieveElementNodesFormingEdge ( int **c2n,
				       int elem_id,
				       int type,
				       int edge_id,
				       int *left,
				       int *right,
				       int order )
{
  int node;                         // Temporary node id if swapping.
  int lnode, rnode;                 // Edge node ids. (May need to be swapped.)
  int num_vert;                     // Number of vertices.

  num_vert = NumberNodesForElement(type);    // Get the number of vertices.

  // Do a switch on the element type and then a switch on the edge id to retreive the proper nodes.
  switch (type)
    {
    case Tri:
      switch (edge_id)                     // Vertex offset in the last column before bracket.
	{                                  //  \/
	case 0:
	  lnode = c2n[type][elem_id*num_vert + 0];
	  rnode = c2n[type][elem_id*num_vert + 1];
	  break;
	case 1:
	  lnode = c2n[type][elem_id*num_vert + 1];
	  rnode = c2n[type][elem_id*num_vert + 2];
	  break;
	case 2:
	  lnode = c2n[type][elem_id*num_vert + 2];
	  rnode = c2n[type][elem_id*num_vert + 0];
	  break;
	default:
	  printf("FATAL ERROR: Edge ID is out of range in 'RetrieveElementNodesFormingEdge'. ID=%d, type=%d\n",edge_id,type);
	  exit(1);
	  break;
	}
      break;
    case Quad:
      switch (edge_id)
	{
	case 0:
	  lnode = c2n[type][elem_id*num_vert + 0];
	  rnode = c2n[type][elem_id*num_vert + 1];
	  break;
	case 1:
	  lnode = c2n[type][elem_id*num_vert + 1];
	  rnode = c2n[type][elem_id*num_vert + 2];
	  break;
	case 2:
	  lnode = c2n[type][elem_id*num_vert + 2];
	  rnode = c2n[type][elem_id*num_vert + 3];
	  break;
	case 3:
	  lnode = c2n[type][elem_id*num_vert + 3];
	  rnode = c2n[type][elem_id*num_vert + 0];
	  break;
	default:
	  printf("FATAL ERROR: Edge ID is out of range in 'RetrieveElementNodesFormingEdge'. ID=%d, type=%d\n",edge_id,type);
	  exit(1);
	  break;
	}
      break;
    default:
      printf("FATAL ERROR: Element tpye is out of range in 'RetrieveElementNodesFormingEdge'. type=%d\n",type);
      exit(1);
      break;
    }

  // Now swap lnode and rnode if needed and desired.
  if ( order )
    {
      if ( lnode > rnode )
	{
	  node = lnode;
	  lnode = rnode;
	  rnode = node;
	}
    }

  *left = lnode;
  *right = rnode;

  return;
}

//=============================================================
// 
//  GlobalToLocalID()
//
//  Given a global element id, returns the element id local to that
//  element type (for searching in c2n) and the type.
//
//  int *num_elem;                    // Number of elements in the grid.
//  int elem_id;                      // Global element id.
//  int *e;                           // Element type.
//  
//=============================================================

int GlobalToLocalID ( int *num_elem,
		      int elem_id,
		      int *e )
{
  int i;                             // Loop counter.

  // Loop over all element types and test if elem_id is within the bounds of the current element type.
  for ( i=Tri; i <= Quad; i++ )
    {
      if ( elem_id <= num_elem[i] )
	{
	  *e = i;        // Set the element type.
	  return elem_id;      // Return local id.
	}
      else
	elem_id -= num_elem[i];     // Subtract this off now since this was previously added when global ids were constructed.
    }

  // At this point the local element id has not been found so kill the program.

  // Add back to elem_id the amount subtracted and display the offender.
  for ( i=Tri; i <= Quad; i++ )  // Compatible with 2D and 3D.
    elem_id += num_elem[i];
  
  printf("FATAL ERROR: 'GlobalToLocalID' has received an invalid Global Element ID: %d.\n",elem_id);
  return 0;

}



//=============================================================
// 
//  LocalToGlobalID()
//
//  Given a local element id and element type, returns the 
//  global element ID.
//
//  int *num_elem;                    // Number of elements in the grid.
//  int elem_id;                      // Local element id.
//  int e;                            // Element type.
//  
//=============================================================

int LocalToGlobalID ( int *num_elem,
		      int elem_id,
		      int e )
{
  int i;                             // Loop counter.

  // Check that element type is valid.
  if ( e < Tri || e > Quad )
    {
      printf("FATAL ERROR: LocalToGlobalID received an invalid element type ( type=%d ).\n",e);
      exit(1);
    }

  // Loop over element types up to the element's type and add to ID the total number of elements
  // of that type.
  for ( i=Tri; i < e; i++ )
    {
      elem_id += num_elem[i];
    }

  return elem_id;

}


//=============================================================
// 
//  Swap()
//
//  Given left and right node, make sure that the left node is
//  lower. If not, swap them.
//
//  int *nodeL;                    // The left node.
//  int *nodeR;                    // The right node.
//=============================================================

void Swap ( int *nodeL, int *nodeR )
{
  int temp;

  if ( (*nodeL) > (*nodeR) )
    {
      temp = (*nodeL);
      (*nodeL) = (*nodeR);
      (*nodeR) = temp;
    }

  return;
}


//=============================================================
// 
//  ConvertConservedToPrimitive()
//
//  Given gamma and the conserved variables, return the associated
//  primitive values.
//
//  double gamma;                    // The value of gamma.
//  double *Q;                       // The conserved variables.
//=============================================================

void ConvertConservedToPrimitive ( double gamma , double *Q )
{
  double rho, u, v, P;

  rho = Q[0];
  u = Q[1] / Q[0];
  v = Q[2] / Q[0];
  P = (gamma - 1.)*(Q[3] - 0.5*rho*(u*u + v*v));

  Q[0] = rho;
  Q[1] = u;
  Q[2] = v;
  Q[3] = P;

  return;
}


//=============================================================
// 
//  ConvertPrimitiveToConserved
//
//  Given gamma and the primitive variables, return the associated
//  conserved values.
//
//  double gamma;                    // The value of gamma.
//  double *Q;                       // The primitive variables.
//=============================================================

void ConvertPrimitiveToConserved ( double gamma , double *Q )
{
  double rho, um, vm, E;

  rho = Q[0];
  um = rho*Q[1];
  vm = rho*Q[2];
  E = Q[3] / (gamma - 1.) + 0.5*Q[0]*(Q[1]*Q[1] + Q[2]*Q[2]);

  Q[0] = rho;
  Q[1] = um;
  Q[2] = vm;
  Q[3] = E;

  return;
}

