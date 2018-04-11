//=============================================================
// 
//  reorder.cpp
//  
//  Applies the Cuthill-McKee Alogrithm to sort the nodes in
//  the mesh.
//
//  Written by - Shane Sawyer
//
//=============================================================

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "grid.h"
#include "maps.h"
#include "mesh_utils.h"
#include "Defined_Variables.h"
#include "reorder.h"

extern "C" int compare_degree ( const void *a, const void *b );

struct node
{
  int id;
  int degree;
};

extern "C" int compare_degree ( const void *a, const void *b )
{
  struct node * aa = (struct node*)a;
  struct node * bb = (struct node*)b; 
  //struct node node1 = (struct node)(*(struct node*))a;
  //struct node node2 = (struct node)(*(struct node*))b;

  if ( aa->degree < bb->degree )
  //if ( node1.degree > node2.degree )
    return 1;
  //else if ( node1.degree < node2.degree )
  else if ( aa->degree > bb->degree )
    return -1;
  else
    return 0;
}


//=============================================================
// 
//  Cuthill_McKee()
//
//  Reorders the vertices in the mesh by applying the Cuthill-
//  McKee Algorithm as outlined in Dr. Hyam's dissertation.
//  
//
//=============================================================

void Cuthill_McKee (GRID *grid)
{
  int i,j,k;                             // Loop counters.
  int num_vert;                          // Number of vertices for an element.
  int seed;                              // Seed node to process.
  int inode;                             // The new node counter - used to renumber the nodes.
  int start_idx;                         // Starting index in an array.
  int end_idx;                           // Ending index in an array.
  int degree;                            // Degree of the node.
  int maxdegree=0;                       // Maximum nodal degree in the mesh.
  int *R;                                // The R queue that contains the nodes connected to the processing node.
  struct node *Rnode;                    // The R queue cast as collection of node objects used for sorting by degree.
  int *S;                                // The queue containing nodes that need processing.
  int *Queue;                            // The queue of unorderd nodes. A 1 or 0 is stored for every node to indicate
                                         // whether or not a node has been reordered.
  int startR=0;                          // Starting position of the R queue.
  int endR=0;                            // Ending position of the R queue.            Initially the queues are empty,
  int startS=0;                          // Starting position of the S queue.          hence the 0's.
  int endS=0;                            // Ending position of the S queue.
  double *tx;                            // Temporary array to move the spatial coordinates into the x array in the 
                                         // new ordering.
  
  // For the queues, they are empty if endX-startX is equal to 0.

  // Allocate space to store the new index of the vertices.
  grid->nn_map = (int*)malloc((grid->nn+1)*sizeof(int));
  if ( grid->nn_map == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'nn_map'.\n"); exit(1); }

  // Allocate space for Queue and initialize all nodes to 0 to indicate unordered.
  Queue = (int*)malloc((grid->nn+1)*sizeof(int));
  if ( Queue == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'Queue'.\n"); exit(1); }

  // To save memory allocations, assume worst case scenario that all nodes will be added to R and
  // S queues at the same time.
  R = (int*)malloc(grid->nn*sizeof(int));
  S = (int*)malloc(grid->nn*sizeof(int));

  if ( R == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'R'.\n"); exit(1); }
  if ( S == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'S'.\n"); exit(1); }

  for ( i=1; i <= grid->nn; i++ )
    {
      Queue[i] = 0;
    }

  // To begin the algorithm, we need to know the nodes surrounding the vertices. This is accomplished by building first
  // the cells surrounding nodes and then the nodes surrounding nodes maps.
  printf("  Cuthill-McKee(): Building list of cells surrounding vertex...");
  Cells_Surrounding_Node(grid->nn,
			 grid->num_elem,
			 grid->c2n,
			 &(grid->csn),
			 &(grid->ncsn));
  printf("Done\n");

  printf("  Cuthill-McKee(): Building list of nodes surrounding vertex...");
  Nodes_Surrounding_Node(grid->nn,
			 grid->num_elem,
			 grid->c2n,
			 grid->csn,
			 grid->ncsn,
			 &(grid->nsn),
			 &(grid->nnsn));
  printf("Done\n");

  // Setting the first new node to 1.
  inode = 1;

  // Choose the initial seed node as the current node '1' even though this is probably not the optimal choice.
  /*
  seed = 1;
  nn_map[seed] = inode;
  inode++;
  */

  // Search through the mesh to find the node of highest degree and choose that as the seed node.
  
  seed = 1;                          // Initialize the seed search.
  maxdegree = grid->nnsn[2]-grid->nnsn[1];
  
  for ( i=1; i <= grid->nn; i++ )
    {
      degree = grid->nnsn[i+1] - grid->nnsn[i];
      
      if ( degree > maxdegree )
	{
	  maxdegree = degree;
	  seed = i;
	}
    }

  // Reorder the seed node.
  grid->nn_map[seed] = inode;
  inode++;

  // Indicate that seed has been reordered.
  Queue[seed] = 1;

  // Start the algorithm in earnest.
  while ( inode <= grid->nn )   // There are still nodes to process.
    {
      // Add the surrounding nodes of the seed vertex to R if they haven't been reordered yet.
      start_idx = grid->nnsn[seed];
      end_idx = grid->nnsn[seed+1];

      for ( i=start_idx; i < end_idx; i++ )
	{
	  if ( Queue[grid->nsn[i]] == 0 )
	    {
	      R[endR] = grid->nsn[i];
	      endR++;

	      // Now 'remove' the current node from Queue.
	      Queue[grid->nsn[i]] = 1;
	    }
	}
      
      // Get the degree of the current seed vertex.                        //
      //degree = end_idx - start_idx;                                      // This is handled above.

      // Track the maximum degree.
      //maxdegree = ( degree > maxdegree ) ? degree : maxdegree ;          //

      // Now sort R so that nodes of maximum degree have lower indices (will be popped off sooner).

      // First we make an array of node struct objects and populate it with the nodes that are in the R queue.
      Rnode = (struct node*) malloc( endR*sizeof(struct node) );
      if ( Rnode == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'Rnode'.\n"); exit(1); }

      for ( i=startR; i < endR; i++ )
	{
	  Rnode[i].id = R[i];
	  Rnode[i].degree = grid->nnsn[R[i]+1] - grid->nnsn[R[i]];
	}

      // Sort the node objects by degree.
      qsort ( Rnode , (endR-startR), sizeof(struct node), compare_degree );

      // Copy the ordered set out to the R queue.
      for ( i=startR; i < endR; i++ )
	{
	  R[i] = Rnode[i].id;
	}

      // Free up the memory of Rnode.
      freenull(Rnode);

      // Renumber the nodes in the R queue and pop them off.
      while ( startR < endR )
	{
	  assert( inode <= grid->nn );
	  grid->nn_map[R[startR]] = inode;
	  
	  // Increment the new node index.
	  inode++;
	  
	  // 'Pop' the node off and send it to the back of the S queue.
	  S[endS] = R[startR];
	  endS++;  // Move the position back.
	  startR++;  // Move the position up to 'delete' the current node from R.
	}

      // Reset the indices in R.
      startR = 0;
      endR = 0;

      // Choose the next seed node.
      seed = S[startS];

      // Pop off the new seed node from S.
      startS++;
    }

  // Now nn_map has the new ordering of the vertices. -> nn_map[i] = j;
  // Here, i is the old node number and j is new number for the vertex.

  // Free up memory used locally in this function and class-member data that will be reallocated by main().
  freenull(Queue);
  freenull(R);
  freenull(S);

  freenull(grid->ncsn);
  freenull(grid->csn);
  freenull(grid->nnsn);
  freenull(grid->nsn);

  // Swap out the old node numbering in the c2n map, the boundary connectivity, and the x array.
  tx =  (double*)malloc( (grid->nn+1)*NDIM * sizeof(double) );
  if ( tx == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'tx'.\n"); exit(1); }

  for ( i=1; i <= grid->nn; i++ )
    {
      for ( j=0; j < NDIM; j++ )
	{
	  tx[ (grid->nn_map[i])*NDIM+j ] = grid->x[i*NDIM+j];
	}
    }

  for ( i=1; i <= grid->nn; i++ )
    {
      for ( j=0; j < NDIM; j++ )
	{
	  grid->x[i*NDIM+j] = tx[i*NDIM+j];
	}
    }

  freenull(tx);

  // Swap out element connectivity.
  for ( i=Tri; i <= Quad; i++ )
    {
      // Get the number of nodes for the current element type.
      num_vert = NumberNodesForElement(i);

      for ( j=1; j <= grid->num_elem[i]; j++ )
	{
	  for ( k=0; k < num_vert; k++ )
	    {
	      grid->c2n[i][j*num_vert+k] = grid->nn_map[ grid->c2n[i][j*num_vert+k] ];
	    }
	}
    }

  // Swap out boundary connectivity.
  for ( i=1; i <= grid->nb; i++ )
    {
      for ( j=1; j <= grid->nbs[i]; j++ )
	{
	  grid->bs[i][j][0] = grid->nn_map[ grid->bs[i][j][0] ];
	  grid->bs[i][j][1] = grid->nn_map[ grid->bs[i][j][1] ];
	}
    }

  printf("Cuthill_McKee(): Finished reordering the vertices. Maximum degree in mesh is <%d>.\n",maxdegree);

  // Debug code to print to screen.
  #if 0
  printf("Old vertex id -> Reordered vertex id.\n");
  for ( i=1; i <= grid->nn; i++ )
    {
      printf("%d  ->  %d.\n",i,grid->nn_map[i]);
    }
  #endif
  
  return;
}
