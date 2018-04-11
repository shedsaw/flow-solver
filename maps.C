//=============================================================
// 
//  maps.C
//  
//  Functions that build additional connectivity maps needed
//  for various opertations.
//
//  Written by - Shane Sawyer
//
//=============================================================

#include <stdio.h>
#include <stdlib.h>
#include "grid.h"
#include "mesh_utils.h"
#include "Defined_Variables.h"
#include "maps.h"


//=============================================================
// 
//  Cells_Surrounding_Node()
//
//  Generates a list of cells that surrounds each node in the mesh.
//  
//
//  const int nnodes              // Number of nodes.
//  const int num_elem[]          // Number of elements array.
//  int **c2n                     // Cell to node table.
//  int **csn                     // Pointer to the cell surrounding node array.
//  int **ncsn                    // Pointer to the access array.
//
//=============================================================

void Cells_Surrounding_Node (const int nnodes,
			     int *num_elem,
			     int **c2n,
			     int **csn,
			     int **ncsn)
{
  int i,j,k;                          // Loop counters.
  int nv;                             // Number of vertices for a given element type.
  int id;                             // Node ID.
  int global_count=0;                 // IDs are stored in csn as global numbers. Since each element type has local numbering
                                      // in the c2n array, a global count adjustment will be made as the cell IDs are entered
                                      // into csn. Triangle elements have the same global and local ids, but other element types
                                      // will differ as global_count is increased by the number of elements of each type greater
                                      // than triangles.
  int index;                          // Index of starting position of a node in the csn array.

  // Pointers.
  int *csn_ = NULL;
  int *ncsn_ = NULL;

  // Allocate memory for ncsn array since we know the size will be fixed at nn+2. One +1 is to account for Fortran
  // index and the other +1 is that ncsn[1] and ncsn[2] reference the first node and ncsn[2] and ncsn[3] are for the
  // second node and so on, such that the for the last node we have ncsn[nn] to ncsn[nn+1].
  (*ncsn) = (int*)calloc( (nnodes+2),sizeof(int) );
  if ( (*ncsn) == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'ncsn'.\n"); exit(1); }
  ncsn_ = (*ncsn);

  // Loop over the element types and count the number of elements surrounding each node. This count will be used to size
  // the csn array and apply offsets to search through it.
  for ( i=Tri; i <= Quad; i++ )  // Element type loop.
    {
      // Get the number of nodes for the current element type.
      nv = NumberNodesForElement(i);

      for ( j=1; j <= num_elem[i]; j++ )  // Element loop.
	{
	  for ( k=0; k < nv; k++ )  // Element vertex loop.
	    {
	      // Retrieve the vertex id.
	      id = c2n[i][j*nv+k];

	      // Account for the element in the vertex's ncsn entry. The information is stored 1 up from the actual ID
	      // to make indexing into csn easier in the next loop.
	      ncsn_[id+1]++;

	    }  // End element vertex loop.
	}  // End element loop.
    }  // End element type loop.

  // Set up indexing into csn for all nodes.
  for ( i=1; i <= (nnodes+1); i++ )
    {
      ncsn_[i] += ncsn_[i-1];
    }

  // Now allocate memory for csn since it is now known how big it needs to be.
  (*csn) = (int*)calloc( ncsn_[nnodes+1],sizeof(int) );
  if ( (*csn) == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'csn'.\n"); exit(1); }
  csn_ = (*csn);

  // Go through the above loops again. This time the elements will be entered into csn and the entries in ncsn
  // will be adjusted such that they will appear to have been shifted up one in memory. This is fixed afterwards.
  for ( i=Tri; i <= Quad; i++ )  // Element type loop.
    {
      // Get the number of nodes for the current element type.
      nv = NumberNodesForElement(i);

      for ( j=1; j<= num_elem[i]; j++ )  // Element loop.
	{
	  for ( k=0; k < nv; k++ )  // Element vertex loop.
	    {
	      // Retrieve the vertex id.
	      id = c2n[i][j*nv+k];

	      // Retrieve the first open position in the csn array for the node. Initially, this is the starting
	      // position of the list but is altered as elements are added to the array.
	      index = ncsn_[id];

	      csn_[index] = j + global_count;      // Record the global element id.

	      // Move to the next open position in the list for node 'id.'
	      ncsn_[id]++;

	    }  // End element vertex loop.
	}  // End element loop.

      // Add on the number of elements of the current type so that quad 1 will be numbered as 1+number triangles. This
      // will generate a global numbering scheme.
      global_count += num_elem[i];

    }  // End element type loop.
  
  // Shift the entries down one position in memory.
  for ( i=(nnodes+1); i >= 1; i-- )
    {
      ncsn_[i] = ncsn_[i-1];
    }

  // Now, for element i, the list of surrounding elements can be retrieved as :
  //    index1 = ncsn[i];   index2 = ncsn[i+1];
  //    for (j=index1; j < index2; j++) ...


  // Debug code to print to screen.
  #if DEBUG
  printf("CSN DEBUG:\n");
  for ( i=1; i <= nnodes; i++ )
    {
      int index1 = ncsn_[i];
      int index2 = ncsn_[i+1];
      printf("%d: ",i);
      for ( j=index1; j < index2; j++ )
	{
	  printf("%d ",csn_[j]);
	}
      printf("\n");
    }
  #endif

  return;
}


//=============================================================
//
//  Nodes_Surrounding_Node()
//
//  Generates the map containing nodes surrounding each node.
//
//  const int nnodes              // Number of nodes.
//  const int num_elem[]          // Number of elements array.
//  int **c2n                     // Cell to node table.
//  int *csn                      // Cells surrounding node array.
//  int *ncsn                     // Access array for csn
//  int **nsn                     // Pointer to the node surrounding node array.
//  int **nnsn                    // Pointer to the access array.
//
//=============================================================

void Nodes_Surrounding_Node ( const int nnodes,
			      int *num_elem,
			      int **c2n,
			      int *csn,
			      int *ncsn,
			      int **nsn,
			      int **nnsn)
{
  int i,j,k,l;                        // Loop counters.
  int *list;                          // List of nodes surrounding a given node in a given element. This list is only allocated once
                                      // to save memory allocations since it will be reused many times. It will be length 4 since
                                      // for the current list of allowable element types a given vertex connects to up to 4 other
                                      // vertices.
  int length;                         // Number of entries placed into list[].
  int index;                          // Index of starting position of a node in the nsn array.
  int dim;                            // Initial guess of average number of nodes surrounding a node.
  int inc;                            // Increment to resize the nsn array by if needed.
  int tdim = 0;                       // Keep track of total dimension of nsn being used to realloc to at the end.
  int start_idx;                      // Starting index of node i in csn array;
  int end_idx;                        // Ending index of node i in the csn array (also start of node i+1).
  int elem;                           // Current working element.

  // Pointers.
  int *nsn_ = NULL;
  int *nnsn_ = NULL;
  int *tmp = NULL;

  // Allocate memory for the list of surrounding element vertices.
  list = (int*)malloc(MAX_NUM_VERT*sizeof(int));
  if ( list == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'list'.\n"); exit(1); }

  // Allocate memory for nnsn array since we know the size will be fixed at nn+2. One +1 is to account for Fortran
  // index and the other +1 is that nnsn[1] and nnsn[2] reference the first node and nnsn[2] and nnsn[3] are for the
  // second node and so on, such that the for the last node we have nnsn[nn] to nnsn[nn+1].
  (*nnsn) = (int*)calloc( (nnodes+2),sizeof(int) );
  if ( (*nnsn) == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'nnsn'.\n"); exit(1); }
  nnsn_ = (*nnsn);

  dim = NDIM*(NDIM-1)*2;   // This formula is used since it returns 4 for 2D and 12 for 3D.
  dim = dim*nnodes;

  // Allocate our initial guess for memory requirements for csn. The memory will be adjusted dynamically.
  (*nsn) = (int*)calloc( dim,sizeof(int) );
  if ( (*nsn) == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'nsn'.\n"); exit(1); }
  nsn_ = (*nsn);
 
  // Set the increment to increase the size of nsn to the number of nodes.
  inc = nnodes;

  // Loop over the nodes and start adding the surrounding nodes to the nsn array.
  for ( i=1; i <= nnodes; i++ )   // Node loop.
    {
      // Set the inner loop indices.
      start_idx = ncsn[i];
      end_idx = ncsn[i+1];

      // Set the indices in nnsn so that start/end point to the same memory location.
      nnsn_[i+1] = nnsn_[i];
      
      for ( j=start_idx; j < end_idx; j++ )   // Surrounding element loop.
	{
	  elem = csn[j];
	  
	  // Get the list of surrounding vertices from the current element.
	  RetrieveElementNodesSurroundingVertex ( c2n, num_elem, elem, i, list, &length );

	  // Make sure there is enough room.
	  if ( (tdim+length) > dim )
	    {
	      printf("...Resizing nsn array.\n");
	      dim += inc;
	      tmp = (int*)realloc( (void*)(*nsn) , dim*sizeof(int) );
	      if ( tmp == NULL ) { printf("MEMORY ERROR: COULD NOT REALLOCATE 'nsn'.\n"); exit(1); }
	      (*nsn) = tmp;
	      nsn_ = (*nsn);
	      tmp = NULL;
	    }

	  // Begin inserting unique surrounding nodes into the list.
	  // If this is from the first element in the list, insert all nodes since they are unique.
	  if ( j == start_idx )
	    {
	      for ( k=0; k < length; k++ )   // Surrounding node loop.
		{
		  index = nnsn_[i+1];   // Get the position to insert the node into. Initially its zero.
		  
		  // Insert the current neighboring node.
		  nsn_[index] = list[k];

		  // Update the next position to insert the next node (also acts as the upper limit if no other nodes are inserted).
		  nnsn_[i+1]++;

		  // Update the total number of nodes inserted.
		  tdim++;
		}    // End surrounding node loop.

	      // Continue to the next element.
	      continue;
	    }

	  // Now insert nodes from the element (not the first in list) while checking that only unique nodes are added.
	  for ( k=0; k < length; k++ )   // Surrounding node loop.
	    {
	      index = nnsn_[i+1];   // Get the position to insert the node into. Initially its zero.

	      // Uniqueness check.
	      for ( l=nnsn_[i]; l < index; l++ )   // Current list loop.
		{
		  if ( list[k] == nsn_[l] )
		    break;                       // Break to the next node in the list if the current node to insert is
                                                 // already in the nsn list for node i.
		}    // End current list loop.

	      if ( l < index )  // i.e., We have not traveserd the entire list in the above loop, so we have a
		continue;       // duplicate node, skip to the next one.

	      // If we make down here, then the node has not appeared in the nsn list.
	      
	      // Insert the current neighboring node.
	      nsn_[index] = list[k];

	      // Update the next position to insert the next node (also acts as the upper limit if no other nodes are inserted).
	      nnsn_[i+1]++;

	      // Update the total number of nodes inserted.
	      tdim++;
	    }   // End surrounding node loop.

	}   // End surrounding element loop.

    }   // End node loop.

  // Reduce the excess memory of nsn array.
  printf("Nodes_Surrounding_Node(): Collapsing nsn array from %d to %d.\n",dim,tdim);
  tmp = (int*)realloc( (void*)(*nsn) , tdim*sizeof(int) );
  (*nsn) = tmp;
  nsn_ = (*nsn);
  tmp = NULL;

  // Make an assertion that these two quantities should be equal.
  if ( tdim != nnsn_[nnodes+1] )
    {
      printf("FATAL ERROR: A discrepancy exists in the length of 'nsn' array.\n");
      exit(1);
    }

  // The total number of nodes surrounding vertices should be exactly twice the number of edges in the mesh.
  //  Check that nnsn[nnodes+1] is therefore divisible by two.
  if ( (nnsn_[nnodes+1] % 2 ) != 0 )
    {
      printf("FATAL ERROR: In Nodes_Surrounding_Node(), the total number of nodes surrounding nodes is invalid.\n");
      exit(1);
    }

  // Clean up memory.
  freenull(list);

  // Debug code to print to screen.
  #if DEBUG
  printf("NNSN DEBUG:\n");
  for ( i=0; i <= (nnodes+1); i++ )
    {
      printf("%d\n",nnsn_[i]);
    }

  printf("NSN DEBUG:\n");
  for ( i=1; i <= nnodes; i++ )
    {
      int index1 = nnsn_[i];
      int index2 = nnsn_[i+1];
      printf("%d: ",i);
      for ( j=index1; j < index2; j++ )
	{
	  printf("%d ",nsn_[j]);
	}
      printf("\n");
    }
  #endif
  
  return;
}



//=============================================================
// 
//  Edges_Surrounding_Node()
//
//  Generates the map containing edges surrounding each node.
//  As a consequence, this will also fill in the edges list
//  for the mesh.
//
//  int nnodes                        // Number of nodes.
//  int *nsn                          // Nodes surrounding nodes array.
//  int *nnsn                         // Access array for nsn.
//  int *nedges                       // Number of edges.
//  int **edges                       // Edges array.
//  int **esn                         // Array of edges surrounding node.
//
//=============================================================

void Edges_Surrounding_Node ( int nnodes,
			      int *nsn,
			      int *nnsn,
			      int *nedges,
			      int **edges,
			      int **esn)
{
  int i,j,k;                          // Loop counters.
  int dim;                            // Initial guess of number of edges in the mesh.
  int start_idx;                      // Starting index of node i in csn array;
  int end_idx;                        // Ending index of node i in the csn array (also start of node i+1).
  int left_node;                      // Left node of edge.
  int right_node;                     // Right node of edge.

  // Pointers.
  int *edges_ = NULL;
  int *esn_ = NULL;

  // Initially there are no edges formed.
  *nedges = 0;

  // Set the dimension of edges has half the total length of the nsn array.
  dim = nnsn[nnodes+1] / 2;

  // Allocate the memory.
  (*edges) = (int*)calloc( (dim+1)*2 , sizeof(int) );
  if ( (*edges) == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'edges'.\n"); exit(1); }
  edges_ = (*edges);

  // Allocate space for the edges surrounding nodes. This memory requirement is fixed as the number of nodes surrounding a node
  // since each surrounding node represents an edge.
  (*esn) = (int*)calloc( nnsn[nnodes+1] , sizeof(int) );
  if ( (*esn) == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'esn'.\n"); exit(1); }
  esn_ = (*esn);

  // Loop over the nodes, and then loop over the surrounding nodes to generate a list of edges and the edges
  // surrounding node map.
  for ( i=1; i <= nnodes; i++ )      // Node loop.
    {
      // Set the left node.
      left_node = i;
      
      // Set up the surrounding node indices.
      start_idx = nnsn[i];
      end_idx = nnsn[i+1];

      for ( j=start_idx; j < end_idx; j++ )        // Surrounding node loop.
	{
	  // To make everything consistent and easy, new edges will be created only if the left node is of lower
	  // index than the right node (the surrounding nodes).
	  
	  // Grab the right node.
	  right_node = nsn[j];

	  // Skip to the next node if the edge would be malformed.
	  if ( left_node > right_node )
	    continue;

	  // Make sure that there is no violation of geometric certaintity, i.e. if we are about to create
	  // a new edge that should not exist. 'dim' is the projected number of edges in the mesh. 
	  if ( ((*nedges)+1) > dim )
	    {
	      printf("FATAL ERROR: Edges_Surrounding_Node() tried to create an excess edge.\n");
	      printf("            Diagnostics -> nedges=%d , dim=%d , nnsn[nn+1]=%d\n",(*nedges),dim,nnsn[nnodes+1]);
	      exit(1);
	    }

	  // Increment edge count.
	  (*nedges)++;

	  // Create the edge.
	  edges_[ (*nedges)*2 ] = left_node;
	  edges_[ (*nedges)*2 + 1 ] = right_node;

	  // Insert the edge into the edge surrounding node map for node i into the position that corresponds to the 
	  // right node.
	  esn_[j] = (*nedges);

	  // Now insert the edge into the map for right_node that corresponds to the index of left_node in its nsn map.
	  
	  // First find the position of i in nsn of right_node.
	  for ( k=nnsn[right_node]; k < nnsn[right_node + 1]; k++ )        // Search loop.
	    {
	      if ( nsn[k] == left_node )
		{
		  break;
		}
	    }                                                              // End search loop.

	  // Insert the edge into the map.
	  esn_[k] = (*nedges);

	}                                    // End surrounding node loop.

    }                               // End node loop.

  printf("  Created %d edges.\n",(*nedges));

  // Debug code to print to screen.
  #if DEBUG
  printf("nedges = %d\n\n",*nedges);
  // Print edges.
  for ( i=1; i <=(*nedges); i++ )
    {
      printf("%d: %d -> %d\n",i,edges_[2*i],edges_[2*i+1]);
    }

  // Print surrounding edges.
  printf("ESN DEBUG:\n");
  for ( i=1; i <= nnodes; i++ )
    {
      int index1 = nnsn[i];
      int index2 = nnsn[i+1];
      printf("%d: ",i);
      for ( j=index1; j < index2; j++ )
	{
	  printf("%d ",esn_[j]);
	}
      printf("\n");
    }
  #endif


  return;
}

//=============================================================
// 
//  Build_Cell_to_Edge_Map()
//
//  Generates the map that contains the edges constituting the
//  elements in the mesh.
//
//  int **c2n;                        // Cell to node map.
//  int *num_elem                     // Number of elements in the grid
//  int *nsn                          // Node surrounding node array.
//  int *nnsn                         // Access array for nsn.
//  int *esn                          // Edges surrounding node array.
//  int ***c2e                        // Cell to edge map.
//
//=============================================================

void Build_Cell_to_Edge_Map ( int **c2n,
			      int *num_elem,
			      int *nsn,
			      int *nnsn,
			      int *esn,
			      int ***c2e)
{
  int i,j,k,l;                        // Loop counters.
  int num_edge;                       // Number of edges of an element.
  int left_node;                      // Left node of edge.
  int right_node;                     // Rigth node of edge.
  int index;                          // Index of right_node in nsn of left_node.
  int start_idx;                      // Start index in nsn.
  int end_idx;                        // Ending index in nsn.

  // Pointer.
  int **c2e_ = NULL;

  index = 0;


  // Make sure the memory has not been allocated. If so, clean it out.
  if ( (*c2e) != 0 )
    {
      for (i=Tri; i <= Quad; i++)
	{
	  freenull( (*c2e)[i] );
	}
      
      freenull( (*c2e) );
    }
  
  // Now allocate memory for element types in c2e.
  (*c2e) = (int**)malloc(NUM_ELEM_TYPE*sizeof(int*));
  if ( (*c2e) == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'c2e'.\n"); exit(1); }
  c2e_ = (*c2e);

  // Initialize the pointers to NULL.
  for ( i=Tri; i <= Quad; i++ )
    {
      c2e_[i] = NULL;
    }

  // Allocate memory for the number of elements of each type.
  for ( i=Tri; i <= Quad; i++ )
    {
      num_edge = NumberEdgesForElement(i);
      (*c2e)[i] = (int*)malloc( ((num_elem[i]+1)*num_edge)*sizeof(int) );
      if ( (*c2e)[i] == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'c2e[%d]'.\n",i); exit(1); }
      c2e_[i] = (*c2e)[i];
    }


  // Loop all elements and determine the edge ID of each forming edge.
  for ( i=Tri; i <= Quad; i++ )     //  Element type loop.
    {
      // Get the number of forming edges for the current element type.
      num_edge = NumberEdgesForElement(i);

      for ( j=1; j <= num_elem[i]; j++ )    //  Element loop.
	{
	  // Loop over the forming edges.
	  for ( k=0; k < num_edge; k++ )    //  Element edge loop.
	    {
	      // Get the edge node ID's and have them sorted so that left_node has a lower index
	      // than right_node.
	      RetrieveElementNodesFormingEdge( c2n, j, i, k, &left_node, &right_node, 1);

	      // Since esn array was generated consistently so that the edge connecting
	      // left_node to right_node has the same index as right_node does in the nsn
	      // map portion for left_node, search for the index of right_node in nsn around
	      // left_node.
	      
	      start_idx = nnsn[left_node];
	      end_idx = nnsn[left_node+1];

	      for ( l=start_idx; l < end_idx; l++ )   // NSN search.
		{
		  if ( nsn[l] == right_node )
		    {
		      index = l;
		      break;
		    }
		}                                   // End NSN search.

	      // Store the edge id.
	      c2e_[i][j*num_edge+k] = esn[index];

	    }                            // End element edge loop.
	}                       // End element loop.
    }                     // End element type loop.
	      

  // Debug code to print to screen.
  #if DEBUG
  printf("C2E DEBUG:\n");
  for ( i=Tri; i <= Quad; i++ )
    {
      num_edge = NumberEdgesForElement(i);
      printf("ELEMENT TYPE = %d\n",i);
      for ( j=1; j <= num_elem[i]; j++ )
	{
	  printf("  ELEM ID %d : ",j);
	  for ( k=0; k < num_edge; k++ )
	    {
	      printf("%d ",c2e_[i][j*num_edge+k]);
	    }
	  printf("\n");
	}
      printf("\n");
    }
  #endif


  return;
}



//=============================================================
// 
//  Elements_Surrounding_Element()
//
//  Generates the map that defines the neighboring elements sharing
//  a face for all elements in the mesh.
//
//  int *num_elem                      // Number of elements in the grid.
//  int **c2n                          // Cell to node table.
//  int *csn                           // Cells surrounding node array.
//  int *ncsn                          // Access array for csn.
//  int **ese                          // Elements surrounding element array.
//  int **nese                         // Access array for ese.
//
//=============================================================

void Elements_Surrounding_Element ( int *num_elem,
				    int **c2n,
				    int *csn,
				    int *ncsn,
				    int **ese,
				    int **nese)
{
  int i, j, k, e;                      // Loop counters.
  int total_faces=0;                   // Total number of faces in the mesh which will correspond to the
                                       // number of expected neighbors to be found.
  int total_elems=0;                   // Total number of elements in the mesh.
  int gelem=0;                         // Global element count offset.
  int idx=0;                           // Current position of index into ese.
  int *degree_list;                    // List of all element degrees relative to a particular element.
  int degree;                          // Degree of a single element relative to a particular element.
  int offset=0;                        // Offset for recording neighboring elements into ese.
  int start_idx;                       // Starting index in csn.
  int end_idx;                         // Ending index in csn.
  int *temp_list;                      // List of all elements surrounding the nodes forming an element (will probably
                                       // contain duplicate entries).
  int ntemp_list;                      // The full length of the temp_list array.
  int tlen;                            // The length of temp_list for a particular element.
  int tl_inc;                          // Increment to increase the size of temp_list if it runs out of space on an element.
  int num_vert;                        // Number of vertices for an element type.
  int node;                            // An element node.
  int ielem;                           // Global element ID of the element being processed.
  int deg1;                            // Degree check to determine if we have an element neighbor.
  int deg2;                            // Second degree check. Two are needed since in 3D we could have a triangle face match (3)
                                       // or a quad face match (4). In 2D, its only 2.
  int num_bs;                          // Number of boundary segments.

  // Pointers.
  int *ese_ = NULL;
  int *nese_ = NULL;
  int *tmp;

  // Count up the elements in the mesh.
  for ( e=Tri; e <= Quad; e++ )
    {
      total_elems += num_elem[e];
    }

  // Count up the total number of neighbors that should exist in the mesh.
  for ( e=Tri; e <= Quad; e++ )
    {
      total_faces += num_elem[e] * NumberNeighborsForElement(e);
    }

  // We now have enough information to allocate memory for the ese and nese arrays.

  if ( (*nese) == NULL )
    {
      (*nese) = (int*)malloc( (total_elems+2)*sizeof(int) );
      if ( (*nese) == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'nese'.\n"); exit(1); }
      nese_ = (*nese);
    }

  if ( (*ese) == NULL )
    {
      (*ese) = (int*)calloc( (total_faces),sizeof(int) );
      if ( (*ese) == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'ese'.\n"); exit(1); }
      ese_ = (*ese);
    }

  // Allocate memory for the temp_list and for degree_list.
  degree_list = (int*)calloc( total_elems+1 , sizeof(int) );
  if ( degree_list == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'degree_list'.\n"); exit(1); }

  ntemp_list = 1000;
  temp_list = (int*)malloc( ntemp_list*sizeof(int) );
  if ( temp_list == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'temp_list'.\n"); exit(1); }

  // 1000 should be enough, but just in case set up an increment to increase the size if needed.
  tl_inc = 500;

  // Set up nese to point into ese (using global element ids).
  for ( e=Tri; e <= Quad; e++ )                              // Element type loop.
    {
      for ( i=1; i <= num_elem[e]; i++ )                     // Element loop.
	{
	  // Set the pointer up.
	  nese_[i+gelem] = idx;

	  // Adjust the idx index position.
	  idx += NumberNeighborsForElement(e);
	}                                                    // End element loop.

      // Adjust the global element count offset.
      gelem += num_elem[e];
    }                                                        // End element type loop.

  // Set up the last index pointer.
  nese_[total_elems+1] = total_faces;

  // Now go through and fill in the ese array with the element neighbors.
  for ( e=Tri; e <= Quad; e++ )                             // Element type loop.
    {
      // Set up the degree match limits for the element type.
      deg1 = 2;  // Always the case in 2D.
      deg2 = 2;

      // ** MAY NEED THIS BUT I THINK IT WORKS FOR 3D!
      //if ( e == Tri )
      //{
      //      deg1 = 3;
      //      deg2 = 3;
      //    }
      //  else if ( e == Quad )
      //    {
      //      deg1 = 4;
      //      deg2 = 4;
      //    }

      // Get the number of element vertices.
      num_vert = NumberNodesForElement(e);

      for ( i=1; i <= num_elem[e]; i++ )                    // Element loop.
	{
	  // Set the temp_list length to 0 for the element.
	  tlen = 0;

	  // Get the global element ID that we are processing.
	  ielem = LocalToGlobalID(num_elem,i,e);

	  for ( j=0; j < num_vert; j++ )                    // Element node loop.
	    {
	      // Grab the next element vertex.
	      node = c2n[e][i*num_vert+j];

	      // Set up the search indices in the csn array.
	      start_idx = ncsn[node];
	      end_idx =   ncsn[node+1];

	      // Loop through the surrounding elements and add them to the temp_list.
	      for ( idx=start_idx; idx < end_idx; idx++ )   // Surrounding element loop.
		{
		  // Grab the element (global ID).
		  gelem = csn[idx];

		  // Record the element and increment the length of the list.
		  temp_list[tlen++] = gelem;

		  // Adjust the degree of the surrounding element relative to the element.
		  degree_list[gelem]++;

		  // Now check that we haven't exceeded our memory limits on temp_list. If so, increase our limits.
		  if ( tlen >= ntemp_list )
		    {
		      ntemp_list += tl_inc;
		      tmp = (int*)realloc( (void*)temp_list , ntemp_list*sizeof(int) );
		      if ( tmp == NULL ) { printf("MEMORY ERROR: COULD NOT REALLOCATE 'temp_list'.\n"); exit(1); }
		      temp_list = tmp;
		      tmp = NULL;
		    }
		}                                           // End surrounding element loop.
	    }                                               // End element node loop.

	  // Set the offset for the element to 0;
	  offset = 0;

	  // Now go through the temp_list and add the neighbors to the ese list for that element.
	  for ( j=0; j < tlen; j++ )                       // temp_list loop.
	    {
	      gelem = temp_list[j];
	      degree = degree_list[gelem];

	      // Skip element that have 0 degree or if its the element we are currently processing.
	      if ( degree == 0 || gelem == ielem )
		continue;

	      if ( degree == deg1 || degree == deg2 )  // We have a face match (hopefully).
		{
		  ese_[nese_[ielem] + offset] = gelem;
		  offset++;
		  degree_list[gelem] = 0;  // Remove it from further consideration.
		}

	      // Sanity check.
	      if ( offset > NumberNeighborsForElement(e) )
		{
		  printf("\nFATAL ERROR: In Elements_Surrounding_Element(): A topological error was encountered in the mesh.\n");
		  printf("             For element type <%d>, element id <%d>, an excess neighbor was added.\n",e,i);
		  exit(1);
		}

	      // At this point, I'm assuming that all elements have a full neighbor listing, so we are done.
	      // I'm assuming here that any 0's are boundary segments! May cause trouble later!!!
	    }                                             // End temp_list loop.

	  // Now go back through and clean up the degree of the elements in temp_list to prevent neighbor corruption.
	  for ( j=0; j < tlen; j++ )
	    {
	      gelem = temp_list[j];
	      degree_list[gelem] = 0;
	    }

	}                                                 // End element loop.

    }                                                     // End element type loop.

  // Do a final sanity check. Make sure that all elements have their full list of neighbors.
  #if 1
  num_bs = 0;
  for ( e=Tri; e <= Quad; e++ )
    {
      for ( i=1; i <= num_elem[e]; i++ )
	{
	  j = LocalToGlobalID(num_elem,i,e);
	  start_idx = nese_[j];
	  end_idx   = nese_[j+1];

	  for ( k=start_idx; k < end_idx; k++ )
	    {
	      if ( ese_[k] == 0 )
		{
		  //printf("\nFATAL ERROR: In Elements_Surrounding_Element(): A topological error was encountered in the mesh.\n");
		  //printf("             For element type <%d>, element id <%d>, position <%d>, a neighbor was not found.\n",e,i,k);
		  //exit(1);
		  num_bs++;
		}
	    }
	}
    }
  printf("Detected %d boundary segments in the grid.\n",num_bs);
  #endif

  // Clean up memory.
  freenull( temp_list );
  freenull( degree_list );
	      
  // Debug code to print to screen.
  #if DEBUG
  printf("ESE DEBUG:\n");
  for ( e=Tri; e <= Quad; e++ )
    {
      printf("ELEMENT TYPE %d:\n",e);
      for ( i=1; i <= num_elem[e]; i++ )
	{
	  j = LocalToGlobalID(num_elem,i,e);
	  int index1 = nese_[j];
	  int index2 = nese_[j+1];
	  printf("%d: ",i);
	  for ( k=index1; k < index2; k++ )
	    {
	      printf("%d ",ese_[k]);
	    }
	  printf("\n");
	}
    }
  #endif

  return;
}



//=============================================================
//
//  Nodes_Surrounding_Node2()
//
//  Generates the map containing nodes surrounding each node
//  including first and second degree neighbors.
//
//  const int nnodes              // Number of nodes.
//  const int num_elem[]          // Number of elements array.
//  int **c2n                     // Cell to node table.
//  int *esn                      // Edges surrounding node array.
//  int *nnsn                     // Access array for nsn and esn.
//  int *edges                    // Edges data array.
//  int **nsn2                    // Pointer to the nodes surrounding node2 array.
//  int **nnsn2                   // Pointer to the access array.
//
//=============================================================

void Nodes_Surrounding_Node2 ( const int nnodes,
			       int *num_elem,
			       int **c2n,
			       int *esn,
			       int *nnsn,
			       int *edges,
			       int **nsn2,
			       int **nnsn2)
{
  int i,j,k;                                    // Loop counters.
  int candidate_node;                           // A potential neighbor.
  int num_neighbors = 0;                        // The number of neighbors for a node.
  int flag = 0;                                 // Logic flag.
  int neighbor;                                 // Neighboring node.
  int mem_alloc = 0;                            // Memory already allocated.
  int node1, node2, iedge;                      // Global edge and its nodes.
  int *list = NULL;                             // List of neighbors.
  int total_mem_needed = 0;                     // Total amount of memory needed.
  int ptr;                                      // Array index pointer.
  int offset;                                   // Memory offset into array.
  int looplimit;                                // Memory bound for array search.

  // Pointers.
  int *nsn2_ = NULL;
  int *nnsn2_ = NULL;
  int *tmp = NULL;

  // The procedure I'm going to use is to loop over all the nodes
  // once to calculate how many neighbors I have for each one by
  // storing all unique neighbors in a junk array. Then allocate 
  // memory and go through the process again this time keeping the
  // unique neighbors in the nsn2 array.

  //Start by allocating a chunk of memory.
  mem_alloc = 100;
  list = (int*)malloc(mem_alloc*sizeof(int));
  if ( list == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'list'.\n"); exit(1); }

  // Loop over all the nodes in the mesh.
  for ( i=1; i <= nnodes; i++ )                                                      // Node loop.
    {
      num_neighbors = 0;

      // Loop over and store all the first degree neighbors.
      for ( j=nnsn[i]; j < nnsn[i+1]; j++ )                                          // First degree neighbor loop.
	{
	  if ( num_neighbors >= mem_alloc )
	    {
	      mem_alloc += 50;
	      tmp = (int*)realloc( (void*)list , mem_alloc*sizeof(int) );
	      if ( tmp == NULL ) { printf("MEMORY ERROR: COULD NOT REALLOCATE 'list'.\n"); exit(1); }
	      list = tmp;
	      tmp = NULL;
	    }

	  // Store the one that's not the node itself.
	  iedge = esn[j];
	  node1 = edges[iedge*2+0];
	  node2 = edges[iedge*2+1];

	  if ( node1 != i )
	    {
	      list[num_neighbors] = node1;
	      num_neighbors++;
	    }
	  else
	    {
	      list[num_neighbors] = node2;
	      num_neighbors++;
	    }
	}                                                                            // End first degree neighbor loop.

      looplimit = num_neighbors;

      // Now we loop over all the first degree neighbors and add their
      // respective neighbors if they haven't already appeared in the list.
      for ( j=0; j < looplimit; j++ )                                                // First degree neighbor search.
	{
	  // Grab a neighbor.
	  candidate_node = list[j];
	  
	  // Loop over its neighbors.
	  for ( ptr=nnsn[candidate_node]; ptr < nnsn[candidate_node+1]; ptr++ )      // Second degree neighbor search.
	    {
	      // Grab a candidate node.
	      iedge = esn[ptr];
	      node1 = edges[iedge*2+0];
	      node2 = edges[iedge*2+1];
	
	      if ( candidate_node != node1 )
		neighbor = node1;
	      else
		neighbor = node2;

	      // See if its in the list.
	      for ( k=0; k < num_neighbors; k++ )
		{
		  if ( neighbor == list[k] || neighbor == i )
		    {
		      flag = 1;
		      break;
		    }
		}
	      if ( flag == 0 )
		{
		  if ( num_neighbors >= mem_alloc )
		    {
		      mem_alloc += 50;
		      tmp = (int*)realloc( (void*)list , mem_alloc*sizeof(int) );
		      if ( tmp == NULL ) { printf("MEMORY ERROR: COULD NOT REALLOCATE 'list'.\n"); exit(1); }
		      list = tmp;
		      tmp = NULL;
		    }

		  list[num_neighbors] = neighbor;
		  num_neighbors++;
		}
	      flag = 0;
	    }                                                                        // End second degree neighbor search.
	}                                                                            // End first degree neighbor search.
      total_mem_needed += num_neighbors;
    }                                                                                // End node loop.

  // Now we know how much space is needed for the map.
  // Do our allocations and now store the neighbors and and set up pointers.
  (*nnsn2) = (int*)calloc( (nnodes+2),sizeof(int) );
  if ( (*nnsn2) == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'nnsn2'.\n"); exit(1); }
  nnsn2_ = (*nnsn2);

  (*nsn2) = (int*)calloc( total_mem_needed,sizeof(int) );
  if ( (*nsn2) == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'nsn2'.\n"); exit(1); }
  nsn2_ = (*nsn2);

  nnsn2_[1] = 0;

  // Loop over all the nodes in the mesh and keep the results this time.
  for ( i=1; i <= nnodes; i++ )
    {
      num_neighbors = 0;
      offset = nnsn2_[i];

      // Loop over and store all the first degree neighbors.
      for ( j=nnsn[i]; j < nnsn[i+1]; j++ )
	{
	  // Store the one that's not the node itself.
	  iedge = esn[j];
	  node1 = edges[iedge*2+0];
	  node2 = edges[iedge*2+1];
	  
	  if ( node1 != i )
	    {
	      nsn2_[offset + num_neighbors] = node1;
	      num_neighbors++;
	    }
	  else
	    {
	      nsn2_[offset + num_neighbors] = node2;
	      num_neighbors++;
	    }
	}
      looplimit = num_neighbors;
      
      // Now we loop over all the first degree neighbors and add their
      // respective neighbors if they haven't already appeared in the list.
      for ( j=0; j < looplimit; j++ )
	{
	  // Grab a neighbor.
	  candidate_node = nsn2_[offset + j];

	  // Loop over its neighbors.
	  for ( ptr=nnsn[candidate_node]; ptr < nnsn[candidate_node+1]; ptr++ )
	    {
	      // Grab a candidate node.
	      iedge = esn[ptr];
	      node1 = edges[iedge*2+0];
	      node2 = edges[iedge*2+1];

	      if ( candidate_node != node1 )
		neighbor = node1;
	      else
		neighbor = node2;

	      // See if its in the list.
	      for ( k=0; k < num_neighbors; k++ )
		{
		  if ( neighbor == nsn2_[offset + k] || neighbor == i )
		    {
		      flag = 1;
		      break;
		    }
		}

	      // If we haven't found a  match then store the node.
	      if ( flag == 0 )
		{
		  nsn2_[offset + num_neighbors] = neighbor;
		  num_neighbors++;
		}
	      flag = 0;
	    }
	}
      // Set the pointer.
      nnsn2_[i+1] = offset + num_neighbors;
    }

  // Quick sanity check.
  if ( nnsn2_[nnodes+1] != total_mem_needed )
    {
      printf("FATAL ERROR: A topological error was encountered in Nodes_Surrounding_Node2().\n");
      printf("  nnsn2_[nnodes+1] = %d , total_mem_needed = %d\n",nnsn2_[nnodes+1],total_mem_needed);
      exit(0);
    }

  // Free.
  freenull(list);

  // Debug code to print to screen.
  #if DEBUG
  printf("NSN2 DEBUG:\n");
  for ( i=1; i <= nnodes; i++ )
    {
      int index1 = nnsn2_[i];
      int index2 = nnsn2_[i+1];
      printf("%d: ",i);
      for ( j=index1; j < index2; j++ )
	{
	  printf("%d ",nsn2_[j]);
	}
      printf("\n");
    }
  #endif

  return;
}


//=============================================================
//
//  Initialize_Generic_Nodes_Surrounding_Node()
//
//  Initializes a node surrounding node list. Fill in with the
//  node itself only.
//
//  const int nnodes              // Number of nodes.
//  const int num_elem[]          // Number of elements array.
//  int **c2n                     // Cell to node table.
//  int *esn                      // Edges surrounding node array.
//  int *nnsn                     // Access array for nsn and esn.
//  int *edges                    // Edges data array.
//  int **gnsn                    // Pointer to the array.
//  int **gnnsn1                  // Pointer array to the array that can access the array. Only stores the start.
//  int **gnnsn                   // Pointer to the access array.
//  int **node_fill_level         // How many layers of neighbors are around node.
//=============================================================

void Initialize_Generic_Nodes_Surrounding_Node ( const int nnodes,
						 int *num_elem,
						 int **c2n,
						 int *esn,
						 int *nnsn,
						 int *edges,
						 int **gnsn,
						 int **gnnsn,
						 int **gnnsn1,
						 int **node_fill_level )
{
  int i,j,k;                                    // Loop counters.
  int ptr;                                      // Array index pointer.
  int sidx, eidx;                               // Array indices.

  // Pointers.
  int *gnsn_ = NULL;
  int *gnnsn_ = NULL;
  int *gnnsn1_ = NULL;
  int *fill_ = NULL;

  // Allocate memory for gnnsn array since we know the size will be fixed at nn+2. One +1 is to account for Fortran
  // index and the other +1 is that gnnsn[1] and gnnsn[2] reference the first node.

  // This array will never be realloced, only the contents will be adjusted.
  (*gnnsn1) = (int*)calloc( (nnodes+2),sizeof(int) );
  if ( (*gnnsn1) == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'gnnsn1'.\n"); exit(1); }
  gnnsn1_ = (*gnnsn1);

  (*gnnsn) = (int*)calloc( (nnodes+2),sizeof(int) );
  if ( (*gnnsn) == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'gnnsn'.\n"); exit(1); }
  gnnsn_ = (*gnnsn);

  // Allocate our initial array for the nodes themselves.
  (*gnsn) = (int*)calloc( (nnodes+1),sizeof(int) );
  if ( (*gnsn) == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'gnsn'.\n"); exit(1); }
  gnsn_ = (*gnsn);

  // Allocate space for the node_fill_level array. It is a fixed size independent of the number of layers.
  (*node_fill_level) = (int*)malloc( (nnodes+1) * sizeof(int) );
  if ( (*node_fill_level) == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'node_fill_level'.\n"); exit(1); }
  fill_ = (*node_fill_level);
  
  for ( i=0; i <= nnodes; i++ )
    {
      fill_[i] = 0;
    }
  
  ptr = 0;  // Set the pointer to the beginning of the array.

  // Put in the nodes.
  for ( i=1; i <= nnodes; i++ )
    {
      gnsn_[ptr] = i;
      gnnsn_[i] = ptr;
      gnnsn1_[i] = i;  // initially it points to the address in gnnsn that == the node number.

      ptr++;
    }

  // last entry. Here i=nnodes+1
  gnnsn_[i] = ptr;
  gnnsn1_[i] = i;

  // Sanity check.
  for ( i=1; i <= nnodes; i++ )
    {
      ptr = gnnsn1_[i];
      ptr = gnnsn_[ptr];

      if ( gnsn_[ptr] != i )
	{
	  printf("FATAL ERROR: In Initialize_Generic_Nodes_Surrounding_Node(), node %d did not equal itself in gns",i);
	  printf(" -> ( %d ).\n",gnsn_[ptr] );
	  fflush(stdout);
	  exit(1);
	}
    }
     
  #if DEBUG
  printf("Initialize_Generic_Nodes_Surrounding_Node(): GNNSN1 DEBUG:\n");
  for ( i=0; i <= (nnodes+1); i++ )
    {
      printf("%d\n",gnnsn1_[i]);
    }

  printf("Initialize_Generic_Nodes_Surrounding_Node(): GNNSN DEBUG:\n");
  for ( i=0; i <= (nnodes+1); i++ )
    {
      printf("%d\n",gnnsn_[i]);
    }

  printf("Initialize_Generic_Nodes_Surrounding_Node(): GNSN DEBUG:\n");
  for ( i=1; i <= nnodes; i++ )
    {
      sidx = gnnsn1_[i];
      eidx = gnnsn1_[i+1];
      sidx = gnnsn_[sidx];
      eidx = gnnsn_[eidx];
      printf("%d: ",i);
      for ( j=sidx; j < eidx; j++ )
	{
	  printf("%d ",gnsn_[j]);
	}
      printf("\n");
    }
  fflush(stdout);
  #endif

  return;
}



//=============================================================
//
//  Add_Layer_to_GNSN()
//
//  Adds a layer of neighbors to the given node into gnsn and adjusts
//  gnnsn and num_fill_level.
//
//  const int nnodes              // Number of nodes.
//  const int num_elem[]          // Number of elements array.
//  int **c2n                     // Cell to node table.
//  int *esn                      // Edges surrounding node array.
//  int *nnsn                     // Access array for nsn and esn.
//  int *edges                    // Edges data array.
//  int **gnsn                    // Pointer to the nodes surrounding node array.
//  int **gnnsn                   // Pointer to the access array.
//  int *gnnsn1                   // Pointer to the array that accesses the array.
//  int *node_fill_level          // Fill level of each node.
//  int the_node                  // The node to which we are adding a layer.
//
//=============================================================

void Add_Layer_to_GNSN ( const int nnodes,
			 int *num_elem,
			 int **c2n,
			 int *esn,
			 int *nnsn,
			 int *edges,
			 int **gnsn,
			 int **gnnsn,
			 int *gnnsn1,
			 int *node_fill_level,
			 int the_node )
{
  int i,j,k;                                    // Loop counters.
  int candidate_node;                           // A potential neighbor.
  int num_neighbors = 0;                        // The number of neighbors for a node.
  int flag = 0;                                 // Logic flag.
  int neighbor;                                 // Neighboring node.
  int mem_alloc = 0;                            // Memory already allocated.
  int node1, node2, iedge;                      // Global edge and its nodes.
  int *list = NULL;                             // List of neighbors to consider adding.
  int *store_list = NULL;                       // List of neighbors we will add.
  int list_ptr;                                 // Pointer to the list that limits the search.
  int looplimit;                                // Memory bound for array search.
  int old_fill_level;                           // The previous fill level.
  int new_fill_level;                           // The new fill level.
  int search_index_start;                       // Array indices.
  int search_index_stop;
  int dim;

  // Pointers. These will be adjusted dynamically.
  int *gnsn_ = NULL;
  int *gnnsn_ = NULL;
  int *tmp = NULL;

  // Set the pointers for easy access.
  gnsn_ = (*gnsn);
  gnnsn_ = (*gnnsn);

  old_fill_level = node_fill_level[the_node];
  new_fill_level = old_fill_level + 1;

  //Start by allocating a chunk of memory.
  mem_alloc = 100;
  list = (int*)malloc(mem_alloc*sizeof(int));
  if ( list == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'list'.\n"); exit(1); }

  list_ptr = 0;

  // The first step is to figure out how many unique neighbors are we going to be adding to the stencil for the node.
  
  // Make a list of the candidate nodes. These are nodes that are neighbors to the_node's largest current fill level.
  // So if the_node already has first degree neighbors, the candidate list will contain first degree neighbors of those first
  // degree neighbors.

  // The original fill level was given above.
  search_index_start = gnnsn1[the_node] + old_fill_level;

  search_index_stop  = gnnsn1[the_node + 1];

  num_neighbors = 0;

  // Populate the list.
  for ( i=gnnsn_[search_index_start]; i < gnnsn_[search_index_stop]; i++ )
    {
      neighbor = gnsn_[i];

      // Loop over the neighbor's first degree neighbors.
      for ( j=nnsn[neighbor]; j < nnsn[neighbor+1]; j++ )
	{
	  if ( num_neighbors >= mem_alloc )
	    {
	      mem_alloc += 50;
	      tmp = (int*)realloc( (void*)list , mem_alloc*sizeof(int) );
	      if ( tmp == NULL ) { printf("MEMORY ERROR: COULD NOT REALLOCATE 'list'.\n"); exit(1); }
	      list = tmp;
	      tmp = NULL;
	    }      
      
	  // Store the one that's not the node itself.
	  iedge = esn[j];
	  node1 = edges[iedge*2+0];
	  node2 = edges[iedge*2+1];
      
	  if ( node1 != neighbor )
	    {
	      list[num_neighbors] = node1;
	      num_neighbors++;
	    }
	  else
	    {
	      list[num_neighbors] = node2;
	      num_neighbors++;
	    }
	  
	}
    }

  looplimit = num_neighbors;

  // Now we have a list of candidates (that is probably going to contain duplicates). Now lets loop over them and keep the ones
  // are actual valid neighbors of the desired fill level.

  store_list = (int*)malloc(num_neighbors*sizeof(int));  // Maximum possible size.
  if ( store_list == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'store_list'.\n"); exit(1); }

  num_neighbors = 0;  // num_neighbors now refers to the number of unique neighbors we will add.

  search_index_start = gnnsn1[the_node];
  search_index_stop  = gnnsn1[the_node+1];

  for ( i=0; i < looplimit; i++ )
    {
      // Grab the candidate. 
      candidate_node = list[i];
      flag = 0;

      // Is this neighbor unique?
      for ( j=gnnsn_[search_index_start]; j < gnnsn_[search_index_stop]; j++ )
	{
	  node1 = gnsn_[j];

	  if ( node1 == candidate_node )
	    {
	      flag = 1;
	      break;
	    }
	}

      for ( j=0; j < num_neighbors; j++ )
	{
	  node1 = store_list[j];

	  if ( node1 == candidate_node )
	    {
	      flag = 1;
	      break;
	    }
	}

      // If the flag is set to one, we know that the neighbor is a duplicate.
      if ( flag == 0 )
	{
	  store_list[num_neighbors] = candidate_node;
	  num_neighbors++;
	}
    }

  // Now we have the list of nodes to add and how many we are adding. Time to start adjusting pointers and
  // reallocating memory.

  // find the size of the gnsn array. That is to say what is the last memory spot addressed.
  dim = gnnsn_[ gnnsn1[nnodes+1] ] + num_neighbors;

  tmp = (int*)realloc( (void*)(*gnsn) , (dim)*sizeof(int) );
  if ( tmp == NULL ) {  printf("MEMORY ERROR: COULD NOT REALLOCATE 'gnsn'.\n"); exit(1); }
  (*gnsn) = tmp;
  gnsn_ = (*gnsn);
  tmp = NULL;
  
  // find the size of the gnnsn array.
  dim = gnnsn1[nnodes+1] + 1 + 1;  // adding a single pointer to the array.
  
  tmp = (int*)realloc( (void*)(*gnnsn) , (dim)*sizeof(int) );                  // had this has gnnsn1[] and then gnnsn_[]
  if ( tmp == NULL ) {  printf("MEMORY ERROR: COULD NOT REALLOCATE 'gnnsn'.\n"); exit(1); }
  (*gnnsn) = tmp;
  gnnsn_ = (*gnnsn);
  tmp = NULL;
  
  // Shift everything in gnsn corresponding to the_node+1 and above down in memory by num_neighbors.
  search_index_start = gnnsn_[ gnnsn1[nnodes+1] ] -  1;                  // had a - 1 at the end.
  search_index_stop  = gnnsn_[ gnnsn1[the_node+1] ];

  for ( i=search_index_start; i >= search_index_stop; i-- )
    {
      gnsn_[i+num_neighbors] = gnsn_[i];
    }

  // Now put in the new neighbors.
  //search_index_start = gnnsn_[ gnnsn1[the_node] ] + 1;  // gnnsn_[ gnnsn1[the_node+1] ];
  search_index_start = gnnsn_[ gnnsn1[the_node+1] ];

  for ( i=0; i < num_neighbors; i++ )
    {
      gnsn_[search_index_start + i] = store_list[i];
    }

  // Now we shift everything in gnnsn down 1 and adjust the pointers along the way.
  search_index_start = gnnsn1[nnodes+1];  // had - 1 here.
  search_index_stop  = gnnsn1[the_node+1];

  for ( i=search_index_start; i >= search_index_stop; i-- )
    {
      gnnsn_[i+1] = gnnsn_[i] + num_neighbors;    // this was gnsn_ on both side of =
    }

  // Last thing to do now is to adjust the pointers in gnnsn1.
  for ( i=(the_node+1); i <= (nnodes+1); i++ )
    {
      gnnsn1[i] += 1;
    }

  node_fill_level[the_node] += 1;
  
  freenull(list);
  freenull(store_list);

  if ( DEBUG )
    {
      printf("FINISHED PROCESSING NODE %d\n",the_node);

      printf("GNNSN1: ");
      for ( i=0; i <= nnodes+1; i++ )
	printf("%d ",gnnsn1[i]);
      printf("\n\n");

      search_index_start = 0;
      search_index_stop  = gnnsn1[nnodes+1];

      printf("GNNSN: ");
      for ( i=search_index_start; i <= search_index_stop; i++ )
	printf("%d ",gnnsn_[i]);
      printf("\n\n");

      search_index_start = gnnsn_[ gnnsn1[1] ];
      search_index_stop  = gnnsn_[ gnnsn1[nnodes+1] ];

      printf("GNSN: ");
      for ( i=search_index_start; i < search_index_stop; i++ )
	printf("%d ",gnsn_[i]);
      printf("\n\n");

      fflush(stdout);
    }

  return;
}
