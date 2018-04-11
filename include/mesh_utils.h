//=============================================================
// 
//  mesh_utils.h
//  
//  Function prototypes for utility functions.
//
//  Written by - Shane Sawyer
//
//=============================================================

#include <stdio.h>
#include <stdlib.h>
#include "Defined_Variables.h"
#include "grid.h"

#ifndef MESH_UTILS_H
#define MESH_UTILS_H


void rewind2D ( GRID *grid );

int NumberNodesForElement ( int e );

int NumberEdgesForElement ( int e );

int NumberNeighborsForElement ( int e );

void RetrieveElementNodesSurroundingVertex ( int **c2n,
					     int *num_elem,
					     int elem_id,
					     int node_id,
					     int *list,
					     int *length );

void RetrieveElementNodesFormingEdge ( int **c2n,
				       int elem_id,
				       int type,
				       int edge_id,
				       int *left,
				       int *right,
				       int order );

int GlobalToLocalID ( int *num_elem,
		      int elem_id,
		      int *e );

int LocalToGlobalID ( int *num_elem,
		      int elem_id,
		      int e );

void Swap ( int *nodeL, int *nodeR );

void ConvertConservedToPrimitive ( double gamma , double *Q );

void ConvertPrimitiveToConserved ( double gamma , double *Q );

#endif
