//=============================================================
// 
//  reconstruct.h
//  
//  Function prototypes for computing the solution reconstruction.
//
//  Written by - Shane Sawyer
//
//=============================================================

#include <stdio.h>
#include <stdlib.h>
#include "Defined_Variables.h"
#include "grid.h"
#include "params.h"
#include "sys_mem.h"

#ifndef RECONSTRUCT_H
#define RECONSTRUCT_H

void Reconstruct ( int isedge, int nodeL, int nodeR, double QL[4], double QR[4], GRID *grid, PARAMS p, int ts );

void generic_reconstruct ( int node, double Q[4], double dx[2], GRID *grid, PARAMS p );

void Reconstruct_Gauss ( int isedge, int nodeL, int nodeR, double GP[2], double QL[4], double QR[4], GRID *grid, PARAMS p, int ts );

void Reconstruct_Gauss_Boundary ( int nodeL, double GP[2], double QL[4], GRID *grid, PARAMS p );

// Functions to reconstruct with no limiter applied.

void Reconstruct_UL ( int isedge, int nodeL, int nodeR, double QL[4], double QR[4], GRID *grid, PARAMS p, int ts );

void generic_reconstruct_UL ( int node, double Q[4], double dx[2], GRID *grid, PARAMS p );

void Reconstruct_Gauss_UL ( int isedge, int nodeL, int nodeR, double GP[2], double QL[4], double QR[4], GRID *grid, PARAMS p, int ts );

void Reconstruct_Gauss_Boundary_UL ( int nodeL, double GP[2], double QL[4], GRID *grid, PARAMS p );

#endif
