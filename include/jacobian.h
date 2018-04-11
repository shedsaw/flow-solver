//=============================================================
// 
//  jacobian.h
//  
//  Function prototypes for computing the flux Jacobian.
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

#ifndef JACOBIAN_H
#define JACOBIAN_H

void compute_numerical_jacobian_roe ( GRID *grid, SYS_MEM *smem, PARAMS p, int nodeL, int nodeR,
				      double QL[4], double QR[4], double flux[4],
				      double nx, double ny, double len );

void compute_boundary_numerical_jacobian_roe ( GRID *grid, SYS_MEM *smem, PARAMS p, int ibedge,
					       double QL[4], double flux[4], double nx, double ny, double len );

void compute_approximate_jacobian_roe ( GRID *grid, SYS_MEM *smem, PARAMS p, int nodeL, int nodeR,
					double QL[4], double QR[4], double flux[4],
					double nx, double ny, double len );

void compute_boundary_approximate_jacobian_roe ( GRID *grid, SYS_MEM *smem, PARAMS p, int ibedge, double QL[4],
						 double QR[4], double flux[4], double nx, double ny, double len );


#endif


