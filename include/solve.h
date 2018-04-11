//=============================================================
// 
//  solve.h
//  
//  Function prototypes for solving the equations.
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

#ifndef SOLVE_H
#define SOLVE_H

double solve ( GRID *grid, SYS_MEM *smem, PARAMS p, int ts );

void Compute_Extrema ( GRID *grid, PARAMS p, int m );

void debug_flux_integral ( GRID *grid, PARAMS p );

void flux_divergence ( double *x, double *q );

void MMS_Compute_Q_Exact ( GRID * grid );

void MMS_Compute_Source_Terms ( GRID * grid );

void MMS_Get_Point_Qe ( double *x, double *q );

void MMS_Get_Flux_Divergence ( double *x, double *q );

void MMS_Compute_Source_Terms_Exp ( GRID *grid );

#endif
