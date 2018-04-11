//=============================================================
// 
//  residual.h
//  
//  Function prototypes for computing the residual.
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

#ifndef RESIDUAL_H
#define RESIDUAL_H

void build_interior_residuals ( GRID *grid, SYS_MEM *smem, PARAMS p, int jupdate, int ts );

void build_boundary_residuals ( GRID *grid, SYS_MEM *smem, PARAMS p, int jupdate );

// HACK
void fix_farfield_boundary ( GRID *grid, SYS_MEM *smem, PARAMS p );
void do_hard_boundary_flux_jac ( GRID *grid, SYS_MEM *smem, PARAMS p );
void do_soft_boundary_flux_jac ( GRID *grid, SYS_MEM *smem, PARAMS p );

#endif
