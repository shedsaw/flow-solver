//=============================================================
// 
//  init_memory.h
//  
//  Function prototypes for initializing the memory needed to
//  solve the problem.
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

#ifndef INIT_MEM_H
#define INIT_MEM_H

void Build_CRS ( GRID *grid, PARAMS p, SYS_MEM *smem );

void build_matrix ( GRID *grid, PARAMS p, SYS_MEM *smem );

#endif
