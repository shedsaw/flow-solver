//=============================================================
// 
//  limiter.h
//  
//  Function prototypes for computing the reconstruction limiter.
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

#ifndef LIMITER_H
#define LIMITER_H

void compute_limiter_barth ( GRID *grid, PARAMS p );

void compute_limiter_venkatakrishnan ( GRID *grid, PARAMS p );

void compute_limiter_venkatakrishnan2 ( GRID *grid, PARAMS p );

void compute_limiter_gooch1 ( GRID *grid, PARAMS p );

void compute_limiter_gooch2 ( GRID *grid, PARAMS p );

#endif
