//=============================================================
// 
//  cvbc.h
//  
//  Function prototypes for computing the characteristic value
//  boundary conditions.
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

#ifndef CVBC_H
#define CVBC_H

void Set_Boundary_Values ( GRID *grid, PARAMS p );

void Get_Boundary_Value ( GRID *grid, PARAMS p, int isubedge, double Qi[4], double Qb[4] );

#endif
