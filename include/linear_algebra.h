#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Defined_Variables.h"
#include "params.h"
#include "grid.h"
#include "sys_mem.h"

void crossProduct ( double *u, double *v, double *w );

double dotProduct ( double *u, double *v );

double dotProductGeneral ( double *u, double *v, int n );

void normalizeVectorGeneral ( double *u, int n );

double explicit_update ( GRID *grid, SYS_MEM *smem, PARAMS p );

double linear_solve ( GRID *grid, SYS_MEM *smem, PARAMS p );

void lu( double *A );

void backstab ( double *A, double *x, double *y );

void forwardstab ( double *A, double *y, double *b );

void MMMult_DtG ( double *A, double *B, double *X, int stride );

void MMMult_GtG ( double *A, double *B, double *X, int stride );

void MVMult ( double *A, double *X, double *B, int stride );

void QR_House ( double *A, double *Q, int m, int n );

void QR_QtB ( double *Q, double *B, double *X, int m );

void QR_Backsub ( double *R, double *B, int n );

#endif
