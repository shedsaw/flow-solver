//=============================================================
// 
//  liblapack_interface.h
//  
//  Define the functions from Lapack that we'll be calling.
//
//  Written by - Shane Sawyer
//
//=============================================================

#ifndef LAPACK_INT_H
#define LAPACK_INT_H

extern "C"
{

  void dgesvd_(char *jobu, char *jobvt, int *m, int *n, double *A, int *lda,
	       double *S, double *U, int *ldu, double *VT, int *ldvt,
	       double *Work, int *lwork, int *info);

}

#endif
