//=============================================================
// 
//  sys_mem.h
//  
//  A simple class that will contain all the data needed for
//  the memory used in solving the problem (Euler equations).
//
//  Written by - Shane Sawyer
//
//=============================================================

#include <stdio.h>
#include <stdlib.h>
#include "Defined_Variables.h"

#ifndef SYS_MEM_H
#define SYS_MEM_H

class SYS_MEM
{
  public:

  // Initialize everything to zero.
  SYS_MEM()
    {
      ia = NULL;
      iau = NULL;
      ja = NULL;
      LHS = NULL;
      RHS = NULL;
    }

  ~SYS_MEM()
    {
      freenull(ia);
      freenull(iau);
      freenull(ja);
      freenull(LHS);
      freenull(RHS);
    }

  // Data members.
  
  int *ia;                      // Starting location in ja for each node.
  int *iau;                     // Location of the diagonal entry in ja for each node.
  int *ja;                      // The compressed array of nodes surrounding node map plus the node itself.
  double *LHS;                  // The left hand side matrix for the compressed rows.
  double *RHS;                  // The right hand side of the system.

};

#endif
