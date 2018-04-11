//=============================================================
// 
//  mesh_io.h
//  
//  Prototype all functions needed for reading/writing a mesh
//  file in Dr. Steve Karman's format.
//
//  Written by - Shane Sawyer
//
//=============================================================

#include <stdio.h>
#include "params.h"
#include "grid.h"

#ifndef MESH_IO_H
#define MESH_IO_H

void read_mesh ( char *file_name , GRID *grid, PARAMS p );

void write_mesh ( char *file_name , GRID *grid );

#endif
