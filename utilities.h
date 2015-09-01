/***************************************************************************
  Name 		  : utilities.h
  Version     : 1.0
  Author(s)   : Panayiotis Danassis (panos_dan@hotmail.com)
  Date        : May 6, 2015
  Description : Auxiliary functions and other useful stuff.
  -----------
  Copyright (C) 2015  Panayiotis Danassis
  School of ECE, National Technical University of Athens.
****************************************************************************/

#ifndef UTILITIES_H
#define UTILITIES_H

/* Constants for a random number generator, for details see numerical recipes in C */
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

/* Type Definitions */

/* Global Variables and matrices. */
extern long int seed;       /* seed for the random number generator */

/* Function Prototypes */
int compareHypernets(const void *pa, const void *pb);
int compareNets(const void *pa, const void *pb);
double rand01(long *idum);
int randomInRange(int minNumber, int maxNumber);
gridPosition_struct findPositionInGrid(int node, int pos);
int subblock0Available(int pos, int **grid, int layer);
double maximum_of_four(double d1, double d2, double d3, double d4);
node_t *list_create_node(int data);
node_t *list_insert_beginning(node_t *list, int data);
#endif