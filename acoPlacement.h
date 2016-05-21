/***************************************************************************
  Name 		  : acoPlacement.h
  Version     : 1.0
  Author(s)   : Panayiotis Danassis (panos_dan@hotmail.com)
  Date        : May 6, 2015
  Description : 3-D FPGA placement algorithm based on Ant Colony Optimization.
  -----------
  Copyright (C) 2015  Panayiotis Danassis
  School of ECE, National Technical University of Athens.
****************************************************************************/

#ifndef ACOPLACEMENT_H
#define ACOPLACEMENT_H

#ifndef TRUE
	#define TRUE	1
#endif
#ifndef FALSE
	#define FALSE	0
#endif

#define restart				INT_MAX
#define printStep			5
#define exportPlacementStep	5

/* Global Variables and matrices. */
extern int grid_size;
extern int numberOfLayers;
extern int numberOfThreads;
extern int maxIterations; /* maximum number of iterations to perform */
extern int iteration;	    /* current iteration */
extern int foundBetter;   /* A better solution haw been found since last call of print_results() */
#endif