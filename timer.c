/***************************************************************************
  Name 		  : timer.c
  Version     : 1.0
  Author(s)   : Panayiotis Danassis (panos_dan@hotmail.com)
  Date        : May 6, 2015
  Description : Functions for a timer implementation.
  -----------
  Copyright (C) 2015  Panayiotis Danassis
  School of ECE, National Technical University of Athens.
****************************************************************************/

#include <stdio.h>
#include <time.h>

#include "timer.h"

/* Variables and matrices definitions. */
static clock_t start_time;
static double elapsed;

/* Function Implementations */
void start_timers(void) {
	/*	FUNCTION:		Starts the timer.
		INPUT:			none
		OUTPUT:			none
		SIDE EFFECTS:	none */
	start_time = clock();
}

double elapsed_time(void) {
	/*	FUNCTION:		Returns the elapsed time since the timers started.
		INPUT:			none
		OUTPUT:			Elapsed time
		SIDE EFFECTS:	none */
	elapsed = clock()- start_time;
	return elapsed / CLOCKS_PER_SEC;
}