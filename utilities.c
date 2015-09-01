/***************************************************************************
  Name 		  : utilities.c
  Version     : 1.0
  Author(s)   : Panayiotis Danassis (panos_dan@hotmail.com)
  Date        : May 6, 2015
  Description : Auxiliary functions and other useful stuff.
  -----------
  Copyright (C) 2015  Panayiotis Danassis
  School of ECE, National Technical University of Athens.
****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "acoPlacement.h"
#include "ants.h"
#include "inOut.h"
#include "utilities.h"

/* Variables and matrices definitions. */
long int seed;				/* seed for the random number generator */

/* Function Implementations */
int compareHypernets(const void *pa, const void *pb) {
	/*	FUNCTION:		Auxiliary function used by qsort  to sort the hypernets in descending order.
		INPUT:			Two pointers that point to the objects being compared.
		OUTPUT:			Returns an integer less than, equal to, or greater than zero if the first argument
						is considered to be respectively less than, equal to, or greater than the second.
		SIDE EFFECTS:	none */
	const hypernet_t *a = pa;
	const hypernet_t *b = pb;
	if ( a->size > b->size )
		return -1;
	else if ( a->size == b->size )
		return 0;
	else
		return 1;
}

int compareNets(const void *pa, const void *pb) {
	/*	FUNCTION:		Auxiliary function used by qsort  to sort the nets in descending order.
		INPUT:			Two pointers that point to the objects being compared.
		OUTPUT:			Returns an integer less than, equal to, or greater than zero if the first argument
						is considered to be respectively less than, equal to, or greater than the second.
		SIDE EFFECTS:	none */
	const net_t *a = pa;
	const net_t *b = pb;
	if ( a->size > b->size )
		return -1;
	else if ( a->size == b->size )
		return 0;
	else
		return 1;
}

double rand01(long *idum) {
	/*	FUNCTION:		Generates a random number that is uniformly distributed in [0,1].
		INPUT:			Pointer to variable with the current seed
		OUTPUT:			Random number uniformly distributed in [0,1]
		SIDE EFFECTS:	Random number seed is modified (important, this has to be done!)
		ORIGIN:			numerical recipes in C */
  long k;
  double ans;

  k =(*idum)/IQ;
  *idum = IA * (*idum - k * IQ) - IR * k;
  if (*idum < 0 ) *idum += IM;
  ans = AM * (*idum);
  return ans;
}

int randomInRange(int minNumber, int maxNumber) {
	/*	FUNCTION:		Creates a random integer in the range [minNumber, maxNumber].
		INPUT:			Two integers defining the range [minNumber, maxNumber].
		OUTPUT:			Random number uniformly distributed in [minNumber, maxNumber]
		SIDE EFFECTS:	None */
	int randomInt;

	randomInt = (int)(rand01(&seed) * (maxNumber + 1)); /* Creates a random integer in the range [minNumber, maxNumber] */
	while (randomInt < minNumber || randomInt > maxNumber) {
		randomInt = (int)(rand01(&seed) * (maxNumber + 1));
	}
	return randomInt;
}

gridPosition_struct findPositionInGrid(int node, int pos) {
	/*	FUNCTION:		Finds the coordinates of the given node in the grid (row and column).
		INPUT:			Node to find its possition and his position in the nodePositions matrix.
		OUTPUT:			gridPosition struct containing the node's position
		SIDE EFFECTS:	None */
	gridPosition_struct position;
	int row, col, exactRow, exactCol, subblock = 0;

	/*	Example of a 3x3 Logic Block Array. Numbers 01-09 are CLBs and 10-33 I/O pads:
		----------------------			I/O	----------------------
		--- -- 13 14 15 -- ---				--- -- 04 05 06 -- ---
		--- -- 10 11 12 -- ---				--- -- 01 02 03 -- ---
		-33 30 01 02 03 16 19-				-24 21 -- -- -- 07 10-
		-32 29 04 05 06 17 20-				-23 20 -- -- -- 08 11-
		-31 28 07 08 09 18 21-				-22 19 -- -- -- 09 12-
		--- -- 24 23 22 -- ---				--- -- 15 14 13 -- ---
		--- -- 27 26 25 -- ---				--- -- 18 17 16 -- ---
		----------------------				----------------------
	*/

	if (!inOuts[node]) {
		/* CLB spot */
		row = pos / grid_size;
		col = pos % grid_size;
		
		if (col != 0 )
			row = row + 1;
		else
			col = grid_size;
	}
	else {
		/* I/O spot */
		pos = pos - (grid_size * grid_size);
		if (1 > pos || pos > 8 * grid_size)
			printf("hi\n");
		assert(1 <= pos && pos <= 8 * grid_size);
		if (pos <= 2 * grid_size) {
			/* Top I/O blocks */
			if (pos > 1 * grid_size) {
				subblock = 1;
				exactRow = -1;
			}
			row = 0;
			col = pos % grid_size;
			if (col == 0)
				col = grid_size;
		}
		else if (pos <= 4 * grid_size) {
			/* Right I/O blocks */
			if (pos > 3 * grid_size) {
				subblock = 1;
				exactCol = grid_size + 2;
			}
			row = pos % grid_size;
			if (row == 0)
				row = grid_size;
			col = grid_size + 1;
		}
		else if (pos <= 6 * grid_size) {
			/* Bottom I/O blocks */
			if (pos > 5 * grid_size) {
				subblock = 1;
				exactRow = grid_size + 2;
			}
			row = grid_size + 1;
			col = pos % grid_size;
			if (col == 0)
				col = grid_size;
		}
		else {
			/* Left I/O blocks */
			if (pos > 7 * grid_size) {
				subblock = 1;
				exactCol = -1;
			}
			row = pos % grid_size;
			if (row == 0)
				row = grid_size;
			col = 0;
		}
	}

	position.row = row;
	position.col = col;
	position.subblock = subblock;

	return position;
}

int subblock0Available(int pos, int **grid, int layer) {
	/*	FUNCTION:		Checks weather the subblock 0 is available in the specified IO pad.
		INPUT:			The position in the logic block array.
		OUTPUT:			TRUE if subblock 0 is available, FALSE otherwise
		SIDE EFFECTS:	None */
	int newpos;

	newpos = pos - grid_size;

	/* Sanity checks */
	pos = pos - (grid_size * grid_size);
	assert(1 <= pos && pos <= 8 * grid_size);
	if ((0 * grid_size + 1 <= pos && pos <= 1 * grid_size) ||
		(2 * grid_size + 1 <= pos && pos <= 3 * grid_size) ||
		(4 * grid_size + 1 <= pos && pos <= 5 * grid_size) ||
		(6 * grid_size + 1 <= pos && pos <= 7 * grid_size) ) {
		/* We shouldn't reach this point */
		printf("ERROR: subblock0Available was called with coordinates that doesn't correspond to some subblock 1 position\n");
		exit(EXIT_FAILURE);
	}

	if (grid[layer][newpos])
		return TRUE;
	else
		return FALSE;
}

double maximum_of_four(double d1, double d2, double d3, double d4) {
	/*	FUNCTION:		Finds the maximun between four doubles.
		INPUT:			4 doubles.
		OUTPUT:			The maximun of the 4 inputs.
		SIDE EFFECTS:	None */
	if ((d1 >= d2) && (d1 >= d3) && (d1 >= d4))
		return d1;
	if ((d2 >= d1) && (d2 >= d3) && (d2 >= d4))
		return d2;
	if ((d3 >= d1) && (d3 >= d2) && (d3 >= d4))
		return d3;
	if ((d4 >= d1) && (d4 >= d2) && (d4 >= d3))
		return d4;
	else {
		/* We should never reach this point */
		printf("maximum_of_four FAILED!\n");
		exit(EXIT_FAILURE);
	}
	return -1;
}

/* Functions for creating and manipulating a linked list */
node_t *list_create_node(int data) {
	/*	FUNCTION:		The memory needed to store the node is allocated, and the pointers are set up.
		INPUT:			Data to store in the node.
		OUTPUT:			Pointer to the node.
		SIDE EFFECTS:	None */
	node_t *node;
	if (!(node = malloc(sizeof(node_t))))
		return NULL;
	node->data = data;
	node->next = NULL;
	return node;
}

node_t *list_insert_beginning(node_t *list, int data) {
	/*	FUNCTION:		Creates and inserts the new node and returns the new head of the list.
		INPUT:			Data to store and the list to store them to.
		OUTPUT:			Pointer to the head of the list.
		SIDE EFFECTS:	None */
	node_t *newnode;
	newnode = list_create_node(data);
	newnode->next = list;
	return newnode;
}