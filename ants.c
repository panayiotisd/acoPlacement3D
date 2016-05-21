/***************************************************************************
  Name 		  : ants.c
  Version     : 1.0
  Author(s)   : Panayiotis Danassis (panos_dan@hotmail.com)
  Date        : May 6, 2015
  Description : Impementation of ant's procedures.
  -----------
  Copyright (C) 2015  Panayiotis Danassis
  School of ECE, National Technical University of Athens.
****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <assert.h>
#include <float.h>
#include <unistd.h>
#include <string.h>

#include "acoPlacement.h"
#include "ants.h"
#include "inOut.h"
#include "utilities.h"

/* Variables and matrices definitions. */
double rho;			/* parameter for evaporation */
double alpha;		/* importance of trail */
double beta;		/* importance of heuristic evaluate */
double q_0;			/* probability of best choice in tour construction */
double xi;			/* local pheromone update rule */
double lambda;		/* cost's function weight parameter */
double tau_0;		/* initial pheromone value */
double tau_min;		/* Minimum pheromone value */
double tau_max;		/* Maximum pheromone value */

int maxIterations;				/* maximum number of iterations to perform */
int updateIterationBestStep;	/* Every "updateIterationBestStep" number of iterations, use iteration_best ant to update pheromone matrix */
double tau_min_divisor;			/* Divisor for the initialization of tau_min */

int logic_block_array_size;

double	**pheromone;		/* Pheromone matrix */
double	**heuristic;		/* Heuristic matrix */
double	**total;			/* total = pheromone * heuristic */
#ifndef PARALLEL
	double	*probabilityOfSelection;
#endif

ant_struct ant[n_ants + 1];
ant_struct best_so_far_ant;
ant_struct iteration_best_ant;

/* Timing analysis matrices and variables */
double d_max, *t_arrival, *t_required, **slack, **delay_table;

double maxCLBHeur, minCLBHeur, maxIOHeur, minIOHeur;

/* Function Implementations */
void set_acs_parameters(double r, double a, double b, double q, double x, double l, int maxIter, int step, double divisor) {
	/*	FUNCTION:		Sets the values for the various parameters for the ACO.
		INPUT:			rho, alpha, beta, q_0, xi, lambda, maxIterations, updateIterationBestStep and tau_min_divisor
		OUTPUT:			none
		SIDE EFFECTS:	the values of the above parameters will be changed to the given ones.
						Also the seed for random number generator and the logic_block_array_size will be set.
						Finally the required memory for the probabilityOfSelection matrix will be allocated */
	
	/* Set seed for random number generator */
	seed = (long int)time(NULL);

	/* Size of Logic Block Array */
	/* Two I/O pads fit into the height of one row or the width of one column of the FPGA's
	   logic block array. As a result we need (grid_size * grid_size + 1) for the clbs and
	   (8 * grid_size) for the I/O */
	logic_block_array_size = (grid_size * grid_size) + (8 * grid_size);

	/* Allocate the required memory for the probabilityOfSelection matrix */
	#ifndef PARALLEL
		probabilityOfSelection = (double *)malloc(sizeof(double) * (logic_block_array_size + 2));
	#endif

	rho = r;
	alpha = a;
	beta = b;
	q_0 = q;
	xi = x;
	lambda = l;
	maxIterations = maxIter;
	updateIterationBestStep = step;
	tau_min_divisor = divisor;

	/* Print parameters */
	printf("Grid Size	= %d\n", grid_size);
	printf("Number of ants	= %d\n", n_ants);
	printf("\n");
	printf("netlist = %s\n", netlist_file);
	printf("hypernets = %s\n", hypernet_file);
	printf("nets = %s\n", nets_file);
	printf("layers = %s\n", layers_file);
	printf("placement = %s\n", placement_file);
	printf("rho	= %f\n", rho);
	printf("alpha	= %f\n", alpha);
	printf("beta	= %f\n", beta);
	printf("q_0	= %f\n", q_0);
	printf("xi	= %f\n", xi);
	printf("lambda	= %f\n", lambda);
	if (wire_timing_flag)
		printf("costFun	= wire + timing\n");
	else if (quadratic_estimate_flag)
		printf("costFun	= quadratic estimate\n");
	else if (hops_flag)
		printf("costFun	= number of hops\n");
	if (acs_flag)
		printf("pherUpd	= ACS\n");
	else if (mmas_flag)
		printf("pherUpd	= MMAS\n");
	printf("maxIter	= %d\n", maxIterations);
	#ifdef PARALLEL
		printf("threads	= %d\n", numberOfThreads);
	#else
		printf("threads	= 1\n");
	#endif
	printf("step	= %d\n", updateIterationBestStep);
	printf("divisor	= %f\n", tau_min_divisor);
	if (heuristic_tour_flag)
		printf("tour	= heuristic\n");
	else
		printf("tour	= random\n");
	if (customized_heuristic_flag)
		printf("customHr= YES\n");
	else
		printf("customHr= NO\n");
	printf("\n");

	return;
}

void init_ants(void) {
	/*	FUNCTION:		Initiates the variables and matrices in the ant_struct.
		INPUT:			none
		OUTPUT:			none
		SIDE EFFECTS:	nodePosition and grid matrices will be initiated for every ant including the best_so_far_ant and iteration_best_ant */
	#ifdef DEBUG
		printf("------------ <init_ants> ------------\n");
	#endif
	
	int i, j, k;
	
	/* Initiate every ant */
	for (i = 1; i <= n_ants; i++) {
		ant[i].nodePosition = (int *)malloc(sizeof(int) * (n + 1));
		for (j = 1; j <= n; j++) {
			ant[i].nodePosition[j] = 0;		/* Set all nodes as unassigned */
		}

		ant[i].grid = (int **)malloc(sizeof(int *) * numberOfLayers);
		for (k = 0; k < numberOfLayers; k++) {
			ant[i].grid[k] = (int *)malloc(sizeof(int) * (logic_block_array_size + 1));
		}

		for (k = 0; k < numberOfLayers; k++) {
			ant[i].grid[k][0] = FALSE;
			for (j = 1; j <= logic_block_array_size; j++) {
				ant[i].grid[k][j] = TRUE;			/* Set all clbs as available */
			}
		}
	}
	
	/* Initiate best_so_far ant */
	best_so_far_ant.nodePosition = (int *)malloc(sizeof(int) * (n + 1));
	for (j = 1; j <= n; j++) {
		best_so_far_ant.nodePosition[j] = 0;	/* Set all nodes as unassigned */
	}
	
	best_so_far_ant.grid = (int **)malloc(sizeof(int *) * numberOfLayers);
	for (k = 0; k < numberOfLayers; k++) {
		best_so_far_ant.grid[k] = (int *)malloc(sizeof(int) * (logic_block_array_size + 1));
	}

	for (k = 0; k < numberOfLayers; k++) {
		best_so_far_ant.grid[k][0] = FALSE;
		for (j = 1; j <= logic_block_array_size; j++) {
			best_so_far_ant.grid[k][j] = TRUE;			/* Set all clbs as available */
		}
	}
	best_so_far_ant.cost = DBL_MAX;

	/* Initiate iteration_best ant */
	iteration_best_ant.nodePosition = (int *)malloc(sizeof(int) * (n + 1));
	for (j = 1; j <= n; j++) {
		iteration_best_ant.nodePosition[j] = 0;	/* Set all nodes as unassigned */
	}
	
	iteration_best_ant.grid = (int **)malloc(sizeof(int *) * numberOfLayers);
	for (k = 0; k < numberOfLayers; k++) {
		iteration_best_ant.grid[k] = (int *)malloc(sizeof(int) * (logic_block_array_size + 1));
	}

	for (k = 0; k < numberOfLayers; k++) {
		iteration_best_ant.grid[k][0] = FALSE;
		for (j = 1; j <= logic_block_array_size; j++) {
			iteration_best_ant.grid[k][j] = TRUE;			/* Set all clbs as available */
		}
	}
	iteration_best_ant.cost = DBL_MAX;
	#ifdef DEBUG
		printf("------------ </init_ants> ------------\n");
	#endif
	return;
}

void init_timing_matrices(void) {
	/*	FUNCTION:		Initiates all the matrices required for timing analysis.
		INPUT:			none
		OUTPUT:			none
		SIDE EFFECTS:	At the end, Tarrival, Trequired, Slack and delay_table will have been initialized */
	int i, j;

	/* Allocate enough memory for the various matrices */
	t_arrival = (double *)malloc(sizeof(double) * (n + 1));
	t_required = (double *)malloc(sizeof(double) * (n + 1));
	slack = (double **)malloc(sizeof(double *) * (n + 1));
	for (i = 1; i <= n; i++)
		slack[i] = (double *)malloc(sizeof(double) * (n + 1));
	delay_table = (double **)malloc(sizeof(double *) * (n + 1));
	for (i = 1; i <= n; i++)
		delay_table[i] = (double *)malloc(sizeof(double) * (n + 1));

	/* Initiate arrays */
	for (i = 1; i <= n; i++) {
		t_arrival[i] = -1.0;
		t_required[i] = -1.0;
		for (j = 1; j <= n; j++) {
			slack[i][j] = -1.0;
			delay_table[i][j] = -1.0;
		}
	}
	return;
}

void init_pheromone_matrix(void) {
	/*	FUNCTION:		Initiates the pheromone matrix.
		INPUT:			none
		OUTPUT:			none
		SIDE EFFECTS:	pheromone matrix initiated with value tau_0.
						tau_0 stores the initial trail */
	#ifdef DEBUG
		printf("------------ <init_pheromone_matrix> ------------\n");
	#endif

	int i, j;
	double cost, fitness;

	/* Allocate enough space for the matrix */
	pheromone = (double **)malloc(sizeof(double *) * (n + 1));
	for (i = 0; i <= n; i++) {
		pheromone[i] = (double *)malloc(sizeof(double) * (logic_block_array_size + 1));
	}
	
	//cost = 1.0 / (n * heuristic_tour());	/* Initial pheromone level = 1 / (Problem Size * Best Tour Cost) */
	if (heuristic_tour_flag)
		cost = heuristic_tour();
		//cost = 3 * rho * heuristic_tour();
	else
		cost = random_tour();
		//cost = 3 * rho * random_tour();

	fitness = 1 / cost;
	
	tau_max = 1 / (rho * 0.35 * cost);
	tau_min = tau_max / tau_min_divisor;

	tau_0 = tau_max; /* Store for later use in local pheromone update routine */
	
	for (i = 1; i <= n; i++) {
		for (j = 1; j <= logic_block_array_size; j++) {
			pheromone[i][j] = tau_max;
		}
	}
	#if defined(DEBUG) || defined(PRINT)
		printf("tau_0 = %f\n", tau_0);
		printf("tau_max = %f\n", tau_max);
		printf("tau_min = %f\n", tau_min);
	#endif
	#ifdef DEBUG
		printf("------------ </init_pheromone_matrix> ------------\n");
	#endif
	return;
}

void reinit_pheromone_matrix(void) {
	/*	FUNCTION:		Reinitiates the pheromone matrix to avoid stagnation.
		INPUT:			none
		OUTPUT:			none
		SIDE EFFECTS:	pheromone matrix reinitiated with value tau_0. */
	#if defined(DEBUG) || defined(PRINT)
		printf("Reinitialize pheromone matrix(%d)\t", iteration);
	#endif
	
	int i, j;
	
	for (i = 1; i <= n; i++) {
		for (j = 1; j <= logic_block_array_size; j++) {
			pheromone[i][j] = tau_0;
		}
	}

	reinit_total_matrix();

	#if defined(DEBUG) || defined(PRINT)
		printf("[DONE]\n");
	#endif
	return;
}

void init_heuristic_matrix(void) {
	/*	FUNCTION:		Initiates the heuristic matrix.
		INPUT:			none
		OUTPUT:			none
		SIDE EFFECTS:	heuristic matrix initiated.
						If there is a file ("{benchmark}.heur") containing the heuristic matrix, then the values are loaded
						from there. Otherwise the values are computed from scratch. A .heur file is then created containing
						a copy of the heuristic matrix.
		DESCRIPTION:	h(i,j) = Cross(i) * Connect(j)
		Where Cross(i) denotes the connective extent between block i and other blocks.
		The straightforward way to evaluate Cross(i) is counting how many blocks are connected to block i.
		Connect(j) is defined for the sum distance degree between location (j) and other locations.
		Intuitively, the greater connection a block has in the circuit netlist, the more centric it should be implemented. */
	#ifdef DEBUG
		printf("------------ <init_heuristic_matrix> ------------\n");
	#endif

	int i, j, tmp;
	double cross_i, connect_j;
	char *heuristics_file, *temp;
	FILE *fp;

	maxCLBHeur = -1;
	minCLBHeur = DBL_MAX;
	maxIOHeur = -1;
	minIOHeur = DBL_MAX;

	/* Allocate enough space for the matrix */
	heuristic = (double **)malloc(sizeof(double *) * (n + 1));
	for (i = 0; i <= n; i++) {
		heuristic[i] = (double *)malloc(sizeof(double) * (logic_block_array_size + 1));
	}
	
	temp = (char *)malloc(strlen(netlist_file) + 1);
	memcpy(temp, netlist_file, strlen(netlist_file));
	heuristics_file = strtok(temp, ".");
	strcat(heuristics_file, ".heur");

	/* Check if the hueristic values have already been computed */
	if(access(heuristics_file, F_OK) != -1) {
		/* Heuristics file exists */
		fp = fopen(heuristics_file, "r"); /* Open a file for reading. The file must exist. */
		if (fp == NULL) {
			printf("fopen failed\n");
			exit(EXIT_FAILURE);
		}

		printf("Heuristic file found: %s\n", heuristics_file);

		if(EOF == fscanf(fp,"%d", &tmp)) {
			printf("ERROR: Reading heuristics file. Could not read grid_size.\n");
			exit(EXIT_FAILURE);
		}
		if (tmp == grid_size) {
			/* We have hueristic information for the correct grid size */
			printf("Reading heuristic values from: %s ", heuristics_file);

			for (i = 1; i <= n; i++) {
				for (j = 1; j <= logic_block_array_size; j++) {
					if(EOF == fscanf(fp,"%lf", &heuristic[i][j])) {
						printf("ERROR: Reading heuristics file. Not enough data\n");
						exit(EXIT_FAILURE);
					}
				}
			}

			if(EOF == fscanf(fp,"%lf", &minCLBHeur)) {
				printf("ERROR: Reading heuristics file. Not enough data\n");
				exit(EXIT_FAILURE);
			}
			if(EOF == fscanf(fp,"%lf", &maxCLBHeur)) {
				printf("ERROR: Reading heuristics file. Not enough data\n");
				exit(EXIT_FAILURE);
			}
			if(EOF == fscanf(fp,"%lf", &minIOHeur)) {
				printf("ERROR: Reading heuristics file. Not enough data\n");
				exit(EXIT_FAILURE);
			}
			if(EOF == fscanf(fp,"%lf", &maxIOHeur)) {
				printf("ERROR: Reading heuristics file. Not enough data\n");
				exit(EXIT_FAILURE);
			}

			/* Sanity check */
			fscanf(fp,"%s", temp);
			if (!feof(fp)) {
				printf("ERROR - init_heuristic_matrix: Heuristics file (%s) contains more information that it should\n", heuristics_file);
				exit(EXIT_FAILURE);
			}
			fclose(fp);
			printf("[ OK ]\n");

			#if defined(DEBUG) || defined(PRINT)
				printf("minCLBHeur	= %f\n", minCLBHeur);
				printf("maxCLBHeur	= %f\n", maxCLBHeur);
				printf("minIOHeur	= %f\n", minIOHeur);
				printf("maxIOHeur	= %f\n", maxIOHeur);
			#endif
			return;
		}
		else {
			printf("                      Wrong grid size! Current grid = %d, File's grid = %d. Recomputing the values...\n", grid_size, tmp);
		}
	}

	/* Heuristics file doesn't exist */
	/* We need to compute the values */
	printf("Computing heuristic values... ");

	for (i = 1; i <= n; i++) {
		cross_i = cross(i);
		assert(cross_i != 0);
		if (!inOuts[i]) {
			/* CLB spot */
			for (j = 1; j <= grid_size * grid_size; j++) {
				connect_j = connect(findPositionInGrid(i, j));
				heuristic[i][j] = 10000 * cross_i * connect_j;
				if (heuristic[i][j] > maxCLBHeur)
					maxCLBHeur = heuristic[i][j];
				else if (heuristic[i][j] < minCLBHeur)
					minCLBHeur = heuristic[i][j];
			}
			for (j = grid_size * grid_size + 1; j <= logic_block_array_size; j++) {
				heuristic[i][j] = 0;
			}
		}
		else {
			/* I/O spot */
			for (j = 1; j <= grid_size * grid_size; j++) {
				heuristic[i][j] = 0;
			}
			for (j = grid_size * grid_size + 1; j <= logic_block_array_size; j++) {
				connect_j = connect(findPositionInGrid(i, j));
				heuristic[i][j] = 10000 * cross_i * connect_j; //FIXME Better initial value for IO heuristic
				if (heuristic[i][j] > maxIOHeur)
					maxIOHeur = heuristic[i][j];
				else if (heuristic[i][j] < minIOHeur)
					minIOHeur = heuristic[i][j];
			}		
		}
	}

	/* Store the matrix in a file for later use */
	fp = fopen(heuristics_file, "w"); /*Opens a text file for writing, if it does not exist then a new file is created. Writing starts from the beginning of the file. */
	if (fp == NULL) {
		printf("fopen failed\n");
		exit(EXIT_FAILURE);
	}

	/* Print grid_size */
	fprintf(fp,"%d\n", grid_size);

	/* Print heuristic values */
	for (i = 1; i <= n; i++)
		for (j = 1; j <= logic_block_array_size; j++)
			fprintf(fp,"%lf\n", heuristic[i][j]);

	fprintf(fp,"%lf\n", minCLBHeur);
	fprintf(fp,"%lf\n", maxCLBHeur);
	fprintf(fp,"%lf\n", minIOHeur);
	fprintf(fp,"%lf\n", maxIOHeur);

	fclose(fp);

	printf("[ OK ]\n");

	#if defined(DEBUG) || defined(PRINT)
		printf("minCLBHeur	= %f\n", minCLBHeur);
		printf("maxCLBHeur	= %f\n", maxCLBHeur);
		printf("minIOHeur	= %f\n", minIOHeur);
		printf("maxIOHeur	= %f\n", maxIOHeur);
	#endif
	#ifdef DEBUG
		printf("------------ </init_heuristic_matrix> ------------\n");
	#endif
	return;
}

double cross(int x) {
	/*	FUNCTION:		Computes the cross (connective extent) of block x.
		INPUT:			Block x
		OUTPUT:			Value of cross
		SIDE EFFECTS:	none */
	int i;
	double cross = 0;
	
	for (i = 1; i <= n; i++) {
		cross += graph[x][i];
	}
	return cross;
}

double connect(gridPosition_struct position) {
	/*	FUNCTION:		Computes the connect of a given grid position.
		INPUT:			Grid position
		OUTPUT:			Value of connect
		SIDE EFFECTS:	none */
	int i, j, row, col;
	double connect = 0;
	
	row = position.row;
	col = position.col;
			
	for (i = 1; i <= grid_size; i++) {
		for (j = 1; j <= grid_size; j++) {
			connect += abs(row - i) + abs(col - j); 
		}
	}
	return (1 / connect);
}

void init_total_matrix(void) {
	/*	FUNCTION:		Initiates the total matrix.
		INPUT:			none
		OUTPUT:			none
		SIDE EFFECTS:	total matrix initiated with value (heuristic * pheromone). */
	#ifdef DEBUG
		printf("------------ <init_total_matrix> ------------\n");
	#endif
	
	int i, j;

	/* Allocate enough space for the matrix */
	total = (double **)malloc(sizeof(double *) * (n + 1));
	for (i = 0; i <= n; i++) {
		total[i] = (double *)malloc(sizeof(double) * (logic_block_array_size + 1));
	}
	
	for (i = 1; i <= n; i++) {
		for (j = 1; j <= logic_block_array_size; j++) {
			total[i][j] = pow(pheromone[i][j], alpha) * pow(heuristic[i][j], beta);
		}
	}
	#ifdef DEBUG
		printf("------------ </init_total_matrix> ------------\n");
	#endif
	return;
}

void reinit_total_matrix(void) {
	/*	FUNCTION:		Reinitiates the total matrix.
		INPUT:			none
		OUTPUT:			none
		SIDE EFFECTS:	total matrix reinitiated with value (heuristic * pheromone). */
	#ifdef DEBUG
		printf("------------ <reinit_total_matrix> ------------\n");
	#endif
	
	int i, j;
	
	for (i = 1; i <= n; i++) {
		for (j = 1; j <= logic_block_array_size; j++) {
			total[i][j] = pow(pheromone[i][j], alpha) * pow(heuristic[i][j], beta);
		}
	}
	#ifdef DEBUG
		printf("------------ </reinit_total_matrix> ------------\n");
	#endif
	return;
}

double heuristic_tour(void) {
	/*	FUNCTION:		Auxiliary routine for finding a "nearest neighbour" tour to initialize pheromone table.
						(for every step we choose the best possible node according only to the heuristic information) 
		INPUT:			none
		OUTPUT:			Returns the cost of the heuristic tour.
		SIDE EFFECTS:	ant[1] will be set as best_so_far_ant.
						memory of ant[1] will be erased */
	#ifdef DEBUG
		printf("------------ <heuristic_tour> ------------\n");
	#endif
	
    int step;
	double cost, currHeuristicValue;

    ant_empty_memory(&ant[1]);

    step = 1; /* counter of the construction steps */
    while (step <= n) {
		currHeuristicValue = choose_heuristicly_best(&ant[1]);
		#if defined(DEBUG) || defined(PRINT)
			if (step == 1)
				printf("choose_heuristicly[max] = %f\n", currHeuristicValue);
			else if (step == n)
				printf("choose_heuristicly[min] = %f\n", currHeuristicValue);
		#endif
		step++;
    }

	compute_placement_quality(&ant[1]);
	set_as_best_ant(&ant[1]);
	
	cost = ant[1].cost;
		
    ant_empty_memory(&ant[1]);
	
	#ifdef DEBUG
		printf("------------ </heuristic_tour> ------------\n");
	#endif
    return cost;
}

double choose_heuristicly_best(ant_struct *a) {
	/*	FUNCTION:		Chooses the next edge for the tour using only the heuristic info.
		INPUT:			ant_struct used to store the computed placement. 
		OUTPUT:			Returns the maximum heuristic value found.
		SIDE EFFECTS:	Assigns the node with the maximum heuristic value to the corresponding clb.
						As a result nodePosition matrix and grid matrix in the ant's struct are updated accordingly. */
	#ifdef DEBUG
		printf("------------ <choose_heuristicly_best> ------------\n");
	#endif
	
    int node, next_node = 0, clb, next_clb = 0;
	double max_heuristic;

    max_heuristic = -1.0;             		/* Search edge with maximum heuristic value */    
    for (node = 1 ; node <= n; node++) {
		if (a->nodePosition[node] == 0) {	/* Circuit node not yet assigned */
			if (!inOuts[node]) {
				/* CLB spot */
				for (clb = 1; clb <= grid_size * grid_size; clb++){
					if (a->grid[nodeLayer[node]][clb] == TRUE) {	/* Spot available */
						if (heuristic[node][clb] > max_heuristic) {
							next_node = node;
							next_clb = clb;
							max_heuristic = heuristic[node][clb];
						}
					}
				}
			}
			else {
				/* I/O spot */
				for (clb = grid_size * grid_size + 1; clb <= logic_block_array_size; clb++){
					if (a->grid[nodeLayer[node]][clb] == TRUE) {	/* Spot available */
						if (heuristic[node][clb] > max_heuristic) {
							next_node = node;
							next_clb = clb;
							max_heuristic = heuristic[node][clb];
						}
					}
				}
			}
		}
    }
    assert(next_node != 0 && next_clb != 0);
    a->nodePosition[next_node] = next_clb;
	a->grid[nodeLayer[next_node]][next_clb] = FALSE;
	
	#ifdef DEBUG
		printf("------------ </choose_heuristicly_best> ------------\n");
	#endif
	return max_heuristic;
}

double random_tour(void) {
	/*	FUNCTION:		Auxiliary routine for finding a random tour to initialize pheromone table.
		INPUT:			none
		OUTPUT:			Returns the cost of the random tour.
		SIDE EFFECTS:	ant[1] will be set as best_so_far_ant.
						memory of ant[1] will be erased */
	#ifdef DEBUG
		printf("------------ <random_tour> ------------\n");
	#endif
	
    int node, clb;
	double cost;

    ant_empty_memory(&ant[1]);

    for (node = 1 ; node <= n; node++) {
    	if (!inOuts[node]) {
			/* CLB spot */
			clb = randomInRange(1, grid_size * grid_size);		/* Choose a CLB spot randomly */
			while (ant[1].grid[nodeLayer[node]][clb] != TRUE) {	/* Until you find one that is available */
				clb = randomInRange(1, grid_size * grid_size);
			}
			ant[1].nodePosition[node] = clb;
			ant[1].grid[nodeLayer[node]][clb] = FALSE;
		}
		else {
			/* I/O spot */
			clb = randomInRange(grid_size * grid_size + 1, logic_block_array_size);	/* Choose an IO spot randomly */
			while (ant[1].grid[nodeLayer[node]][clb] != TRUE) {						/* Until you find one that is available */
				clb = randomInRange(grid_size * grid_size + 1, logic_block_array_size);
			}
			ant[1].nodePosition[node] = clb;
			ant[1].grid[nodeLayer[node]][clb] = FALSE;
		}
    }

	compute_placement_quality(&ant[1]);
	set_as_best_ant(&ant[1]);
	
	cost = ant[1].cost;
		
    ant_empty_memory(&ant[1]);
	
	#ifdef DEBUG
		printf("------------ </random_tour> ------------\n");
	#endif
    return cost;
}

int termination_condition(void){
	/*	FUNCTION:		Checks whether a termination condition is met.
		INPUT:			none 
		OUTPUT:			TRUE if condition is met, FALSE otherwise 
		SIDE EFFECTS:	none */
	#ifdef DEBUG
		printf("------------ <termination_condition> ------------\n");
	#endif
	if (iteration > maxIterations) {
		#ifdef DEBUG
			printf("------------ </termination_condition: TRUE maxIterations> ------------\n");
		#endif
		return TRUE;
	}
	else {
		#ifdef DEBUG
			printf("------------ </termination_condition: FALSE> ------------\n");
		#endif
		return FALSE;
	}
}

void change_parameters_according_to_schedule(void) {
	/*	FUNCTION:		Changes the Ant Colony's parameters according to a cooling schedule, to modify the exploration/expoitation ratio.
		INPUT:			none 
		OUTPUT:			none
		SIDE EFFECTS:	Parameters q_0 and beta will be changed. */
	int probability, it;
	double step1Percentage = 0.5, step2Percentage = 0.9;

	if (iteration <= step1Percentage * maxIterations) {
		/* Step 1 - For the first 70% of the iterations */
		probability = 0.9 / (step1Percentage * maxIterations);
		q_0 = (iteration - 1) * probability;
	}
	else if (iteration <= step2Percentage * maxIterations) {
		/* Step 2 - For the next 20% of the iterations */
		probability = 0.99 / ((step2Percentage - step1Percentage) * maxIterations);
		it = (iteration - 1) - (step1Percentage * maxIterations);
		q_0 = 0.9 + it * probability;
		beta = beta / 2;
	}
	else {
		/* Step 3 - For the rest 10% of the iterations */
		q_0 = 0.99;
		beta = beta / 2;
	}
	return;
}

void compute_placement_quality(ant_struct *a) {
	/*	FUNCTION:		Computes the current placement's quality.
		INPUT:			The ant_struct for which we want to compute it's placement's quality 
		OUTPUT:			none
		SIDE EFFECTS:	ant's wiringCost, timingCost, numberOfHops, quadraticEstimate and cost values will be set */
	#ifdef DEBUG
		printf("------------ <compute_placement_quality> ------------\n");
	#endif
	
	int i, j, l, node, pos, row, col, currentBBX, currentBBY, x1, x2, y1, y2, z1, z2;
	int minX = INT_MAX, maxX = 0, minY = INT_MAX, maxY = 0, minZ = INT_MAX, maxZ = 0;
	double k, q, sumBB = 0.0, timingCost = 0.0, wiringCost = 0.0, quadraticEstimate = 0.0, numberOfHops = 0.0, cost = 0.0;
	gridPosition_struct position;
	
	#ifdef DEBUG
		/* Print current node assignment */
		printf("------------ <current node assignment> ------------\n");
		printf("nodePosition:\n");
		for (i = 1; i <= n; i++) {
			printf("node = %d,	poss = %d,	layer=%d\n", i, a->nodePosition[i], nodeLayer[i]);
		}
		print_placement(a->nodePosition);
		printf("------------ <current node assignment> ------------\n");
	#endif

	/* Reinitialize timing tables */
	if (((lambda > 0) && (wire_timing_flag)) || (hops_flag))
		reinit_timing_tables();

	if (wire_timing_flag) {
		sumBB = 0.0;
		k = 0.025;
		/* Compute the sum of the Bounding Boxes for the k% of the largest nets (Wiring_Cost) */
		for (i = 1; i <= (int)ceil(k * numberOfNets); ++i) {
			minX = INT_MAX;
			minY = INT_MAX;
			minZ = INT_MAX;
			maxX = 0;
			maxY = 0;
			maxZ = 0;
			for (j = 1; j <= nets[i].size; j++) {
				node = nets[i].blocks[j];
				pos = a->nodePosition[node];
				position = findPositionInGrid(node, pos);
				row = position.row;
				col = position.col;
				/* Compute current bounding box */
				if (row < minX) {
					minX = row;
				}
				if (row > maxX) {
					maxX = row;
				}
				if (col < minY) {
					minY = col;
				}
				if (col > maxY) {
					maxY = col;
				}
				if (nodeLayer[node] < minZ) {
					minZ = nodeLayer[node];
				}
				if (nodeLayer[node] > maxZ) {
					maxZ = nodeLayer[node];
				}
			}

			/**************/
			
			if (nets[i].size <= 3)
				q = 1;
			else if (nets[i].size <= 50)
				q = 0.0358 * nets[i].size + 1;	//FIXME
			else
				q = 2.79 + 0.02616 * (nets[i].size - 50);
			
			currentBBX = (maxX - minX) + 1;
			currentBBY = (maxY - minY) + 1;
			//currentBBZ = (maxZ - minZ) + 1; /* We don't need bbz since each node goes to a prespecified layer */
			//currentBB  = currentBBX * currentBBY;
			assert(currentBBX > 0 && currentBBY > 0);
			assert(currentBBX <= grid_size + 2 && currentBBY <= grid_size + 2);
			sumBB += q * (currentBBX + currentBBY);
		}
		wiringCost = sumBB;

		/* Perform timing analysis */
		if (lambda > 0) {
			perform_timing_analysis(a->nodePosition);

			/* Compute the Timing_Cost */
			timingCost = 0.0;
			for (i = 1; i <= n; i++) {
				for (j = 1; j <= n; j++) {
					if ((blocks[j].input1Number == i) ||
						(blocks[j].input2Number == i) ||
						(blocks[j].input3Number == i) ||
						(blocks[j].input4Number == i)) {
						timingCost += delay(i, j, a->nodePosition) * pow(criticality(i, j), criticality_exponent);
					}
				}
			}
		}

		cost = lambda * timingCost + (1 - lambda) * wiringCost;
	}
	else if (quadratic_estimate_flag) {
		/* Compute the quadratic estimate for the k-point net using k(k-1)/2 2-point nets */
		quadraticEstimate = 0.0;
		k = 0.025;
		for (i = 1; i <= (int)ceil(k * numberOfNets); ++i) {
			for (j = 1; j <= nets[i].size; j++) {
				node = nets[i].blocks[j];
				pos = a->nodePosition[node];
				position = findPositionInGrid(node, pos);
				x1 = position.col;
				y1 = position.row;
				z1 = nodeLayer[node];
				for (l = j + 1; l <= nets[i].size; l++) {
					node = nets[i].blocks[l];
					pos = a->nodePosition[node];
					position = findPositionInGrid(node, pos);
					x2 = position.col;
					y2 = position.row;
					z2 = nodeLayer[node];
					quadraticEstimate += (1.0 / (nets[i].size - 1.0)) * (pow((x1 - x2), 2.0) + pow((y1 - y2), 2.0) + pow((z1 - z2), 2.0));
				}
			}
		}

		cost = quadraticEstimate;
	}
	else {
		/* Compute number of hops */
		numberOfHops = 0.0;
		for (i = 1; i <= n; i++) {
			for (j = 1; j <= n; j++) {
				if ((blocks[j].input1Number == i) ||
					(blocks[j].input2Number == i) ||
					(blocks[j].input3Number == i) ||
					(blocks[j].input4Number == i)) {
					numberOfHops += delay(i, j, a->nodePosition);
				}
			}
		}

		cost = numberOfHops;
	}

	a->wiringCost = wiringCost;
	a->timingCost = timingCost;
	a->numberOfHops = numberOfHops;
	a->quadraticEstimate = quadraticEstimate;
	a->cost = cost;
	
	#ifdef DEBUG
	/* Print placement's quality info for debugging purposes. */
		printf("------------ <placement_quality_info> ------------\n");
		printf("wiringCost			= %lf\n", a->wiringCost);
		printf("previousWiringCost	= %lf\n", best_so_far_ant.wiringCost);
		printf("timingCost			= %lf\n", a->timingCost);
		printf("previousTimingCost	= %lf\n", best_so_far_ant.timingCost);
		printf("numberOfHops		= %lf\n", a->numberOfHops);
		printf("previousNumberOfHops= %lf\n", best_so_far_ant.numberOfHops);
		printf("quadraticEstimate	= %lf\n", a->quadraticEstimate);
		printf("previousquadraticEst= %lf\n", best_so_far_ant.quadraticEstimate);
		printf("cost				= %lf\n", a->cost);
		printf("------------ </placement_quality_info> ------------\n");

		printf("------------ </compute_placement_quality> ------------\n");
	#endif
	return;
}

void compute_placement_quality_parallel(ant_struct *a) {
	/*	FUNCTION:		Computes the current placement's quality.
						Only Bounding Box cost function works in the parallel implementation because of the shared matrices in delay tables.
		INPUT:			The ant_struct for which we want to compute it's placement's quality 
		OUTPUT:			none
		SIDE EFFECTS:	ant's wiringCost, timingCost, numberOfHops, quadraticEstimate and cost values will be set */
	#ifdef DEBUG
		printf("------------ <compute_placement_quality> ------------\n");
	#endif
	
	int i, j, node, pos, row, col, currentBBX, currentBBY;
	int minX = INT_MAX, maxX = 0, minY = INT_MAX, maxY = 0, minZ = INT_MAX, maxZ = 0;
	double k, q, sumBB = 0.0, timingCost = 0.0, wiringCost = 0.0, quadraticEstimate = 0.0, numberOfHops = 0.0, cost = 0.0;
	gridPosition_struct position;

	sumBB = 0.0;
	k = 0.025;
	/* Compute the sum of the Bounding Boxes for the k% of the largest nets (Wiring_Cost) */
	for (i = 1; i <= (int)ceil(k * numberOfNets); ++i) {
		minX = INT_MAX;
		minY = INT_MAX;
		minZ = INT_MAX;
		maxX = 0;
		maxY = 0;
		maxZ = 0;
		for (j = 1; j <= nets[i].size; j++) {
			node = nets[i].blocks[j];
			pos = a->nodePosition[node];
			position = findPositionInGrid(node, pos);
			row = position.row;
			col = position.col;
			/* Compute current bounding box */
			if (row < minX) {
				minX = row;
			}
			if (row > maxX) {
				maxX = row;
			}
			if (col < minY) {
				minY = col;
			}
			if (col > maxY) {
				maxY = col;
			}
			if (nodeLayer[node] < minZ) {
				minZ = nodeLayer[node];
			}
			if (nodeLayer[node] > maxZ) {
				maxZ = nodeLayer[node];
			}
		}

		/**************/
		
		if (nets[i].size <= 3)
			q = 1;
		else if (nets[i].size <= 50)
			q = 0.0358 * nets[i].size + 1;	//FIXME
		else
			q = 2.79 + 0.02616 * (nets[i].size - 50);
		
		currentBBX = (maxX - minX) + 1;
		currentBBY = (maxY - minY) + 1;
		//currentBBZ = (maxZ - minZ) + 1; /* We don't need bbz since each node goes to a prespecified layer */
		//currentBB  = currentBBX * currentBBY;
		assert(currentBBX > 0 && currentBBY > 0);
		assert(currentBBX <= grid_size + 2 && currentBBY <= grid_size + 2);
		sumBB += q * (currentBBX + currentBBY);
	}
	wiringCost = sumBB;

	cost = wiringCost;

	a->wiringCost = wiringCost;
	a->timingCost = timingCost;
	a->numberOfHops = numberOfHops;
	a->quadraticEstimate = quadraticEstimate;
	a->cost = cost;

	#ifdef DEBUG
		printf("------------ </compute_placement_quality> ------------\n");
	#endif
	return;
}

void reinit_timing_tables(void) {
	/*	FUNCTION:		Reinitializes the timing tables
		INPUT:			None
		OUTPUT:			None
		SIDE EFFECTS:	Matrices t_arrival, t_required, slack, delay_table and variable d_max will be reinitialized to value -1.0 */
	#ifdef DEBUG
		printf("------------ <reinit_timing_tables> ------------\n");
	#endif
	int i, j;
	d_max = -1.0;
	for (i = 1; i <= n; i++) {
		t_arrival[i] = -1.0;
		t_required[i] = -1.0;
		for (j = 1; j <= n; j++) {
			slack[i][j] = -1.0;
			delay_table[i][j] = -1.0;
		}
	}
	#ifdef DEBUG
		printf("------------ </reinit_timing_tables> ------------\n");
	#endif
	return;
}

double delay(int i, int j, int *nodePosition) {
	/*	FUNCTION:		Computes the delay of the connection between the two given blocks.
		INPUT:			The two nodes for which we want to compute the delay and the
						nodePosition matrix with their positions in the grid.
		OUTPUT:			The delay between the two nodes
		SIDE EFFECTS:	None */
	double delay;
	int xi, yi, zi, xj, yj, zj;
	gridPosition_struct position;

	/* The value has allready been computed */
	if (delay_table[i][j] != -1.0)
		return delay_table[i][j];

	position = findPositionInGrid(i, nodePosition[i]);
	xi = position.col;
	yi = position.row;
	zi = nodeLayer[i];
	position = findPositionInGrid(j, nodePosition[j]);
	xj = position.col;
	yj = position.row;
	zj = nodeLayer[j];

	delay = abs(xi - xj) + abs(yi - yj) + (abs(zi - zj) * cross_layer_delay_weight);

	/* Store the value for further use */
	delay_table[i][j] = delay;
	delay_table[j][i] = delay;

	return delay;
}

void perform_timing_analysis(int *nodePosition) {
	/*	FUNCTION:		Performs timing analysis to the circuit.
		INPUT:			nodePosition matrix containing the positions of each node in the grid.
		OUTPUT:			none
		SIDE EFFECTS:	At the end, Tarrival and Trequired will have been computed for every node
						and also the Slack(i, j) for every connection (i, j). Finally we store the
						critical path delay Dmax (maximun arrival time of all sinks of the circuit). */
	int i, j;

	/* Fill t_arrival for INPUTS */
	for(i = 1; i <= n; i++)
		if (blocks[i].type == INPUT)
			t_arrival[i] = 0.0;
	for(i = 1; i <= n; i++)
		if (blocks[i].type == OUTPUT)
			compute_arrival_time(i, nodePosition);

	/* Compute d_max */
	d_max = -1.0;
	for(i = 1; i <= n; i++)
		if (blocks[i].type == OUTPUT)
			if (t_arrival[i] > d_max)
				d_max = t_arrival[i];

	/* Sanity check */
	for(i = 1; i <= n; i++)
		assert((t_arrival[i] <= d_max) && (t_arrival[i] > -1.0)); //DELETE MEEE!!

	/* Fill t_required */
	for(i = 1; i <= n; i++)
		if (blocks[i].type == OUTPUT)
			t_required[i] = d_max;
	for(i = 1; i <= n; i++)
		if (blocks[i].type == CLB)
			compute_required_time(i, nodePosition);

	/* Sanity check */
	for(i = 1; i <= n; i++)
		if (blocks[i].type != INPUT)
			assert((t_required[i] <= d_max) && (t_required[i] > -1.0)); //DELETE MEEE!!!

	/* Fill slack */
	for (i = 1; i <= n; i++) {
		for (j = 1; j <= n; j++)
			slack[i][j] = t_required[j] - t_arrival[i] - delay(i, j, nodePosition);
	}

	return;
}

double compute_arrival_time(int i, int *nodePosition) {
	/*	FUNCTION:		Computes the arival time at the node i.
		INPUT:			The node we want to compute the arrival time and nodePosition matrix containing the positions of each node in the grid.
		OUTPUT:			Arrival time
		SIDE EFFECTS:	none */
	double t1 = 0.0, t2 = 0.0, t3 = 0.0, t4 = 0.0;
	int j;

	if(t_arrival[i] >= 0.0)
		return t_arrival[i];

	j = blocks[i].input1Number;
	assert(j != -1);
	t1 = compute_arrival_time(j, nodePosition) + delay(j, i, nodePosition);
	j = blocks[i].input2Number;
	if (j != -1)
		t2 = compute_arrival_time(j, nodePosition) + delay(j, i, nodePosition);
	j = blocks[i].input3Number;
	if (j != -1)
		t3 = compute_arrival_time(j, nodePosition) + delay(j, i, nodePosition);
	j = blocks[i].input4Number;
	if (j != -1)
		t4 = compute_arrival_time(j, nodePosition) + delay(j, i, nodePosition);

	t_arrival[i] = maximum_of_four(t1, t2, t3, t4);
	return t_arrival[i];
}

double compute_required_time(int i, int *nodePosition) {
	/*	FUNCTION:		Computes the required time at the node i.
		INPUT:			The node we want to compute the required time and nodePosition matrix containing the positions of each node in the grid.
		OUTPUT:			Required time
		SIDE EFFECTS:	none */
	int j;
	double t, t_req = DBL_MAX;
	node_t *node;

	if(t_required[i] >= 0.0)
		return t_required[i];

	node = fanout[i];
	while(node) {
		j = node->data;
		t = compute_required_time(j, nodePosition) - delay(i, j, nodePosition);
		if (t_req > t)
			t_req = t; /* Keep minimum t_req forall j belonging fanout[i] */
		node = node->next;
	}

	t_required[i] = t_req;

	return t_required[i];
}

double criticality(int i, int j) {
	/*	FUNCTION:		Computes the criticality of the connection between the two given blocks.
						Criticality(i, j) = 1 - (Slack(i,j) / Dmax)
		INPUT:			The two nodes for which we want to compute the criticality and the
						nodePosition matrix with their positions in the grid.
		OUTPUT:			The criticality between the two nodes
		SIDE EFFECTS:	None */
	double criticality;

	if ((slack[i][j] <= -1.0) || (d_max <= -1.0))
		printf("i %d j %d slack %lf d_max %lf\n", i, j, slack[i][j], d_max);
	assert((slack[i][j] > -1.0) && (d_max > -1.0));
	criticality = 1 - (slack[i][j] / d_max);

	return criticality;
}

void ant_empty_memory(ant_struct *a) {
	/*	FUNCTION:		Initializes ant's memory
		INPUT:			ant to initialize 
		OUTPUT:			none
		SIDE EFFECTS:	ant's nodePosition and grid matrix will be initialized */
	int i, k;
	for (i = 1; i <= n; i++) {
		a->nodePosition[i] = 0;	/* Set all nodes as unassigned */
	}
	for (k = 0; k < numberOfLayers; k++) {
		for (i = 1; i <= logic_block_array_size; i++) {
			a->grid[k][i] = TRUE;		/* Set all clbs as available */
		}
	}
	return;
}

void construct_solution(ant_struct *a) {
	/*	FUNCTION:		Manages the solution construction phase for each ant
		INPUT:			ant for which to construct a solution
		OUTPUT:			none
		SIDE EFFECTS:	ant's nodePosition and grid matrix will be updated accordingly */
	#ifdef DEBUG
		printf("------------ <construct_solution> ------------\n");
	#endif
	int i, j, step = 1, node, prevNode;
	
	ant_empty_memory(a);

	double k = 0.025;
	if (customized_heuristic_flag) {
		/* First place the nodes for the k% of the largest nets */
		for (i = 1; i <= (int)ceil(k * numberOfNets); ++i) {
			for (j = 1; j <= nets[i].size; j++) {
				node = nets[i].blocks[j];
				if (j == 1)
					prevNode = -1;
				else
					prevNode = nets[i].blocks[j - 1];
				if (a->nodePosition[node] == 0) {	/* Node hasn't yet been assigned */
					place_specific_node(node, prevNode, a);
					step++;
				}
			}
		}
	}
	
	/* Place the rest of the nodes the old-fashioned way */
	for (step; step <= n; step++) {
		place_ant(a);
	}
	
	#ifdef DEBUG
		printf("------------ </construct_solution> ------------\n");
	#endif
	return;
}

void place_ant(ant_struct *a) {
	/*	FUNCTION:		Place an ant on a randomly chosen node of the search space
						using the pseudorandom proportional rule of ACS
						(i.e. place a netlist node onto a location on the FPGA)
		INPUT:			ant to place 
		OUTPUT:			none
		SIDE EFFECTS:	ant's nodePosition and grid matrix will be updated accordingly */
	#ifdef DEBUG
		printf("------------ <place_ant> ------------\n");
	#endif
	
	int position, node, minNumber, maxNumber, j, maxAttractiveness;
	double random, partial_sum = 0.0, sum_prob = 0.0;

	#ifdef PARALLEL
		/* In parallel implementation, every ant (thread) has its own probabilityOfSelection matrix */
		double *probabilityOfSelection;

		/* Allocate the required memory for the probabilityOfSelection matrix */
		probabilityOfSelection = (double *)malloc(sizeof(double) * (logic_block_array_size + 2));
	#endif
	
	minNumber = 1;
	maxNumber = n;
	node = randomInRange(minNumber, maxNumber); /* Creates a random integer in the range [minNumber, maxNumber] */
	while (a->nodePosition[node] != 0) {	/* Find a node that hasn't yet been assigned */
		node = randomInRange(minNumber, maxNumber);
	}
	
	#ifdef DEBUG
		printf("Node = %d\n", node);
	#endif
	
	/*	pseudorandom proportional rule */
	/*	with probability q_0 make the best possible choice according to pheromone trails and heuristic information (exploitation)
		with probability 1 - q_0 use the random proportional rule to perform biased exploration (exploration)
		we first check whether q_0 > 0.0, to avoid having to compute a random number (which is computationally expensive) in the
		very common case of q_0 = 0.0 */
	if ((q_0 > 0.0) && (rand01(&seed) <= q_0)) {
		/* exploitation */
		#ifdef DEBUG
			printf("exploitation\n");
		#endif
		position = 0;
		maxAttractiveness = -1;
		if (!inOuts[node]) {
			/* CLB */
			for (j = 1; j <= grid_size * grid_size; j++) {
				if (a->grid[nodeLayer[node]][j] == TRUE) {	/* Spot available */
					if (total[node][j] > maxAttractiveness) {
						maxAttractiveness = total[node][j];
						position = j;
					}
				}
			}
		}
		else {
			/* I/O */
			for (j = grid_size * grid_size + 1; j <= logic_block_array_size; j++) {
				if (a->grid[nodeLayer[node]][j] == TRUE) {	/* Spot available */
					if (total[node][j] > maxAttractiveness) {
						maxAttractiveness = total[node][j];
						position = j;
					}
				}
			}
		}
		#ifdef DEBUG
			printf("position = %d\n", position);
		#endif
    }
	else {
		/* exploration */
		/* Roulette wheel selection */
		#ifdef DEBUG
			printf("exploration\n");
		#endif

		partial_sum = 0.0;
		sum_prob = 0.0;
		if (!inOuts[node]) {
			/* CLB */
			for (j = 1; j <= grid_size * grid_size; j++) {
				if (a->grid[nodeLayer[node]][j] == FALSE) {
					probabilityOfSelection[j] = 0.0; /* CLB already assigned */
				}
				else {
					probabilityOfSelection[j] = total[node][j];
					sum_prob += probabilityOfSelection[j];
				}
			}
			for (j = grid_size * grid_size + 1; j <= logic_block_array_size; j++) {
				probabilityOfSelection[j] = 0.0; /* Can't place a CLB in an I/O spot */
			}
		}
		else {
			/* I/O */
			for (j = 1; j <= grid_size * grid_size; j++) {
				probabilityOfSelection[j] = 0.0; /* Can't place an I/O in a CLB spot */
			}
			for (j = grid_size * grid_size + 1; j <= logic_block_array_size; j++) {
				if (a->grid[nodeLayer[node]][j] == FALSE) {
					probabilityOfSelection[j] = 0.0; /* CLB already assigned */
				}
				else {
					probabilityOfSelection[j] = total[node][j];
					sum_prob += probabilityOfSelection[j];
				}
			}
		}
		
		random = rand01(&seed); /* random number uniformly distributed in [0, 1] */
		random *= sum_prob;
		probabilityOfSelection[logic_block_array_size + 1] = HUGE_VAL;
		j = 1;
		partial_sum = probabilityOfSelection[j];
		/* This loop always stops because probabilityOfSelection[grid_size * grid_size +1 ] == HUGE_VAL  */
		while (partial_sum <= random) {
			j++;
			partial_sum += probabilityOfSelection[j];
		}
		assert(j > 0 && j <= logic_block_array_size); /* For debuging purposes: make sure we are in range */
		if (!inOuts[node])
			assert(j >= 1 && j <= grid_size * grid_size);
		else
			assert(j >= grid_size * grid_size + 1 && j <= logic_block_array_size);
		position = j;
		#ifdef DEBUG
			printf("position = %d, layer = %d\n", position, nodeLayer[node]);
		#endif		
	}

	a->nodePosition[node] = position;
	a->grid[nodeLayer[node]][position] = FALSE;

	#ifdef PARALLEL
		free(probabilityOfSelection);
	#endif
	
	#ifdef DEBUG
		printf("------------ </place_ant> ------------\n");
	#endif
	return;
}

void place_specific_node(int node, int prevNode, ant_struct *a) {
	/*	FUNCTION:		Place a specific netlist node onto a location on the FPGA.
		INPUT:			node to place, the node that was previously placed (in order to use the distance between them as heuristic info),
						ant that constructs the solution
		OUTPUT:			none
		SIDE EFFECTS:	ant's nodePosition and grid matrix will be updated accordingly */
	#ifdef DEBUG
		printf("------------ <place_specific_node> ------------\n");
	#endif
	
	int j, position, maxAttractiveness, xi, yi, zi, xj, yj, zj;
	double random, partial_sum = 0.0, sum_prob = 0.0, distance = 0.0;
	gridPosition_struct gridPosition;

	#ifdef PARALLEL
		/* In parallel implementation, every ant (thread) has its own probabilityOfSelection matrix */
		double *probabilityOfSelection;

		/* Allocate the required memory for the probabilityOfSelection matrix */
		probabilityOfSelection = (double *)malloc(sizeof(double) * (logic_block_array_size + 2));
	#endif

	assert(a->nodePosition[node] == 0); /* Node unassigned */

	#ifdef DEBUG
		printf("Node = %d, Name = %s\n", node, blocks[node].name);
	#endif

	/* Recompute heuristic values according to the placement of the previous node of the net */
	if (prevNode != -1) {
		gridPosition = findPositionInGrid(prevNode, a->nodePosition[prevNode]);
		xi = gridPosition.col;
		yi = gridPosition.row;
		zi = nodeLayer[prevNode];

		if (!inOuts[node]) {
			/* CLB */
			for (j = 1; j <= grid_size * grid_size; j++) {	/* For every possible possition, compute a heuristic value */
				gridPosition = findPositionInGrid(node, j);
				xj = gridPosition.col;
				yj = gridPosition.row;
				zj = nodeLayer[node];

				distance = abs(xi - xj) + abs(yi - yj) + (abs(zi - zj) * cross_layer_delay_weight);

				if (distance > 0)
					heuristic[node][j] = maxCLBHeur / distance;
				else
					heuristic[node][j] = maxCLBHeur;

				total[node][j] = pow(pheromone[node][j], alpha) * pow(heuristic[node][j], beta);
			}			
		}
		else {
			/* I/O */
			for (j = grid_size * grid_size + 1; j <= logic_block_array_size; j++) {	/* For every possible possition, compute a heuristic value */
				gridPosition = findPositionInGrid(node, j);
				xj = gridPosition.col;
				yj = gridPosition.row;
				zj = nodeLayer[node];

				distance = abs(xi - xj) + abs(yi - yj) + (abs(zi - zj) * cross_layer_delay_weight);

				if (distance > 0)
					heuristic[node][j] = maxIOHeur / distance;
				else
					heuristic[node][j] = maxIOHeur;

				total[node][j] = pow(pheromone[node][j], alpha) * pow(heuristic[node][j], beta);
			}
		}
	}
	
	/*	pseudorandom proportional rule */
	/*	with probability q_0 make the best possible choice according to pheromone trails and heuristic information (exploitation)
		with probability 1 - q_0 use the random proportional rule to perform biased exploration (exploration)
		we first check whether q_0 > 0.0, to avoid having to compute a random number (which is computationally expensive) in the
		very common case of q_0 = 0.0 */
	if ((q_0 > 0.0) && (rand01(&seed) <= q_0)) {
		/* exploitation */
		#ifdef DEBUG
			printf("exploitation\n");
		#endif
		position = 0;
		maxAttractiveness = -1;
		if (!inOuts[node]) {
			/* CLB */
			for (j = 1; j <= grid_size * grid_size; j++) {
				if (a->grid[nodeLayer[node]][j] == TRUE) {	/* Spot available */
					if (total[node][j] > maxAttractiveness) {
						maxAttractiveness = total[node][j];
						position = j;
					}
				}
			}
		}
		else {
			/* I/O */
			for (j = grid_size * grid_size + 1; j <= logic_block_array_size; j++) {
				if (a->grid[nodeLayer[node]][j] == TRUE) {	/* Spot available */
					if (total[node][j] > maxAttractiveness) {
						maxAttractiveness = total[node][j];
						position = j;
					}
				}
			}
		}
		#ifdef DEBUG
			printf("position = %d\n", position);
		#endif
    }
	else {
		/* exploration */
		/* Roulette wheel selection */
		#ifdef DEBUG
			printf("exploration\n");
		#endif

		partial_sum = 0.0;
		sum_prob = 0.0;
		if (!inOuts[node]) {
			/* CLB */
			for (j = 1; j <= grid_size * grid_size; j++) {
				if (a->grid[nodeLayer[node]][j] == FALSE) {
					probabilityOfSelection[j] = 0.0; /* CLB already assigned */
				}
				else {
					probabilityOfSelection[j] = total[node][j];
					sum_prob += probabilityOfSelection[j];
				}
			}
			for (j = grid_size * grid_size + 1; j <= logic_block_array_size; j++) {
				probabilityOfSelection[j] = 0.0; /* Can't place a CLB in an I/O spot */
			}
		}
		else {
			/* I/O */
			for (j = 1; j <= grid_size * grid_size; j++) {
				probabilityOfSelection[j] = 0.0; /* Can't place an I/O in a CLB spot */
			}
			for (j = grid_size * grid_size + 1; j <= logic_block_array_size; j++) {
				if (a->grid[nodeLayer[node]][j] == FALSE) {
					probabilityOfSelection[j] = 0.0; /* CLB already assigned */
				}
				else {
					probabilityOfSelection[j] = total[node][j];
					sum_prob += probabilityOfSelection[j];
				}
			}
		}
		
		random = rand01(&seed); /* random number uniformly distributed in [0, 1] */
		random *= sum_prob;
		probabilityOfSelection[logic_block_array_size + 1] = HUGE_VAL;
		j = 1;
		partial_sum = probabilityOfSelection[j];
		/* This loop always stops because probabilityOfSelection[grid_size * grid_size +1 ] == HUGE_VAL  */
		while (partial_sum <= random) {
			j++;
			partial_sum += probabilityOfSelection[j];
		}

		assert(j > 0 && j <= logic_block_array_size); /* For debuging purposes: make sure we are in range */
		if (!inOuts[node])
			assert(j >= 1 && j <= grid_size * grid_size);
		else
			assert(j >= grid_size * grid_size + 1 && j <= logic_block_array_size);
		position = j;
		#ifdef DEBUG
			printf("position = %d, layer = %d\n", position, nodeLayer[node]);
		#endif		
	}

	a->nodePosition[node] = position;
	a->grid[nodeLayer[node]][position] = FALSE;

	#ifdef PARALLEL
		free(probabilityOfSelection);
	#endif
	
	#ifdef DEBUG
		printf("------------ </place_specific_node> ------------\n");
	#endif
	return;
}

void compare_best_so_far_ant(ant_struct *a) {
	/*	FUNCTION:		Compares current solution with best so far solution.
						If current is better, best so far becomes current.
		INPUT:			ant to compare with best_so_far_ant
		OUTPUT:			none
		SIDE EFFECTS:	best_so_far_ant might be replaced with input ant */
	#ifdef DEBUG
		printf("------------ <compare_best_so_far_ant> ------------\n");
		printf("Best so far:\n");
		printf("	wiringCost	= %lf\n", best_so_far_ant.wiringCost);
		printf("	timingCost	= %lf\n", best_so_far_ant.timingCost);
		printf("	numberOfHops= %lf\n", best_so_far_ant.numberOfHops);
		printf("	quadraticEst= %lf\n", best_so_far_ant.quadraticEstimate);
		printf("	cost		= %lf\n", best_so_far_ant.cost);
		printf("Current ant:\n");
		printf("	wiringCost	= %lf\n", a->wiringCost);
		printf("	timingCost	= %lf\n", a->timingCost);
		printf("	numberOfHops= %lf\n", a->numberOfHops);
		printf("	quadraticEst= %lf\n", a->quadraticEstimate);
		printf("	cost		= %lf\n", a->cost);
	#endif
	
	if (a->cost < best_so_far_ant.cost) {
		set_as_best_ant(a);
	}
	return;
	#ifdef DEBUG
		printf("------------ </compare_best_so_far_ant> ------------\n");
	#endif
}

void compare_iteration_best_ant(ant_struct *a) {
	/*	FUNCTION:		Compares current solution with iteration's best solution.
						If current is better, iteration best becomes current.
		INPUT:			ant to compare with iteration_best_ant
		OUTPUT:			none
		SIDE EFFECTS:	iteration_best_ant might be replaced with input ant */
	#ifdef DEBUG
		printf("------------ <compare_iteration_best_ant> ------------\n");
		printf("Best so far:\n");
		printf("	wiringCost	= %lf\n", iteration_best_ant.wiringCost);
		printf("	timingCost	= %lf\n", iteration_best_ant.timingCost);
		printf("	numberOfHops= %lf\n", iteration_best_ant.numberOfHops);
		printf("	quadraticEst= %lf\n", iteration_best_ant.quadraticEstimate);
		printf("	cost		= %lf\n", iteration_best_ant.cost);
		printf("Current ant:\n");
		printf("	wiringCost	= %lf\n", a->wiringCost);
		printf("	timingCost	= %lf\n", a->timingCost);
		printf("	numberOfHops= %lf\n", a->numberOfHops);
		printf("	quadraticEst= %lf\n", a->quadraticEstimate);
		printf("	cost		= %lf\n", a->cost);
	#endif
	
	if (a->cost < iteration_best_ant.cost) {
		set_as_iteration_best_ant(a);
	}
	return;
	#ifdef DEBUG
		printf("------------ </compare_iteration_best_ant> ------------\n");
	#endif
}

void set_as_best_ant(ant_struct *a) {
	/*	FUNCTION:		Sets given ant as best.
		INPUT:			ant to set as best_so_far_ant
		OUTPUT:			none
		SIDE EFFECTS:	best_so_far_ant is replaced with input ant */
    int   i, k;
	#if defined(DEBUG) || defined(PRINT)
		printf("Replacing "GREEN"best"RESET" so far ant with current ant...\n");
	#endif
	foundBetter = TRUE;
    for (i = 1; i <= n; i++) {
		best_so_far_ant.nodePosition[i] = a->nodePosition[i];
    }
    for (k = 0; k < numberOfLayers; k++) {
		for (i = 1; i <= logic_block_array_size; i++) {
			best_so_far_ant.grid[k][i] = a->grid[k][i];
		}
	}
	best_so_far_ant.wiringCost = a->wiringCost;
	best_so_far_ant.timingCost = a->timingCost;
	best_so_far_ant.numberOfHops = a->numberOfHops;
	best_so_far_ant.quadraticEstimate = a->quadraticEstimate;
	best_so_far_ant.cost = a->cost;
	return;
}

void set_as_iteration_best_ant(ant_struct *a) {
	/*	FUNCTION:		Sets given ant as iteration best ant.
		INPUT:			ant to set as iteration_best_ant
		OUTPUT:			none
		SIDE EFFECTS:	iteration_best_ant is replaced with input ant */
    int   i, k;
	#if defined DEBUG
		printf("Replacing iteration best ant with current ant...\n");
	#endif
    for (i = 1; i <= n; i++) {
		iteration_best_ant.nodePosition[i] = a->nodePosition[i];
    }
    for (k = 0; k < numberOfLayers; k++) {
		for (i = 1; i <= logic_block_array_size; i++) {
			iteration_best_ant.grid[k][i] = a->grid[k][i];
		}
	}
	iteration_best_ant.wiringCost = a->wiringCost;
	iteration_best_ant.timingCost = a->timingCost;
	iteration_best_ant.numberOfHops = a->numberOfHops;
	iteration_best_ant.quadraticEstimate = a->quadraticEstimate;
	iteration_best_ant.cost = a->cost;
	return;
}

void reinit_iteration_best_ant(void) {
	/*	FUNCTION:		Reinitializes iteration_best_ant for the next iteration.
		INPUT:			none
		OUTPUT:			none
		SIDE EFFECTS:	iteration_best_ant's cost will be set to DBL_MAX */
		iteration_best_ant.cost = DBL_MAX;
	return;
}

void local_pheromone_trail_update(ant_struct *a) {
	/*	FUNCTION:		Removes some pheromone on the edges of ant's solution.
		INPUT:			ant whose solution will be used to update the pheromone matrix
		OUTPUT:			none
		SIDE EFFECTS:	pheromone and total matrices will be updated */
	int i, node, clb;

    for (i = 1; i <= n; i++) {
		node = i;
		clb = a->nodePosition[node];
		
		if (acs_flag) {
			/* Use Ants Colony System's (ACS) local pheromone update rule */
			pheromone[node][clb] = (1.0 - xi) * pheromone[node][clb] + xi * tau_0;
		}
		else {
			pheromone[node][clb] = (1.0 - xi) * pheromone[node][clb];
			if (pheromone[node][clb] < tau_min)
				pheromone[node][clb] = tau_min;
		}
		total[node][clb] = pow(pheromone[node][clb], alpha) * pow(heuristic[node][clb], beta);
    }
	return;
}

void global_pheromone_trail_update(void) {
	/*	FUNCTION:		Reinforces the edges used in best ant's solution.
		INPUT:			none
		OUTPUT:			none
		SIDE EFFECTS:	pheromone and total matrices will be updated */
	#ifdef DEBUG
		printf("------------ <global_pheromone_trail_update> ------------\n");
	#endif
	
	int i, j, node, clb;
    double d_tau;

	d_tau = 1.0 / best_so_far_ant.cost;
	//d_tau = best_so_far_ant.cost;

	if (acs_flag) {
		/* Use Ants Colony System's (ACS) global pheromone update rule */
		for (i = 1; i <= n; i++) {
			node = i;
			clb = best_so_far_ant.nodePosition[node];

			pheromone[node][clb] = (1.0 - rho) * pheromone[node][clb] + rho * d_tau;
			total[node][clb] = pow(pheromone[node][clb], alpha) * pow(heuristic[node][clb], beta);
		}
	}
	else {
		/* Use MAX-MIN Ant System's (MMAS) global pheromone update rule */
		for (i = 1; i <= n; i++) {
			for (j = 1; j <= logic_block_array_size; j++) {
				pheromone[i][j] = (1 - rho) * pheromone[i][j];	/* Evaporation */

				if (iteration % updateIterationBestStep == 0) {
					if (j == iteration_best_ant.nodePosition[i])
						pheromone[i][j] = pheromone[i][j] + d_tau;
				}
				else {
					if (j == best_so_far_ant.nodePosition[i])
						pheromone[i][j] = pheromone[i][j] + d_tau;				
				}

				if (pheromone[i][j] < tau_min) {
					pheromone[i][j] = tau_min;
					#ifdef DEBUG
						//printf("tau_min reached! (%d, %d)\n", i, j);
					#endif
				}
				if (pheromone[i][j] > tau_max) {
					pheromone[i][j] = tau_max;
					#ifdef DEBUG
						//printf("tau_max reached! (%d, %d)\n", i, j);
					#endif
				}

				total[i][j] = pow(pheromone[i][j], alpha) * pow(heuristic[i][j], beta);
			}
		}
	}
	
	#ifdef DEBUG
		printf("------------ </global_pheromone_trail_update> ------------\n");
	#endif
	return;
}

void free_allocated_memory(void) {
	/*	FUNCTION:		Frees all the allocated memory.
		INPUT:			none
		OUTPUT:			none
		SIDE EFFECTS:	All the allocated memory will be returned to the OS */

	int i, k;

	#ifndef PARALLEL
		free(probabilityOfSelection);
	#endif

	/* Free ants */
	for (i = 1; i <= n_ants; i++) {
		free(ant[i].nodePosition);
		for (k = 0; k < numberOfLayers; k++)
			free(ant[i].grid[k]);
		free(ant[i].grid);
	}

	/* Free best_so_far_ant */
	free(best_so_far_ant.nodePosition);
	for (k = 0; k < numberOfLayers; k++)
		free(best_so_far_ant.grid[k]);
	free(best_so_far_ant.grid);

	/* Free iteration_best_ant */
	free(iteration_best_ant.nodePosition);
	for (k = 0; k < numberOfLayers; k++)
		free(iteration_best_ant.grid[k]);
	free(iteration_best_ant.grid);

	if (wire_timing_flag) {
		/* Free timing matrices */
			free(t_arrival);
			free(t_required);
			for (i = 1; i <= n; i++)
				free(slack[i]);
			free(slack);
			for (i = 1; i <= n; i++)
				free(delay_table[i]);
			free(delay_table);
	}

	/* Free pheromone matrix */
	for (i = 0; i <= n; i++)
		free(pheromone[i]);
	free(pheromone);

	/* Free heuristics matrix */
	for (i = 0; i <= n; i++)
		free(heuristic[i]);
	free(heuristic);

	/* Free total matrix */
	for (i = 0; i <= n; i++)
		free(total[i]);
	free(total);

	/* Free input data structs */
	for (i = 0; i <= numberOfBlocks; i++)
		free(graph[i]);
	free(graph);
	free(blocks);
	free(globalSignals);
	free(inOuts);
	free(nodeLayer);
	if (hypernets_flag)
		free(hypernets);
	if (nets_flag)
		free(nets);

	return;
}