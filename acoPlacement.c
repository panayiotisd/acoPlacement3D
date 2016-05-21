/***************************************************************************
  Name 		  : acoPlacement.c
  Version     : 1.0
  Author(s)   : Panayiotis Danassis (panos_dan@hotmail.com)
  Date        : May 6, 2015
  Description : 3-D FPGA placement algorithm based on Ant Colony Optimization.
  -----------
  Copyright (C) 2015  Panayiotis Danassis
  School of ECE, National Technical University of Athens.
****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include "acoPlacement.h"
#include "ants.h"
#include "inOut.h"
#include "timer.h"

int grid_size;
int numberOfLayers;
int numberOfThreads = 1;
int iteration = 0;			/* current iteration */
int foundBetter = FALSE;	/* A better solution has been found since last call of print_results() */

double time_used;

/* Main routine for running the ACO algorithm */
int main(int argc, char *argv[]) {
	int i;
	printf(BOLDWHITE"\nacoPlacement - A 3-D FPGA placement algorithm based on Ant Colony Optimization.\n");
	printf("Created by Panayiotis Danassis (panos_dan@hotmail.com) with the valuable contribution of Kostas Siozios (ksiop@microlab.ntua.gr)\n");
	printf("This software is licensed under the MIT License (see README.md).\n\n"RESET);

	#ifdef DEBUG
		printf("------------ <main> ------------\n");
		printf("DEBUG MODE ON\n");
	#endif
	#ifdef PARALLEL
		printf("PARALLEL MODE ON\n");
	#endif

	start_timers();
	parse_parameters(argc, argv);	/* Set ACO parameters and I/O filenames */
	parse_netlist();				/* Read input netlist */
	parse_hypernets();				/* Read hypernets */
	parse_nets();					/* Read nets */
	parse_layers();					/* Read nade's layer info */
	time_used = elapsed_time();
	printf("Reading Input Files [ OK ]\n");
	printf("Reading took %.10f seconds\n", time_used);
	
	start_timers();
	init_ants();
	printf("Initializing ants [ OK ]\n");
	init_timing_matrices();
	printf("Initializing timing matrices [ OK ]\n");
	time_used = elapsed_time();
	printf("Initialization took %.10f seconds\n", time_used);

	start_timers();
	init_heuristic_matrix();
	time_used = elapsed_time();
	printf("Initializing heuristic matrix [ OK ]\n");
	printf("Initialization took %.10f seconds\n", time_used);
	
	start_timers();
	init_pheromone_matrix();
	time_used = elapsed_time();
	printf("Initializing pheromone matrix [ OK ]\n");
	printf("Initialization took %.10f seconds\n", time_used);
	
	start_timers();
    init_total_matrix();
    time_used = elapsed_time();
    printf("Initializing total matrix [ OK ]\n");
	printf("Initialization took %.10f seconds\n", time_used);
	printf("\n");
	
	print_results();
	foundBetter = FALSE;

	#ifdef PARALLEL
		omp_set_num_threads(numberOfThreads);
	#endif

	iteration = 1;
	while (!termination_condition()) {
		printf("Iteration = %d\n", iteration);
		//change_parameters_according_to_schedule();

		#ifdef PARALLEL
			#pragma omp parallel for private(i) shared(iteration, foundBetter, best_so_far_ant)
		#endif
		for (i = 1; i <= n_ants; i++) {
			#ifdef DEBUG
				printf("Ant = %d\n", i);
			#endif
			construct_solution(&ant[i]);
			#ifdef PARALLEL
				compute_placement_quality_parallel(&ant[i]);
				compare_best_so_far_ant(&ant[i]);
			#else
				compute_placement_quality(&ant[i]);
				compare_best_so_far_ant(&ant[i]);
				compare_iteration_best_ant(&ant[i]);
				local_pheromone_trail_update(&ant[i]);
			#endif
		}


		if (foundBetter == TRUE) {
			if (iteration % printStep == 0) {
				print_results(); /* Print results periodically */
				foundBetter = FALSE;
			}
		}
		
		if (iteration % restart == 0) {
			reinit_pheromone_matrix(); /* Reinitialize pheromone matrix to avoid stagnation */
		}

		if (iteration % exportPlacementStep == 0) {
			export_placement(best_so_far_ant.nodePosition, best_so_far_ant.grid);	/* Export placement file */
		}

		global_pheromone_trail_update();
		#ifndef PARALLEL
			reinit_iteration_best_ant();
		#endif
		iteration++;
	}

	export_placement(best_so_far_ant.nodePosition, best_so_far_ant.grid);
	printf("Done!\n");

	/* Free all allocated memory */
	free_allocated_memory();
	#ifdef DEBUG
		printf("------------ </main> ------------\n");
	#endif
	return 0;
}