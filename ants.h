/***************************************************************************
  Name 		  : ants.h
  Version     : 1.0
  Author(s)   : Panayiotis Danassis (panos_dan@hotmail.com)
  Date        : May 6, 2015
  Description : Impementation of ant's procedures.
  -----------
  Copyright (C) 2015  Panayiotis Danassis
  School of ECE, National Technical University of Athens.
****************************************************************************/

#ifndef ANTS_H
#define ANTS_H

#ifndef TRUE
	#define TRUE	1
#endif
#ifndef FALSE
	#define FALSE	0
#endif

/* Number of ants */
#define n_ants	256

/* Vertical connections are much faster so we multiply the delay with a small factor < 1 */
#define cross_layer_delay_weight 0.1

/* The purpose of the exponent is to heavily weight connections that are critical, while giving
   less weight to connections that are non-critical. */
#define criticality_exponent 3

/* Size of Logic Block Array */
/* Two I/O pads fit into the height of one row or the width of one column of the FPGA's
   logic block array. As a result we need (grid_size * grid_size + 1) for the clbs and
   (8 * grid_size) for the I/O */
extern int logic_block_array_size;

/* Type Definitions */
typedef struct {
	int		*nodePosition;
	int		**grid;
	double	wiringCost;
	double	timingCost;
	double	numberOfHops;
	double	quadraticEstimate;
	double	cost;
} ant_struct;

typedef struct {
	int	row;
	int	col;
	int subblock;
	int	exactRow;
	int	exactCol;
} gridPosition_struct;

/* Global Variables and matrices. */
extern ant_struct ant[n_ants + 1];
extern ant_struct best_so_far_ant;
extern ant_struct iteration_best_ant;

extern double maxCLBHeur, minCLBHeur, maxIOHeur, minIOHeur;

/* Function Prototypes */
void set_acs_parameters(double r, double a, double b, double q, double x, double l, int maxIter, int step, double divisor);
void init_ants(void);
void init_timing_matrices(void);
void init_pheromone_matrix(void);
void reinit_pheromone_matrix(void);
void init_heuristic_matrix(void);
double cross(int x);
double connect(gridPosition_struct position);
void init_total_matrix(void);
void reinit_total_matrix(void);
double heuristic_tour(void);
double random_tour(void);
double choose_heuristicly_best(ant_struct *a);
int termination_condition(void);
void change_parameters_according_to_schedule(void);
void compute_placement_quality(ant_struct *a);
void compute_placement_quality_parallel(ant_struct *a);
double delay(int i, int j, int *nodePosition);
void perform_timing_analysis(int *nodePosition);
double compute_arrival_time(int i, int *nodePosition);
double compute_required_time(int i, int *nodePosition);
double criticality(int i, int j);
void reinit_timing_tables(void);
void ant_empty_memory(ant_struct *a);
void construct_solution(ant_struct *a);
void place_ant(ant_struct *a);
void place_specific_node(int node, int prevNode, ant_struct *a);
void compare_best_so_far_ant(ant_struct *a);
void compare_iteration_best_ant(ant_struct *a);
void set_as_best_ant(ant_struct *a);
void set_as_iteration_best_ant(ant_struct *a);
void reinit_iteration_best_ant(void);
void local_pheromone_trail_update(ant_struct *a);
void global_pheromone_trail_update(void);
void free_allocated_memory(void);
#endif