/***************************************************************************
  Name 		  : inOut.h
  Version     : 1.0
  Author(s)   : Panayiotis Danassis (panos_dan@hotmail.com)
  Date        : May 6, 2015
  Description : Functions for the required input/output operations.
  -----------
  Copyright (C) 2015  Panayiotis Danassis
  School of ECE, National Technical University of Athens.
****************************************************************************/

#ifndef INOUT_H
#define INOUT_H

/* Color definitions for printf */
#define RESET		"\x1B[0m"
#define RED			"\x1B[31m"
#define GREEN		"\x1B[32m"
#define YELLOW		"\x1B[33m"
#define BLUE		"\x1B[34m"
#define MAGENTA		"\x1B[35m"
#define CYAN		"\x1B[36m"
#define WHITE		"\x1B[37m"
#define BOLDWHITE	"\033[1m\033[37m"

#ifndef TRUE
	#define TRUE	1
#endif
#ifndef FALSE
	#define FALSE	0
#endif

/* Block type */
#define	INPUT	1
#define	OUTPUT	2
#define CLB		3

/* Maximum value for block name */
#define	MAX_BLOCK_NAME	64

/* Initial array capacity */
#define INIT_ARRAY_CAPACITY 1000

/* Type Definitions */
typedef struct block {
	short int	type;
	char		name[MAX_BLOCK_NAME];
	char		input1[MAX_BLOCK_NAME];
	char		input2[MAX_BLOCK_NAME];
	char		input3[MAX_BLOCK_NAME];
	char		input4[MAX_BLOCK_NAME];
	char		output[MAX_BLOCK_NAME];
	char		clock[MAX_BLOCK_NAME];
	int			number;
	int			input1Number;
	int			input2Number;
	int			input3Number;
	int			input4Number;
	int			outputNumber;
} block_t;

typedef struct global_signal {
	short int	type;
	char		name[MAX_BLOCK_NAME];
	int			netNumber;
} global_signal_t;

typedef struct hypernet {
	int size;
	char name[MAX_BLOCK_NAME];
	int *blocks;
} hypernet_t;

typedef struct net {
	int size;
	char name[MAX_BLOCK_NAME];
	int *blocks;
} net_t;

typedef struct node {
	int data;
	struct node *next;
} node_t;

/* Global Variables and matrices. */
extern int				**graph;
extern block_t			*blocks;
extern global_signal_t	*globalSignals;
extern int				*inOuts;
extern int				*nodeLayer;
extern hypernet_t		*hypernets;
extern net_t			*nets;

extern node_t			**fanout;
extern node_t			**fanin;

extern int n;						/* number of nodes */
extern int numberOfBlocks;			/* number of blocks in the input netlist */
extern int numberOfGlobalSignals;	/* number of global signals in the input netlist */
extern int numberOfHypernets;		/* number of hypernets */
extern int numberOfNets;			/* number of nets */

extern int hypernets_flag, nets_flag;	/* Flags for the existance of hypernets/nets input files */

extern int wire_timing_flag;		/* Use bounding box and timing analysis as cost function */
extern int hops_flag;				/* Use number of hops as cost function - DEFAULT */
extern int quadratic_estimate_flag;	/* Use the quadratic estimate as cost function */

extern int acs_flag;			/* Use Ant Colony System's (ACS) pheromone update rules */
extern int mmas_flag;			/* Use MAX-MIN Ant System's (MMAS)  pheromone update rules - DEFAULT */

extern int customized_heuristic_flag; /* Use customized heuristic for the largests nets/paths for better results */

extern int heuristic_tour_flag;	/* Initiate pheromone matrix using a heuristic tour */
								/* If not set, a random tour is used */

extern char	*optarg;			/* getopt's variables */
extern int	optind;

extern char *netlist_file;		/* Input Files */
extern char *hypernet_file;
extern char *nets_file;
extern char *layers_file;
extern char	*placement_file;	/* Output file */

/* Function Prototypes */
void parse_netlist(void);
void parse_hypernets(void);
void parse_nets(void);
void parse_layers(void);
void parse_parameters(long int argc, char *argv[]);
void print_results(void);
void print_placement(int nodePositions[]);
void print_layers_placement(int nodePositions[], int layer);
void export_placement(int nodePositions[], int **grid);
#endif