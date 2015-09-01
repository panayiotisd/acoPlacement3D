/***************************************************************************
  Name 		  : inOut.c
  Version     : 1.0
  Author(s)   : Panayiotis Danassis (panos_dan@hotmail.com)
  Date        : May 6, 2015
  Description : Functions for the required input/output operations.
  -----------
  Copyright (C) 2015  Panayiotis Danassis
  School of ECE, National Technical University of Athens.
****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <getopt.h>

#include "acoPlacement.h"
#include "ants.h"
#include "inOut.h"
#include "utilities.h"

/* Variables and matrices definitions. */
int				**graph;
block_t			*blocks;
global_signal_t	*globalSignals;
int				*inOuts;
int				*nodeLayer;
int				*isGlobalSignal;
hypernet_t		*hypernets;
net_t			*nets;

node_t			**fanout;
node_t			**fanin;

int n;						/* number of nodes */
int numberOfBlocks;			/* number of blocks in the input netlist */
int numberOfGlobalSignals;	/* number of global signals in the input netlist */
int numberOfHypernets;		/* number of hypernets */
int numberOfNets;			/* number of nets */

int hypernets_flag, nets_flag; 	/* Flags for the existance of hypernets/nets input files */

int wire_timing_flag = FALSE;			/* Use bounding box and timing analysis as cost function */
int hops_flag = TRUE;					/* Use number of hops as cost function - DEFAULT */
int quadratic_estimate_flag = FALSE;	/* Use the quadratic estimate as cost function */

int acs_flag = FALSE;					/* Use Ant Colony System's (ACS) pheromone update rules */
int mmas_flag = TRUE;					/* Use MAX-MIN Ant System's (MMAS)  pheromone update rules - DEFAULT */

int customized_heuristic_flag = TRUE;	/* Use customized heuristic for the largests nets/paths for better results */

int heuristic_tour_flag		 = FALSE;	/* Initiate pheromone matrix using a heuristic tour */
										/* If not set, a random tour is used */

int print_block_numbers_flag = FALSE;	/* print block's number in print_placement() */

char *netlist_file, *hypernet_file, *nets_file, *layers_file, *placement_file;	/* File names */

/* Function Implementations */
void parse_netlist(void) {
	/*	FUNCTION:		Parses the input netlist.
		INPUT:			none
		OUTPUT:			none
		SIDE EFFECTS:	The following variables and matrices will be set:
						Number of nodes: n
						Adjacency matrix: graph
						Number of blocks in the netlist: numberOfBlocks
						Matrix that stores all the blocks: blocks
						Matrix that determine if a block is an IO block: inOuts */
	#ifdef DEBUG
		printf("------------ <parse_netlist> ------------\n");
	#endif

	int i, j, cnt, counter = -1, globalSignalsCounter = 0, currentCapacity = INIT_ARRAY_CAPACITY, globalSignalsCurrentCapacity = INIT_ARRAY_CAPACITY, n_clb, n_io, found = FALSE;
	char buf[MAX_BLOCK_NAME];
	block_t *nets = (block_t*)malloc(sizeof(block_t) * INIT_ARRAY_CAPACITY);
	globalSignals = (global_signal_t*)malloc(sizeof(global_signal_t) * INIT_ARRAY_CAPACITY);
	
	FILE *fp = fopen(netlist_file, "r"); /* Open a file for reading. The file must exist. */
	if (fp == NULL) {
		printf("fopen failed\n");
		exit(EXIT_FAILURE);
	}
	
	/* Count number of blocks */
	while (fscanf(fp,"%s", buf) != EOF) {
		if (strcmp(".global", buf) == 0) {
			/* Found GLOBAL signal */
			#ifdef DEBUG
				printf("Global signal found! (%d)\n", globalSignalsCounter + 1);
			#endif
			/* Store global signal */
			globalSignalsCounter++;
			/* If array is full, reallocate more space */
			if (globalSignalsCounter >= currentCapacity) {
				#ifdef DEBUG
					printf("Array full (global signals)! Resizing...\n");
				#endif
				globalSignalsCurrentCapacity *= 2;
				globalSignals = (global_signal_t*)realloc(globalSignals, sizeof(global_signal_t) * globalSignalsCurrentCapacity);
			}
			fscanf(fp,"%s", buf);				/* Scan name */
			strcpy(globalSignals[globalSignalsCounter].name, buf);
			globalSignals[globalSignalsCounter].netNumber = -1;
		}
		else if (strcmp(".input", buf) == 0) {
			/* Found INPUT block */
			#ifdef DEBUG
				printf("Input found! (%d)\n", counter + 1);
			#endif
			/* Store block */
			counter++;
			/* If array is full, reallocate more space */
			if (counter >= currentCapacity) {
				#ifdef DEBUG
					printf("Array full! Resizing...\n");
				#endif
				currentCapacity *= 2;
				nets = (block_t*)realloc(nets, sizeof(block_t) * currentCapacity);
			}
			nets[counter].type = INPUT;
			fscanf(fp,"%s", buf);				/* Scan name */
			strcpy(nets[counter].name, buf);
			nets[counter].number = counter;
			/* Read pinlist */
			fscanf(fp,"%s", buf);
			if (strcmp("pinlist:", buf) == 0) {
				fscanf(fp,"%s", buf);				/* Scan input1 */
				strcpy(nets[counter].input1, buf);
			}
		}
		else if (strcmp(".output", buf) == 0) {
			/* Found OUTPUT block */
			#ifdef DEBUG
				printf("Output found! (%d)\n", counter + 1);
			#endif
			/* Store block */
			counter++;
			/* If array is full, reallocate more space */
			if (counter >= currentCapacity) {
				#ifdef DEBUG
					printf("Array full! Resizing...\n");
				#endif
				currentCapacity *= 2;
				nets = (block_t*)realloc(nets, sizeof(block_t) * currentCapacity);
			}
			nets[counter].type = OUTPUT;
			fscanf(fp,"%s", buf);				/* Scan name */
			strcpy(nets[counter].name, buf);
			nets[counter].number = counter;
			/* Read pinlist */
			fscanf(fp,"%s", buf);
			if (strcmp("pinlist:", buf) == 0) {
				fscanf(fp,"%s", buf);				/* Scan input1 */
				strcpy(nets[counter].input1, buf);
			}
		}
		else if (strcmp(".clb", buf) == 0) {
			/* Found CLB block */
			#ifdef DEBUG
				printf("Clb found! (%d)\n", counter + 1);
			#endif
			/* Store block */
			counter++;
			/* If array is full, reallocate more space */
			if (counter >= currentCapacity) {
				#ifdef DEBUG
					printf("Array full! Resizing...\n");
				#endif
				currentCapacity *= 2;
				nets = (block_t*)realloc(nets, sizeof(block_t) * currentCapacity);
			}
			nets[counter].type = CLB;
			fscanf(fp,"%s", buf);				/* Scan name */
			strcpy(nets[counter].name, buf);
			nets[counter].number = counter;
			/* Read pinlist */
			fscanf(fp,"%s", buf);
			if (strcmp("pinlist:", buf) == 0) {
				fscanf(fp,"%s", buf);				/* Scan input1 */
				strcpy(nets[counter].input1, buf);
				fscanf(fp,"%s", buf);				/* Scan input2 */
				strcpy(nets[counter].input2, buf);
				fscanf(fp,"%s", buf);				/* Scan input3 */
				strcpy(nets[counter].input3, buf);
				fscanf(fp,"%s", buf);				/* Scan input4 */
				strcpy(nets[counter].input4, buf);
				fscanf(fp,"%s", buf);				/* Scan output */
				strcpy(nets[counter].output, buf);
				fscanf(fp,"%s", buf);				/* Scan clock */
				strcpy(nets[counter].clock, buf);
			}
		}
	}
	/* Close file */
	fclose(fp);

	/* Number of global signals in the input netlist */
	numberOfGlobalSignals = globalSignalsCounter;

	/* Number of blocks in the input netlist */
	numberOfBlocks = counter + 1 - numberOfGlobalSignals;

	/* Reallocate Globals Signals' Table to the exact size */
	globalSignals = (global_signal_t*)realloc(globalSignals, sizeof(global_signal_t) * (numberOfGlobalSignals + 1));

	/* Fill Global Signals' Table information */
	for(i = 0; i <= counter; i++) {
		for (j = 1; j <= numberOfGlobalSignals; j++) {
			if (strcmp(nets[i].name, globalSignals[j].name) == 0) {
				globalSignals[j].type = nets[i].type;
				globalSignals[j].netNumber = i;
			}
		}
	}
	for (j = 1; j <= numberOfGlobalSignals; j++) {
		assert(globalSignals[j].netNumber != -1); /* For debugging purposes. We have found all the global signals in the netlist file */
	}

	/* Allocate enough space for the isGlobalSignal table */
	/* isGlobalSignal table is used to quickly determine if a block is a global signal or not */
	isGlobalSignal = (int *)malloc(sizeof(int) * (counter + 1));
	for(i = 0; i <= counter; i++) {
		for (j = 1; j <= numberOfGlobalSignals; j++) {
			if (strcmp(nets[i].name, globalSignals[j].name) == 0) {
				isGlobalSignal[i] = TRUE;
			}
			else {
				isGlobalSignal[i] = FALSE;
			}
		}
	}


	/* Allocate enough space for the global blocks table */
	blocks = (block_t*)malloc(sizeof(block_t) * (numberOfBlocks + 1));
	cnt = 0;
	/* First store the CLBs then the INPUTs and finally the OUTPUTs */
	for(i = 0; i <= counter; i++) {
		if (isGlobalSignal[i])
			continue; /* Skip if global signal */

		if (nets[i].type == CLB) {
			cnt++;
			blocks[cnt].type = nets[i].type;
			strcpy(blocks[cnt].name, nets[i].name);
			strcpy(blocks[cnt].input1, nets[i].input1);
			strcpy(blocks[cnt].input2, nets[i].input2);
			strcpy(blocks[cnt].input3, nets[i].input3);
			strcpy(blocks[cnt].input4, nets[i].input4);
			strcpy(blocks[cnt].output, nets[i].output);
			strcpy(blocks[cnt].clock, nets[i].clock);
			blocks[cnt].number = cnt;
		}
	}
	for(i = 0; i <= counter; i++) {
		if (isGlobalSignal[i])
			continue; /* Skip if global signal */

		if (nets[i].type == INPUT) {
			cnt++;
			blocks[cnt].type = nets[i].type;
			strcpy(blocks[cnt].name, nets[i].name);
			strcpy(blocks[cnt].input1, nets[i].input1);
			strcpy(blocks[cnt].input2, nets[i].input2);
			strcpy(blocks[cnt].input3, nets[i].input3);
			strcpy(blocks[cnt].input4, nets[i].input4);
			strcpy(blocks[cnt].output, nets[i].output);
			strcpy(blocks[cnt].clock, nets[i].clock);
			blocks[cnt].number = cnt;
		}
	}
	for(i = 0; i <= counter; i++) {
		if (isGlobalSignal[i])
			continue; /* Skip if global signal */

		if (nets[i].type == OUTPUT) {
			cnt++;
			blocks[cnt].type = nets[i].type;
			strcpy(blocks[cnt].name, nets[i].name);
			strcpy(blocks[cnt].input1, nets[i].input1);
			strcpy(blocks[cnt].input2, nets[i].input2);
			strcpy(blocks[cnt].input3, nets[i].input3);
			strcpy(blocks[cnt].input4, nets[i].input4);
			strcpy(blocks[cnt].output, nets[i].output);
			strcpy(blocks[cnt].clock, nets[i].clock);
			blocks[cnt].number = cnt;
		}
	}

	assert(cnt == numberOfBlocks);
	
	/* Allocate enough space for the inOuts table */
	/* inOuts table is used to quickly determine if a block is an IO block in order to map it on IO pins */
	inOuts = (int *)malloc(sizeof(int) * (numberOfBlocks + 1));
	for (i = 0; i <= numberOfBlocks; i++) {
		inOuts[i] = FALSE;
	}
	
	/* Fill inOuts table */
	n_clb = 0;
	n_io = 0;
	for(i = 1; i <= numberOfBlocks; i++) {
		if (blocks[i].type == INPUT) {
			inOuts[i] = TRUE;
			n_io++;
		}
		else if (blocks[i].type == OUTPUT) {
			inOuts[i] = TRUE;
			n_io++;
		}
		else {
			n_clb++;
		}
	}
	
	n = numberOfBlocks;
	
	/* Check if the Logic Block Array is large enough */
	if (numberOfLayers == 1 && ((n_clb > grid_size * grid_size) || (n_io > 8 * grid_size))) {
		printf("Error: Logic Block Array is too small!!\n");
		printf("       Grid size=%d(%dx%d), I/O spots=%d, #Nodes=%d, #I/Os=%d\n", grid_size * grid_size, grid_size, grid_size, 8 * grid_size, n_clb, n_io);
		exit(EXIT_FAILURE);
	}
	#if defined(DEBUG) || defined(PRINT)
		if (numberOfLayers == 1)
			printf("Grid size = %d(%dx%d), I/O spots = %d, #Nodes = %d, #I/Os = %d\n", grid_size * grid_size, grid_size, grid_size, 8 * grid_size, n_clb, n_io);
	#endif
	
	/* Allocate enough space for the graph's adjacency matrix */
	graph = (int**)malloc(sizeof(int *) * (numberOfBlocks + 1));
	for (i = 0; i <= numberOfBlocks; i++) {
		graph[i] = (int*)malloc(sizeof(int) * (numberOfBlocks + 1));
	}
	/* Initialize graph with zeros */
	for (i = 0; i <= numberOfBlocks; i++) {
		for (j = 0; j <= numberOfBlocks; j++) {
			graph[i][j] = 0;
		}
	}
	/* Fill the adjacency matrix */
	for(i = 1; i <= numberOfBlocks; i++) {
		found = FALSE;
		blocks[i].input1Number = -1;
		blocks[i].input2Number = -1;
		blocks[i].input3Number = -1;
		blocks[i].input4Number = -1;
		blocks[i].outputNumber = -1;
		if (blocks[i].type == CLB) {
			for (j = 1; j <= numberOfBlocks; j++) {
				if (strcmp(blocks[j].name, blocks[i].input1) == 0) {
					graph[i][j] = 1;
					graph[j][i] = 1;
					found = TRUE;
					blocks[i].input1Number = blocks[j].number;
				}
				else if (strcmp(blocks[j].name, blocks[i].input2) == 0) {
					graph[i][j] = 1;
					graph[j][i] = 1;
					found = TRUE;
					blocks[i].input2Number = blocks[j].number;
				}
				else if (strcmp(blocks[j].name, blocks[i].input3) == 0) {
					graph[i][j] = 1;
					graph[j][i] = 1;
					found = TRUE;
					blocks[i].input3Number = blocks[j].number;
				}
				else if (strcmp(blocks[j].name, blocks[i].input4) == 0) {
					graph[i][j] = 1;
					graph[j][i] = 1;
					found = TRUE;
					blocks[i].input4Number = blocks[j].number;
				}
				else if (strcmp(blocks[j].name, blocks[i].output) == 0) {
					graph[i][j] = 1;
					graph[j][i] = 1;
					found = TRUE;
					blocks[i].outputNumber = blocks[j].number;
				}
				/*else {
					// Add "out:" to output string 
					strcpy(buf, "out:");
					strcat(buf, blocks[i].output);
					if (strcmp(blocks[j].name, buf) == 0) {
						graph[i][j] = 1;
						graph[j][i] = 1;
						found = TRUE;
						blocks[i].outputNumber = blocks[j].number;
					}
				}*/
			}
			if (!found) {
				printf("Unable to find block's [%d] (%s) pins\n", i, blocks[i].name);
				exit(EXIT_FAILURE);
			}
		}
		if (blocks[i].type == OUTPUT) {
			for (j = 1; j <= numberOfBlocks; j++) {
				if ((strcmp(blocks[j].name, blocks[i].input1) == 0)) {
					graph[i][j] = 1;
					graph[j][i] = 1;
					found = TRUE;
					blocks[i].input1Number = blocks[j].number;
				}
			}
			if (!found) {
				printf("Unable to find block's [%d] (%s) pins\n", i, blocks[i].name);
				exit(EXIT_FAILURE);
			}
		}
	}
	for(i = 1; i <= numberOfBlocks; i++) {
		graph[i][i] = 0;
	}


	/* Allocate enough space for the fanin matrix */
	fanin = (node_t**)malloc(sizeof(node_t*) * (numberOfBlocks + 1));	
	/* Store the fanin of each node */
	for(i = 1; i <= numberOfBlocks; i++) {
		fanin[i] = NULL;
		if (blocks[i].input1Number != -1)
			fanin[i] = list_insert_beginning(fanin[i], blocks[i].input1Number);
		if (blocks[i].input2Number != -1)
			fanin[i] = list_insert_beginning(fanin[i], blocks[i].input2Number);
		if (blocks[i].input3Number != -1)
			fanin[i] = list_insert_beginning(fanin[i], blocks[i].input3Number);
		if (blocks[i].input4Number != -1)
			fanin[i] = list_insert_beginning(fanin[i], blocks[i].input4Number);
	}

	/* Allocate enough space for the fanout matrix */
	fanout = (node_t**)malloc(sizeof(node_t*) * (numberOfBlocks + 1));	
	/* Store the fanout of each node */
	for(i = 1; i <= numberOfBlocks; i++) {
		fanout[i] = NULL;
		for (j = 1; j <= numberOfBlocks; j++) {
			if ((blocks[j].input1Number == i) ||
				(blocks[j].input2Number == i) ||
				(blocks[j].input3Number == i) ||
				(blocks[j].input4Number == i)){
				fanout[i] = list_insert_beginning(fanout[i], j);
			}
		}
	}
	
	
	#ifdef DEBUG
		/* print netlist */
		printf("------------ <print_nets_table> ------------\n");
		for(i = 0; i <= counter; i++) {
			if (nets[i].type == INPUT) {
				printf("type=INPUT, name=%s, netnumber=%d\n\n", nets[i].name, nets[i].number);
			}
			if (nets[i].type == OUTPUT) {
				printf("type=OUTPUT, name=%s, netnumber=%d\n\n", nets[i].name, nets[i].number);
			}
			if (nets[i].type == CLB) {
				printf("type=CLB, name=%s, netnumber=%d\n", nets[i].name, nets[i].number);
				printf("pinlist: %s %s %s %s %s %s\n\n", nets[i].input1, nets[i].input2, nets[i].input3, nets[i].input4, nets[i].output, nets[i].clock);
			}
		}
		printf("------------ </print_nets_table> ------------\n");
		/* print blocks */
		printf("------------ <print_blocks_table> ------------\n");
		for(i = 1; i <= numberOfBlocks; i++) {
			if (blocks[i].type == INPUT) {
				printf("type=INPUT, name=%s, number=%d\n\n", blocks[i].name, blocks[i].number);
			}
			if (blocks[i].type == OUTPUT) {
				printf("type=OUTPUT, name=%s, number=%d\n\n", blocks[i].name, blocks[i].number);
			}
			if (blocks[i].type == CLB) {
				printf("type=CLB, name=%s, number=%d\n", blocks[i].name, blocks[i].number);
				printf("pinlist: %s %s %s %s %s %s\n\n", blocks[i].input1, blocks[i].input2, blocks[i].input3, blocks[i].input4, blocks[i].output, blocks[i].clock);
			}
		}
		printf("------------ </print_blocks_table> ------------\n");
		/* print inOut table */
		printf("------------ <print_inOuts_table> ------------\n");
		for(i = 1; i <= n; i++) {
			if (inOuts[i]) {
				printf("inOuts[%d]=TRUE\n", i);
			}
			else {
				printf("inOuts[%d]=FALSE\n", i);
			}
		}
		printf("------------ </print_inOuts_table> ------------\n");
		/* print graph table 
		printf("------------ <print_graph> ------------\n");
		for(i = 1; i <= n; i++) {
			for(j = 1; j <= n; j++) {
				printf("graph[%d][%d]=%d\n", i, j, graph[i][j]);
			}
		}
		printf("------------ </print_graph> ------------\n");*/
		/* print inOut table */
		printf("------------ <print_globalSignals_table> ------------\n");
		for (i = 1; i <= numberOfGlobalSignals; i++) {
			if (globalSignals[i].type == INPUT) {
				printf("global: %s, type=INPUT, netNumber=%d\n\n", globalSignals[i].name, globalSignals[i].netNumber);
			}
			if (globalSignals[i].type == OUTPUT) {
				printf("global: %s, type=OUTPUT, netNumber=%d\n\n", globalSignals[i].name, globalSignals[i].netNumber);
			}
			if (globalSignals[i].type == CLB) {
				printf("global: %s, type=CLB, netNumber=%d\n\n", globalSignals[i].name, globalSignals[i].netNumber);
			}
		}
		printf("------------ </print_globalSignals_table> ------------\n");
		/* print fanin table */
		printf("------------ <print_fanin_table> ------------\n");
		node_t *node;
		for (i = 1; i <= numberOfBlocks; i++) {
			node = fanin[i];
			printf("fanin[%d]=", i);
			while(node) {
				printf("%d ", node->data);
				node = node->next;
			}
			printf("\n");
		}
		printf("------------ </print_fanin_table> ------------\n");
		/* print fanout table */
		printf("------------ <print_fanout_table> ------------\n");
		for (i = 1; i <= numberOfBlocks; i++) {
			node = fanout[i];
			printf("fanout[%d]=", i);
			while(node) {
				printf("%d ", node->data);
				node = node->next;
			}
			printf("\n");
		}
		printf("------------ </print_fanout_table> ------------\n");
	#endif

	#if defined(DEBUG)
		/* print netlist with the same format as the input file */
		/* Use diff to see if you parsed the file correctly */
		
		fp = fopen("parse_netlist.out", "w"); /*Opens a text file for writing, if it does not exist then a new file is created. Writing starts from the beginning of the file. */
		char sb0[MAX_BLOCK_NAME], sb1[MAX_BLOCK_NAME], sb2[MAX_BLOCK_NAME], sb3[MAX_BLOCK_NAME], sb4[MAX_BLOCK_NAME], sb5[MAX_BLOCK_NAME];
		for(i = 0; i <= counter; i++) {
			if (nets[i].type == INPUT) {
				fprintf(fp, ".input %s\n", nets[i].name);
				fprintf(fp, "pinlist: %s \n", nets[i].input1);
				fprintf(fp, "\n");
			}
			if (nets[i].type == OUTPUT) {
				fprintf(fp, ".output %s\n", nets[i].name);
				fprintf(fp, "pinlist: %s \n", nets[i].input1);
				fprintf(fp, "\n");
			}
			if (nets[i].type == CLB) {
				fprintf(fp, ".clb %s\n", nets[i].name);
				fprintf(fp, "pinlist: %s %s %s %s %s %s \n", nets[i].input1, nets[i].input2, nets[i].input3, nets[i].input4, nets[i].output, nets[i].clock);
				strcpy(sb0, "0");
				strcpy(sb1, "1");
				strcpy(sb2, "2");
				strcpy(sb3, "3");
				strcpy(sb4, "4");
				strcpy(sb5, "5");
				if (strcmp(nets[i].input1, "open") == 0) {
					strcpy(sb0, "open");
				}
				if (strcmp(nets[i].input2, "open") == 0) {
					strcpy(sb1, "open");
				}
				if (strcmp(nets[i].input3, "open") == 0) {
					strcpy(sb2, "open");
				}
				if (strcmp(nets[i].input4, "open") == 0) {
					strcpy(sb3, "open");
				}
				if (strcmp(nets[i].output, "open") == 0) {
					strcpy(sb4, "open");
				}
				if (strcmp(nets[i].clock, "open") == 0) {
					strcpy(sb5, "open");
				}
				fprintf(fp, "subblock: %s %s %s %s %s %s %s \n", nets[i].name, sb0, sb1, sb2, sb3, sb4, sb5);
				fprintf(fp, "\n");
			}
		}
		fclose(fp);
	#endif
	
	#ifdef DEBUG
		printf("------------ </parse_netlist> ------------\n");
	#endif
	return;
}

void parse_hypernets(void) {
	/*	FUNCTION:		Parses the input hypernets.
		INPUT:			none
		OUTPUT:			none
		SIDE EFFECTS:	The following variables and matrices will be set:
						Hypernets matrix: hypernets
						Number of hypernets: numberOfHypernets */
	#ifdef DEBUG
		printf("------------ <parse_hypernets> ------------\n");
	#endif

	int i, j, size, block, counter, currentCapacity = INIT_ARRAY_CAPACITY, found = FALSE;
	char buf[MAX_BLOCK_NAME * MAX_BLOCK_NAME], tmpbuf[MAX_BLOCK_NAME * MAX_BLOCK_NAME], *string;
	
	if (!hypernets_flag)
		return;

	hypernets = (hypernet_t*)malloc(sizeof(hypernet_t) * (INIT_ARRAY_CAPACITY + 1));

	FILE *fp = fopen(hypernet_file, "r"); /* Open a file for reading. The file must exist. */
	if (fp == NULL) {
		printf("fopen failed\n");
		exit(EXIT_FAILURE);
	}

	while (fscanf(fp,"%s", buf) != EOF) {
		if (strcmp("DEBUG:", buf) == 0) {
			/* Found first hypernet */
			break;
		}
	}

	hypernets[0].size = size;
	strcpy(hypernets[0].name, "");
	hypernets[0].blocks = NULL;
	/* Read first hypernet */
	fscanf(fp,"%d", &size);
	hypernets[1].size = size;
	hypernets[1].blocks = (int *)malloc(sizeof(int) * (size + 1)); /* Allocate enough space */
	fscanf(fp,"%s %s %s %s %s", buf, buf, buf, buf, buf);	/* Go through useless info */
	strtok(buf, ":");										/* Remove colon fron hypernet's node name */
	if (strcmp("out", buf) == 0) {							/* Output block */
		string = strtok(NULL, ":");
		strcpy(tmpbuf, string);
		strcat(buf, ":");
		strcat(buf, tmpbuf);
	}
	strcpy(hypernets[1].name, buf);
	for (i = 1; i <= size; i++) {
		fscanf(fp,"%s", buf);				/* Read net */
		strtok(buf, "(");					/* Remove useless info */
		string = strtok(NULL, "(");
		strtok(string, ")");
		/* Find the number of the net we just read */
		found = FALSE;
		for(j = 1; j <= numberOfBlocks; j++) {
			if (strcmp(blocks[j].name, string) == 0) {
				block = blocks[j].number;
				found = TRUE;
				break;
			}
		}
		if (!found) {
			printf("Unable to find hypernet's [1] block: %s\n", buf);
			exit(EXIT_FAILURE);
		}
		hypernets[1].blocks[i] = block;
	}

	counter = 1;
	while (fscanf(fp,"%s", buf) != EOF) {
		if (strcmp("DEBUG:", buf) == 0) {
			/* Found hypernet */
			counter++;
			/* If array is full, reallocate more space */
			if (counter >= currentCapacity) {
				#ifdef DEBUG
					printf("Array full! Resizing...\n");
				#endif
				currentCapacity *= 2;
				hypernets = (hypernet_t*)realloc(hypernets, sizeof(hypernet_t) * (currentCapacity + 1));
			}
			/* Read hypernet's info */
			fscanf(fp,"%d", &size);
			hypernets[counter].size = size;
			hypernets[counter].blocks = (int *)malloc(sizeof(int) * (size + 1)); /* Allocate enough space */
			fscanf(fp,"%s %s %s %s %s", buf, buf, buf, buf, buf);	/* Go through useless info */
			strtok(buf, ":");										/* Remove colon fron hypernet's node name */
			if (strcmp("out", buf) == 0) {							/* Output block */
				string = strtok(NULL, ":");
				strcpy(tmpbuf, string);
				strcat(buf, ":");
				strcat(buf, tmpbuf);
			}
			strcpy(hypernets[counter].name, buf);
			for (i = 1; i <= size; i++) {
				fscanf(fp,"%s", buf);				/* Read net */
				strtok(buf, "(");					/* Remove useless info */
				string = strtok(NULL, "(");
				strtok(string, ")");
				/* Find the number of the net we just read */
				found = FALSE;
				for(j = 1; j <= numberOfBlocks; j++) {
					if (strcmp(blocks[j].name, string) == 0) {
						block = blocks[j].number;
						found = TRUE;
						break;
					}
				}
				if (!found) {
					printf("Unable to find hypernet's [%d] block: %s\n", counter, buf);
					exit(EXIT_FAILURE);
				}
				hypernets[counter].blocks[i] = block;
			}
		}
		if (strcmp("INFO:", buf) == 0) {
			/* Finished reading the hypernets */
			break;
		}
	}

	numberOfHypernets = counter;
	/* Reallocate to the exact size */
	hypernets = (hypernet_t*)realloc(hypernets, sizeof(hypernet_t) * (numberOfHypernets + 1));

	/* Sort by size in descending order */
	qsort(hypernets + 1, numberOfHypernets, sizeof(hypernet_t), compareHypernets);
	
	#ifdef DEBUG
		/* print hypernets */
		for (i = 1; i <= numberOfHypernets; ++i) {
			printf("Hypernet(%d) size=%d %s: ", i, hypernets[i].size, hypernets[i].name);
			for(j = 1; j <= hypernets[i].size; j++) {
				printf("%d ", hypernets[i].blocks[j]);
			}
			printf("\n");
		}
		printf("------------ </parse_hypernets> ------------\n");
	#endif

	fclose(fp);
	return;
}

void parse_nets(void) {
	/*	FUNCTION:		Parses the input nets.
		INPUT:			none
		OUTPUT:			none
		SIDE EFFECTS:	The following variables and matrices will be set:
						Nets matrix: nets
						Number of nets: numberOfNets */
	#ifdef DEBUG
		printf("------------ <parse_nets> ------------\n");
	#endif

	int i, j, k, l, size, block, type, counter, currentCapacity = INIT_ARRAY_CAPACITY, found = FALSE, maxBlockNumber = -1;
	char buf[MAX_BLOCK_NAME * MAX_BLOCK_NAME];
	block_t *blockList;
	
	nets = (net_t*)malloc(sizeof(net_t) * (INIT_ARRAY_CAPACITY + 1));

	FILE *fp = fopen(nets_file, "r"); /* Open a file for reading. The file must exist. */
	if (fp == NULL) {
		printf("fopen failed\n");
		exit(EXIT_FAILURE);
	}

	while (fscanf(fp,"%s", buf) != EOF) {
		if ((strcmp("Net", buf) == 0)		&&
			(fscanf(fp,"%s", buf) != EOF)	&&
			(strcmp("Name", buf) == 0)		&&
			(fscanf(fp,"%s", buf) != EOF)	&&
			(strcmp("#Pins", buf) == 0)		&&
			(fscanf(fp,"%s", buf) != EOF)	&&
			(strcmp("Driver", buf) == 0)	&&
			(fscanf(fp,"%s", buf) != EOF)	&&
			(strcmp("Recvs.", buf) == 0)	&&
			(fscanf(fp,"%s", buf) != EOF)	&&
			(strcmp("(block,", buf) == 0)	&&
			(fscanf(fp,"%s", buf) != EOF)	&&
			(strcmp("pin)", buf) == 0)) {
			/* Found the start of nets listing */
			break;
		}
	}

	nets[0].size = -1;
	strcpy(nets[0].name, "");
	nets[0].blocks = NULL;
	/* Read nets */
	counter = 0;
	while (fscanf(fp,"%s", buf) != EOF && (strcmp("Block", buf) != 0)) {
		counter++;
		/* If array is full, reallocate more space */
		if (counter >= currentCapacity) {
			#ifdef DEBUG
				printf("Array full! Resizing...\n");
			#endif
			currentCapacity *= 2;
			nets = (net_t*)realloc(nets, sizeof(net_t) * (currentCapacity + 1));
		}
		/* Read net's info */
		fscanf(fp,"%d", &size); 	/* Read and discard net's id */
		fscanf(fp,"%s", buf);
		/* Check if its a global signal */
		found = FALSE;
		for (i = 1; i <= numberOfGlobalSignals; i++) {
			if (strcmp(globalSignals[i].name, buf) == 0) {
				found = TRUE;
			}
		}
		/* Read net's info */
		strcpy(nets[counter].name, buf);
		fscanf(fp,"%d", &size);
		nets[counter].size = size;
		nets[counter].blocks = (int *)malloc(sizeof(int) * (size + 1)); /* Allocate enough space */
		for (i = 1; i <= size; i++) {
			fscanf(fp,"%s", buf);		/* Remove parenthesis */
			if (strcmp("(", buf) != 0) {
				strtok(buf, ",");
				block = atoi(buf + 1);
			}
			else {
				fscanf(fp,"%d", &block);	/* Read net */
				fscanf(fp,"%s", buf);		/* Discard comma */
			}
			nets[counter].blocks[i] = block;
			if (block > maxBlockNumber)
				maxBlockNumber = block;		/* Keep a reference to how many blocks there are */
			fscanf(fp,"%d", &block);		/* Discard useless info */
			fscanf(fp,"%s", buf);
		}

		/* Disregard it, if its a net of a global signal */
		if (found) {
			free(nets[counter].blocks);
			counter--;
		}
		found = FALSE;
	}

	numberOfNets = counter;
	/* Reallocate to the exact size */
	nets = (net_t*)realloc(nets, sizeof(net_t) * (numberOfNets + 1));

	fscanf(fp,"%s", buf); /* Discard useless info */
	while (strcmp("Connections", buf) != 0) {
		fscanf(fp,"%s", buf); /* Discard useless info */
	}
	
	/* Read block list */
	blockList = (block_t*)malloc(sizeof(block_t) * maxBlockNumber + 1);
	for (i = 0; i <= maxBlockNumber; ++i) {
		fscanf(fp,"%d", &block);
		blockList[i].number = block;
		assert(i == block);
		fscanf(fp,"%s", buf); /* Scan name */
		strcpy(blockList[i].name, buf);
		fscanf(fp,"%d", &type);
		if (type == 2) {
			blockList[i].type = INPUT;
		}
		else if (type == 1) {
			blockList[i].type = OUTPUT;
		}
		else if (type == 0) {
			blockList[i].type = CLB;
		}
		fscanf(fp,"%s", buf); /* Discard useless info */
		if (type == 0) {
			fscanf(fp,"%s %s %s %s %s", buf, buf, buf, buf, buf);
		}
	}

	/* For every net, change the net's block numbers from the ones in the net file to ours */
	for (i = 1; i <= numberOfNets; ++i) {
		for(j = 1; j <= nets[i].size; j++) {
			for (k = 0; k <= maxBlockNumber; ++k) {
				if (nets[i].blocks[j] == blockList[k].number) {
					for(l = 1; l <= numberOfBlocks; l++) {
						if (strcmp(blockList[k].name, blocks[l].name) == 0) {
							assert(blockList[k].type == blocks[l].type);
							nets[i].blocks[j] = blocks[l].number;
							break;
						}
					}
					break;
				}
			}
		}
	}

	/* Sort by size in descending order */
	qsort(nets + 1, numberOfNets, sizeof(net_t), compareNets);
	
	#ifdef DEBUG
		/* print nets */
		for (i = 1; i <= numberOfNets; ++i) {
			printf("Net(%d) size=%d %s: ", i, nets[i].size, nets[i].name);
			for(j = 1; j <= nets[i].size; j++) {
				printf("%d ", nets[i].blocks[j]);
			}
			printf("\n");
		}
		printf("------------ </parse_nets> ------------\n");
	#endif

	fclose(fp);
	return;
}

void parse_layers(void) {
	/*	FUNCTION:		Parses the file with the information used to determine in which layer we should place every node.
		INPUT:			none
		OUTPUT:			none
		SIDE EFFECTS:	The following variables and matrices will be set:
						Node's layer: nodeLayer */

	#ifdef DEBUG
		printf("------------ <parse_layers> ------------\n");
	#endif

	int i, k, layer, countCLBs, countIOs, maxCLBs = 0, maxIOs = 0, *clbPerLayer, *ioPerLayer;
	char buf[MAX_BLOCK_NAME];

	/* Allocate enough space for the nodeLayer table */
	/* nodeLayer table is used to find the layer that each node should be placed */
	nodeLayer = (int *)malloc(sizeof(int) * (numberOfBlocks + 1));
	for (i = 0; i <= numberOfBlocks; i++) {
		nodeLayer[i] = -1;
	}

	if (numberOfLayers == 1) {
		/* We don't need to read from file. Only one layer available. */
		for (i = 0; i <= numberOfBlocks; i++) {
			nodeLayer[i] = 0;
		}
		return;
	}

	FILE *fp = fopen(layers_file, "r"); /* Open a file for reading. The file must exist. */
	if (fp == NULL) {
		printf("fopen failed\n");
		exit(EXIT_FAILURE);
	}

	/* Read each node's layer */
	while (fscanf(fp, "%s", buf) != EOF) {
		for (i = 1; i <= numberOfBlocks; i++) {
			if (strcmp(blocks[i].name, buf) == 0) {
				fscanf(fp, "%s", buf); /* disregard x coordinate */
				fscanf(fp, "%s", buf); /* disregard y coordinate */
				fscanf(fp, "%s", buf); /* layer = z coordinate */
				layer = atoi(buf);
				#ifdef DEBUG
					printf("Found block %s,	layer = %d\n", blocks[i].name, layer);
				#endif
				assert(blocks[i].number == i);
				nodeLayer[i] = layer;
			}
		}
	}

	/* Sanity check */
	for (i = 1; i <= numberOfBlocks; i++) {
		assert(nodeLayer[i] != -1);
	}

	/* Check if the Logic Block Array is large enough */
	clbPerLayer = (int *)malloc(sizeof(int) * numberOfLayers + 1);
	ioPerLayer = (int *)malloc(sizeof(int) * numberOfLayers + 1);
	for (k = 0; k < numberOfLayers; k++) {
		countCLBs = 0;
		countIOs = 0;
		for (i = 1; i <= numberOfBlocks; i++) {
			if (nodeLayer[i] == k) {
				if (inOuts[i])
					countIOs++;
				else
					countCLBs++;
			}
		}
		clbPerLayer[k] = countCLBs;
		ioPerLayer[k] = countIOs;
		if (countCLBs > maxCLBs)
			maxCLBs = countCLBs;
		if (countIOs > maxIOs)
			maxIOs = countIOs;
	}
	
	if ((maxCLBs > grid_size * grid_size) || (maxIOs > 8 * grid_size)) {
		printf("Error: Logic Block Array is too small!!\n");
		printf("       Grid size=%d(%dx%d), I/O spots=%d, Layers=%d, #maxNodes/Layer=%d, #maxIO/Layer=%d\n", grid_size * grid_size, grid_size, grid_size, 8 * grid_size, numberOfLayers, maxCLBs, maxIOs);
		for (k = 0; k < numberOfLayers; k++)
			printf("       Layer %d: #Nodes=%d, #I/Os=%d\n", k, clbPerLayer[k], ioPerLayer[k]);
		exit(EXIT_FAILURE);
	}
	#if defined(DEBUG) || defined(PRINT)
		printf("Grid size = %d(%dx%d), I/O spots = %d, Layers=%d, #Nodes = %d, #I/Os = %d\n", grid_size * grid_size, grid_size, grid_size, 8 * grid_size, numberOfLayers, maxCLBs, maxIOs);
		for (k = 0; k < numberOfLayers; k++)
			printf("       Layer %d: #Nodes=%d, #I/Os=%d\n", k, clbPerLayer[k], ioPerLayer[k]);
	#endif

	#ifdef DEBUG
		printf("------------ </parse_layers> ------------\n");
	#endif
	return;
}

void parse_parameters(long int argc, char *argv[]) {
	/*	FUNCTION:		Parses the input parameters.
		INPUT:			argc (argument count) and argv (argument vector)
		OUTPUT:			none
		SIDE EFFECTS:	At the end, the set_acs_parameters will be called so the values of the ACS parameters will be
						changed to the given ones and the following global variebles and the names for the input/output will be set:
						Grid Size: grid_size
						Number of layers: numberOfLayers

						Netlist file: netlist_file
						Hypernets file: hypernet_file
						Nets file: nets_file
						Layers file: layers_file
						Output file: placement_file

						Flag which determines the type of the initial tour: heuristic_tour_flag */
	#if defined(DEBUG) || defined(PRINT)
		printf("------------ <parse_parameters> ------------\n");
	#endif
	
	static char usage[] =	"Usage: %s -g grid_size -c numberOfLayers -i input_netlist -h hypernets -n nets -e layers -p placement\n \
			[-r rho -a alpha -b beta -q q_0 -x xi -l lambda -f costFunction -w pheromoneUpd (mmas/acs) -z customized_heuristic (y/n)]\n \
			[-m maxIterations -t numberOfThreads -u initial_placement(random/heuristic) -s iteration_best_step -d tau_min_divisor]\n";
	int grid_size_flag = 0, netlist_flag = 0, layers_flag = 0, placement_flag = 0, error_flag = 0, zero_flag = 0, tour_flag = 0, opt;
	
	hypernets_flag = 0;
	nets_flag = 0;
	
	/* Default values for the ACO's parameters */
	double rho		= 0.5;		/* parameter for evaporation */
	double alpha	= 1;		/* importance of trail */
	double beta		= 2;		/* importance of heuristic evaluate */
	double q_0		= 0.95;		/* probability of best choice in tour construction */
	double xi		= 0.1;		/* local pheromone update rule */
	double lambda	= 0;		/* cost's function weight parameter */

	int maxIterations = 10;				/* maximum number of iterations to perform */
	int updateIterationBestStep = 3;	/* Every "updateIterationBestStep" number of iterations, use iteration_best ant to update pheromone matrix */
	double tau_min_divisor = 15.0;		/* Divisor for the initialization of tau_min */

	/* Default value for number of layers in the FPGA */
	numberOfLayers = 1;

	while ((opt = getopt(argc, argv, "g:c:i:h:n:e:p:r:a:b:q:x:l:t:m:s:d:f:w:u:z:")) != -1) {
		switch (opt) {
		case 'g':
			grid_size = atoi(optarg);
			if (grid_size == 0)
				grid_size_flag = 1;
			break;
		case 'c':
			numberOfLayers = atoi(optarg);
			if (numberOfLayers == 0)
				zero_flag = 1;
			break;
		case 'i':
			netlist_flag = 1;
			netlist_file = optarg;
			break;
		case 'h':
			hypernets_flag = 1;
			hypernet_file = optarg;
			break;
		case 'n':
			nets_flag = 1;
			nets_file = optarg;
			break;
		case 'e':
			layers_flag = 1;
			layers_file = optarg;
			break;
		case 'p':
			placement_flag = 1;
			placement_file = optarg;
			break;
		case 'r':
			rho = atof(optarg);
			if (rho == 0)
				zero_flag = 1;
			break;
		case 'a':
			alpha = atof(optarg);
			if (alpha == 0)
				zero_flag = 1;
			break;
		case 'b':
			beta = atof(optarg);
			break;
		case 'q':
			q_0 = atof(optarg);
			break;
		case 'x':
			xi = atof(optarg);
			break;
		case 'l':
			lambda = atof(optarg);
			break;
		case 'u':
			if (strcmp("heuristic", optarg) == 0)
				heuristic_tour_flag = TRUE;
			else if (strcmp("random", optarg) == 0)
				heuristic_tour_flag = FALSE;
			else
				tour_flag = 1;
			break;
		case 'm':
			maxIterations = atoi(optarg);
			if (maxIterations == 0)
				zero_flag = 1;
			break;
		case 's':
			updateIterationBestStep = atoi(optarg);
			break;
		case 'd':
			tau_min_divisor = atof(optarg);
			break;
		case 'f':
			if (strcmp("wire_timing", optarg) == 0)
				wire_timing_flag = TRUE;
			else if (strcmp("quadratic_estimate", optarg) == 0)
				quadratic_estimate_flag = TRUE;
			else
				hops_flag = TRUE;
			break;
		case 'w':
			if (strcmp("mmas", optarg) == 0)
				mmas_flag = TRUE;
			else if (strcmp("acs", optarg) == 0)
				acs_flag = TRUE;
			break;
		case 't':
			numberOfThreads = atoi(optarg);
			break;
		case 'z':
			if (strcmp("y", optarg) == 0)
				customized_heuristic_flag = TRUE;
			else if (strcmp("n", optarg) == 0)
				customized_heuristic_flag = FALSE;
			break;
		case '?':
			error_flag = 1;
			break;
		}
	}

	if (grid_size_flag == 1) {
		/* Grid size has not been defined */
		fprintf(stderr, "%s: missing mandatory -g option (grid_size)\n", argv[0]);
		fprintf(stderr, usage, argv[0]);
		exit(EXIT_FAILURE);
	}
	if (netlist_flag == 0) {
		/* Input netlist file has not been defined */
		fprintf(stderr, "%s: missing mandatory -i option (input netilist file)\n", argv[0]);
		fprintf(stderr, usage, argv[0]);
		exit(EXIT_FAILURE);
	}
	else if (nets_flag == 0) {
		/* Input nets file has not been defined */
		fprintf(stderr, "%s: missing mandatory -n option (nets file)\n", argv[0]);
		fprintf(stderr, usage, argv[0]);
		exit(EXIT_FAILURE);
	}
	else if (layers_flag == 0) {
		/* Layers file has not been defined */
		if (numberOfLayers < 1) {
			fprintf(stderr, "%s: Invalid number of layers (%d)\n", argv[0], numberOfLayers);
			exit(EXIT_FAILURE);
		}
		if (numberOfLayers != 1) {
			fprintf(stderr, "%s: missing mandatory -e option (layers file)\n", argv[0]);
			fprintf(stderr, usage, argv[0]);
			exit(EXIT_FAILURE);
		}
	}
	else if (placement_flag == 0) {
		/* Output placement file has not been defined */
		fprintf(stderr, "%s: missing mandatory -p option (placement filename)\n", argv[0]);
		fprintf(stderr, usage, argv[0]);
		exit(EXIT_FAILURE);
	}
	else if (zero_flag) {
		fprintf(stderr, "%s: Invalid value! (Parameters -c, -r, -a and -m can not be zero)\n", argv[0]);
		fprintf(stderr, usage, argv[0]);
		exit(EXIT_FAILURE);
	}
	else if (tour_flag) {
		fprintf(stderr, "%s: Invalid value for initial tour!\n", argv[0]);
		fprintf(stderr, usage, argv[0]);
		exit(EXIT_FAILURE);
	}
	else if (error_flag) {
		fprintf(stderr, usage, argv[0]);
		exit(EXIT_FAILURE);
	}

	set_acs_parameters(rho, alpha, beta, q_0, xi, lambda, maxIterations, updateIterationBestStep, tau_min_divisor);

	#if defined(DEBUG) || defined(PRINT)
		if (optind < argc) {/* these are the arguments after the command-line options */
			printf("The following arguments haven't been processed:\n");
			for (; optind < argc; optind++)
				printf("argument: \"%s\"\n", argv[optind]);
		}
		else {
			printf("no arguments left to process\n");
		}
		printf("------------ </parse_parameters> ------------\n");
	#endif
	return;
}

void print_results(void) {
	/*	FUNCTION:		Prints the final results.
		INPUT:			none
		OUTPUT:			none
		SIDE EFFECTS:	none */
	#ifdef DEBUG
		printf("------------ <print_results> ------------\n");
	#endif
	printf("Placement's quality (iteration = %d):\n", iteration);
	printf("wiringCost	= %lf\n", best_so_far_ant.wiringCost);
	printf("timingCost	= %lf\n", best_so_far_ant.timingCost);
	printf("numberOfHops	= %lf\n", best_so_far_ant.numberOfHops);
	printf("quadraticEst	= %lf\n", best_so_far_ant.quadraticEstimate);
	printf("cost		= %lf\n", best_so_far_ant.cost);
	#if defined(DEBUG) || defined(PRINT)
		/* Print best ant's placement */
		printf("------------ <best_ant_placement> ------------\n");
		print_placement(best_so_far_ant.nodePosition);
		printf("------------ </best_ant_placement> ------------\n");
	#endif
	#ifdef DEBUG
		printf("------------ </print_results> ------------\n");
	#endif
	return;
}

void print_placement(int nodePositions[]) {
	/*	FUNCTION:		Prints the given node placement for all the layers of the FPGA.
		INPUT:			nodePositions matrix
		OUTPUT:			none
		SIDE EFFECTS:	none */
	int k;

	for (k = 0; k < numberOfLayers; k++) {
		printf(BOLDWHITE "LAYER %d:" RESET "\n", k);
		print_layers_placement(nodePositions, k);
		printf("\n");
	}
	return;
}

void print_layers_placement(int nodePositions[], int layer) {
	/*	FUNCTION:		Prints the given node placement for the specified layer of the FPGA.
		INPUT:			nodePositions matrix and layer to print
		OUTPUT:			none
		SIDE EFFECTS:	none */
	#ifdef DEBUG
		printf("------------ <print_placement> ------------\n");
	#endif

	int i, j, node, pos, row, col, **grid;
	//int grid[grid_size + 4 + 1][grid_size + 4 + 1];

	grid = (int **)malloc(sizeof(int *) * (grid_size + 4 + 1));
	for (i = 0; i <= grid_size + 4; i++)
		grid[i] = (int *)malloc(sizeof(int) * (grid_size + 4 + 1));

	for (i = 0; i <= grid_size + 4; i++)
		for (j = 0; j <= grid_size + 4; j++)
			grid[i][j] = 0;

	for (i = 1; i <= n; i++) {
		node = i;
		if (nodeLayer[node] != layer)
			continue;
		if (!inOuts[node]) {
			/* CLB */
			pos = nodePositions[node];
			row = pos / grid_size;
			col = pos % grid_size;
			
			if (col != 0 )
				row = row + 1;
			else
				col = grid_size;
				
			grid[row + 2][col + 2] = node;
		}
		else {
			/* I/O */
			pos = nodePositions[node];
			pos = pos - (grid_size * grid_size);
			assert(1 <= pos && pos <= 8 * grid_size);
			if ( pos <= 1 * grid_size ) {
				/* Top I/O blocks */
				row = 0;
				col = pos % grid_size;
				if (col == 0)
					col = grid_size;
			}
			else if ( pos <= 2 * grid_size ) {
				/* Top I/O blocks */
				row = -1;
				col = pos % grid_size;
				if (col == 0)
					col = grid_size;
			}
			else if ( pos <= 3 * grid_size ) {
				/* Right I/O blocks */
				row = pos % grid_size;
				if (row == 0)
					row = grid_size;
				col = grid_size + 1;
			}
			else if ( pos <= 4 * grid_size ) {
				/* Right I/O blocks */
				row = pos % grid_size;
				if (row == 0)
					row = grid_size;
				col = grid_size + 2;
			}
			else if ( pos <= 5 * grid_size ) {
				/* Bottom I/O blocks */
				row = grid_size + 1;
				col = pos % grid_size;
				if (col == 0)
					col = grid_size;
			}
			else if ( pos <= 6 * grid_size ) {
				/* Bottom I/O blocks */
				row = grid_size + 2;
				col = pos % grid_size;
				if (col == 0)
					col = grid_size;
			}
			else if ( pos <= 7 * grid_size ) {
				/* Left I/O blocks */
				row = pos % grid_size;
				if (row == 0)
					row = grid_size;
				col = 0;
			}
			else if ( pos <= 8 * grid_size ) {
				/* Left I/O blocks */
				row = pos % grid_size;
				if (row == 0)
					row = grid_size;
				col = -1;
			}
			
			grid[row + 2][col + 2] = node;
		}
	}
	
	/* Top I/O pads */
	if (print_block_numbers_flag) {
		printf("---- ");
	}
	else {
		printf("-- ");
	}
	for (i = 3; i <= grid_size + 2; i++) {
		if (grid[1][i] && grid[2][i]) {
			if (print_block_numbers_flag) {
				printf(RED "%2d" RESET CYAN "%2d " RESET, grid[1][i], grid[2][i]);
			}
			else {
				printf(RED "X" RESET CYAN "X " RESET);
			}
		}
		else if (grid[1][i]) {
			if (print_block_numbers_flag) {
				printf(RED "%2d" RESET "%2d ", grid[1][i], grid[2][i]);
			}
			else {
				printf(RED "X" RESET "* ");
			}
		}
		else if (grid[2][i]) {
			if (print_block_numbers_flag) {
				printf("%2d" CYAN "%2d " RESET, grid[1][i], grid[2][i]);
			}
			else {
				printf("*" CYAN "X " RESET);
			}
		}
		else {
			if (print_block_numbers_flag) {
				printf("%2d%2d ", grid[1][i], grid[2][i]);
			}
			else {
				printf("** ");
			}
		}
	}
	if (print_block_numbers_flag) {
		printf("---- \n");
	}
	else {
		printf("-- \n");
	}

	for (i = 3; i <= grid_size + 2; i++) {
		/* Left I/O pads */
		if (grid[i][1] && grid[i][2]) {
			if (print_block_numbers_flag) {
				printf(RED "%2d" RESET CYAN "%2d " RESET, grid[i][1], grid[i][2]);
			}
			else {
				printf(RED "X" RESET CYAN "X " RESET);
			}
		}
		else if (grid[i][1]) {
			if (print_block_numbers_flag) {
				printf(RED "%2d" RESET "%2d ", grid[i][1], grid[i][2]);
			}
			else {
				printf(RED "X" RESET "* ");
			}
		}
		else if (grid[i][2]) {
			if (print_block_numbers_flag) {
				printf("%2d" CYAN "%2d " RESET, grid[i][1], grid[i][2]);
			}
			else {
				printf("*" CYAN "X " RESET);
			}
		}
		else {
			if (print_block_numbers_flag) {
				printf("%2d%2d ", grid[i][1], grid[i][2]);
			}
			else {
				printf("** ");
			}
		}

		/* CLBs */
		for (j = 3; j <= grid_size + 2; j++) {
			if (grid[i][j]) {
				if (print_block_numbers_flag) {
					printf(GREEN "%4d " RESET, grid[i][j]);
				}
				else {
					printf(GREEN "XX " RESET);
				}
			}
			else {
				if (print_block_numbers_flag) {
					printf("%4d ", grid[i][j]);
				}
				else {
					printf("** ");
				}
			}
		}

		/* Right I/O pads */
		if (grid[i][grid_size + 3] && grid[i][grid_size + 4]) {
			if (print_block_numbers_flag) {
				printf(RED "%2d" RESET CYAN "%2d " RESET, grid[i][grid_size + 3], grid[i][grid_size + 4]);
			}
			else {
				printf(RED "X" RESET CYAN "X " RESET);
			}
		}
		else if (grid[i][grid_size + 3]) {
			if (print_block_numbers_flag) {
				printf(RED "%2d" RESET "%2d ", grid[i][grid_size + 3], grid[i][grid_size + 4]);
			}
			else {
				printf(RED "X" RESET "* ");
			}
		}
		else if (grid[i][grid_size + 4]) {
			if (print_block_numbers_flag) {
				printf("%2d" CYAN "%2d " RESET, grid[i][grid_size + 3], grid[i][grid_size + 4]);
			}
			else {
				printf("*" CYAN "X " RESET);
			}
		}
		else {
			if (print_block_numbers_flag) {
				printf("%2d%2d ", grid[i][grid_size + 3], grid[i][grid_size + 4]);
			}
			else {
				printf("** ");
			}
		}
		printf("\n");
	}

	/* Bottom I/O pads */
	if (print_block_numbers_flag) {
		printf("---- ");
	}
	else {
		printf("-- ");
	}
	for (i = 3; i <= grid_size + 2; i++) {
		if (grid[grid_size + 3][i] && grid[grid_size + 4][i]) {
			if (print_block_numbers_flag) {
				printf(RED "%2d" RESET CYAN "%2d " RESET, grid[grid_size + 3][i], grid[grid_size + 4][i]);
			}
			else {
				printf(RED "X" RESET CYAN "X " RESET);
			}
		}
		else if (grid[grid_size + 3][i]) {
			if (print_block_numbers_flag) {
				printf(RED "%2d" RESET "%2d ", grid[grid_size + 3][i], grid[grid_size + 4][i]);
			}
			else {
				printf(RED "X" RESET "* ");
			}
		}
		else if (grid[grid_size + 4][i]) {
			if (print_block_numbers_flag) {
				printf("%2d" CYAN "%2d " RESET, grid[grid_size + 3][i], grid[grid_size + 4][i]);
			}
			else {
				printf("*" CYAN "X " RESET);
			}
		}
		else {
			if (print_block_numbers_flag) {
				printf("%2d%2d ", grid[grid_size + 3][i], grid[grid_size + 4][i]);
			}
			else {
				printf("** ");
			}
		}
	}
	if (print_block_numbers_flag) {
		printf("---- \n");
	}
	else {
		printf("-- \n");
	}

	for (i = 0; i <= grid_size + 4; i++)
		free(grid[i]);
	free(grid);

	#ifdef DEBUG
		printf("------------ </print_placement> ------------\n");
	#endif
	return;
}

void export_placement(int nodePositions[], int **grid) {
	/*	FUNCTION:		Exports the given node placement to the specified file.
		INPUT:			nodePositions matrix to export and grid matrix which holds the clb's availability status.
		OUTPUT:			none
		SIDE EFFECTS:	Opens a text file (whose name was specified by the -p parameter)
						for writing, if it does not exist then a new file is created, and
						stores the given placement following the format given here:
						http://www.eecg.toronto.edu/~vaughn/challenge/netlist.html */
	int i, j, k, node, pos, row, col, subblock, layer, found;
	char name[MAX_BLOCK_NAME];
	gridPosition_struct position;

	printf("Exporting Placement(%d)\t", iteration);

	FILE *fp = fopen(placement_file, "w"); /*Opens a text file for writing, if it does not exist then a new file is created. Writing starts from the beginning of the file. */
	if (numberOfLayers == 1) {
		fprintf(fp, "Netlist file: %s   Architecture file: 4lut_sanitized.arch\n", netlist_file);
		fprintf(fp, "Array size: %d x %d logic blocks\n", grid_size, grid_size);
		fprintf(fp, "\n");
		fprintf(fp, "#block name\t\t\t\tx\ty\tsubblk\t\tblock number\n");
		fprintf(fp, "#----------\t\t\t\t--\t--\t------\t\t------------\n");
	}
	else {
		fprintf(fp, "Netlist file: %s   Architecture file: 4lut.arch\n", netlist_file);
		fprintf(fp, "Array size: %d x %d x %d logic blocks\n", grid_size, grid_size, numberOfLayers);
		fprintf(fp, "\n");
		fprintf(fp, "#block name\t\t\t\tx\ty\tz\tsubblk\t\tblock number\n");
		fprintf(fp, "#----------\t\t\t\t--\t--\t--\t------\t\t------------\n");
	}

	/* Place global signal in random locations */
	for (k = 1; k <= numberOfGlobalSignals; ++k) {
		strcpy(name, globalSignals[k].name);
		if (globalSignals[k].type == CLB) {
			/* CLB spot */
			pos = randomInRange(1, grid_size * grid_size);	/* Choose a CLB spot randomly */
			layer = randomInRange(0, numberOfLayers - 1);	/* Choose a layer randomly */
			while (grid[layer][pos] != TRUE) {				/* Until you find one that is available */
				pos = randomInRange(1, grid_size * grid_size);
				layer = randomInRange(0, numberOfLayers - 1);
			}
			/* Choose a random CLB just to give as first parameter in findPositionInGrid */
			node = -1;
			for (i = 1; i <= n; i++) {
				if (!inOuts[i]) {
					node = i;
					position = findPositionInGrid(i, pos);
					break;
				}
			}
			assert(node != -1);
			grid[layer][pos] = FALSE; /* Occupy that spot */
		}
		else {
			/* Choose a random IO just to give as first parameter in findPositionInGrid */
			node = -1;
			for (i = 1; i <= n; i++) {
				if (inOuts[i]) {
					node = i;
					break;
				}
			}
			assert(node != -1);
			/* I/O spot */
			pos = randomInRange(grid_size * grid_size + 1, logic_block_array_size);	/* Choose an IO spot randomly */
			layer = randomInRange(0, numberOfLayers - 1);							/* Choose a layer randomly */
			position = findPositionInGrid(node, pos);
			subblock = position.subblock;
			while (grid[layer][pos] != TRUE || subblock == 1) {							/* Until you find one that is available */
				pos = randomInRange(grid_size * grid_size + 1, logic_block_array_size);
				layer = randomInRange(0, numberOfLayers - 1);
				position = findPositionInGrid(node, pos);
				subblock = position.subblock;
			}
			grid[layer][pos] = FALSE; /* Occupy that spot */
			
		}
		row = position.row;
		col = position.col;

		if (numberOfLayers == 1)
			fprintf(fp, "%s\t\t\t\t%d\t%d\t\t%d\t\t\t#\n", name, row, col, subblock);
		else
			fprintf(fp, "%s\t\t\t\t%d\t%d\t%d\t\t%d\t\t\t#\n", name, row, col, layer, subblock);
	}

	/* Export placement of the rest of the nodes */
	for (i = 1; i <= n; i++) {
		node = i;
		/* Find the name of the block */
		found = FALSE;
		for(j = 1; j <= numberOfBlocks; j++) {
			if (blocks[j].number == node) {
				strcpy(name, blocks[j].name);
				found = TRUE;
				break;
			}
		}
		if (!found) {
			printf("Unable to find nodes's [%d] name!\n", node);
			exit(EXIT_FAILURE);
		}
		/* 	Find block's coordinates */
		pos = nodePositions[node];
		position = findPositionInGrid(node, pos);
		row = position.row;
		col = position.col;
		subblock = position.subblock;
		if (subblock == 1)
			if (subblock0Available(pos, grid, nodeLayer[node]))
				subblock = 0;
		if (numberOfLayers == 1)
			fprintf(fp, "%s\t\t\t\t%d\t%d\t\t%d\t\t\t#%d\n", name, row, col, subblock, node);
		else
			fprintf(fp, "%s\t\t\t\t%d\t%d\t%d\t\t%d\t\t\t#%d\n", name, row, col, nodeLayer[node], subblock, node);
	}


	fclose(fp);
	printf("[DONE]\n");
	return;
}