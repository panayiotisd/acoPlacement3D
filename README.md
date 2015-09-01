ACO-Placement3D
===========

This is the README file for ACO-Placement3D v1.0, a novel algorithm for 3-D FPGA
placement based on Ant Colony Optimization.

Copyright & Licensing
------------------------

If you use ACO-Placement3D in your research, I would appreciate a citation in your
publication(s). Please cite it as:  

Panayiotis Danassis. ACO-Placement3D, Version 1.0. Available from
<http://FIXME>, 2015.

The software is licensed under the MIT License:

Copyright (c) 2015 Panayiotis Danassis

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

Authors
--------

This software was developed by [Panayiotis Danassis](http://panosd.eu)
with the valuable contribution of [Kostas Siozios](http://proteas.microlab.ntua.gr/ksiop/)

Manifest
--------

Main control routines:

- acoPlacement.c
- acoPlacement.h

Implementation of ants' procedures:

- ants.c
- ants.h

Input / output routines:

- inOut.c
- inOut.h

Time measurement:

- timer.c
- timer.h

Auxiliary routines:

- utilities.c
- utilities.h

Other:

- Makefile
- README.md

Description
------------

ACO-Placement3D is a novel placer, based on ant colony optimization (ACO), targeting
3-D FPGAs. It is a distributed algorithm, based on indirect communication (stigmergy)
between artificial ants which work asynchronously to build a feasible placement. As so it is
characterized by inherent parallelism and can greatly benefit from multi-core processors.
It combines a positive feedback mechanism with stochastic decision policy and swarm
intelligence.

The software was developed as part of my diploma thesis for the National Technical
University of Athens (NTUA), school of Electrical and Computer Engineering. For more
information please visit the following link:
<http://FIXME>

Installation
------------

Use the provided Makefile to compile the program under Linux. Executable
"acoPlacement3D" is produced. There are four different compilation modes, which are
presented in detail as follows:

__Mode 1 (Default installation):__  
This is the default installation mode. Optimization flag '-O3' is used.  
To compile in Default Mode  just type:

	make

__Mode 2 (Print Mode):__  
In "Print Mode" the program displays, using a simple text based graphical representation,
the best placement found so far, in order to visualize and give a better perspective of the
solution found. It also informs you every time a better solution is found. Optimization flag
'-O3' is used.   
To compile in Print Mode type:

	make print

__Mode 3 (Debug Mode):__  
In "Debug Mode" the software prints a ton of debugging info. More specifically:

- The values of various key variables.
- The input netlist.
- The input/output pins.
- The global signals.
- The fanin and fanout table for every block.
- Exports in a file called 'parse_netlist.out' the input netlist with the same format as the
input file (with the exception of global signals).  
- The hypernets (paths), if given.
- The nets.
- States explicitly almost every function call and return in order to track the progress of
the program.

No optimization flags are used. Also the software is compiled with the '-g' flag in order to
generate debugging information to be used by a debugger such as GDB.  
To compile in Debug Mode type:

	make debug

__Mode 4 (Parallel Mode):__  
In this mode the appropriate compiler flag is used ('-fopenmp') to "turn on" OpenMP, thus
effectively parallelizing the application. Other than that, the "Parallel Mode" is similar to the
"Print Mode". Optimization flag '-O3' is used here as well.  
To compile in Parallel Mode type:

	make parallel
	
__Cleanup:__  
Use the following command to clean the directory from the executable and all the already
compiled object files:

	make clean

Usage
-------

To display the help text with the usage of the executable, type:

	./acoPlacement3D
	
The following text will be displayed:  
Usage: ./acoPlacement3D -g grid_size -c numberOfLayers -i input_netlist -h hypernets -n nets -p placement -e layers -v TSVs  
	[-r rho -a alpha -b beta -q q_0 -x xi -l lambda -f costFunction -w pheromoneUpd (mmas/acs) -z dynamic_heuristic (y/n)]  
	[-m maxIterations -t numberOfThreads -u initial_placement(random/heuristic) -s iteration_best_step -d tau_min_divisor]  

Breaking down the above list, the command line options that the executable
'acoPlacement3D' provides, are the following:

__Mandatory Options:__  

* -g		grid size.
* -i		(input) netlist file name.
* -n		(input) nets file name.
* -p		(output) placement file name.
* -e		(input) layers file (file that contains the layer that each block will be placed).
			Mandatory only if number of layers > 1. The file must be in the same format as the
			placement file. (see 'Output' section for the format of the placement file)

__Other Options:__  

* -m		maximum number of iterations to perform (termination condition).
* -c		number of layers (for 3D placement).
* -h		(input) hypernets (paths) file name.
* -v		(input) TSVs locations (file containing the location os the TSVs in a heterogeneous fabric.
			(see "TSV File Format" section for the format of the file)
* -t		number of threads (for parallel - OMP version).
* -u		initial placement('random', 'heuristic'). Chooses between a random initial
			placement or one based only on heuristic information.
* -f		changes cost function (available choices: 'wire_timing', 'quadratic_estimate', 'hops').
* -l		lambda, changes the relative importance between wire cost (bounding box) and
			timing cost in the 'wire_timing' cost function.
			cost = lambda * timing_cost + (1 - lambda) * wire_cost;
* -w		pheromoneUpd, changes between the two implemented pheromone update routines
			('mmas', 'acs'). Type 'mmas' to use the pheromone update rule of the MAX-MIN
			Ant System (MMAS) or 'acs' to use the rule of Ant Colony System (ACS). 
* -z		dynamic heuristic ('y'/'n'). Apply a dynamic heuristic criterion for the largest nets
			in order to improve the quality of the solution. A drawback is that runtime increases
			drastically.

__ACO Related Options:__  

* -r		rho, pheromone evaporation rate.
* -a		alpha, pheromone trail influence.
* -b		beta, heuristic information influence.
* -q		q_0, Ant Colony System's (ACS) pseudorandom proportional rule's parameter.
* -x		xi, Ant Colony System's (ACS) local pheromone update rule's evaporation rate.
* -s		iteration_best_step, defines the relative frequency with which we choose to update
			the trails based on either the best so far solution', or the 'iteration best solution'.
			Must be a positive non-zero value. E.g. if -s 1, then we only use 'iteration best
			solution', while if -s 999999 > max number of iterations, then we only use 'best so
			far solution'. If -s 3, then we use 'iteration best solution' every 3 iterations.
* -d		tau_min_divisor, defines the lower pheromone trail limit for the MAX-MIN Ant
			System (MMAS), and as a result changes the stagnation behavior of the algorithm.
			tau_min = tau_max / tau_min_divisor.

__Default Values:__  

As default, the pheromone update rules of MAX-MIN Ant System and the pseudorandom
proportional rule of Ant Colony System (ACS) are used. The cost function that the algorithm
tries to minimize is the total number of hops. Specifically the default values for the above
parameters are the following:

* -m		10
* -c		1
* -h		(null)
* -t		1
* -u		'random'
* -f		'hops'
* -l		0
* -w		'mmas'
* -z		'y'
* -r		0.1
* -a		1
* -b		2
* -q		0.95
* -x		0
* -s		3
* -d		15

__Pre-Compilation Options:__  

There are some parameters of the algorithm that are defined as constants. Those can be
found at the .h files and are the following:

* n_ants						number of ants in the colony (value: 256 ants)
* restart						reinitialize pheromone matrix to avoid stagnation
									(value: INT_MAX iterations = disabled)
* printStep					print results periodically (value: 5 iterations)
* exportPlacementStep	export placement file (value: 5 iterations)

Note that options -c, -r, -a and -m can not take zero value.

Examples for running a benchmark:  

./acoPlacement3D -i apex4.net -n apex4_net.echo -g 36 -p placement.p -m 10

or  

./acoPlacement3D -i s38417.net -n s38417_net.echo -g 81 -p placement.p -t random -q 0.9 -b 2 -m 10 -z n

Output
--------

Every run of the algorithm produces the following two files:  
placement.p
'input netlist's name'.heur  

The first one is the output placement file, which has the following format:  

__Placement File Format__  
The first line of the placement file lists the netlist file and the architecture description file.
The second line of the placement file gives the size of the logic block array (e.g. 20 x 20
logic blocks).

All the following lines have the format:

block_name    x        y   subblock_number

The block name is the name of this block, as given in the input netlist. x and y are the row
and column in which the block is placed, respectively. The subblock number is meaningful
only for pads. Since we can place two pads in a row or column the subblock number
specifies which of the possible pad locations (either location 0 or location 1) in row x and
column y contains this pad. Note that the first pad occupied at some (x, y) location is
always that with subblock number 0. For logic blocks (.clbs), the subblock number is always
zero.

The placement files also include a fifth field as a comment. You can ignore this field.

The second one is an auxiliary file with the values of the heuristic information. You can
ignore this file. It's only used in subsequent runs of the algorithm to avoid the
recomputation of the heuristic information and thus save some time.

__TSVs File Format__  
This file is used in heterogeneous fabric architectures and contains the locations for all
the TSVs (Through-Silicon Via). The file is just an array of (grid_size x grid_size) values
of {0, 1}. An example of such file is presented below:  

0 0 0 0 0 0 0 0 0 0  
0 0 0 0 0 0 0 0 0 0  
0 0 0 0 0 0 0 0 0 0  
0 0 0 1 1 1 1 0 0 0  
0 0 0 1 1 1 1 0 0 0  
0 0 0 1 1 1 1 0 0 0  
0 0 0 1 1 1 1 0 0 0  
0 0 0 0 0 0 0 0 0 0  
0 0 0 0 0 0 0 0 0 0  
0 0 0 0 0 0 0 0 0 0  

Known Bugs (& Future Work)
--------------------------------

- Our netlist parser does not take into account comments. __If your netlist file includes comments, the program is going to crash!__
  To solve this, use under Linux the following command to get rid of the comments:  
	sed -e 's/#.*$//' bench.net > bench_noComments.net  
- For some reason OpenMP implementation doesn't work in some systems. Specifically the
software works as it should if compiled under gcc version 4.7.2 (Debian 4.7.2-5) but if
compiled under gcc version 4.8.2 20140120 (Red Hat 4.8.2-16) the program crashes.
- Parallel implementation works only with 'Bounding Box' as cost function, due to
dependencies in the usage of global timing matrices.

Contribute
----------

- Issue Tracker: github.com/$project/$project/issues
- Source Code: github.com/$project/$project

Support
-------

If you are having issues, contact us at: