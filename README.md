# README:

### Emergent multi*level* and multi*scale* selection in a computational model of the evolution of altruism

# Description of the simulation code, 1D and 2D

Rutger Hermsen, *r.hermsen@uu.nl* (2021)

## Introduction

This document contains a brief description of the code provided in this repository. 

The code implements a simulation of an agent-based model of the evolution of altruism and analyzes the results in various ways. This code was used to obtain simulation results presented in the following two manuscripts:

* Doekes, Hilje .M. and Rutger Hermsen  *Multiscale selection in spatially structured population*, in preparation, 2021.
* Hermsen, Rutger.  *Emergent multilevel selection in a simple model of the evolution of altruism*, PLOS Computational Biology 18, no. 10 (October 25, 2022): e1010612. https://doi.org/10.1371/journal.pcbi.1010612.

Detailed descriptions of the model and methods can be found in these texts. If despite these sources and the current document any questions remain, never hesitate to send us an email.

In the text below and in the code, the terms "trait value" and "level of altruism" are used interchangeably. In the context of multilevel selection theory, we frequently refer to colonies of individuals as "collectives", following Okasha (2006). 

## Two versions of the model

Two versions of the program are included: one simulates the dynamics of organisms in a two-dimensional habitat (folder `2D`), the other in a one-dimensional one (folder `1D`). 

### The 2D version

The folder `2D` contains code for the 2D model. In this version, organisms live in a two-dimensional space with periodic boundary conditions. Under the right parameter regime, once a minimal level of altruism has evolved the organisms spontaneously form a hexagonal pattern of colonies that themselves act as replicators (Hermsen, 2022).

### The 1D version

The folder `1D`  contains the code for the 1D model. In this version, organisms live in a one-dimensional space, again with periodic boundary conditions. The behavior in 1D is similar as in 2D case: colonies form at regular distances. Simulations of the 1D model are faster than those of the 2D model and in 1D the colonies can be identified more easily. This version of the model was therefore used for detailed analysis (Hermsen, 2022) and hence its code is more elaborate.

## Compilation

For each version  of the model, a `makefile` is provided. To compile under Linux, run `make` on the command line; this produces an executable file `altr.out`. The makefile assumes that you are using the GNU compiler `gfortran`, but other compilers should also be able to compile the code. 

## Structure of the code

The programs are written in Fortran. For convenience, they are split over several `*.f90` files.  The main programs are found in `altruism_default_1D.f90`  (1D version) and `altruism_default_2D.f90` (2D version),  which each rely on the other files (modules).  Input parameters and settings are specified in `parameters.f90`. The other modules contain functions and subroutines organized by function.

## Acknowledgements

All code was written by us, except for two modules: 

1. The module `fftpack5.1d.f90` for Fast Fourier Transforms was written by P. Swarztrauber, R. Valent, and J. Burkardt (Swarztrauber, 1984); further references are included in the file.
2. The module `ziggurat.f90` implements the fast Ziggurat method for normally distributed random numbers; it was originally written in C by G. Marsaglia & W.W. Tsang (2000) and translated to Fortran by Alan Miller. Details are in the file.

## What does the program do?

The program simulates the basic evolutionary dynamics as presented in the papers Hermsen *et al* (2022) and Doekes et al (2021). In addition, it performs analyses on the resulting dynamics, writing the results to comma-separated `*.txt` files (see section [Output](#output)). There is no graphical interface.

### Analyses implemented in both the 1D and the 2D version

The following analyses have been implemented in both the 1D and the 2D version.

#### Price equation: quantifying the effects of selection , drift, and transmission

The Price equation quantifies the effect of several evolutionary forces on a trait $\phi$.  More precisely, it describes the change in the mean trait value $\Delta\bar{\phi}$ over a certain time window $\Delta t = t_2 - t_1$ as a result of selection, transmission effects, and possibly other forces (drift, migration, …). See Marshall (2015) for an accessible introduction.

After each timestep of the simulation, all terms of the Price equation are calculated. The computational timestep is used as the time window $\Delta t$ of the Price equation. 

* In the standard version of the Price equation, the fitness of a "parent" organism at time $t_1$ is defined as the actual number of offspring it has at time $t_2$ (including itself, if it survives until $t_2$), sometimes called the realized fitness. In this case, the right-hand side of the equation has two term: 

  $$\Delta \bar{\phi} = S + T,$$

where $S$ is the selection differential and $T$ the transmission term. With this definition of fitness, random drift does not exist / is an integral part of the selection differential. This behavior can be chosen by setting `REALIZED_FITNESS = .true.` in `parameters.f90`. 

  Alternatively, one can define fitness of an ancestor at time $t_1$ as the *expected* number of offspring it has at time $t_2$. Any deviation from that expectation is then attributed to random drift $D$, in which case we can write

  $$\Delta \phi = S + D + T.$$

(See Okasha (2006) and Hermsen *et al* (2021) for details.) This behavior is obtained by setting `REALIZED_FITNESS = .true`. 

* Because $S$, $D$, and $T$ vary from one timestep to the next, the code also calculates a moving average of these quantities, using a bin width `MEAN_INTERVAL` that is currently set to 1000 timesteps.

*  Lastly, cumulative values of $S$, $D$, and $T$ are stored; if the simulation is starts with $\bar{\phi}=0$, the value of $\bar{\phi}$ at time $t$ is precisely the sum of the cumulative values of $S$, $D$, and $T$ at time $t$.

#### Price equation: quantifying purifying/stabilizing selection

Each time step, also a value is calculated for the *purifying* or *stabilizing* selection acting during the previous time step. Here, the purifying selection is quantified as the effect of selection on the variance in phenotype – the mean squared deviation of $\phi$​​ from the mean of $\bar{\phi}$ (Rice, 2004)). Again, a moving average and a cumulative value are stored as well. 

#### Radial distribution function $g(r)$

In both the 1D and the 2D version, the population spontaneously organizes into colonies/collectives/clusters. To quantify this clustering and to analyze the scales involved in the spatial patterns, the code can calculate the *radial distribution function*, often denoted as $g(r)$.  This analysis is performed every `REPORT_INTERVAL` timesteps, provided the parameter `OUTPUT_G_OF_R` is set to `.true.`.

#### Multi*scale* selection

The program can additionally perform the multi*scale* selection analysis as described by Doekes *et al* (2021). The scales $\sigma$ used for this analysis are hard-coded in the module `kernels.f90`, function `choose_s_curve_vals()`. The total number of scales defined there has to be consistent with the parameter `NR_SCALES` defined in `parameters.f90`. The code calculates $S_\mathrm{local}$ and $S_\mathrm{interlocal}$ for each scale; we here refer to the results as $S$-curves. 

For our purposes, it was sufficient to calculate $S$-curves only after the evolutionary dynamics had equilibrated and the mean phenotype no longer changed systematically. The timepoints at which $S$-curves are calculated are determined by two parameters in `parameters.f90`: calculation of $S$-curves starts no earlier than `NR_AV` timesteps before the end of the simulation, and is then performed once every `SC_IN` timesteps, but in such a way that the last analysis is performed right at the end of the simulation. The program also calculates the mean and standard deviation of all $S$-curves.

The multiscale selection analysis relies on the Local Selection Differential (LSD), which in turn depends on a kernel function (see Doekes *et al* (2021). The program supports a choice between two kernel functions: a step function or a normal distribution. To choose the normal distribution, set `STEP_FUNCTION_KERNEL` to `.false.` in `parameters.f90`. 

### Analyses implemented exclusively in the 1D version

In the 1D version of the code, we also implemented two methods to quantify multilevel selection, referred to as Multilevel Selection (MLS) 1 and 2.  (See Damuth and Heisler (1988), Okasha (2006), and Hermsen (2021).)

Because both methods are based on versions of the Price equation, they again are involve a timestep $\Delta t$. Because the demographics – births and deaths – of colonies/collectives take place at a longer time scale than the demographics of individuals, we choose a larger timestep of `ANALYSIS_INTERVAL` timesteps, currently set to 1000. Because of this large timestep, the expected number of offspring of an individual or a colony can no longer be calculated analytically. Therefore we use the *realized* number of offspring of an individual or colony as fitness (rather than the expectation value).

#### Multilevel Selection 1

MLS 1 is based on the selection differential $S$ of the Price equation for the evolution of the mean phenotype of individuals. It splits the selection differential $S$ into two parts, acting within and among groups:

$$S = S_\mathrm{within} + S_\mathrm{among}.$$

Throughout the simulation, these terms are calculated after each `ANALYSIS_INTERVAL` timesteps.  The method relies on the automated recognition of colonies/collectives, which is done as detailed in Hermsen et al (2021).

#### Multilevel Selection 2

MLS 2 is based on the selection differential $S$ and the transmission term $T$ of the Price equation describing the evolution of the mean trait value of colonies/collectives, which is in turn defined as the mean of the trait values of the organisms that belong to the collective. As such, the relevant fitness is the fitness of the colonies, *i.e.*, the number of offspring colonies each ancestor colony has after the timestep $\Delta t$ set by `ANALYSIS_INTERVAL`. Now, $S$ should be interpreted as the effect of selection at the colony level, whereas $T$ represents the internal evolution of colonies.   Again, throughout the simulation, these terms are calculated after each `ANALYSIS_INTERVAL` timesteps.

Like MLS 1, MLS 2 requires the automated recognition of colonies/collectives.  In addition, to determine the fitness of colonies, the program needs to track lineages of colonies, and hence birth and death events of colonies.  This is done as described in Hermsen *et al* (2021).  

## Relevant parameters

We here discuss several parameters that can be set in `parameters.f90` but are not mentioned above or below:

* `RESOLUTION`, `DT`, and `N`: The simulation uses discrete space and time (a grid), but is designed to approximate a (stochastic) process in continuous space-time. The approximation is better if a finer grid is used. The grid is controlled by parameters `RESOLUTION`, specifying the number of grid cells that fit in the range of altruism, and `DT`, the timestep in units such that the death rate of individuals is 1. The linear size of the grid (in grid cells) is `N`. The Fast Fourier method used repeatedly by the code is faster if `N` is a product of small primes; we therefore recommend to choose a power of 2.
* `T_MAX`: The number of time steps to be simulated.
* `RATE_DEATH`: The death rate; is set to 1 by definition of the unit of time and therefore should not be altered.
* `RATE_GROWTH`: The basal growth rate, in the absence of competition and altruism.
* `RANGE_SOCIAL`: The range of the social interaction (altruism); is set to 1 by definition of the unit of length and therefore should not be altered.
* `RANGE_COMPETITION`: The range of competitive interactions. 
* `P_MUT`: Probability of a mutation in the phenotype for each newborn individual.
* `MEAN_SIZE_MUT`: Mutations have a random effect on phenotype, drawn from an exponential distribution with mean `MEAN_SIZE_MUT`. The sign is subsequently chosen at random.
* `K_DIFF`: Individuals move by diffusion with diffusion constant `K_DIFF`.
* `K_COMPETITION`: To simulate resource competition, the reproduction rate of individuals reduces linearly with the local density. `K_COMPETITION` is the $x$-intercept: the local density at which reproduction is fully suppressed. It scales the carrying capacity of the habitat.
* `COST`: The reproductive deficit due to altruism; is set to 1 by definition of the scale of the trait $\phi$ and therefore should not be altered.
* `MAX_BENEFIT`: The reproduction-rate benefit gleaned from being in the presence of altruists saturates at this value.
* `BENEFIT`: If the level of altruism experienced by an individual is in the low limit, the benefit grows linearly with the level of altruism. Parameter `BENEFIT` is the slope of that linear increase.
* `WELL_MIXED`: To demonstrate that the emergent spatial structure is essential for the results, one can run simulations in which the position of individuals is randomized each time step. To do so, set `WELL_MIXED` to `.true.`.  The result should be similar to setting a very large value of `K_DIFF`. 
* `seed`: Random seed used to initialize the random-number generator. Set explicitly for reproducibility; should be varied to produce independent replicates.

## Initial conditions

Initially, an estimate is made of the total number of organisms that can be supported by the simulation field given the parameters. This number of organisms is then placed at uniformly random positions on the grid. Initially, each organism has trait value $\phi = 0$.

## Output

### Output at initiation

When the program is run,  several preparations are made before the actual simulation starts. The program reports the progress of these steps to the standard output.

### Progress reports and intermediate results at regular intervals

Next, at regular intervals (set by parameter `REPORT_INTERVAL`, currently set to 5000 simulation timesteps) the simulation reports its progress by printing the following information to the standard output: the current *timestep*, *population size*, and *mean level of altruism*. 

At these times, several quantities are also calculated and printed to files. The desired output can be selected by switching the corresponding Boolean parameters in `parameters.f90`. 

* Various values quantifying the evolutionary process are written to the file `selection.txt`. Among others, its columns include *timestep*, *time*, *population size*, *mean phenotype* (level of altruism), *mean absolute fitness* of the organisms, *selection*, *drift*, and *transmission* terms of the Price equation (calculated over the previous computational time step), *purifying selection*, *rolling means* of these quantities (averaged over the last `MEAN_INTERVAL` timesteps, currently set to 1000), and *cumulative* values for these quantities. 
* If `OUTPUT_COUNT_GLOBAL_LIST` is selected in `parameters.f90`, the program writes a simple tab-separated file with only the variables *time*, *population size*, and *mean phenotype*, as `countglobal.txt`.
* If `OUTPUT_O_LIST` is selected in `parameters.f90`, the program writes a list with key properties of all organisms to `o_list*.txt`: *position* in real units, position specified as the *grid-point*, *phenotype*, and *relative fitness*.
* If `OUTPUT_OCC_FIELD` is selected in `parameters.f90`, the program writes a vector (1D version) or matrix (2D version) with the number of organisms that currently reside on each grid cell (the "occupancy") to  `occ_field*.txt`. 
* If `OUTPUT_P_FIELD` is selected in `parameters.f90`, the program writes a vector (1D version) or matrix (2D version) with, for each grid cell, the sum of the trait values of organisms that reside on that grid cell:  `p_field*.txt`. 
* If `OUTPUT_ALTR_FIELD`is selected in `parameters.f90`, the program writes a vector (1D version) or matrix (2D version) to  `altr*.txt` with, for each grid cell, the level of altruism (perhaps the amount of a common good) that would be experienced by an organism at that site. 
* If `OUTPUT_HISTOGRAMS` is selected in `parameters.f90`, the program writes a binned frequency table of the phenotypes of organisms to the file `histograms.txt`. Frequency tables at different time points are added as new lines in the same file.
* If `OUTPUT_G_OF_R` is selected in `parameters.f90`, the program writes the pair-distribution function $g(r)$ to the file `g_or_r*.txt`.

### Results of the multi*scale* selection analysis

In addition, if `OUTPUT_ALL_S_CURVES` is selected in `parameters.f90`, the program outputs $S_{\mathrm{local}}$ and $S_{\mathrm{interlocal}}$ for the selected scales $\sigma$ at selected times (see section Multiscale selection above); all values are written to `S_curves.txt`. The *averages* of these values over all time points at which they were calculated are written to `S_curves_mean.txt`.

### Results of the multi*level* selection analyses (1D version only)

Additional output can be selected in the 1D version to report results for Multilevel Selection 1 and 2.

* To identify colonies, a Kernel Density Estimate (smoothed occupancy field) is calculated every time the analyses of MLS 1 and 2 are performed. If `OUTPUT_SMOOTHED_FIELD` is selected, these densities are written to `smoothed_field.txt`. The resolution (number of "pixels") is set by `RES_SMOOTH`.
* The values of the various selection terms calculated using the analyses MLS 1 and 2 are written to `selection_coll.txt`
* If `OUTPUT_COLL` is selected, several properties of colonies/collectives are written to file `coll_list_time.txt` each time the analyses of MLS 1 and 2 are performed. In particular: *timestep*, *time*, *left border* of the colony, *right border* of the colony, *center of mass* of the colony, *mean altruism level* of the colony.
* If `OUTPUT_COLL_DEMOGRAPHICS` is selected, all birth and death events of colonies are written to the file `coll_demographics.txt`. In column `type`, "b" indicates birth, "d" death of a colony.
* If `OUTPUT_COLL_LIST_TIME_LINES` is selected, a file with the name `coll_list_time_lines.txt` is created. It specifies line segments in space-time that connect ancestor colonies with their offspring, thus creating a phylogenetic tree in space-time at the level of colonies.

Note that, if many of the output options are switched on (`.true.`) and/or the `REPORT_INTERVAL` is set to a small value, the program produces a large number of files taking up many gigabytes of disk space.

 ## References

* Damuth, J., and I. L. Heisler  *Alternative Formulations of Multilevel Selection*, Biology and Philosophy 3, no. 4 (October 1, 1988): 407–30.
* Doekes, H.M. & R. Hermsen  *Multiscale selection in spatially structured population*, in preparation, 2021.
* Hermsen, R. *Emergent multilevel selection in a simple model of the evolution of altruism*, PLOS Computational Biology 18, no. 10, 2022.
* Marsaglia, G. & W.W. Tsang. *The ziggurat method for generating random variables* , J. Statist. Software, v5(8), 2000.
* Marshall, J.A.R. *Social Evolution and Inclusive Fitness Theory: An Introduction*, Princeton University Press, 2015. 
* Rice, S.H.  *Evolutionary Theory: Mathematical and Conceptual Foundations*, Sinauer Associates, 2004.
* Okasha, S. *Evolution and the Levels of Selection*, Oxford University Press, 2006.
* Swarztrauber, P.  *Fast Fourier Transform Algorithms for Vector Computers*, Parallel Computing, pages 45-63, 1984.
