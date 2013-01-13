This port was made from the FORTRAN code into the NEURON enviroment based on 

	Traub RD, Buhl EH, Gloveli T, Whittington MA. 
	Fast Rhythmic Bursting Can Be Induced in Layer 2/3 Cortical Neurons by 
	Enhancing Persistent Na(+) Conductance or by Blocking BK Channels.
	J Neurophysiol. 2003 Feb; 89(2):909-21

This port was made by Roger D Traub and Maciej Lazarewicz (mlazarew@seas.upenn.edu)

Thanks to Ashlen P Reid for help with porting a morphology of the cell.

This port is reproducing figures 2, 4, 5, 6 and 7.

Under unix systems:
-------------------
to compile the mod files use the command

 nrnivmod

and run the simulation hoc file with the command

 nrngui pyr3.hoc


Under Windows systems:
----------------------

to compile the mod files use the "mknrndll" command.
A double click on the simulation file

 pyr3.hoc

will open the simulation window.

Questions on how to use this port should be directed to mlazarew@seas.upenn.edu