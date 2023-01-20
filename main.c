#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"


/*! \file main.c
 *  \brief start of the program
 */

/*!
 *  This function initializes the MPI communication packages, and sets
 *  cpu-time counters to 0.  Then begrun() is called, which sets up
 *  the simulation either from IC's or from restart files.  Finally,
 *  run() is started, the main simulation loop, which iterates over
 *  the timesteps.
 */
int main(int argc, char **argv)
{
  double t0, t1;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &NTask);

  if(NTask <= 1)
    {
      if(ThisTask == 0)
	printf
	  ("Note: This is a massively parallel code, but you are running with 1 processor only.\nCompared to an equivalent serial code, there is some unnecessary overhead.\n");
    }

  for(PTask = 0; NTask > (1 << PTask); PTask++);

  if(argc < 2)
    {
      if(ThisTask == 0)
	{
	  printf("Parameters are missing.\n");
	  printf("Call with <ParameterFile> [<RestartFlag>]\n");
	}
      endrun(0);
    }

  strcpy(ParameterFile, argv[1]);

  if(argc >= 3)
    {RestartFlag = atoi(argv[2]);}
  else
    {RestartFlag = 0;}

  All.CPU_TreeConstruction = All.CPU_TreeWalk = All.CPU_Gravity = All.CPU_Potential = All.CPU_Domain =
  All.CPU_Snapshot = All.CPU_Total = All.CPU_CommSum = All.CPU_Imbalance = All.CPU_Hydro =
  All.CPU_HydCompWalk = All.CPU_HydCommSumm = All.CPU_HydImbalance =
  All.CPU_EnsureNgb = All.CPU_Predict = All.CPU_TimeLine = All.CPU_PM = All.CPU_Peano = 0;

  CPUThisRun = 0;

  t0 = second();

#ifdef SINK
SinkFlag = 0 ;
#endif 


  begrun();			/* set-up run  */

  t1 = second();
  CPUThisRun += timediff(t0, t1);
  All.CPU_Total += timediff(t0, t1);

  run();			/* main simulation loop */

  MPI_Finalize();		/* clean up & finalize MPI */

  return 0;
}




