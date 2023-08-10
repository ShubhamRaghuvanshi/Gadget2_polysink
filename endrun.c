#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <signal.h>
#include <unistd.h>

#include "allvars.h"
#include "proto.h"


/*! \file endrun.c 
 *  \brief Termination of simulation
 *
 *  This file contains routines for termination of the simulation.
 */

/*!  This function aborts the simulations. If a single processors wants an
 *   immediate termination, the function needs to be called with ierr>0. A
 *   bunch of MPI-error messages may also appear in this case.  For ierr=0,
 *   MPI is gracefully cleaned up, but this requires that all processors
 *   call endrun().
 */
void endrun(int ierr)
{
  if(ierr)
    {
      printf("task %d: endrun called with an error level of %d\n\n\n", ThisTask, ierr);
      fflush(stdout);

// 			#ifdef SINK
//    	char filename[100];
//    	sprintf(filename, "gaslist_%d.txt", ThisTask);
//    	FILE *gaslist;
//    	gaslist = fopen(filename,"w");
//    	for(int i=0; i< NumPart + N_accrete; i++){
//    		fprintf(gaslist, "%d	%d 	%d 	%d 	%d\n", i, NumPart, N_accrete, P[i].ID, P[i].Type);
//    	}
//    	fclose(gaslist);
//    	#endif


#ifdef DEBUG
      terminate_processes();
      raise(SIGABRT);
      sleep(60);
#else
      MPI_Abort(MPI_COMM_WORLD, ierr);
#endif
      exit(0);
    }

  MPI_Finalize();
  exit(0);
}








