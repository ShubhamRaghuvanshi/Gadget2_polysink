#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"

/*! \file accel.c
 *  \brief driver routine to carry out force computation
 */


/*! This routine computes the accelerations for all active particles.
 *  First, the long-range PM force is computed if the TreePM algorithm is
 *  used and a "big" PM step is done.  Next, the gravitational tree forces
 *  are computed. This also constructs the tree, if needed.
 *
 *  If gas particles are present, the density-loop for active SPH particles
 *  is carried out. This includes an iteration on the correct number of
 *  neighbours.  Finally, the hydrodynamical forces are added.
 */
void compute_accelerations(int mode)
{
  double tstart, tend;
  FLOAT hacc_abs, gacc_abs;

  if(ThisTask == 0)
    {
      printf("Start force computation...\n");
      fflush(stdout);
    }

#ifdef PMGRID
  if(All.PM_Ti_endstep == All.Ti_Current)
    {
      tstart = second();
      long_range_force();
      tend = second();
      All.CPU_PM += timediff(tstart, tend);
    }
#endif

  tstart = second();		/* measure the time for the full force computation */

  gravity_tree();		/* computes gravity accel. */

  if(All.TypeOfOpeningCriterion == 1 && All.Ti_Current == 0)
    gravity_tree();		/* For the first timestep, we redo it
				 * to allow usage of relative opening
				 * criterion for consistent accuracy.
				 */
  tend = second();
  All.CPU_Gravity += timediff(tstart, tend);

#ifdef FORCETEST
  gravity_forcetest();
#endif

  if(All.TotN_gas > 0)
    {
      if(ThisTask == 0)
	{
	  printf("Start density computation...\n");
	  fflush(stdout);
	}

      tstart = second();


      density();		/* computes density, and pressure */


#ifdef SINK
#ifdef SETSINKBND
      if(All.TotN_sink > 0 ){
        MPI_Barrier(MPI_COMM_WORLD);           
        setdens();
        //endrun(901); 
      }
#endif
#endif
      for(int igas=0; igas<N_gas; igas++){
        if(isinf(SphP[igas].DhsmlDensityFactor) || isnan(SphP[igas].DhsmlDensityFactor) ){
          printf("********************************************* After density correction\n") ;
          printf("ThisTask : %d, Dhsmlrho is inf, pid: %d, rho: %g, numngb: %g, hsml: %g, divv: %g, curl: %g, dhsmlrho: %g\n", ThisTask, P[igas].ID, SphP[igas].Density, SphP[igas].NumNgb, SphP[igas].Hsml ,SphP[igas].DivVel, SphP[igas].CurlVel, SphP[igas].DhsmlDensityFactor);
          printf("posx: %g, posy: %g, posz: %g, type: %d \n", P[igas].Pos[0], P[igas].Pos[1], P[igas].Pos[2], P[igas].Type);
          printf("********************************************* After density correction\n") ; 
          endrun(999);
        }    
      }

      tend = second();
      All.CPU_Hydro += timediff(tstart, tend);

      tstart = second();
      force_update_hmax();      /* tell the tree nodes the new SPH smoothing length such that they are guaranteed to hold the correct max(Hsml) */
      tend = second();

      All.CPU_Predict += timediff(tstart, tend);

      if(ThisTask == 0)
	{
	  printf("Start hydro-force computation...\n");
	  fflush(stdout);
	}


      tstart = second();
      hydro_force();		/* adds hydrodynamical accelerations and computes viscous entropy injection  */        
      tend = second();
      All.CPU_Hydro += timediff(tstart, tend);



//#ifdef SINK
//#ifdef SETSINKBND
 //     if(All.TotN_sink > 0 ){
 //     savepositions(All.SnapshotFileCount++); 
 //     if(ThisTask == 0 ){
 //       printf("snapshot count : %d \n", All.SnapshotFileCount) ; 
 //     }  
 //     if(All.SnapshotFileCount > 100 ) {endrun(777);}
 //     }
//#endif 
//#endif 

    }

  if(ThisTask == 0)
    {
      printf("force computation done.\n");
      fflush(stdout);
    }
}
