#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <unistd.h>

#include "allvars.h"
#include "proto.h"

/*! \file run.c
 *  \brief  iterates over timesteps, main loop
 */

/*! This routine contains the main simulation loop that iterates over single
 *  timesteps. The loop terminates when the cpu-time limit is reached, when a
 *  `stop' file is found in the output directory, or when the simulation ends
 *  because we arrived at TimeMax.
 */
void run(void)
{
  FILE *fd;
  int stopflag = 0;
  char stopfname[200], contfname[200];
  double t0, t1;


  sprintf(stopfname, "%sstop", All.OutputDir);
  sprintf(contfname, "%scont", All.OutputDir);
  unlink(contfname);

  do				/* main loop */
    {
     t0 = second();



      find_next_sync_point_and_drift();	/* find next synchronization point and drift particles to this time.
					 * If needed, this function will also write an output file
					 * at the desired time.
					 */

      every_timestep_stuff();	/* write some info to log-files */


	
      domain_Decomposition();	/* do domain decomposition if needed */

      compute_accelerations(0);  /* compute accelerations for 
                                 * the particles that are to be advanced  
                                 *                                  *                                  */

/*

			double n_crit_cgs = 5e8, mh =  2.408544e-24;
			double rho =  n_crit_cgs * mh/All.UnitDensity_in_cgs;           
   		if(All.NumCurrentTiStep >=2000){
				for(int i=0; i<N_gas; i++){
					if( SphP[i].Density > rho){ 
					//	SphP[i].Density = rho;
					//	printf("Maximum density reached\n")	;			
						endrun(10945);			
					}
				}
   		}
*/

	//	 if(All.NumCurrentTiStep%1000 ==0 ){cutoff();}


	
#ifdef SINK
     // only check for accreted particles every 10 timesteps. this is sort of arbitrary.
     // don't destroy particles as soon as they're marked.

      if(All.NumCurrentTiStep >= All.CriticalNumstep){

          
     //   if(All.NumCurrentTiStep%100 ==0 ){cutoff();}
 
          
				if(All.NumCurrentTiStep%10 == 0){
					create_sink();
					MPI_Barrier(MPI_COMM_WORLD);
				}  

				if(All.TotN_sink > 0 ){
				//	All.ExternalPressure =0;	
//					RemoveParticles();
		      MPI_Barrier(MPI_COMM_WORLD);
		      identify_doomed_particles();			
		//			move_doomed_particles();
		     	destroy_doomed_particles();
		      MPI_Barrier(MPI_COMM_WORLD);


				}

//				if(All.TotN_sink > 1 && All.NumCurrentTiStep%100 == 0){
//					sinkmerger();
//				}


			N_sink = NumPart-N_gas;
      MPI_Barrier(MPI_COMM_WORLD);
			
			if(N_sink > 0){
				char filename[100];
				FLOAT sr, sv, sa, L ;
				
				for(int i=N_gas; i<NumPart; i++){
			  	FILE *sinklist;
					sprintf(filename, "%s/sink_%d.txt", All.OutputDir, P[i].ID);
					sinklist = fopen(filename,"a");

    			sr = (P[i].Pos[0]-SysState.CenterOfMass[0]) * (P[i].Pos[0]-SysState.CenterOfMass[0]) +
               (P[i].Pos[1]-SysState.CenterOfMass[1]) * (P[i].Pos[1]-SysState.CenterOfMass[1]) +
               (P[i].Pos[2]-SysState.CenterOfMass[2]) * (P[i].Pos[2]-SysState.CenterOfMass[2]);

					sr = sqrt(sr);
					sv = sqrt( P[i].Vel[0]*P[i].Vel[0] + P[i].Vel[1]*P[i].Vel[1] + P[i].Vel[2]*P[i].Vel[2]  );
					sa = sqrt( P[i].GravAccel[0]*P[i].GravAccel[0] +  P[i].GravAccel[1]*P[i].GravAccel[1] + P[i].GravAccel[2]*P[i].GravAccel[2]); 
			 		L = sqrt(P[i].Spin[0]*P[i].Spin[0] + P[i].Spin[1]*P[i].Spin[1] + P[i].Spin[2]*P[i].Spin[2]) ; 
			 	
					fprintf(sinklist, "%d		%0.6f		%0.4f		%0.4f		%0.4f		%0.4f		%d		%0.4f		%0.4f		%0.4f\n", All.NumCurrentTiStep, All.Time, sr, sv, sa, P[i].Mass, P[i].NAccreted, P[i].Spin[2], L, TotMassInSinks );
					fclose(sinklist);
				}				
			}

			

//			if(All.TotN_sink > 1){
//				sinkmerger();
//			}

        //This will cause domain decomposition to be performed next time this loop is iterated
        //MPI_Allreduce(&AccNum, &AccNumTot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);		
        if(All.TotN_accrete > 0) All.NumForcesSinceLastDomainDecomp =  All.TotNumPart * All.TreeDomainUpdateFrequency + 1;			 

			 if(ThisTask == 0 ){
					
					TotMassInSinks  = AccNumAll*P[0].Mass; 
					
					printf("ThisTask: %d, NumPart:%d, N_gas:%d, N_sinks:%d, AccNumAll: %d, TotN_gas:%d, TotN_accrete:%d, All.TotN_sink:%d, TotNumPart:%d, TotMassInSinks: %g, AccretionRadius: %g \n", ThisTask, NumPart, N_gas, N_sink, AccNumAll, All.TotN_gas,   All.TotN_accrete , All.TotN_sink,  All.TotN_gas+All.TotN_sink, TotMassInSinks, All.AccretionRadius);				
			 
			 }
       MPI_Barrier(MPI_COMM_WORLD);
			 MPI_Bcast(&TotMassInSinks , 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

					if(TotMassInSinks > 50.0){
						if(ThisTask == 0 ){
							printf("ThisTask : %d, Total mass in sinks has reached 50 , you've defeated the final boss \n", ThisTask);	
							endrun(420);
						}
					}

      }
      
      
#endif




      /* check whether we want a full energy statistics */
      if((All.Time - All.TimeLastStatistics) >= All.TimeBetStatistics)
	{
#ifdef COMPUTE_POTENTIAL_ENERGY
	  compute_potential();
#endif
	  energy_statistics();	/* compute and output energy statistics */
	  All.TimeLastStatistics += All.TimeBetStatistics;
	}

 
  
      advance_and_find_timesteps();	/* 'kick' active particles in
					 * momentum space and compute new
					 * timesteps for them
					 */
	

      All.NumCurrentTiStep++;

      /* Check whether we need to interrupt the run */
      if(ThisTask == 0)
	{
	  /* Is the stop-file present? If yes, interrupt the run. */
	  if((fd = fopen(stopfname, "r")))
	    {
	      fclose(fd);
	      stopflag = 1;
	      unlink(stopfname);
	    }

	  /* are we running out of CPU-time ? If yes, interrupt run. */
	  if(CPUThisRun > 0.85 * All.TimeLimitCPU)
	    {
	      printf("reaching time-limit. stopping.\n");
	      stopflag = 2;
	    }
	}

      MPI_Bcast(&stopflag, 1, MPI_INT, 0, MPI_COMM_WORLD);

      if(stopflag)
	{
	  restart(0);		/* write restart file */
	  MPI_Barrier(MPI_COMM_WORLD);

	  if(stopflag == 2 && ThisTask == 0)
	    {
	      if((fd = fopen(contfname, "w")))
		fclose(fd);
	    }

	  if(stopflag == 2 && All.ResubmitOn && ThisTask == 0)
	    {
	      close_outputfiles();
	      system(All.ResubmitCommand);
	    }
	  return;
	}

      /* is it time to write a regular restart-file? (for security) */
      if(ThisTask == 0)	
	{
	  if((CPUThisRun - All.TimeLastRestartFile) >= All.CpuTimeBetRestartFile)
	    {
	      All.TimeLastRestartFile = CPUThisRun;
	      stopflag = 3;
	    }
	  else
	    stopflag = 0;
	}

      MPI_Bcast(&stopflag, 1, MPI_INT, 0, MPI_COMM_WORLD);

      if(stopflag == 3)
	{
	  restart(0);		/* write an occasional restart file */
	  stopflag = 0;
	}

      t1 = second();

      All.CPU_Total += timediff(t0, t1);
      CPUThisRun += timediff(t0, t1);
    }
  while(All.Ti_Current < TIMEBASE && All.Time <= All.TimeMax);

  restart(0);



 savepositions(All.SnapshotFileCount++);	
 /* write a last snapshot
						 * file at final time (will
						 * be overwritten if
						 * All.TimeMax is increased
						 * and the run is continued)
						 */
}


/*! This function finds the next synchronization point of the system (i.e. the
 *  earliest point of time any of the particles needs a force computation),
 *  and drifts the system to this point of time.  If the system drifts over
 *  the desired time of a snapshot file, the function will drift to this
 *  moment, generate an output, and then resume the drift.
 */
void find_next_sync_point_and_drift(void)
{
  int n, min, min_glob, flag, *temp;
  double timeold;
  double t0, t1;

  t0 = second();

  timeold = All.Time;

  for(n = 1, min = P[0].Ti_endstep; n < NumPart; n++)
    if(min > P[n].Ti_endstep)
      min = P[n].Ti_endstep;

  MPI_Allreduce(&min, &min_glob, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

  /* We check whether this is a full step where all particles are synchronized */
  flag = 1;
  for(n = 0; n < NumPart; n++)
    if(P[n].Ti_endstep > min_glob)
      flag = 0;

  MPI_Allreduce(&flag, &Flag_FullStep, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

#ifdef PMGRID
  if(min_glob >= All.PM_Ti_endstep)
    {
      min_glob = All.PM_Ti_endstep;
      Flag_FullStep = 1;
    }
#endif

  /* Determine 'NumForceUpdate', i.e. the number of particles on this processor that are going to be active */
  for(n = 0, NumForceUpdate = 0; n < NumPart; n++)
    {
      if(P[n].Ti_endstep == min_glob)
#ifdef SELECTIVE_NO_GRAVITY
        if(!((1 << P[n].Type) & (SELECTIVE_NO_GRAVITY)))
#endif
          NumForceUpdate++;
    }

  /* note: NumForcesSinceLastDomainDecomp has type "long long" */
  temp = malloc(NTask * sizeof(int));
  MPI_Allgather(&NumForceUpdate, 1, MPI_INT, temp, 1, MPI_INT, MPI_COMM_WORLD);
  for(n = 0; n < NTask; n++)
    All.NumForcesSinceLastDomainDecomp += temp[n];
  free(temp);



  t1 = second();

  All.CPU_Predict += timediff(t0, t1);

  while(min_glob >= All.Ti_nextoutput && All.Ti_nextoutput >= 0)
    {
      move_particles(All.Ti_Current, All.Ti_nextoutput);

      All.Ti_Current = All.Ti_nextoutput;

      if(All.ComovingIntegrationOn)
	All.Time = All.TimeBegin * exp(All.Ti_Current * All.Timebase_interval);
      else
	All.Time = All.TimeBegin + All.Ti_Current * All.Timebase_interval;

#ifdef OUTPUTPOTENTIAL
      All.NumForcesSinceLastDomainDecomp = 1 + All.TotNumPart * All.TreeDomainUpdateFrequency;
      domain_Decomposition();
      compute_potential();
#endif
      savepositions(All.SnapshotFileCount++);	/* write snapshot file */

      All.Ti_nextoutput = find_next_outputtime(All.Ti_nextoutput + 1);
    }

  move_particles(All.Ti_Current, min_glob);

  All.Ti_Current = min_glob;

  if(All.ComovingIntegrationOn)
    All.Time = All.TimeBegin * exp(All.Ti_Current * All.Timebase_interval);
  else
    All.Time = All.TimeBegin + All.Ti_Current * All.Timebase_interval;

  All.TimeStep = All.Time - timeold;
}



/*! this function returns the next output time that is equal or larger to
 *  ti_curr
 */
int find_next_outputtime(int ti_curr)
{
  int i, ti, ti_next, iter = 0;
  double next, time;

  ti_next = -1;

  if(All.OutputListOn)
    {
      for(i = 0; i < All.OutputListLength; i++)
	{
	  time = All.OutputListTimes[i];

	  if(time >= All.TimeBegin && time <= All.TimeMax)
	    {
	      if(All.ComovingIntegrationOn)
		ti = log(time / All.TimeBegin) / All.Timebase_interval;
	      else
		ti = (time - All.TimeBegin) / All.Timebase_interval;

	      if(ti >= ti_curr)
		{
		  if(ti_next == -1)
		    ti_next = ti;

		  if(ti_next > ti)
		    ti_next = ti;
		}
	    }
	}
    }
  else
    {
      if(All.ComovingIntegrationOn)
	{
	  if(All.TimeBetSnapshot <= 1.0)
	    {
	      printf("TimeBetSnapshot > 1.0 required for your simulation.\n");
	      endrun(13123);
	    }
	}
      else
	{
	  if(All.TimeBetSnapshot <= 0.0)
	    {
	      printf("TimeBetSnapshot > 0.0 required for your simulation.\n");
	      endrun(13123);
	    }
	}

      time = All.TimeOfFirstSnapshot;

      iter = 0;

      while(time < All.TimeBegin)
	{
	  if(All.ComovingIntegrationOn)
	    time *= All.TimeBetSnapshot;
	  else
	    time += All.TimeBetSnapshot;

	  iter++;

	  if(iter > 1000000)
	    {
	      printf("Can't determine next output time.\n");
	      endrun(110);
	    }
	}

      while(time <= All.TimeMax)
	{
	  if(All.ComovingIntegrationOn)
	    ti = log(time / All.TimeBegin) / All.Timebase_interval;
	  else
	    ti = (time - All.TimeBegin) / All.Timebase_interval;

	  if(ti >= ti_curr)
	    {
	      ti_next = ti;
	      break;
	    }

	  if(All.ComovingIntegrationOn)
	    time *= All.TimeBetSnapshot;
	  else
	    time += All.TimeBetSnapshot;

	  iter++;

	  if(iter > 1000000)
	    {
	      printf("Can't determine next output time.\n");
	      endrun(111);
	    }
	}
    }

  if(ti_next == -1)
    {
      ti_next = 2 * TIMEBASE;	/* this will prevent any further output */

      if(ThisTask == 0)
	printf("\nThere is no valid time for a further snapshot file.\n");
    }
  else
    {
      if(All.ComovingIntegrationOn)
	next = All.TimeBegin * exp(ti_next * All.Timebase_interval);
      else
	next = All.TimeBegin + ti_next * All.Timebase_interval;

      if(ThisTask == 0)
	printf("\nSetting next time for snapshot file to Time_next= %g\n\n", next);
    }

  return ti_next;
}




/*! This routine writes one line for every timestep to two log-files.  In
 *  FdInfo, we just list the timesteps that have been done, while in FdCPU the
 *  cumulative cpu-time consumption in various parts of the code is stored.
 */
void every_timestep_stuff(void)
{
  double z;
  double mh = (2.408544e-24), rho_m=0 ;
  int thistask=0, index=0;	

  if(All.NumCurrentTiStep%2 == 0){
    find_max_dens(&thistask, &index);
    if(thistask == ThisTask){
      rho_m = SphP[index].Density;
      rho_n = rho_m*All.UnitDensity_in_cgs/mh;
      smthl = SphP[index].Hsml;
    }

  MPI_Bcast(&rho_n, 1, MPI_DOUBLE, thistask, MPI_COMM_WORLD);
  MPI_Bcast(&rho_m, 1, MPI_DOUBLE, thistask, MPI_COMM_WORLD);
  MPI_Bcast(&smthl , 1, MPI_DOUBLE, thistask, MPI_COMM_WORLD);
  }	


  if(ThisTask == 0)
    {

//	printf("rho_m: %g, unitrgo: %g, mh: %g\n", rho_m,All.UnitDensity_in_cgs,mh);		
      if(All.ComovingIntegrationOn)
	{
	  z = 1.0 / (All.Time) - 1;
	  fprintf(FdInfo, "\nBegin Step %d, Time: %g, Redshift: %g, Systemstep: %g, Dloga: %g, maxdens: %g, hsml: %g, TotNumPart: %d\n",
		  All.NumCurrentTiStep, All.Time, z, All.TimeStep,
		  log(All.Time) - log(All.Time - All.TimeStep), rho_n, smthl, All.TotNumPart);
	  printf("\nBegin Step %d, Time: %g, Redshift: %g, Systemstep: %g, Dloga: %g, maxdens: %g, hsml: %g, TotNumPart: %d\n", All.NumCurrentTiStep,
		 All.Time, z, All.TimeStep, log(All.Time) - log(All.Time - All.TimeStep), rho_n, smthl, All.TotNumPart);
	  fflush(FdInfo);
	}
      else
	{
	  fprintf(FdInfo, "\nBegin Step %d, Time: %g, Systemstep: %g, maxdens: %g, hsml: %g, TotNumPart: %d\n", All.NumCurrentTiStep, All.Time,
		  All.TimeStep, rho_n, smthl, All.TotNumPart);
	  printf("\nBegin Step %d, Time: %g, Systemstep: %g, maxdens: %g, hsml: %g, 	TotNumPart: %d\n", All.NumCurrentTiStep, All.Time, All.TimeStep,  rho_n, smthl,All.TotNumPart);
	  fflush(FdInfo);
	}

      fprintf(FdCPU, "Step %d, Time: %g, CPUs: %d\n", All.NumCurrentTiStep, All.Time, NTask);

      fprintf(FdCPU,
	      "%10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n",
	      All.CPU_Total, All.CPU_Gravity, All.CPU_Hydro, All.CPU_Domain, All.CPU_Potential,
	      All.CPU_Predict, All.CPU_TimeLine, All.CPU_Snapshot, All.CPU_TreeWalk, All.CPU_TreeConstruction,
	      All.CPU_CommSum, All.CPU_Imbalance, All.CPU_HydCompWalk, All.CPU_HydCommSumm,
	      All.CPU_HydImbalance, All.CPU_EnsureNgb, All.CPU_PM, All.CPU_Peano);
      fflush(FdCPU);
    }

  set_random_numbers();
}


/*! This routine first calls a computation of various global quantities of the
 *  particle distribution, and then writes some statistics about the energies
 *  in the various particle components to the file FdEnergy.
 */
void energy_statistics(void)
{
  compute_global_quantities_of_system();

  if(ThisTask == 0)
    {
      fprintf(FdEnergy,
	      "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	      All.Time, SysState.EnergyInt, SysState.EnergyPot, SysState.EnergyKin, SysState.EnergyIntComp[0],
	      SysState.EnergyPotComp[0], SysState.EnergyKinComp[0], SysState.EnergyIntComp[1],
	      SysState.EnergyPotComp[1], SysState.EnergyKinComp[1], SysState.EnergyIntComp[2],
	      SysState.EnergyPotComp[2], SysState.EnergyKinComp[2], SysState.EnergyIntComp[3],
	      SysState.EnergyPotComp[3], SysState.EnergyKinComp[3], SysState.EnergyIntComp[4],
	      SysState.EnergyPotComp[4], SysState.EnergyKinComp[4], SysState.EnergyIntComp[5],
	      SysState.EnergyPotComp[5], SysState.EnergyKinComp[5], SysState.MassComp[0],
	      SysState.MassComp[1], SysState.MassComp[2], SysState.MassComp[3], SysState.MassComp[4],
	      SysState.MassComp[5]);

      fflush(FdEnergy);
    }
}


void find_max_dens(int *thstsk, int *indx){

	double max_dens=0;
	int tempindex=0;                    
	
	if(ThisTask == 0){
		double max_glob=0;
		*thstsk =0;
		*indx =0;
	
		for(int i = 0; i < N_gas; i++){		 		
			if(SphP[i].Density > max_glob){		
					#ifdef SINK
						if(SphP[i].DivVel <= 0 && P[i].Ti_endstep == All.Ti_Current){
							max_glob = SphP[i].Density;	
							*indx  = i;					
						}
					#else
						max_glob = SphP[i].Density;	
						*indx  = i;				
					#endif					 
			}
		}
//		printf("Ngas = %d, Local maxima  = %.5f by %d for particle id %d  at %d\n", N_gas, max_glob, ThisTask, P[*indx].ID ,*indx);							

		for(int rank=1; rank<NTask; rank++){
			MPI_Recv(&max_dens   , 1, MPI_DOUBLE, rank, 101, MPI_COMM_WORLD, MPI_STATUS_IGNORE);		
			MPI_Recv(&tempindex  , 1, MPI_DOUBLE, rank, 101, MPI_COMM_WORLD, MPI_STATUS_IGNORE);		

//			printf("i recieved max_dens = %.5f from %d at %d\n", max_dens, rank, tempindex);
			if( max_glob < max_dens ){
				max_glob   = max_dens;
				*thstsk    = rank;
				*indx      = tempindex;
			}
		}
//		printf("     ##########################   All processess done, Global maxima is %0.5f by %d at %d   ###########################   \n",max_glob, *thstsk, *indx);

	} //root 	
	else {	
		for(int i = 0; i < N_gas; i++){
			if(SphP[i].Density > max_dens){		
					#ifdef SINK
						if(SphP[i].DivVel <= 0  && P[i].Ti_endstep == All.Ti_Current ){
							max_dens = SphP[i].Density;	
							tempindex       = i;
						}
					#else
						max_dens = SphP[i].Density;	
						tempindex       = i;
					#endif					 
			}
		}
		MPI_Send(&max_dens   , 1, MPI_DOUBLE, 0, 101, MPI_COMM_WORLD);
		MPI_Send(&tempindex  , 1, MPI_INT   , 0, 101, MPI_COMM_WORLD);
	}
	MPI_Barrier(MPI_COMM_WORLD); 	
		
	if(ThisTask == 0){	
		for(int rank=1; rank<NTask; rank++){
			MPI_Send(thstsk, 1, MPI_INT, rank, 111, MPI_COMM_WORLD);
    	MPI_Send(indx  , 1, MPI_INT, rank, 111, MPI_COMM_WORLD);
		}
	}
	else{
		MPI_Recv(thstsk, 1, MPI_INT, 0, 111, MPI_COMM_WORLD, MPI_STATUS_IGNORE);		
		MPI_Recv(indx  , 1, MPI_INT, 0, 111, MPI_COMM_WORLD, MPI_STATUS_IGNORE);				
	}
	MPI_Barrier(MPI_COMM_WORLD); 
}


//		
//		MPI_Bcast(vel  , 3, MPI_INT, thistask, MPI_COMM_WORLD);
//		relvel = (vel[0]-P[i].Vel[0])*(vel[0]-P[i].Vel[0]) + (vel[1]-P[i].Vel[1])*(vel[1]-P[i].Vel[1]) + (vel[2]-P[i].Vel[2])*(vel[2]-P[i].Vel[2]);
//			printf("  Particle id: %d at %d in prcoesssor: %d tagged as sink at pos:(%f, %f, %f), vel:(%f, %f, %f) \n",P[index].ID, index, thistask, P[index].Pos[0], P[index].Pos[1], P[index].Pos[2], P[index].Vel[0], P[index].Vel[1], P[index].Vel[2]  );

#ifdef SINK
void create_sink(){

	int thistask, index;
	int numsinks, numsinkstot;

  FLOAT *local_sink_posx, *local_sink_posy, *local_sink_posz; 
  int *local_sink_ID;  
  
  FLOAT *list_sink_posx, *list_sink_posy, *list_sink_posz;
  int *list_sink_ID;  
		
  numsinks = NumPart - N_gas;
  MPI_Allreduce(&numsinks, &numsinkstot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); 

  local_sink_posx = malloc(sizeof(double) * numsinkstot);
  local_sink_posy = malloc(sizeof(double) * numsinkstot);  
  local_sink_posz = malloc(sizeof(double) * numsinkstot); 
  local_sink_ID   = malloc(sizeof(int) * numsinkstot);   

  list_sink_posx = malloc(sizeof(double) * numsinkstot * NTask);
  list_sink_posy = malloc(sizeof(double) * numsinkstot * NTask);  
  list_sink_posz = malloc(sizeof(double) * numsinkstot * NTask); 
  list_sink_ID   = malloc(sizeof(int) * numsinkstot * NTask);    

  for(int i = 0; i < numsinkstot; i++) {
		local_sink_ID[i]   = -1;
		local_sink_posx[i] = 0;
		local_sink_posy[i] = 0;
		local_sink_posz[i] = 0; 
  }  
  for(int i = 0; i < numsinks; i++){
    local_sink_posx[i] = P[i+N_gas].Pos[0];
    local_sink_posy[i] = P[i+N_gas].Pos[1];    
    local_sink_posz[i] = P[i+N_gas].Pos[2];
    local_sink_ID[i]   = P[i+N_gas].ID;    
   // printf("ThisTask: %d, sinknum: %d, sinkid: %d\n\n", ThisTask, i+1, local_sink_ID[i]);    
  }  

  MPI_Barrier(MPI_COMM_WORLD);   
  MPI_Allgather(local_sink_posx, numsinkstot, MPI_DOUBLE, list_sink_posx, numsinkstot, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Allgather(local_sink_posy, numsinkstot, MPI_DOUBLE, list_sink_posy, numsinkstot, MPI_DOUBLE, MPI_COMM_WORLD);  
  MPI_Allgather(local_sink_posz, numsinkstot, MPI_DOUBLE, list_sink_posz, numsinkstot, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Allgather(local_sink_ID  , numsinkstot, MPI_INT   , list_sink_ID  , numsinkstot, MPI_INT   , MPI_COMM_WORLD);           
  MPI_Barrier(MPI_COMM_WORLD); 

	// collected positions of all the sink particles

//	FLOAT hsml; 

	find_max_dens(&thistask, &index);	

	int   issink, sinkid, i_ngb, num, startnode;
  FLOAT sp_pos[3], sp_vel[3], t_pos[3], t_vel[3]; 
	FLOAT distsq, dist,  twohsml, hsml, u, Temp;

	struct particle_data      temp_p;
	struct sph_particle_data  temp_sph;

	FLOAT mh =  2.408544e-24;
	FLOAT CriticalMassDensity = All.CriticalNumberDensity*mh / All.UnitDensity_in_cgs;

	AccNum=0;
	 for(int itag=0; itag<N_gas; itag++){
		 SphP[itag].AccretionTarget =0;
	 }  


	if(ThisTask == thistask){	

		u = SphP[index].Entropy / ETA_MINUS1 * pow(SphP[index].Density, GAMMA_MINUS1);
		
		Temp = u*TempFac;  
		
		printf("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Temp = %g \n", Temp);
		
		if(SphP[index].Density >= CriticalMassDensity && Temp >= All.CriticalTemperature){	
			
		//		hsml = SphP[index].Hsml;
		//	All.AccretionRadius = 4.0*hsml;			//comment this to set accretion radius manually	
			FLOAT tworaccsq = (2.0*All.AccretionRadius)*(2.0*All.AccretionRadius);
		
			sinkid = P[index].ID;
			hsml  = SphP[index].Hsml;	
			issink = 1;
			for(int i=0; i<3; i++){
				sp_pos[i] = P[index].Pos[i];
				sp_vel[i] = P[index].Vel[i];
			}

		  for(int i = 0; i < numsinkstot * NTask; i++){	 /* go through all the sink particles (From all processors) */
  	  	if(list_sink_ID[i] >= 0 && list_sink_ID[i] !=  sinkid ){
					distsq = ( list_sink_posx[i] - sp_pos[0] ) * ( list_sink_posx[i] - sp_pos[0] ) + 	
									 ( list_sink_posy[i] - sp_pos[1] ) * ( list_sink_posy[i] - sp_pos[1] ) +
									 ( list_sink_posz[i] - sp_pos[2] ) * ( list_sink_posz[i] - sp_pos[2] ) ;	

					if( distsq <= tworaccsq  ){
						printf("ThisTask; %d wanted to form sink but there is another sink within %g  at dist=%g \n", ThisTask, sqrt(tworaccsq), sqrt(distsq));
						issink = 0;
						break;
					}
				}
			}
		} //crit dens

		if(issink){

			P[index].Type = 4;

			temp_p   =  P[index];
			temp_sph =  SphP[index];

	    for(int k = index; k<=NumPart-2- numsinks ; k++){ 
  	  	P[k]    = P[k+1];
  	   	SphP[k] = SphP[k+1]; 
  	  } 
 	    P[NumPart -1 - numsinks]   = temp_p;
  	  SphP[NumPart-1 - numsinks] = temp_sph;      

			printf("ThisTask = %d,  moving sink from %d to %d, sink ID: %d, numsinklocal: %d, numsinkstot: %d, N_gas: %d, NumPart: %d\n", ThisTask, index, NumPart -1 - numsinks, P[NumPart -1 - numsinks].ID, numsinks, numsinkstot, N_gas, NumPart);
  	  N_gas--;

		}
	} //thistask
	MPI_Bcast(&issink , 1, MPI_INT	 , thistask, MPI_COMM_WORLD);

	if(issink){
		MPI_Bcast(&sinkid , 1, MPI_INT	 , thistask, MPI_COMM_WORLD);
		MPI_Bcast(&sp_pos , 1, MPI_DOUBLE, thistask, MPI_COMM_WORLD);
		MPI_Bcast(&sp_vel , 1, MPI_DOUBLE, thistask, MPI_COMM_WORLD);
		MPI_Bcast(&hsml   , 1, MPI_DOUBLE, thistask, MPI_COMM_WORLD);


//		twohsml = 0.00004;
		twohsml = 2.0*hsml;
		startnode = All.MaxPart;
		do{
			num = ngb_treefind_variable(sp_pos, twohsml, &startnode);
			for(int n = 0; n < num; n++){
				i_ngb = Ngblist[n];
				if(P[i_ngb].Type == 0  && P[i_ngb].Ti_endstep == All.Ti_Current ){  
					dist = 0;
				  for(int j = 0; j < 3; j++){ 
				  	dist += (P[i_ngb].Pos[j]-sp_pos[j]) * (P[i_ngb].Pos[j]-sp_pos[j]);
				  }   
				  dist = sqrt(dist);     
				  if(dist <= twohsml){
						SphP[i_ngb].AccretionTarget = sinkid;
					}
				}
			} 
 		}while(startnode>=0);


		FLOAT dposx=0, dposy=0, dposz=0, dvelx=0, dvely=0, dvelz=0, dmass=0, spinx=0, spiny=0, spinz=0;
		FLOAT dposxtot=0, dposytot=0, dposztot=0, dvelxtot=0, dvelytot=0, dvelztot=0, dmasstot=0, spinxtot=0, spinytot=0, spinztot=0;      

		for(int i=0; i<N_gas; i++){
			if(SphP[i].AccretionTarget == sinkid && P[i].Type == 0 ){
				AccreteList[AccNum] = i;
				AccNum++;	

        dvelx += P[i].Mass*P[i].Vel[0];
        dvely += P[i].Mass*P[i].Vel[1];	    
        dvelz += P[i].Mass*P[i].Vel[2];	
        dposx += P[i].Mass*P[i].Pos[0];
        dposy += P[i].Mass*P[i].Pos[1];	    
        dposz += P[i].Mass*P[i].Pos[2];   
        dmass += P[i].Mass;  

				t_pos[0] = P[i].Pos[0] - sp_pos[0];
				t_pos[1] = P[i].Pos[1] - sp_pos[1];
				t_pos[2] = P[i].Pos[2] - sp_pos[2];

				t_vel[0] = P[i].Vel[0] - sp_vel[0];
				t_vel[1] = P[i].Vel[1] - sp_vel[1];
				t_vel[2] = P[i].Vel[2] - sp_vel[2];

				spinx += t_pos[1]*t_vel[2] -  t_pos[2]*t_vel[1];   
				spiny += t_pos[2]*t_vel[0] -  t_pos[0]*t_vel[2];   
				spinz += t_pos[0]*t_vel[1] -  t_pos[1]*t_vel[0];   
								
			}
		}
		   
    MPI_Barrier(MPI_COMM_WORLD);       
    MPI_Allreduce(&dvelx, &dvelxtot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
    MPI_Allreduce(&dvely, &dvelytot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);  
    MPI_Allreduce(&dvelz, &dvelztot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
    MPI_Allreduce(&dposx, &dposxtot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
    MPI_Allreduce(&dposy, &dposytot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);  
    MPI_Allreduce(&dposz, &dposztot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);     
    MPI_Allreduce(&dmass, &dmasstot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);        
    MPI_Allreduce(&spinx, &spinxtot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
    MPI_Allreduce(&spiny, &spinytot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
    MPI_Allreduce(&spinz, &spinztot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
    MPI_Barrier(MPI_COMM_WORLD);     
	
		if(ThisTask == thistask){
			dposxtot += P[index].Pos[0] * P[index].Mass;
			dposytot += P[index].Pos[1] * P[index].Mass;	
			dposztot += P[index].Pos[2] * P[index].Mass;
			dmasstot += P[index].Mass;

			//Move position to centre of mass
			P[index].Pos[0] = dposxtot / dmasstot;
			P[index].Pos[1] = dposytot / dmasstot;	  	  
			P[index].Pos[2] = dposztot / dmasstot;	  

			//Should this be predicted velocity?
			dvelxtot += P[index].Mass * P[index].Vel[0] ;	  
			dvelytot += P[index].Mass * P[index].Vel[1] ;	  
			dvelztot += P[index].Mass * P[index].Vel[2] ;

			//Add momentum to the sink
			P[index].Vel[0] = dvelxtot / dmasstot ;
			P[index].Vel[1] = dvelytot / dmasstot ;
			P[index].Vel[2] = dvelztot / dmasstot ;	  	  

			//Add the mass to the sink
			P[index].Mass = dmasstot;

			P[index].Spin[0] += spinxtot;
			P[index].Spin[1] += spinytot;
			P[index].Spin[2] += spinztot;

			P[index].Ti_endstep = All.Ti_Current;		
		}
		MPI_Barrier(MPI_COMM_WORLD); 
	
		if(AccNum > 1) qsort(AccreteList, AccNum, sizeof(int), index_compare_key);
  	int acc_counter = 0;
  	for(int n = 0;n < AccNum;n++){   
    	i_ngb = AccreteList[n] - acc_counter;
    	if(SphP[i_ngb].AccretionTarget > -2){      
    	  if(P[i_ngb].Ti_endstep == All.Ti_Current){
    	    NumForceUpdate--;
    	    NumSphUpdate--;
    	  }
 
 	  	  for(int k = i_ngb; k<=NumPart-2 ; k++){ 
	 		  	P[k]    = P[k+1];
  			  SphP[k] = SphP[k+1]; 
  		  }     
 
      	NumPart--;   // decrement the local countrs of particles and gas particles
      	N_gas--;
	      acc_counter++; 
    	}
  	}
	
		N_sink = NumPart - N_gas;
		MPI_Barrier(MPI_COMM_WORLD); 
		MPI_Allreduce(&N_gas	   	, &All.TotN_gas		, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&NumPart, &All.TotNumPart	 , 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&N_sink , &All.TotN_sink	 , 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&AccNum , &All.TotN_accrete, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD); 


		header.npart[0] 		= All.TotN_gas;
		header.npart[4] 		= All.TotN_sink;
  	AccNumAll += All.TotN_accrete;

		if(ThisTask == thistask){
			printf("ThisTask: %d, Created sink with %d neigbours within the radius %g\n", ThisTask, All.TotN_accrete, twohsml);
		}

	}	 //issink
	AccNum=0;
  free(list_sink_posx);
  free(list_sink_posy);
  free(list_sink_posz);  
  free(list_sink_ID);
  free(local_sink_posx);
  free(local_sink_posy);
  free(local_sink_posz);  
  free(local_sink_ID);      
}



void sinkmerger(){

  int numsinks, numsinkstot;
	
  FLOAT *local_sink_posx, *local_sink_posy, *local_sink_posz;
  FLOAT *local_sink_velx, *local_sink_vely, *local_sink_velz;
  FLOAT *local_sink_mass;  
  int *local_sink_ID;  
 
  FLOAT *list_sink_posx, *list_sink_posy, *list_sink_posz;
  FLOAT *list_sink_velx, *list_sink_vely, *list_sink_velz;
  FLOAT *list_sink_mass;  
  int *list_sink_ID;  

  numsinks = NumPart - N_gas; 
  MPI_Allreduce(&numsinks, &numsinkstot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); 
  MPI_Barrier(MPI_COMM_WORLD);

  local_sink_posx = malloc(sizeof(double) * numsinkstot);
  local_sink_posy = malloc(sizeof(double) * numsinkstot);  
  local_sink_posz = malloc(sizeof(double) * numsinkstot); 
  local_sink_velx = malloc(sizeof(double) * numsinkstot);
  local_sink_vely = malloc(sizeof(double) * numsinkstot);  
  local_sink_velz = malloc(sizeof(double) * numsinkstot);  
  local_sink_ID = malloc(sizeof(int) * numsinkstot);   
  local_sink_mass = malloc(sizeof(double) * numsinkstot);     
  list_sink_posx = malloc(sizeof(double) * numsinkstot * NTask);
  list_sink_posy = malloc(sizeof(double) * numsinkstot * NTask);  
  list_sink_posz = malloc(sizeof(double) * numsinkstot * NTask); 
  list_sink_velx = malloc(sizeof(double) * numsinkstot * NTask);
  list_sink_vely = malloc(sizeof(double) * numsinkstot * NTask);  
  list_sink_velz = malloc(sizeof(double) * numsinkstot * NTask); 
  list_sink_ID = malloc(sizeof(int) * numsinkstot * NTask);    
  list_sink_mass = malloc(sizeof(double) * numsinkstot * NTask);  
  
	for(int is = N_gas; is < NumPart; is++){ 
 		P[is].AccreteOntoSInk = -1;
 	}
  
  for(int i = 0; i < numsinkstot; i++) local_sink_mass[i] = -1;
  
  for(int i = 0; i < numsinks; i++){
    local_sink_posx[i] = P[i+N_gas].Pos[0];
    local_sink_posy[i] = P[i+N_gas].Pos[1];    
    local_sink_posz[i] = P[i+N_gas].Pos[2];
    local_sink_velx[i] = P[i+N_gas].Vel[0];
    local_sink_vely[i] = P[i+N_gas].Vel[1];    
    local_sink_velz[i] = P[i+N_gas].Vel[2];
    local_sink_ID[i] = P[i+N_gas].ID;    
    local_sink_mass[i] = P[i+N_gas].Mass;     
  }  
  
  MPI_Barrier(MPI_COMM_WORLD);   
  MPI_Allgather(local_sink_posx, numsinkstot, MPI_DOUBLE, list_sink_posx, numsinkstot, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Allgather(local_sink_posy, numsinkstot, MPI_DOUBLE, list_sink_posy, numsinkstot, MPI_DOUBLE, MPI_COMM_WORLD);  
  MPI_Allgather(local_sink_posz, numsinkstot, MPI_DOUBLE, list_sink_posz, numsinkstot, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Allgather(local_sink_velx, numsinkstot, MPI_DOUBLE, list_sink_velx, numsinkstot, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Allgather(local_sink_vely, numsinkstot, MPI_DOUBLE, list_sink_vely, numsinkstot, MPI_DOUBLE, MPI_COMM_WORLD);  
  MPI_Allgather(local_sink_velz, numsinkstot, MPI_DOUBLE, list_sink_velz, numsinkstot, MPI_DOUBLE, MPI_COMM_WORLD);  
  MPI_Allgather(local_sink_mass, numsinkstot, MPI_DOUBLE, list_sink_mass, numsinkstot, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Allgather(local_sink_ID, numsinkstot, MPI_INT, list_sink_ID, numsinkstot, MPI_INT, MPI_COMM_WORLD);           
  MPI_Barrier(MPI_COMM_WORLD); 

	FLOAT x_js, y_js, z_js, vx_js, vy_js, vz_js, m_js, ID_js ; 
	FLOAT x_is, y_is, z_is, vx_is, vy_is, vz_is, m_is, ID_is; 
  FLOAT seperation, relvel, relenergy, little_L, KeplerL2, sinkrad;
	int notestflag=0;
	int Nacc=0, NaccTot=0, AccList[100]={0};

	sinkrad = All.AccretionRadius;
  for(int js = 0; js < numsinkstot * NTask; js++){	 
    if(list_sink_mass[js] > 0){
    	x_js  = list_sink_posx[js];
    	y_js  = list_sink_posy[js];
    	z_js  = list_sink_posz[js];
    	vx_js = list_sink_velx[js];
    	vy_js = list_sink_vely[js];
    	vz_js = list_sink_velz[js];
    	m_js  = list_sink_mass[js]; 
			ID_js = list_sink_ID[js];    	
    	
      KeplerL2 = All.G * m_js * sinkrad;
		  for(int is = N_gas; is < NumPart; is++){	 
    		if(P[is].ID != ID_js &&  P[is].Ti_endstep == All.Ti_Current ){
					x_is  = P[is].Pos[0];
					y_is  = P[is].Pos[1];
					z_is  = P[is].Pos[2];
					vx_is = P[is].Vel[0];
					vy_is = P[is].Vel[1];
					vz_is = P[is].Vel[2];
					m_is  = P[is].Mass; 
					ID_is = P[is].ID;    	

    			seperation = (x_js - x_is)*(x_js - x_is) + (y_js - y_is)*(y_js - y_is) + (z_js - z_is)*(z_js - z_is) ;      			
 					seperation = sqrt(seperation);
 					
    			if(seperation < sinkrad/2.0){
    				relvel = (vx_js - vx_is)*(vx_js - vx_is) + (vy_js - vy_is)*(vy_js - vy_is) + (vz_js - vz_is)*(vz_js - vz_is) ; 
    				relenergy = 0.5*relvel - All.G*m_js/seperation;	

            if(notestflag) {relenergy = -1;}  
            if(relenergy < 0){                          
              little_L  = pow( (y_js-y_is)*(vz_js - vz_is) - (z_js - z_is)*(vy_js - vy_is) , 2  )  ;
              little_L += pow( (z_js-z_is)*(vx_js - vx_is) - (x_js - x_is)*(vz_js - vz_is) , 2  )  ;
              little_L += pow( (x_js-x_is)*(vy_js - vy_is) - (y_js - y_is)*(vx_js - vx_is) , 2  )  ;
            	if(notestflag) little_L = KeplerL2 /2.0;								  						
              if(little_L < KeplerL2){  
              		if( m_is < m_js && P[is].AccreteOntoSInk == -1 ){   
              			P[is].AccreteOntoSInk = ID_js; 
              	//		AccList[Nacc] = is;  
              			Nacc++;
              			printf(" ##################################################################### Bound Merger event .....\n");
              		}
              	            		                
              }           
  					} //relenergy
    			} //dist	
    			else if (seperation < sinkrad/4.0){
		    		if( m_is < m_js && P[is].AccreteOntoSInk == -1 ){   
		    			P[is].AccreteOntoSInk = ID_js; 
		    	//		AccList[Nacc] = is;  
		    			Nacc++;
		    			printf(" ##################################################################### Dist Merger event .....\n");
		    		}    			
    			}
				} //m_is>0
			} //is
		} //m_js >0 
	} //js
		
	MPI_Allreduce(&Nacc	   	, &NaccTot		, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	if(NaccTot > 0 ){
	Nacc=0;	
  for(int is = N_gas; is < NumPart; is++){
    if(P[is].AccreteOntoSInk > -1){
      AccList[Nacc] = is;
      Nacc++;
    }
  } 
	printf("Nacc = %d, NaccTot = %d \n",Nacc, NaccTot);
	
		FLOAT dposx=0, dposy=0, dposz=0, dvelx=0, dvely=0, dvelz=0, dmass=0, spinx=0, spiny=0, spinz=0;
		FLOAT dposxtot=0, dposytot=0, dposztot=0, dvelxtot=0, dvelytot=0, dvelztot=0, dmasstot=0, spinxtot=0, spinytot=0, spinztot=0;   
		
		   
		int target;
		for(int js = 0; js < numsinkstot * NTask; js++){	 
		  if(list_sink_mass[js] > 0){
				target = 	list_sink_ID[js];	
				for(int is = N_gas; is < NumPart; is++){	 
					if( P[is].AccreteOntoSInk == target ){
						dvelx += P[is].Mass*P[is].Vel[0];
          	dvely += P[is].Mass*P[is].Vel[1];	    
          	dvelz += P[is].Mass*P[is].Vel[2];	
          	dposx += P[is].Mass*P[is].Pos[0];
          	dposy += P[is].Mass*P[is].Pos[1];	    
          	dposz += P[is].Mass*P[is].Pos[2]; 
	
						spinx += P[is].Spin[0];   
						spiny += P[is].Spin[1];   
						spinz += P[is].Spin[2];   
										
	          dmass += P[is].Mass;           	
          	  
					}
				}

				MPI_Barrier(MPI_COMM_WORLD);       
				MPI_Allreduce(&dvelx, &dvelxtot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
				MPI_Allreduce(&dvely, &dvelytot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);  
				MPI_Allreduce(&dvelz, &dvelztot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
				MPI_Allreduce(&dposx, &dposxtot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
				MPI_Allreduce(&dposy, &dposytot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);  
				MPI_Allreduce(&dposz, &dposztot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);     
				MPI_Allreduce(&dmass, &dmasstot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);        
				MPI_Allreduce(&spinx, &spinxtot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
				MPI_Allreduce(&spiny, &spinytot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
				MPI_Allreduce(&spinz, &spinztot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
				MPI_Barrier(MPI_COMM_WORLD);     

				for(int j = N_gas; j < NumPart; j++){
					if(P[j].ID == target){
    	      dposxtot += P[j].Pos[0] * P[j].Mass;
    	      dposytot += P[j].Pos[1] * P[j].Mass;	
    	      dposztot += P[j].Pos[2] * P[j].Mass;

  	        dvelxtot += P[j].Mass * ( P[j].Vel[0] );	  
  	        dvelytot += P[j].Mass * ( P[j].Vel[1] );	  
  	        dvelztot += P[j].Mass * ( P[j].Vel[2] );

    	      dmasstot += P[j].Mass;

         
	          //Move position to centre of mass
  	        P[j].Pos[0] = dposxtot / dmasstot;
  	        P[j].Pos[1] = dposytot / dmasstot;	  	  
  	        P[j].Pos[2] = dposztot / dmasstot;	  
          
          
	          //Add momentum to the sink
  	        P[j].Vel[0] = dvelxtot / dmasstot ;
  	        P[j].Vel[1] = dvelytot / dmasstot ;
  	        P[j].Vel[2] = dvelztot / dmasstot ;	  	  
          
	          //Add the mass to the sink
  	        P[j].Mass = dmasstot;

  	        P[j].Spin[0] += spinxtot;
  	        P[j].Spin[1] += spinytot;
  	        P[j].Spin[2] += spinztot;

            P[j].Ti_endstep = All.Ti_Current;
        	}   
      	}
  			 
			} //mass_js
		}			



//		printf("I have made it till here 0\n");
		MPI_Barrier(MPI_COMM_WORLD);
		if(Nacc > 1) qsort(AccList, Nacc, sizeof(int), index_compare_key);

		int acc_counter = 0, indx;
		for(int n = 0;n < Nacc ;n++){
  		indx = AccList[n] - acc_counter;
  		if(P[indx].Ti_endstep == All.Ti_Current){
  			NumForceUpdate--;
    		NumSphUpdate--;
  		}

  		for(int k = indx; k<=NumPart-2 ; k++){ 
  			P[k]    = P[k+1];
   			SphP[k] = SphP[k+1]; 
  		} 

  		NumPart--;   // decrement the local countrs of particles and gas particles
	    acc_counter++;     
		}

		N_sink = NumPart-N_gas;
		MPI_Barrier(MPI_COMM_WORLD); 
		MPI_Allreduce(&NumPart, &All.TotNumPart	 , 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&N_sink , &All.TotN_sink	 , 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		header.npart[4] 		= All.TotN_sink;
		All.NumForcesSinceLastDomainDecomp =  All.TotNumPart * All.TreeDomainUpdateFrequency + 1;
	} //NaccTot

  free(list_sink_posx);
  free(list_sink_posy);
  free(list_sink_posz);  
  free(list_sink_velx);
  free(list_sink_vely);
  free(list_sink_velz);    
  free(list_sink_ID);
  free(list_sink_mass);  
  free(local_sink_posx);
  free(local_sink_posy);
  free(local_sink_posz);  
  free(local_sink_velx);
  free(local_sink_vely);
  free(local_sink_velz); 
  free(local_sink_ID);      
  free(local_sink_mass);      
}


#endif 














































