#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include "allvars.h"
#include "proto.h"


/*! \file density.c 
 *  \brief SPH density computation and smoothing length determination
 *
 *  This file contains the "first SPH loop", where the SPH densities and
 *  some auxiliary quantities are computed.  If the number of neighbours
 *  obtained falls outside the target range, the correct smoothing
 *  length is determined iteratively, if needed.
 */


#ifdef PERIODIC
static double boxSize, boxHalf;

#ifdef LONG_X
static double boxSize_X, boxHalf_X;
#else
#define boxSize_X boxSize
#define boxHalf_X boxHalf
#endif
#ifdef LONG_Y
static double boxSize_Y, boxHalf_Y;
#else
#define boxSize_Y boxSize
#define boxHalf_Y boxHalf
#endif
#ifdef LONG_Z
static double boxSize_Z, boxHalf_Z;
#else
#define boxSize_Z boxSize
#define boxHalf_Z boxHalf
#endif
#endif


/*! This function computes the local density for each active SPH particle,
 *  the number of neighbours in the current smoothing radius, and the
 *  divergence and curl of the velocity field.  The pressure is updated as
 *  well.  If a particle with its smoothing region is fully inside the
 *  local domain, it is not exported to the other processors. The function
 *  also detects particles that have a number of neighbours outside the
 *  allowed tolerance range. For these particles, the smoothing length is
 *  adjusted accordingly, and the density computation is executed again.
 *  Note that the smoothing length is not allowed to fall below the lower
 *  bound set by MinGasHsml.
 */
void density(void)
{
  long long ntot, ntotleft;
  int *noffset, *nbuffer, *nsend, *nsend_local, *numlist, *ndonelist;
  int i, j, n, ndone, npleft, maxfill, source, iter = 0;
  int level, ngrp, sendTask, recvTask, place, nexport;
  double dt_entr, tstart, tend, tstart_ngb = 0, tend_ngb = 0;
  double sumt, sumcomm, timengb, sumtimengb;
  double timecomp = 0, timeimbalance = 0, timecommsumm = 0, sumimbalance;
  MPI_Status status;

  //double dhsml_bef, divv_bef, curlv_bef, numngb_bef, rho_bef, hsml_bef; 
  //float dhsml_bef_f, divv_bef_f, curlv_bef_f, numngb_bef_f, rho_bef_f, hsml_bef_f; 

#ifdef PERIODIC
  boxSize = All.BoxSize;
  boxHalf = 0.5 * All.BoxSize;
#ifdef LONG_X
  boxHalf_X = boxHalf * LONG_X;
  boxSize_X = boxSize * LONG_X;
#endif
#ifdef LONG_Y
  boxHalf_Y = boxHalf * LONG_Y;
  boxSize_Y = boxSize * LONG_Y;
#endif
#ifdef LONG_Z
  boxHalf_Z = boxHalf * LONG_Z;
  boxSize_Z = boxSize * LONG_Z;
#endif
#endif


  noffset = malloc(sizeof(int) * NTask);	/* offsets of bunches in common list */
  nbuffer = malloc(sizeof(int) * NTask);
  nsend_local = malloc(sizeof(int) * NTask);
  nsend = malloc(sizeof(int) * NTask * NTask);
  ndonelist = malloc(sizeof(int) * NTask);

  for(n = 0, NumSphUpdate = 0; n < N_gas; n++)
    {
      SphP[n].Left = SphP[n].Right = 0;

      if(P[n].Ti_endstep == All.Ti_Current)
	NumSphUpdate++;
    }

  numlist = malloc(NTask * sizeof(int) * NTask);
  MPI_Allgather(&NumSphUpdate, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
  for(i = 0, ntot = 0; i < NTask; i++)
    ntot += numlist[i];
  free(numlist);



  /* we will repeat the whole thing for those particles where we didn't
   * find enough neighbours
   */
  do
    {
      i = 0;			/* beginn with this index */
      ntotleft = ntot;		/* particles left for all tasks together */

      while(ntotleft > 0)
	{
	  for(j = 0; j < NTask; j++)
	    nsend_local[j] = 0;

	  /* do local particles and prepare export list */
	  tstart = second();
	  for(nexport = 0, ndone = 0; i < N_gas && nexport < All.BunchSizeDensity - NTask; i++)
	    if(P[i].Ti_endstep == All.Ti_Current)
	      {
		ndone++;

		for(j = 0; j < NTask; j++)
		  Exportflag[j] = 0;

		density_evaluate(i, 0);

		for(j = 0; j < NTask; j++)
		  {
		    if(Exportflag[j])
		      {
			DensDataIn[nexport].Pos[0] = P[i].Pos[0];
			DensDataIn[nexport].Pos[1] = P[i].Pos[1];
			DensDataIn[nexport].Pos[2] = P[i].Pos[2];
			DensDataIn[nexport].Vel[0] = SphP[i].VelPred[0];
			DensDataIn[nexport].Vel[1] = SphP[i].VelPred[1];
			DensDataIn[nexport].Vel[2] = SphP[i].VelPred[2];
			DensDataIn[nexport].Hsml = SphP[i].Hsml;
			DensDataIn[nexport].Index = i;
			DensDataIn[nexport].Task = j;
			nexport++;
			nsend_local[j]++;
		      }
		  }
	      }
	  tend = second();
	  timecomp += timediff(tstart, tend);

	  qsort(DensDataIn, nexport, sizeof(struct densdata_in), dens_compare_key);

	  for(j = 1, noffset[0] = 0; j < NTask; j++)
	    noffset[j] = noffset[j - 1] + nsend_local[j - 1];

	  tstart = second();

	  MPI_Allgather(nsend_local, NTask, MPI_INT, nsend, NTask, MPI_INT, MPI_COMM_WORLD);

	  tend = second();
	  timeimbalance += timediff(tstart, tend);


	  /* now do the particles that need to be exported */

	  for(level = 1; level < (1 << PTask); level++)
	    {
	      tstart = second();
	      for(j = 0; j < NTask; j++)
		nbuffer[j] = 0;
	      for(ngrp = level; ngrp < (1 << PTask); ngrp++)
		{
		  maxfill = 0;
		  for(j = 0; j < NTask; j++)
		    {
		      if((j ^ ngrp) < NTask)
			if(maxfill < nbuffer[j] + nsend[(j ^ ngrp) * NTask + j])
			  maxfill = nbuffer[j] + nsend[(j ^ ngrp) * NTask + j];
		    }
		  if(maxfill >= All.BunchSizeDensity)
		    break;

		  sendTask = ThisTask;
		  recvTask = ThisTask ^ ngrp;

		  if(recvTask < NTask)
		    {
		      if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
			{
			  /* get the particles */
			  MPI_Sendrecv(&DensDataIn[noffset[recvTask]],
				       nsend_local[recvTask] * sizeof(struct densdata_in), MPI_BYTE,
				       recvTask, TAG_DENS_A,
				       &DensDataGet[nbuffer[ThisTask]],
				       nsend[recvTask * NTask + ThisTask] * sizeof(struct densdata_in),
				       MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, &status);
			}
		    }

		  for(j = 0; j < NTask; j++)
		    if((j ^ ngrp) < NTask)
		      nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
		}
	      tend = second();
	      timecommsumm += timediff(tstart, tend);


	      tstart = second();
	      for(j = 0; j < nbuffer[ThisTask]; j++)
		density_evaluate(j, 1);
	      tend = second();
	      timecomp += timediff(tstart, tend);

	      /* do a block to explicitly measure imbalance */
	      tstart = second();
	      MPI_Barrier(MPI_COMM_WORLD);
	      tend = second();
	      timeimbalance += timediff(tstart, tend);

	      /* get the result */
	      tstart = second();
	      for(j = 0; j < NTask; j++)
		nbuffer[j] = 0;
	      for(ngrp = level; ngrp < (1 << PTask); ngrp++)
		{
		  maxfill = 0;
		  for(j = 0; j < NTask; j++)
		    {
		      if((j ^ ngrp) < NTask)
			if(maxfill < nbuffer[j] + nsend[(j ^ ngrp) * NTask + j])
			  maxfill = nbuffer[j] + nsend[(j ^ ngrp) * NTask + j];
		    }
		  if(maxfill >= All.BunchSizeDensity)
		    break;

		  sendTask = ThisTask;
		  recvTask = ThisTask ^ ngrp;

		  if(recvTask < NTask)
		    {
		      if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
			{
			  /* send the results */
			  MPI_Sendrecv(&DensDataResult[nbuffer[ThisTask]],
				       nsend[recvTask * NTask + ThisTask] * sizeof(struct densdata_out),
				       MPI_BYTE, recvTask, TAG_DENS_B,
				       &DensDataPartialResult[noffset[recvTask]],
				       nsend_local[recvTask] * sizeof(struct densdata_out),
				       MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, &status);

			  /* add the result to the particles */
			  for(j = 0; j < nsend_local[recvTask]; j++)
			    {
			      source = j + noffset[recvTask];
			      place = DensDataIn[source].Index;

			      SphP[place].NumNgb += DensDataPartialResult[source].Ngb;
			      SphP[place].Density += DensDataPartialResult[source].Rho;
			      SphP[place].DivVel += DensDataPartialResult[source].Div;

			      SphP[place].DhsmlDensityFactor += DensDataPartialResult[source].DhsmlDensity;

			      SphP[place].Rot[0] += DensDataPartialResult[source].Rot[0];
			      SphP[place].Rot[1] += DensDataPartialResult[source].Rot[1];
			      SphP[place].Rot[2] += DensDataPartialResult[source].Rot[2];
			    }
			}
		    }

		  for(j = 0; j < NTask; j++)
		    if((j ^ ngrp) < NTask)
		      nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
		}
	      tend = second();
	      timecommsumm += timediff(tstart, tend);

	      level = ngrp - 1;
	    }

	  MPI_Allgather(&ndone, 1, MPI_INT, ndonelist, 1, MPI_INT, MPI_COMM_WORLD);
	  for(j = 0; j < NTask; j++)
	    ntotleft -= ndonelist[j];
	}



    /* do final operations on results */
    tstart = second();
    for(i = 0, npleft = 0; i < N_gas; i++){
      
	    if(P[i].Ti_endstep == All.Ti_Current){
	    {
     //   if(isinf(SphP[i].DhsmlDensityFactor) || isnan(SphP[i].DhsmlDensityFactor) ){
     //     printf("********************************************* Dhsmlrho is inf Before final dens operation\n") ;
     //     printf("ThisTask : %d, Step: %d, , pid: %d, rho: %g, numngb: %g, hsml: %g, divv: %g, curl: %g, dhsmlrho: %g\n", ThisTask, All.NumCurrentTiStep ,P[i].ID, SphP[i].Density, SphP[i].NumNgb, SphP[i].Hsml ,SphP[i].DivVel, SphP[i].CurlVel, SphP[i].DhsmlDensityFactor);
     //     printf("********************************************* Dhsmlrho is inf Before final dens operation\n") ;
     //   }    
     //   dhsml_bef = SphP[i].DhsmlDensityFactor;
     //   divv_bef  = SphP[i].DivVel;
     //   curlv_bef = SphP[i].CurlVel;
     //   rho_bef   = SphP[i].Density;
     //   numngb_bef = SphP[i].NumNgb;
     //   hsml_bef = SphP[i].Hsml;

     //   dhsml_bef_f = SphP[i].DhsmlDensityFactor;
     //   divv_bef_f  = SphP[i].DivVel;
     //   curlv_bef_f = SphP[i].CurlVel;
     //   rho_bef_f   = SphP[i].Density;
     //   numngb_bef_f = SphP[i].NumNgb;
     //   hsml_bef_f = SphP[i].Hsml;


        SphP[i].DhsmlDensityFactor =
		      1.0 / (1.0 + SphP[i].Hsml * SphP[i].DhsmlDensityFactor / (3.0 * SphP[i].Density));

		    SphP[i].CurlVel = sqrt(SphP[i].Rot[0] * SphP[i].Rot[0] +
				  SphP[i].Rot[1] * SphP[i].Rot[1] +
				  SphP[i].Rot[2] * SphP[i].Rot[2]) / SphP[i].Density;

    		SphP[i].DivVel /= SphP[i].Density;
		    dt_entr = (All.Ti_Current - (P[i].Ti_begstep + P[i].Ti_endstep) / 2) * All.Timebase_interval;

        #ifdef VARPOLYTROPE
		    SphP[i].Pressure = SphP[i].Entropy * pow(SphP[i].Density, SphP[i].Gamma);
        #else /* CHEMCOOL */
		    SphP[i].Pressure =
		      (SphP[i].Entropy + SphP[i].DtEntropy * dt_entr) * pow(SphP[i].Density, GAMMA) ;
        #endif

      //  if(isinf(SphP[i].DhsmlDensityFactor) || isnan(SphP[i].DhsmlDensityFactor) ){
      //    printf("********************************************* Dhsmlrho is inf After final dens operation\n") ;
      //    printf("ThisTask : %d, Step: %d, pid: %d, rho: %g, numngb: %g, hsml: %g, divv: %g, curl: %g, dhsmlrho: %g\n", ThisTask, All.NumCurrentTiStep ,P[i].ID, SphP[i].Density, SphP[i].NumNgb, SphP[i].Hsml ,SphP[i].DivVel, SphP[i].CurlVel, SphP[i].DhsmlDensityFactor);
      //    printf("ThisTask : %d, Step: %d, pid: %d, posx: %g, posy: %g, posz: %g\n", ThisTask, All.NumCurrentTiStep ,P[i].ID, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);
      //    printf("ThisTask : %d, Step: %d, pid: %d, velx: %g, vely: %g, velz: %g\n", ThisTask, All.NumCurrentTiStep ,P[i].ID, SphP[i].VelPred[0], SphP[i].VelPred[1], SphP[i].VelPred[2]);
      //    printf("ThisTask : %d, Step: %d, pid: %d, accx: %g, accy: %g, accz: %g\n", ThisTask, All.NumCurrentTiStep ,P[i].ID, SphP[i].HydroAccel[0], SphP[i].HydroAccel[1], SphP[i].HydroAccel[2]);
      //    printf("dhsml_bef: %g, divv_bef: %g, curlv_bef: %g, numngb_bef: %g, rho_bef: %g, hsml_bef: %g\n", dhsml_bef, divv_bef, curlv_bef, numngb_bef, rho_bef, hsml_bef);
      //    printf("dhsml_bef_f: %g, divv_bef_f: %g, curlv_bef_f: %g, numngb_bef_f: %g, rho_bef_f: %g, hsml_bef_f: %g\n", dhsml_bef_f, divv_bef_f, curlv_bef_f, numngb_bef_f, rho_bef_f, hsml_bef_f);
      //    printf("Step: %d, pid: %d, new_dhsml: %g , %g\n", All.NumCurrentTiStep, P[i].ID, hsml_bef * dhsml_bef / (3.0 * rho_bef), 1.0 / (1.0 + hsml_bef * dhsml_bef / (3.0 * rho_bef) ) );
      //    printf("Step: %d, pid: %d, new_dhsml_f: %g , %g\n", All.NumCurrentTiStep, P[i].ID,  hsml_bef_f * dhsml_bef_f / (3.0 * rho_bef_f), 1.0 / (1.0 + hsml_bef_f * dhsml_bef_f / (3.0 * rho_bef_f) ) );
      //    printf("********************************************* Dhsmlrho is inf After final dens operation\n") ;
      //    endrun(101010);
      //  }    
	    }
      /* now check whether we had enough neighbours */
      if(SphP[i].NumNgb < (All.DesNumNgb - All.MaxNumNgbDeviation) ||
		    (SphP[i].NumNgb > (All.DesNumNgb + All.MaxNumNgbDeviation)
		    && SphP[i].Hsml > (1.01 * All.MinGasHsml))){
		    /* need to redo this particle */
		    npleft++;
		    if(SphP[i].Left > 0 && SphP[i].Right > 0)
		      if((SphP[i].Right - SphP[i].Left) < 1.0e-3 * SphP[i].Left){
			      /* this one should be ok */
			      npleft--;
			      P[i].Ti_endstep = -P[i].Ti_endstep - 1;	/* Mark as inactive */
			      continue;
		      }

		      if(SphP[i].NumNgb < (All.DesNumNgb - All.MaxNumNgbDeviation))
		        SphP[i].Left = dmax(SphP[i].Hsml, SphP[i].Left);
		      else{
		        if(SphP[i].Right != 0){
			        if(SphP[i].Hsml < SphP[i].Right)
			          SphP[i].Right = SphP[i].Hsml;
			      }
		        else
			        SphP[i].Right = SphP[i].Hsml;
		      }

		      if(iter >= MAXITER - 10){
		        printf("i=%d task=%d ID=%d Hsml=%g Left=%g Right=%g Ngbs=%g Right-Left=%g\n   pos=(%g|%g|%g)\n",
			      i, ThisTask, (int) P[i].ID, SphP[i].Hsml, SphP[i].Left, SphP[i].Right,
			      (float) SphP[i].NumNgb, SphP[i].Right - SphP[i].Left, P[i].Pos[0], P[i].Pos[1],P[i].Pos[2]);
		        fflush(stdout);
		      }

		      if(SphP[i].Right > 0 && SphP[i].Left > 0)
		        SphP[i].Hsml = pow(0.5 * (pow(SphP[i].Left, 3) + pow(SphP[i].Right, 3)), 1.0 / 3);
		      else{
		        if(SphP[i].Right == 0 && SphP[i].Left == 0)
			        endrun(8188);	/* can't occur */

		      if(SphP[i].Right == 0 && SphP[i].Left > 0)
			{
			  if(P[i].Type == 0 && fabs(SphP[i].NumNgb - All.DesNumNgb) < 0.5 * All.DesNumNgb)
			    {
			      SphP[i].Hsml *=
				1 - (SphP[i].NumNgb -
				     All.DesNumNgb) / (NUMDIMS * SphP[i].NumNgb) * SphP[i].DhsmlDensityFactor;
			    }
			  else
			    SphP[i].Hsml *= 1.26;
			}

		      if(SphP[i].Right > 0 && SphP[i].Left == 0)
			{
			  if(P[i].Type == 0 && fabs(SphP[i].NumNgb - All.DesNumNgb) < 0.5 * All.DesNumNgb)
			    {
			      SphP[i].Hsml *=
				1 - (SphP[i].NumNgb -
				     All.DesNumNgb) / (NUMDIMS * SphP[i].NumNgb) * SphP[i].DhsmlDensityFactor;
			    }
			  else
			    SphP[i].Hsml /= 1.26;
			}
		    }

		  if(SphP[i].Hsml < All.MinGasHsml)
		    SphP[i].Hsml = All.MinGasHsml;
		}
	      else
		P[i].Ti_endstep = -P[i].Ti_endstep - 1;	/* Mark as inactive */
	    }
	}
      tend = second();
      timecomp += timediff(tstart, tend);


      numlist = malloc(NTask * sizeof(int) * NTask);
      MPI_Allgather(&npleft, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
      for(i = 0, ntot = 0; i < NTask; i++)
	ntot += numlist[i];
      free(numlist);

      if(ntot > 0)
	{
	  if(iter == 0)
	    tstart_ngb = second();

	  iter++;

	  if(iter > 0 && ThisTask == 0)
	    {
	      printf("ngb iteration %d: need to repeat for %d%09d particles.\n", iter,
		     (int) (ntot / 1000000000), (int) (ntot % 1000000000));
	      fflush(stdout);
	    }

	  if(iter > MAXITER)
	    {
	      printf("failed to converge in neighbour iteration in density()\n");
	      fflush(stdout);
	      endrun(1155);
	    }
	}
      else
	tend_ngb = second();
    }
  while(ntot > 0);


  /* mark as active again */
  for(i = 0; i < NumPart; i++)
    if(P[i].Ti_endstep < 0)
      P[i].Ti_endstep = -P[i].Ti_endstep - 1;

  free(ndonelist);
  free(nsend);
  free(nsend_local);
  free(nbuffer);
  free(noffset);


  /* collect some timing information */
  if(iter > 0)
    timengb = timediff(tstart_ngb, tend_ngb);
  else
    timengb = 0;

  MPI_Reduce(&timengb, &sumtimengb, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timecomp, &sumt, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timecommsumm, &sumcomm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timeimbalance, &sumimbalance, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      All.CPU_HydCompWalk += sumt / NTask;
      All.CPU_HydCommSumm += sumcomm / NTask;
      All.CPU_HydImbalance += sumimbalance / NTask;
      All.CPU_EnsureNgb += sumtimengb / NTask;
    }
}



/*! This function represents the core of the SPH density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 */
void density_evaluate(int target, int mode)
{
  int j, n, startnode, numngb, numngb_inbox;
  double h, h2, fac, hinv, hinv3, hinv4;
  double rho, divv, wk, dwk;
  double dx, dy, dz, r, r2, u, mass_j;
  double dvx, dvy, dvz, rotv[3];
  double weighted_numngb, dhsmlrho;
  FLOAT *pos, *vel;


  if(mode == 0)
    {
      pos = P[target].Pos;
      vel = SphP[target].VelPred;
      h = SphP[target].Hsml;
    }
  else
    {
      pos = DensDataGet[target].Pos;
      vel = DensDataGet[target].Vel;
      h = DensDataGet[target].Hsml;
    }

  h2 = h * h;
  hinv = 1.0 / h;
#ifndef  TWODIMS
  hinv3 = hinv * hinv * hinv;
#else
  hinv3 = hinv * hinv / boxSize_Z;
#endif
  hinv4 = hinv3 * hinv;

  rho = divv = rotv[0] = rotv[1] = rotv[2] = 0;
  weighted_numngb = 0;
  dhsmlrho = 0;

  startnode = All.MaxPart;
  numngb = 0;
  do
    {
      numngb_inbox = ngb_treefind_variable(&pos[0], h, &startnode);

      for(n = 0; n < numngb_inbox; n++)
	{
	  j = Ngblist[n];

	  dx = pos[0] - P[j].Pos[0];
	  dy = pos[1] - P[j].Pos[1];
	  dz = pos[2] - P[j].Pos[2];

#ifdef PERIODIC			/*  now find the closest image in the given box size  */
	  if(dx > boxHalf_X)
	    dx -= boxSize_X;
	  if(dx < -boxHalf_X)
	    dx += boxSize_X;
	  if(dy > boxHalf_Y)
	    dy -= boxSize_Y;
	  if(dy < -boxHalf_Y)
	    dy += boxSize_Y;
	  if(dz > boxHalf_Z)
	    dz -= boxSize_Z;
	  if(dz < -boxHalf_Z)
	    dz += boxSize_Z;
#endif
	  r2 = dx * dx + dy * dy + dz * dz;

	  if(r2 < h2)
	    {
	      numngb++;

	      r = sqrt(r2);

	      u = r * hinv;

	      if(u < 0.5)
		{
		  wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
		  dwk = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
		}
	      else
		{
		  wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
		  dwk = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
		}

	      mass_j = P[j].Mass;

	      rho += mass_j * wk;

	      weighted_numngb += NORM_COEFF * wk / hinv3;

	      dhsmlrho += -mass_j * (NUMDIMS * hinv * wk + u * dwk);
      //  if(isinf(dhsmlrho) || isnan(dhsmlrho) ){
      //    printf("*********************************************\n") ;
      //    printf("ThisTask : %d, Dhsmlrho is inf, pid: %d\n", ThisTask, P[target].ID);
      //    printf("mass_j: %g, wk: %g, dwk: %g, u: %g, hinv: %g\n ", mass_j, wk, dwk, u, hinv);
      //    printf("*********************************************\n") ; 
      //  }

	      if(r > 0)
		{
		  fac = mass_j * dwk / r;

		  dvx = vel[0] - SphP[j].VelPred[0];
		  dvy = vel[1] - SphP[j].VelPred[1];
		  dvz = vel[2] - SphP[j].VelPred[2];

		  divv -= fac * (dx * dvx + dy * dvy + dz * dvz);

		  rotv[0] += fac * (dz * dvy - dy * dvz);
		  rotv[1] += fac * (dx * dvz - dz * dvx);
		  rotv[2] += fac * (dy * dvx - dx * dvy);
		}
	    }
	}
    }
  while(startnode >= 0);

  if(mode == 0)
    {    
      SphP[target].NumNgb = weighted_numngb;
      SphP[target].Density = rho;
            
      SphP[target].DivVel = divv;
      SphP[target].DhsmlDensityFactor = dhsmlrho;
      SphP[target].Rot[0] = rotv[0];
      SphP[target].Rot[1] = rotv[1];
      SphP[target].Rot[2] = rotv[2];
    }
  else
    {
      DensDataResult[target].Rho = rho;

      DensDataResult[target].Div = divv;
      DensDataResult[target].Ngb = weighted_numngb;
      DensDataResult[target].DhsmlDensity = dhsmlrho;
      DensDataResult[target].Rot[0] = rotv[0];
      DensDataResult[target].Rot[1] = rotv[1];
      DensDataResult[target].Rot[2] = rotv[2];
    }
}




/*! This routine is a comparison kernel used in a sort routine to group
 *  particles that are exported to the same processor.
 */
int dens_compare_key(const void *a, const void *b)
{
  if(((struct densdata_in *) a)->Task < (((struct densdata_in *) b)->Task))
    return -1;

  if(((struct densdata_in *) a)->Task > (((struct densdata_in *) b)->Task))
    return +1;

  return 0;
}

#ifdef SINK

void setdens(){

  int startnode;
  int numsinks, numsinkstot, num, numbnd=0, numbndtot=0;
  FLOAT *local_sink_posx, *local_sink_posy, *local_sink_posz;
  FLOAT *local_sink_velx, *local_sink_vely, *local_sink_velz;
  FLOAT *local_sink_spinx, *local_sink_spiny, *local_sink_spinz;
  FLOAT *local_sink_mass;
  FLOAT *list_sink_posx, *list_sink_posy, *list_sink_posz;
  FLOAT *list_sink_velx, *list_sink_vely, *list_sink_velz;
  FLOAT *list_sink_spinx, *list_sink_spiny, *list_sink_spinz;
  FLOAT *list_sink_mass;
  FLOAT *pos, *vel;
  FLOAT sinkrad, seperation;	
  FLOAT SinkPos[3], SinkVel[3], SinkSpin[3];
  int i,j,k,n, igas, jgas;

  //physical variables
  FLOAT r1,r2, zmin, zmax, ymax;
  FLOAT v1, v2;

  FLOAT d, zo, zod, vol, th, phi, thm, zcut;
  FLOAT dt_entr, rho;
  FLOAT r,h,hinv,hinv3,hinv4,wk,dwk,u;
  FLOAT dx,dy,dz,dvx,dvy,dvz,divv,fac,dhsmlrho, dhsmlrho_old;
  int missing_ngb ;
  float costheta;

  FLOAT m_sink, m_gas, rs;
  FLOAT ngb_vabs, ngb_lsq ,vesc, KeplerL, ngb_L[3];
  FLOAT sys_pos[3], sys_vel[3], vec_d[3], rotv[3], vec_rs[3], vec_vs[3], temp_r;

  FLOAT vabs_min, lsq_min, vsx_min, vsy_min, vsz_min;

  srand(time(0));
  
  if(ThisTask == 0 ){
    printf("Starting density correction for boundary particles .... \n");
  }

  numsinks = NumPart - N_gas; 
  MPI_Allreduce(&numsinks, &numsinkstot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); 
  MPI_Barrier(MPI_COMM_WORLD);
 	  
  local_sink_posx = malloc(sizeof(double) * numsinkstot);
  local_sink_posy = malloc(sizeof(double) * numsinkstot);  
  local_sink_posz = malloc(sizeof(double) * numsinkstot); 
  local_sink_velx = malloc(sizeof(double) * numsinkstot);
  local_sink_vely = malloc(sizeof(double) * numsinkstot);  
  local_sink_velz = malloc(sizeof(double) * numsinkstot);

  local_sink_spinx = malloc(sizeof(double) * numsinkstot);
  local_sink_spiny = malloc(sizeof(double) * numsinkstot);
  local_sink_spinz = malloc(sizeof(double) * numsinkstot);


  local_sink_mass = malloc(sizeof(double) * numsinkstot);     
  list_sink_posx = malloc(sizeof(double) * numsinkstot * NTask);
  list_sink_posy = malloc(sizeof(double) * numsinkstot * NTask);  
  list_sink_posz = malloc(sizeof(double) * numsinkstot * NTask); 
  list_sink_velx = malloc(sizeof(double) * numsinkstot * NTask);
  list_sink_vely = malloc(sizeof(double) * numsinkstot * NTask);  
  list_sink_velz = malloc(sizeof(double) * numsinkstot * NTask);

  list_sink_spinx = malloc(sizeof(double) * numsinkstot * NTask);
  list_sink_spiny = malloc(sizeof(double) * numsinkstot * NTask);
  list_sink_spinz = malloc(sizeof(double) * numsinkstot * NTask);
 
  list_sink_mass = malloc(sizeof(double) * numsinkstot * NTask);  

  for(int igas=0; igas<N_gas; igas++){
    SphP[igas].Hsml_old = SphP[igas].Hsml;
    SphP[igas].DhsmlDensityFactor_old = SphP[igas].DhsmlDensityFactor;
    SphP[igas].Density_old = SphP[igas].Density;
  }  

     
  for(i = 0; i < numsinkstot; i++) local_sink_mass[i] = -1;
  
  for(i = 0; i < numsinks; i++){
    local_sink_posx[i] = P[i+N_gas].Pos[0];
    local_sink_posy[i] = P[i+N_gas].Pos[1];    
    local_sink_posz[i] = P[i+N_gas].Pos[2];
    local_sink_velx[i] = P[i+N_gas].Vel[0];
    local_sink_vely[i] = P[i+N_gas].Vel[1];    
    local_sink_velz[i] = P[i+N_gas].Vel[2];
    local_sink_spinx[i] = P[i+N_gas].Spin[0];
    local_sink_spiny[i] = P[i+N_gas].Spin[1];
    local_sink_spinz[i] = P[i+N_gas].Spin[2];
    local_sink_mass[i] = P[i+N_gas].Mass; 
  }  

  MPI_Allgather(local_sink_posx, numsinkstot, MPI_DOUBLE, list_sink_posx, numsinkstot, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Allgather(local_sink_posy, numsinkstot, MPI_DOUBLE, list_sink_posy, numsinkstot, MPI_DOUBLE, MPI_COMM_WORLD);  
  MPI_Allgather(local_sink_posz, numsinkstot, MPI_DOUBLE, list_sink_posz, numsinkstot, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Allgather(local_sink_velx, numsinkstot, MPI_DOUBLE, list_sink_velx, numsinkstot, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Allgather(local_sink_vely, numsinkstot, MPI_DOUBLE, list_sink_vely, numsinkstot, MPI_DOUBLE, MPI_COMM_WORLD);  
  MPI_Allgather(local_sink_velz, numsinkstot, MPI_DOUBLE, list_sink_velz, numsinkstot, MPI_DOUBLE, MPI_COMM_WORLD);  
  MPI_Allgather(local_sink_spinx, numsinkstot, MPI_DOUBLE, list_sink_spinx, numsinkstot, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Allgather(local_sink_spiny, numsinkstot, MPI_DOUBLE, list_sink_spiny, numsinkstot, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Allgather(local_sink_spinz, numsinkstot, MPI_DOUBLE, list_sink_spinz, numsinkstot, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Allgather(local_sink_mass, numsinkstot, MPI_DOUBLE, list_sink_mass, numsinkstot, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD); 
  
  for(i = 0; i < numsinkstot * NTask; i++){	 /* go through all the sink particles (From all processors) and find boundary gas */
    if(list_sink_mass[i] > 0){
      SinkPos[0] = list_sink_posx[i];
      SinkPos[1] = list_sink_posy[i];
      SinkPos[2] = list_sink_posz[i];
      SinkVel[0] = list_sink_velx[i];
      SinkVel[1] = list_sink_vely[i];
      SinkVel[2] = list_sink_velz[i];            
      SinkSpin[0] = list_sink_spinx[i];
      SinkSpin[1] = list_sink_spiny[i];
      SinkSpin[2] = list_sink_spinz[i];

      pos = SinkPos;
      vel = SinkVel;
      
      startnode = All.MaxPart;
      sinkrad = All.AccretionRadius;       
	
      do
      {
        num = ngb_treefind_variable(&pos[0], 2.0*sinkrad,&startnode); 
        for(n = 0; n < num; n++){
          igas = Ngblist[n];

          if(P[igas].Type == 0 && P[igas].Ti_endstep == All.Ti_Current && igas < N_gas ){   
          	for(seperation = 0,j = 0; j < 3; j++) seperation += (P[igas].Pos[j]-pos[j]) * (P[igas].Pos[j]-pos[j]);  
					  seperation = sqrt(seperation);   
            if(seperation <= sinkrad + SphP[igas].Hsml){
              numbnd++;
             //printf("ThisTask: %d, Step: %d,  correcting for %d, pid: %d\n", ThisTask, All.NumCurrentTiStep, igas, P[igas].ID);
              r1 = sinkrad;
              r2 = SphP[igas].Hsml;
              d  = seperation;
	            m_sink = list_sink_mass[i]; 
	            m_gas = P[igas].Mass;
	
	            v1 = 4.0*r1*r1*r1/3.0;
	            v2 = 4.0*r2*r2*r2/3.0;        
	            rho=0;        
	            zo = ( (r1*r1 - r2*r2)/d + d)/2.0;
	            zod = zo-d;
	            vol = 2.0*r1*r1*r1/3.0 - r1*r1*zo + zo*zo*zo/3.0  + 2.0*r2*r2*r2/3.0 + r2*r2*zod - zod*zod*zod/3.0 ;
              ymax = sqrt(r1*r1 - zo*zo);

	            if( (r1>=r2 && d>= r1-r2) || (r1<=r2 && d>=r2-r1)  ){
                rho =  SphP[igas].Density*( 1.0 + vol/v2 );
		            missing_ngb =  round( SphP[igas].NumNgb  * vol/v2 );
		            //printf("ThisTask: %d, Step: %d, Cond. 1(Overlap)  for %d, sep: %g, rs+h: %g, MSNGB: %d, pid: %d\n", ThisTask, All.NumCurrentTiStep, igas, seperation, sinkrad + SphP[igas].Hsml, missing_ngb,P[igas].ID);
                //printf("ThisTask: %d, Step: %d, Cond. 1(Overlap)  for %d\n: ", ThisTask, All.NumCurrentTiStep, P[igas].ID)  ;
	              for(int in=0; in<num; in++){
	                jgas = Ngblist[in];
		              for(rs = 0,j = 0; j < 3; j++) {
	                  rs += (P[jgas].Pos[j]-pos[j]) * (P[jgas].Pos[j]-pos[j]);
                  }  
                  rs = sqrt(rs);
		              zcut = ( (rs*rs - r2*r2)/d + d)/2.0;
		              thm = fabs(acos(zcut/rs)) ;
		              th = fabs(acos( (P[jgas].Pos[2] - pos[2]) /rs )); 
                  if(rs > d-r2 && rs < r1 && th <= thm  ){
                    missing_ngb--;			    
                  }  
	              }	  
              }
              else if( r1 < r2 && d< r2-r1  ){
	              vol = v1;
                //rho =  SphP[igas].Density*( 1.0 + vol/v2 );
                missing_ngb =  round( SphP[igas].NumNgb * vol/v2 );
		            //printf("ThisTask: %d, Step: %d, Cond. 2(sink inside)  for %d, sep: %g, rs+h: %g, MSNGB: %d, pid: %d\n", ThisTask, All.NumCurrentTiStep, igas, seperation, sinkrad + SphP[igas].Hsml, missing_ngb, P[igas].ID);
                //printf("ThisTask: %d, Step: %d, Cond. 2(sink indside)  for %d\n: ", ThisTask, All.NumCurrentTiStep, P[igas].ID)  ;
                for(int in=0; in<num; in++){
                  jgas = Ngblist[in];
                  for(rs = 0,j = 0; j < 3; j++){
                    rs += (P[jgas].Pos[j]-pos[j]) * (P[jgas].Pos[j]-pos[j]);
                  }  
                  rs = sqrt(rs);
                  if(rs<r1 ){
                    missing_ngb--;
                  }  
	              }
              }
	            else if ( r1>r2 && d<r1-r2 ){
                vol = v2;
		            missing_ngb = round( SphP[igas].NumNgb);
                //printf("ThisTask: %d, Step: %d, Cond. 3(gas inside)  for %d, sep: %g, rs+h: %g, MSNGB: %d, pid: %d\n", ThisTask, All.NumCurrentTiStep, igas, seperation, sinkrad + SphP[igas].Hsml, missing_ngb,P[igas].ID);
                //printf("ThisTask: %d, Step: %d, Cond. 3(gas indside)  for %d\n: ", ThisTask, All.NumCurrentTiStep, P[igas].ID)  ;
		            vec_d[0] = P[igas].Pos[0] - SinkPos[0];
                vec_d[1] = P[igas].Pos[1] - SinkPos[1];
                vec_d[2] = P[igas].Pos[2] - SinkPos[2];

		            for(int in=0; in<num; in++){
                  jgas = Ngblist[in];
                  vec_rs[0] = P[jgas].Pos[0] - SinkPos[0];
                  vec_rs[1] = P[jgas].Pos[1] - SinkPos[0];
                  vec_rs[2] = P[jgas].Pos[2] - SinkPos[0];
                  rs = sqrt( vec_rs[0]*vec_rs[0] + vec_rs[1]*vec_rs[1] + vec_rs[2]*vec_rs[2] );
                  costheta = (vec_rs[0]*vec_d[0]+vec_rs[1]*vec_d[1]+vec_rs[2]*vec_d[2])/(rs*d);   
                  if( rs*rs + d*d - 2.0*rs*d*costheta  <r2*r2 )
                    missing_ngb--;
                }
	            }
              else {
                //printf("ThisTask: %d,  should not happen %d %g, %g, %g, %g, NumNgb: %g \n", ThisTask, igas, r1, r2, d , d-r1+r2, SphP[igas].NumNgb );
                endrun(899);
              }

	            if(missing_ngb > 0){ 
                h = pow( 3.0*SphP[igas].NumNgb*m_gas/(4.0*M_PI*rho), 1.0/3.0 );    
                //printf(" ThisTask: %d, igas: %d, NumNgb: %g, corrected msngb: %d, m_sink: %g, m_gas: %g, vol/v2: %g\n", ThisTask, igas, SphP[igas].NumNgb, missing_ngb, m_sink, m_gas, vol/v2);
                // printf( "ThisTask : %d, igas: %d, r1: %g, r2: %g, r1-r2: %g, d: %g \n", ThisTask, igas, r1, r2, r1-r2,d ); 
	              h =  SphP[igas].Hsml;
                hinv = 1.0/h;
                hinv3 = hinv*hinv*hinv;
                hinv4 = hinv3*hinv;
									
     	          rho = 0;
      	        dhsmlrho =0;
                divv=0;
                rotv[0]=0;
                rotv[1]=0;
                rotv[2]=0;  
									
                for(int imiss=0; imiss<missing_ngb; imiss++ ){
                  if( (r1>=r2 && d>= r1-r2) || (r1<=r2 && d>=r2-r1)  ){
                    zmin = d-r2;
		                if (zmin < 0){
		                  zmin = -zmin;
		                }	  
	                  zmax = r1;
		                rs =  ( (double)rand() / (double)RAND_MAX )*(zmax - zmin)  + zmin ;
		                zo = ( (rs*rs -r2*r2)/d + d)/2.0;
		                thm = acos(zo/rs);
		                th =  ( (double)rand() / (double)RAND_MAX ) * thm ;
		                phi = ( (double)rand() / (double)RAND_MAX )*2.0*M_PI;
		                //printf( "Condition 1, ThisTask : %d, igas: %d, r1: %g, r2: %g, r1-r2: %g, d-r2: %g, thm: %g, th: %g, zmin: %g, zmax: %g \n", ThisTask, igas, r1, r2, r1-r2,d-r2, thm, th, zmin, zmax );
                  }                      
		              else if(r1 < r2 && d< r2-r1)  {
		                zmin = 0 ;
	                  zmax = r1;	  
		                thm = M_PI;
	                  rs =  ( (double)rand() / (double)RAND_MAX )*r1;
	                  th = 2.0* ( ( (double)rand() / (double)RAND_MAX ) - 1.0) * M_PI ;
		                phi = ( (double)rand() / (double)RAND_MAX )*2.0*M_PI;
		                //printf( "Condition 2, ThisTask : %d, igas: %d, r1: %g, r2: %g, r1-r2: %g, d: %g, thm: %g, th: %g, zmin: %g, zmax: %g \n", ThisTask, igas, r1, r2, r1-r2,d, thm, th, zmin, zmax );
	                }
	                else if ( r1>r2 && d<r1-r2 ){
                    zmin = d-r2;
                    if (zmin < 0){
		                  zmin = -zmin;
		                }	  
		                zmax = d+r2;
                    rs =  ( (double)rand() / (double)RAND_MAX )*(zmax - zmin)  + zmin ;
                    zo = ( (rs*rs -r2*r2)/d + d)/2.0;
		                thm = acos(zo/rs);
                    th =  ( (double)rand() / (double)RAND_MAX ) * thm ;
		                phi = ( (double)rand() / (double)RAND_MAX )*2.0*M_PI;
		              }	
		              else{
	                  //printf("ThisTask: %d this is impossible, id: %d\n", ThisTask, igas);
		                endrun(999);
		              } 

                  vec_rs[0] = rs*sin(th)*cos(phi);
		              vec_rs[1] = rs*sin(th)*sin(phi);
                  vec_rs[2] = rs*cos(th);
                  r = sqrt( rs*rs + d*d - 2.0*rs*d*cos(th) ); 
                  if(isnan(r)) {
                    printf("ThisTask: %d, r is nan,igas: %d, rs: %g, d: %g, r: %g, costh: %g, th: %g, rssq: %g\n", ThisTask, igas, rs, d,r, cos(th), th, rs*rs + d*d - 2.0*rs*d*cos(th));
                    endrun(1001); 
                  }
		 
		              //printf("ThisTask: %d, IGAS: %d, vec_rs[0]:, %g, vec_rs[1]: %g, vec_rs[2]: %g, rs: %g, th: %g, phi: %g, zmin: %g, zmax: %g, u: %g\n ", ThisTask, igas, vec_rs[0], vec_rs[1], vec_rs[2], rs, th, phi, zmin, zmax, r/r2); 

                  //transform to system coordinates 
                  vec_d[0] = P[igas].Pos[0] - SinkPos[0];
                  vec_d[1] = P[igas].Pos[1] - SinkPos[1];
                  vec_d[2] = P[igas].Pos[2] - SinkPos[2];

                  th = acos( vec_d[2]/( sqrt(vec_d[0]*vec_d[0] + vec_d[1]*vec_d[1] + vec_d[2]*vec_d[2] )  ) );
                  phi = atan2( vec_d[1],vec_d[0] );

                  sys_pos[0] = vec_rs[0]*cos(th)*cos(phi) - vec_rs[1]*sin(phi) + vec_rs[2]*sin(th)*cos(phi);
                  sys_pos[1] = vec_rs[0]*cos(th)*sin(phi) + vec_rs[1]*cos(phi) + vec_rs[2]*sin(th)*sin(phi);
                  sys_pos[2] = -vec_rs[0]*sin(th) + vec_rs[2]*cos(th);

		              temp_r = sys_pos[0]*sys_pos[0] + sys_pos[1]*sys_pos[1] + sys_pos[2]*sys_pos[2];

                  //printf("ThisTask: %d, IGAS: %d, s_pos[0]:, %g, s_pos[1]: %g, s_pos[2]: %g, t_r: %g, r1: %g, d-r2: %g\n ", ThisTask, igas, sys_pos[0], sys_pos[1], sys_pos[2], temp_r, r1*r1, (d-r2)*(d-r2));

		              if( temp_r > r1*r1 || temp_r < (d-r2)*(d-r2) ){
		                //printf("ThisTask: %d Position exeption 1 for %d, temp_r: %g, r2: %g, d-r1: %g \n", ThisTask, igas, temp_r, r1*r1, (d-r2)*(d-r2) );  
		                endrun(101);
		              }

		              sys_pos[0] = sys_pos[0] + SinkPos[0];
		              sys_pos[1] = sys_pos[1] + SinkPos[1];
		              sys_pos[2] = sys_pos[2] + SinkPos[2];

                  vesc = sqrt(2.0*All.G*m_sink/rs);
		              KeplerL = All.G * m_sink * rs;
		              ngb_vabs = 1e10;
		              ngb_lsq = 1e10;

		              vabs_min = ngb_vabs;
		              lsq_min = ngb_lsq;
                  int itr=0;
		              float vdotr = 100.0;
		              float vth, vphi;
		              costheta=-10.0;
      		        while ( ngb_vabs >= vesc || ngb_lsq >= KeplerL || costheta < 0 || vdotr > 0 ){
	                  ngb_vabs = ( (double)rand() / (double)RAND_MAX ) * vesc;
                    vth = ( (double)rand() / (double)RAND_MAX )*M_PI;
                    vphi = ( (double)rand() / (double)RAND_MAX )*2.0*M_PI;   		  
		                vec_vs[0] = ngb_vabs*sin(vth)*cos(vphi);
		                vec_vs[1] = ngb_vabs*sin(vth)*sin(vphi);
		                vec_vs[2] = ngb_vabs*cos(vth);
		                ngb_L[0] = vec_rs[1]*vec_vs[2] - vec_rs[2]*vec_vs[1];
                    ngb_L[1] = vec_rs[2]*vec_vs[0] - vec_rs[0]*vec_vs[2];
		                ngb_L[2] = vec_rs[0]*vec_vs[1] - vec_rs[1]*vec_vs[0];
            
		                for(ngb_lsq = 0, ngb_vabs=0, j = 0;j < 3; j++) {
		                  ngb_lsq  += ngb_L[j]*ngb_L[j];
		                  //  ngb_vabs += vec_vs[j]*vec_vs[j]; 
		                }

		                //ngb_vabs = sqrt(ngb_vabs);
		                vdotr = vec_rs[0]*vec_vs[0] + vec_rs[1]*vec_vs[1] + vec_rs[2]*vec_vs[2];
	                  costheta = (SinkSpin[0]*ngb_L[0] + SinkSpin[1]*ngb_L[1] + SinkSpin[2]*ngb_L[2])/sqrt(ngb_lsq);
                    costheta = costheta/sqrt(SinkSpin[0]*SinkSpin[0] + SinkSpin[1]*SinkSpin[1] + SinkSpin[2]*SinkSpin[2]  );
		                if( vabs_min > ngb_vabs && lsq_min > ngb_lsq && costheta > 0 && vdotr < 0 ){ 
		                  vabs_min = ngb_vabs;
		                  lsq_min = ngb_lsq;
		                  vsx_min = vec_vs[0];
                      vsy_min = vec_vs[1];
		                  vsz_min = vec_vs[2];
		                }   
     		            //printf("ThisTask: %d, itr: %d, Correcting vel for %d, vesc:, %g, KL: %g, ngb_vabs: %g, ngb_lsq: %g \n", ThisTask, itr, igas, vesc, KeplerL, ngb_vabs, ngb_lsq );
		                itr++;
		                if(itr >= 2000 )
	                    break;		   
      		        } //while 
		              if ( ngb_vabs >= vesc || ngb_lsq >= KeplerL || costheta < 0 || vdotr > 0 ){
		                //printf("ThisTask: %d, igas: %d, Velocity not converged, v-vesc: %g, l-lkep: %g, costheta: %g, vdotr: %g \n", ThisTask, igas, vabs_min-vesc, lsq_min-KeplerL, costheta, vdotr);
		                vec_vs[0] = vsx_min;
		                vec_vs[1] = vsy_min;
		                vec_vs[2] = vsz_min;
		                endrun(100);
		              }
		              else{
		                //  printf("ThisTask:%d, igas: %d, Velocity assigned in itr: %d, costheta: %g, vabs-vesc: %g, l-KepL: %g \n", ThisTask, igas, itr, costheta, ngb_vabs-vesc, ngb_lsq-KeplerL );
		              }
		
		              sys_vel[0] = vec_vs[0]*cos(th)*cos(phi) - vec_vs[1]*sin(phi) + vec_vs[2]*sin(th)*cos(phi) + SinkVel[0];
                  sys_vel[1] = vec_vs[0]*cos(th)*sin(phi) + vec_vs[1]*cos(phi) + vec_vs[2]*sin(th)*sin(phi) + SinkVel[1];
                  sys_vel[2] = -vec_vs[0]*sin(th) + vec_vs[2]*cos(th)  + SinkVel[2];
		
		              u = r*hinv;
		              if(u < 0.5){
		                wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
		                dwk = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
		              }
		              else{
		                wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
		                dwk = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
		              } //check this
		              if( wk < 0 || r<0 || u>1.0 || isnan(wk) || isnan(dwk) || isnan(u)){
		                printf("LOL YOU SUCK, wk: %g, dwk: %g, r: %g, u: %g, hinv3: %g, hinv4: %g \n", wk, dwk, r, u, hinv3, hinv4);
		                endrun(454);
		              }
		              rho += m_gas * wk;
		              dhsmlrho += -m_gas * (NUMDIMS * hinv * wk + u * dwk);
		              fac = m_gas * dwk / r;
	
		              dx = P[igas].Pos[0] - sys_pos[0];
	                dy = P[igas].Pos[0] - sys_pos[1];
	                dz = P[igas].Pos[2] - sys_pos[2];	
			
		              dvx = P[igas].Vel[0] - sys_vel[0];
		              dvy = P[igas].Vel[1] - sys_vel[1];
		              dvz = P[igas].Vel[2] - sys_vel[2];

	                divv -= fac * (dx * dvx + dy * dvy + dz * dvz);	
		
		              rotv[0] += fac * (dz * dvy - dy * dvz);
		              rotv[1] += fac * (dx * dvz - dz * dvx);
		              rotv[2] += fac * (dy * dvx - dx * dvy);        
    	          } //msngb		      

	              //printf("\n ThisTask: %d, igas: %d, Oldvalues: div: %g, curl: %g, dhsml: %g, hsml: %g, dens: %g, pres: %g, rot[2]: %g\n", ThisTask, igas, SphP[igas].DivVel, SphP[igas].CurlVel, SphP[igas].DhsmlDensityFactor, SphP[igas].Hsml, SphP[igas].Density, SphP[igas].Pressure, SphP[igas].Rot[2] );

	              //undo the final operations
                //dhsmlrho_old = SphP[igas].DhsmlDensityFactor;
	              SphP[igas].DivVel  = SphP[igas].DivVel*SphP[igas].Density;
	              SphP[igas].CurlVel = SphP[igas].CurlVel*SphP[igas].Density;							
	              if( SphP[igas].DhsmlDensityFactor == 0 ){
	                printf("Houston, we got a problem\n");
		              endrun(212);
	              }
	              else{
	                SphP[igas].DhsmlDensityFactor = (1.0/SphP[igas].DhsmlDensityFactor - 1.0)*NUMDIMS*SphP[igas].Density/SphP[igas].Hsml;   
	              }

                if(isnan(rho)){
                  printf("NAN error occured, 1066") ;
                }
	              //correction
	              SphP[igas].Density += rho;	   
	              SphP[igas].Hsml = pow( 3.0*SphP[igas].NumNgb*m_gas/(4.0*M_PI*SphP[igas].Density), 1.0/3.0 );						
	              SphP[igas].DivVel += divv;
	              SphP[igas].DhsmlDensityFactor += dhsmlrho;
	              SphP[igas].Rot[0] += rotv[0];
	              SphP[igas].Rot[1] += rotv[1];
	              SphP[igas].Rot[2] += rotv[2];
						
	              //redo the final operations
	              SphP[igas].DivVel /= SphP[igas].Density;
	              SphP[igas].DhsmlDensityFactor = 1.0 / (1.0 + SphP[igas].Hsml * SphP[igas].DhsmlDensityFactor / (NUMDIMS * SphP[igas].Density) );
	              SphP[igas].CurlVel = sqrt(SphP[igas].Rot[0] * SphP[igas].Rot[0] +
				                                  SphP[igas].Rot[1] * SphP[igas].Rot[1] +
					                                SphP[igas].Rot[2] * SphP[igas].Rot[2]) / SphP[igas].Density;
	              dt_entr = (All.Ti_Current - (P[igas].Ti_begstep + P[igas].Ti_endstep) / 2) * All.Timebase_interval;
                #ifdef VARPOLYTROPE
		            SphP[i].Pressure = SphP[i].Entropy * pow(SphP[i].Density, SphP[i].Gamma);
                #else /* CHEMCOOL */
		            SphP[i].Pressure = (SphP[i].Entropy + SphP[i].DtEntropy * dt_entr) * pow(SphP[i].Density, GAMMA) ;
                #endif
	              //SphP[igas].Pressure = (SphP[igas].Entropy + SphP[igas].DtEntropy * dt_entr) * pow(SphP[igas].Density, GAMMA);
	      
                //printf("ThisTask: %d, Step: %d,  corrected for %d, pid: %d, dhsml_rho: %g   %g \n", ThisTask, All.NumCurrentTiStep, igas, P[igas].ID, dhsmlrho_old, SphP[igas].DhsmlDensityFactor);
	              //printf(" ThisTask: %d, igas: %d, Correction density: %g, divv: %g, rot[0]: %g, rot[1]: %g, rot[2]: %g\n", ThisTask, igas, rho, divv/SphP[igas].Density, rotv[0]/SphP[igas].Density, rotv[1]/SphP[igas].Density, rotv[2]/SphP[igas].Density ); 
                //printf(" ThisTask: %d, igas: %d, Correction density: %g, divv: %g, rot[0]: %g, rot[1]: %g, rot[2]: %g\n", ThisTask, igas, rho, divv, rotv[0], rotv[1], rotv[2] ); 
	              //printf("ThisTask: %d, igas: %d, Newvalues: div: %g, curl: %g, dhsml: %g, hsml: %g, dens: %g, pres: %g, rot[2]: %g, costheta: %g\n\n", ThisTask, igas, SphP[igas].DivVel, SphP[igas].CurlVel, SphP[igas].DhsmlDensityFactor, SphP[igas].Hsml, SphP[igas].Density, SphP[igas].Pressure, SphP[igas].Rot[2],costheta );
              } //msngb>0
  	        } //seperation
	        } //type=0
	      } //num          
      }while(startnode>=0);
    } //m_sink > 0 		
  }// numsinktot	
  MPI_Allreduce(&numbnd, &numbndtot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(ThisTask ==0){
    printf("Corrected density for %d particles \n", numbndtot);
  }
			
  //  printf("ThisTask: %d, Reached end of density correction\n", ThisTask);
  MPI_Barrier(MPI_COMM_WORLD);
}
#endif


















