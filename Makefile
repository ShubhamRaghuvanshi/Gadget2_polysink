
#--------------------------------------- Basic operation mode of code
#OPT   +=  -DPERIODIC 
OPT   +=  -DUNEQUALSOFTENINGS
OPT   +=  -DADAPTIVE_GRAVSOFT_FORGAS


#--------------------------------------- Things that are always recommended
OPT   +=  -DPEANOHILBERT
OPT   +=  -DWALLCLOCK   


#--------------------------------------- TreePM Options
#OPT   +=  -DPMGRID=128
#OPT   +=  -DPLACEHIGHRESREGION=1
#OPT   +=  -DENLARGEREGION=1.2
#OPT   +=  -DASMTH=1.25
#OPT   +=  -DRCUT=4.5


#--------------------------------------- Single/Double Precision
OPT   +=  -DDOUBLEPRECISION      
#OPT   +=  -DDOUBLEPRECISION_FFTW      


#--------------------------------------- Time integration options
OPT   +=  -DSYNCHRONIZATION
#OPT   +=  -DFLEXSTEPS
#OPT   +=  -DPSEUDOSYMMETRIC
OPT   +=  -DNOSTOP_WHEN_BELOW_MINTIMESTEP
#OPT   +=  -DNOPMSTEPADJUSTMENT


#--------------------------------------- Output options
OPT   +=  -DHAVE_HDF5  
OPT   +=  -DH5_USE_16_API
OPT   +=  -DOUTPUTPOTENTIAL
OPT   +=  -DOUTPUTACCELERATION
#OPT   +=  -DOUTPUTCHANGEOFENTROPY
#OPT   +=  -DOUTPUTTIMESTEP


#--------------------------------------- Things for special behaviour
#OPT   +=  -DNOGRAVITY     
#OPT   +=  -DNOTREERND 
#OPT   +=  -DNOTYPEPREFIX_FFTW        
#OPT   +=  -DLONG_X=60
#OPT   +=  -DLONG_Y=5
#OPT   +=  -DLONG_Z=0.2
#OPT   +=  -DTWODIMS
#OPT   +=  -DSPH_BND_PARTICLES
#OPT   +=  -DNOVISCOSITYLIMITER
#OPT   +=  -DCOMPUTE_POTENTIAL_ENERGY
#OPT   +=  -DLONGIDS
#OPT   +=  -DSELECTIVE_NO_GRAVITY=2+4+8+16

#--------------------------------------- Testing and Debugging options
#OPT   +=  -DFORCETEST=0.1


#--------------------------------------- Glass making
#OPT   +=  -DMAKEGLASS=262144

# sink particle options
OPT += -DSINK
#OPT += -DSINKBND

#OPT +=  -DISOTHERM_EQS
#OPT += -DADIABATIC
#OPT += -DPOLYTROPE
OPT += -DVARPOLYTROPE

#----------------------------------------------------------------------
# Here, select compile environment for the target machine. This may need 
# adjustment, depending on your local system. Follow the examples to add
# additional target platforms, and to get things properly compiled.
#----------------------------------------------------------------------


#--------------------------------------- Adjust settings for target computer

CC       =  mpiicc
OPTIMIZE =  -O3 -Wall
GSL_INCL =  -I/u/sraghu/software/GSL_install/include
GSL_LIBS =  -L/u/sraghu/software/GSL_install/lib
FFTW_INCL=  -I/u/sraghu/software/FFTW2_install/include
FFTW_LIBS=  -L/u/sraghu/software/FFTW2_install/lib
MPICHLIB =  -L/mpcdf/soft/SLE_15/packages/x86_64/intel_oneapi/2022.3/mpi/2021.7.1/lib
HDF5INCL =  -I/u/sraghu/software/HDF5_install/include
HDF5LIB  =  -L/u/sraghu/software/HDF5_install/lib -lhdf5 -lz

OPTIONS =  $(OPTIMIZE) $(OPT) -std=gnu99

EXEC   = Gadget2

OBJS   = main.o  run.o  predict.o begrun.o endrun.o global.o  \
	 timestep.o  init.o restart.o  io.o    \
	 accel.o   read_ic.o  ngb.o  \
	 system.o  allocate.o  density.o  \
	 gravtree.o hydra.o  driftfac.o  \
	 domain.o  allvars.o potential.o  \
         forcetree.o   peano.o gravtree_forcetest.o \
	 pm_periodic.o pm_nonperiodic.o longrange.o 

INCL   = allvars.h  proto.h  tags.h Makefile

CFLAGS = $(OPTIONS) $(GSL_INCL) $(FFTW_INCL) $(HDF5INCL)


ifeq (NOTYPEPREFIX_FFTW,$(findstring NOTYPEPREFIX_FFTW,$(OPT)))    # fftw installed with type prefix?
  FFTW_LIB = $(FFTW_LIBS) -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw
else
ifeq (DOUBLEPRECISION_FFTW,$(findstring DOUBLEPRECISION_FFTW,$(OPT)))
  FFTW_LIB = $(FFTW_LIBS) -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw
else
  FFTW_LIB = $(FFTW_LIBS) -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw
endif
endif


LIBS   =   $(HDF5LIB) -g  $(MPICHLIB)  $(GSL_LIBS) -lgsl -lgslcblas -lm $(FFTW_LIB) 

$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(LIBS)   -o  $(EXEC)  

$(OBJS): $(INCL) 


clean:
	rm -f $(OBJS) $(EXEC)



