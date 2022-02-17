#----------------------------------------------------------------------
# SUPPORTED COMPILE: pgi-cuda, pgi-cuda-sync, pgi, intel, gnu, ibmxl, cray-cuda
# SUPPORTED MODES  : opt, debug [default=debug]
# SUPPORTED PRECISIONS (PREC) : single, double [default=double]
#----------------------------------------------------------------------
COMPILE = "nvfortran"
MODE    = "opt"
PREC    = "double"
#----------------------------------------------------------------------
ifeq (${PREC},"single") # mpif90 shipped with PGI compiler
    PRECFLAGS = -DSINGLE_PRECISION 
endif
#----------------------------------------------------------------------

# singularity build nvhpc_21.1_devel.sif docker://nvcr.io/nvidia/nvhpc:21.1-devel-cuda_multi-ubuntu20.04
# singularity shell --nv nvhpc_21.1_devel.sif
ifeq (${COMPILE},"nvfortran") # mpif90 shipped with PGI compiler
    FC = mpifort # precompiled mpifort included with nvidia hpc
    ifeq (${MODE},"opt") 
        OLEVEL = -O3 -Mcuda=ccall,cuda11.2 -DUSE_CUDA -DCUDA_ASYNC #-DNOCUDAAWAREMPI
    else
        OLEVEL = -Mcuda=ccall -DUSE_CUDA -DCUDA_ASYNC -g -DNOCUDAAWAREMPI
    endif
endif

ifeq (${COMPILE},"gnu") # mpif90 from MPI library
    FC = mpif90
    ifeq (${MODE},"opt") 
        OLEVEL = -O3 #-fno-inline-functions -fno-underscoring
    else
        OLEVEL = -O2 -fbounds-check -g -fbacktrace -fdump-core -ffpe-trap=zero,invalid,overflow -Wuninitialized -Wall
        #OLEVEL = -O2 -fbounds-check -g -fbacktrace -fdump-core
    endif
endif

RM = rm -f

#----------------------------------------------------------------------
# Base code
#----------------------------------------------------------------------
MAIN = main.o
OBJ_FILES = alloc.o bcdf.o bcextr.o bcfree.o bc.o bcrecyc.o bcrelax.o bcshk.o bcswap.o \
    bcwall.o bcwall_pl.o bcwall_staggered.o check_input.o computeav.o computemetrics.o \
    computeresidual.o compute_coeff_xyz_midpi.o constants.o df_par.o digital_filtering_util.o \
    finalize.o generategrid.o generateinflowrand.o generatewmean_channel.o \
    generatewmean.o heatflux_compute.o init_channel.o init.o initurb.o init_windtunnel.o manage_solver.o osw.o pgrad.o \
    prims.o printres.o readdf.o readgrid.o readinp.o readrst.o readstat.o recyc.o rk.o sensor.o \
    setup.o solver.o startmpi.o stats.o step.o target_reystress.o tbforce.o updateghost.o utility.o utyibm.o \
    visflx.o writedf.o writefield.o writefieldvtk.o writegridplot3d.o writerst.o \
    writestatbl.o writestatchann.o writestat.o writestatzbl.o write_wallpressure.o visflx_stag.o

OBJ_MODS = mod_streams.o mod_sys.o euler.o
OBJS = $(OBJ_MODS) $(OBJ_FILES) $(MAIN)
EXEC = streams.exe
LDFLAGS =
#----------------------------------------------------------------------
# Additional package
#----------------------------------------------------------------------
-include makefile.inc
#----------------------------------------------------------------------

$(EXEC): $(OBJS) Makefile $(wildcard makefile.inc)
	$(FC) $(PRECFLAGS) $(OLEVEL) -o $@ $(OBJS) $(FCLIBS) $(LDFLAGS)

%.o:%.F90 Makefile $(wildcard makefile.inc)
	$(FC) $(PRECFLAGS) $(INCLUDE) $(OLEVEL) -c $< -o $@

%.o:%.cpp Makefile $(wildcard makefile.inc)
	$(CXX) $(PRECFLAGS) $(INCLUDE) -O2 -c $< -o $@

%.o:%.cxx Makefile $(wildcard makefile.inc)
	$(CXX) $(PRECFLAGS) $(INCLUDE) -O2 -c $< -o $@

.PHONY: clean
clean:
	$(RM) $(EXEC) $(MAIN) $(OBJS) $(OBJ_MODS) $(GRID) *.mod
