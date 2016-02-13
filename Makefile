#Makefile for OptQC
#Variables
progname = OptQC
compiler = mpif90
lpflags = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm
#switch = -O2 -warn all -nogen-interfaces -xHost -ipo 
switch = -O2
objects = OptQC_Bridge.o OptQC_Common.o OptQC_Common_Module.o OptQC_CSD.o OptQC_Main.o OptQC_Perm.o OptQC_WKVar.o OptQC_CPLX_Main.o

#Makefile
$(progname): $(objects)
	$(compiler) -o $(progname) $(switch) *.o $(lpflags)

%.o: %.f90
	$(compiler) -c $(switch) $<

common_module.mod: OptQC_Common_Module.o
OptQC_Bridge.o: common_module.mod
OptQC_CSD.o: arrays_real.mod arrays_cplx.mod common_module.mod
OptQC_Main.o: common_module.mod csd_tools.mod
OptQC_Perm.o: csd_tools.mod
csd_tools.mod: OptQC_CSD.o
arrays_real.mod: OptQC_WKVar.o
OptQC_CPLX_Main.o: common_module.mod csd_tools.mod
arrays_cplx.mod: OptQC_CPLX_WKVar.o

#Cleaning files
clean:
	rm -rf *.o
	rm -rf *.mod
	rm -rf *__genmod.f90
	rm -rf $(progname)
