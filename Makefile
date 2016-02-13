#Makefile for OptQC
#Variables
progname = OptQC
compiler = mpif90
lpflags = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm
switch = -O2 -i8 -r8
objects = OptQC_Bridge.o OptQC_Common.o OptQC_Common_Module.o OptQC_CSD.o OptQC_Main.o OptQC_Perm.o OptQC_WKVar.o

#Makefile
$(progname): $(objects)
	$(compiler) -o $(progname) $(switch) *.o $(lpflags)

%.o: %.f90
	$(compiler) -c $(switch) $<

OptQC_CSD.o: arrays_real.mod common_module.mod
OptQC_Main.o: csd_real.mod
OptQC_Perm.o: csd_real.mod
arrays_real.mod: OptQC_WKVar.o
common_module.mod: OptQC_Common_Module.o
csd_real.mod: OptQC_CSD.o

#Cleaning files
clean:
	rm -rf *.o
	rm -rf *.mod
	rm -rf $(progname)
