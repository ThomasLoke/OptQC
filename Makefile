#Makefile for OptQC
#Variables
progname = OptQC
compiler = mpif90
lpflags = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm
switch = -O2 -i8 -r8
objects = OptQC_Bridge.o OptQC_Common.o CYG_INDEXTABLE.o OptQC_Main.o CYGR_BLKCSD.o CYGR_CSD.o CYGR_CUTGATE.o \
		CYGR_WRITEF.o OptQC_GateCount.o OptQC_Perm.o OptQC_Output.o \
		OptQC_Common_Module.o OptQC_WKVar.o OptQC_Mem.o

#Makefile
$(progname): $(objects)
	$(compiler) -o $(progname) $(switch) $(objects) $(lpflags)

%.o: %.f90
	$(compiler) -c $(switch) $<

OptQC_Main.o: memwork_real.mod
CYGR_CSD.o: memwork_real.mod
common_module.mod: OptQC_Common_Module.o
arrays_real.mod: OptQC_WKVar.o
memwork_real.mod: OptQC_Mem.o
OptQC_Mem.o: arrays_real.mod common_module.mod

#Cleaning files
clean:
	rm -rf *.o
	rm -rf *.mod
	rm -rf $(progname)
