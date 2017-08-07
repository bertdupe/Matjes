include ../make.sys

#export LD=$(LD)
#export $FC
#export $(CC) $(FFLAGS) $(CPPS)


.SUFFIXES: .F .f .f90 .F90

SUBDIRS = \
	Energy \
	dynamics \
	MonteCarlo \
	parallelization \
	io-inputs \
	io-ouputs \
	random-numbers \
	setup \
	utilities \
	GNEB
	
BIN_DIR= ../bin	

toto:
	@(cd $(SUBDIRS) && make $*)

Matjes : $(OBJS) $(MPI) $(DIP)
	$(LD) -o $@ $(LDFLAGS) $(FIRST) $(OBJS) $(DIP) $(MPI) $(LIBRARIES)

.f.o:
	$(FC) -cpp $(CPPS) -c $< $(FFLAGS) -o $@

.F.o:
	$(FC) -c $< $(FFLAGS) -o $@

.f77.o:
	$(FC) -c $< $(FFLAGS) -o $@

.f90.o:
	$(FC) -cpp $(CPPS) -c $< $(FFLAGS) -o $@

.F90.o:
	$(FC) -cpp $(CPPS) -c $< $(FFLAGS) $(INCLUDES) -o $@

.f.a:
	$(FC) -c $(FFLAGS) $<
	$(AR) $(ARFLAGS) $@ $*.o
	$(RM) $*.o

.F.a:
	$(FC) -c $(FFLAGS) $<
	$(AR) $(ARFLAGS) $@ $*.o
	$(RM) $*.o

.f77.a:
	$(FC) -c $(FFLAGS) $<
	$(AR) $(ARFLAGS) $@ $*.o
	$(RM) $*.o

clean:
	$(RM) $(FIRST) $(OBJS) $(MPI) $(DIP) *.mod core *.o

tests:
	cd ../TESTS && $(make)
