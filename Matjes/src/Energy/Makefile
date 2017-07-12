include ../../make.sys

.SUFFIXES: .F .f .f90 .F90

OBJS = cal_local_energy.o CalculateAverages.o demag-matrix.o demag-tensor.o DeriveValue.o Energy_density.o Energy.o

toto: $(OBJS)
	$(LD) -o $@ $(LDFLAGS) $(OBJS) $(LIBRARIES)
	
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
	$(RM) $(OBJS) *.mod core *.o

Energy.o:
cal_local_energy.o: Energy.o