
.SUFFIXES: .F .f .f90 .F90

OBJS = Beff_sd.o field-eff_sd.o error_correction_SD.o Efield_sd.o spindynamics.o av_sd.o topocharge_sd.o

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

field-eff_sd.o:
Beff_sd.o: field-eff_sd.o