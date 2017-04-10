include= ../make.sys

OBJ= vectors.o

clean:
	@(RM) *.o *.mod

.SUFFIXES: .o .f90 .c

%.o: %.f90
	$(F90) $(FFLAGS) $(INCS) -c $.f90
%.o: %.c
	$(CC) $(CFLAGS) $(INCS) -c $.c

include dependencies