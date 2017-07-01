include ../../make.sys

.SUFFIXES: .o .f90 .c

OFOLDER= OBJECT-FILES

lib: 
	@( cd $(OFOLDER) && make -f ../SRC/Makefile)

clean:
	@( cd $(OFOLDER) && $(RM) *.o *.mod *.a)