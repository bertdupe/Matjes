.PHONY: serial clean cleanall toto


toto:
	@(cd src && make $*)
	
serial: 
	+@(cd src && make $*)

install:  serial
	+@(cd src && mv Matjes ../bin/.)
clean: 
	+@(cd src && make clean)

cleanall: clean
	+@(cd src && make clean)
	+@(cd bin && rm *)
