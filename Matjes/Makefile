.PHONY: serial clean cleanall 

serial: 
	+@(cd src && make $*)

install:  serial
	+@(cd src && mv Matjes ../bin/.)
clean: 
	+@(cd src && make clean)

cleanall: clean
	+@(cd src && make clean)
	+@(cd bin && rm *)
