SHELL = /bin/bash

.PHONY: serial clean cleanall

serial: copy
	cd build && make Matjes

copy:
	cp src/*.f90 src/*/*.f* src/*/*.F* build/.

install:  clean serial
	+@(cd build && mv Matjes ../bin/.)
	
clean: 
	+@(cd build && make clean)

cleanall: clean
	+@(cd build && make clean)
	+@(cd bin && rm *)
