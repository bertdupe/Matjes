#this file tries to find and check libraries to set appropriate environment variables automatically

#check for CUDA and set parameters if used
message("\n\n Search for CUDA using environment CUDA_PATH")
if(( USE_CUDA OR NOT DEFINED USE_CUDA) AND DEFINED ENV{CUDA_PATH})
    message("FOUND CUDA and compiling with it.\n\n" )
    set(USE_CUDA TRUE)
    enable_language (CUDA)
    add_compile_definitions(CPP_CUDA)
    set(CUDA_PATH $ENV{CUDA_PATH})
    link_directories( ${CUDA_PATH}/lib64 )
    include_directories( ${CUDA_PATH}/include )
    if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        set(CMAKE_CUDA_ARCHITECTURES 75)
        message("CMake variable 'CMAKE_CUDA_ARCHITECTURES' not set, using default value: ${CMAKE_CUDA_ARCHITECTURES}")
    endif()
elseif(USE_CUDA AND NOT DEFINED ENV{CUDA_PATH})
    message( FATAL_ERROR "Enforcing to use cuda, but CUDA_PATH is not defined.\n Export the CUDA_PATH environment variable or remove USE_CUDA in config.cmake\n\n" )
else()
    message( "Skipping CUDA-part of compilation\n\n" )
endif()


#check for EIGEN and set parameters if used
message("\n\n Search for EIGEN")
find_package (Eigen3 3.3 NO_MODULE)
if((USE_EIGEN OR NOT DEFINED USE_EIGEN) AND Eigen3_FOUND)
    message("FOUND Eigen library and compiling with it.\n\n" )
    add_compile_definitions(CPP_EIGEN)
    set(USE_EIGEN TRUE)
elseif(USE_EIGEN AND NOT Eigen3_FOUND)
	message( FATAL_ERROR "Enforcing to use EIGEN, but library can not be found.\n\n" )
else()
    message( "Skipping Eigen-part of compilation\n\n" )
endif()


#check for MKL
message("\n\n Searching for MKL (sparse)")
if(NOT DEFINED USE_MKL OR (DEFINED USE_MKL AND USE_MKL))
    if(DEFINED MKL_library_path)
	    message(" Using manually set MKL library path: ${MKL_library_path}")
    else()
        set(MKL_library_path $ENV{MKLROOT}/lib/intel64)
	message(" Using default MKL library path: ${MKL_library_path}")
    endif()

    if(DEFINED MKL_include_path)
	    message(" Using manually set MKL include path: ${MKL_include_path}")
    else()
        set(MKL_include_path $ENV{MKLROOT}/include)
	    message(" Using default MKL include path: ${MKL_include_path}")
    endif()

    if(DEFINED MKL_linker)
    	message(" Using manually set MKL linker: ${MKL_linker}")
    else()
        if(USE_OPENMP)
            set(MKL_linker "-Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl") 
        else()
            set(MKL_linker "-Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl") 
        endif()
	    message(" Using default MKL linker: ${MKL_linker}")
    endif()

#check sparse MKL
    try_compile(FOUND_mkl "${CMAKE_BINARY_DIR}/temp" "${CMAKE_SOURCE_DIR}/cmake/tests/spblas_mkl.f90"
        CMAKE_FLAGS
                  "-DINCLUDE_DIRECTORIES=${MKL_include_path}"
                  "-DLINK_DIRECTORIES=${MKL_library_path}"
        LINK_LIBRARIES "${MKL_linker}"
        OUTPUT_VARIABLE MKL_test_output
        )

    if(FOUND_mkl)
        message(" Success testing MKL (sparse).\n Compiling with MKL_Sparse")
        add_compile_definitions(CPP_MKL)
        set(USE_MKL TRUE)
    elseif(USE_MKL AND NOT FOUND_mkl)
        message("Error information testing MKL (sparse):")
        message("${MKL_test_output}")
        message("\n")
        message( FATAL_ERROR "Unsuccessfull MKL (sparse) test.\n USE_MKL has been set to TRUE, thus aborting.")
    else()
        message("Unsuccessfull MKL (sparse) test.\n Compiling without MKL (sparse).")
        set(USE_MKL FALSE)
    endif()

#check FFTW3 from mkl
    if(USE_MKL AND (NOT DEFINED USE_FFTW OR USE_FFTW))
        try_compile(FOUND_FFTW "${CMAKE_BINARY_DIR}/temp" "${CMAKE_SOURCE_DIR}/cmake/tests/fftw.f90"
            CMAKE_FLAGS
                      "-DINCLUDE_DIRECTORIES=${MKL_include_path}"
                      "-DLINK_DIRECTORIES=${MKL_library_path}"
            LINK_LIBRARIES "${MKL_linker}"
            OUTPUT_VARIABLE FFTW_test_output
            )
        if(FOUND_FFTW)
            message(" Success testing FFTW with mkl.\n Disable other search for fftw3 library implementation.\n Compiling with CPP_FFTW3")
            add_compile_definitions(CPP_FFTW3)
            set(USE_FFTW FALSE)
        endif()

        try_compile(FFTW_threaded "${CMAKE_BINARY_DIR}/temp" "${CMAKE_SOURCE_DIR}/cmake/tests/fftw_thread.f90"
            CMAKE_FLAGS
                      "-DINCLUDE_DIRECTORIES=${MKL_include_path}"
                      "-DLINK_DIRECTORIES=${MKL_library_path}"
            LINK_LIBRARIES "${MKL_linker}"
            )
        if(FFTW_threaded AND FOUND_FFTW)
            message(" Success testing FFTW with threading .\n Compiling with CPP_FFTW3_THREAD")
            add_compile_definitions(CPP_FFTW3_THREAD)
        elseif(NOT FFTW_threaded AND FOUND_FFTW)
            message(" Failure to compile FFTW with threading, but normal fftw should work")
        endif()
    endif()

else()
    message("Found USE_MKL = FALSE, skipping MKL (sparse) test.")
endif()

#check for FFTW3
message("\n\n Searching for FFTW3")
if(NOT DEFINED USE_FFTW OR (DEFINED USE_FFTW AND USE_FFTW))

    if(DEFINED FFTW_library_path)
	    message(" Using manually set FFTW library path: ${FFTW_library_path}")
    else()
        set(FFTW_library_path "/usr/lib")
	     message(" Using default FFTW library path: ${FFTW_library_path}")
    endif()

    if(DEFINED FFTW_include_path)
	    message(" Using manually set FFTW include path: ${FFTW_include_path}")
    else()
        set(FFTW_include_path "/usr/include")
	    message(" Using default FFTW include path: ${FFTW_include_path}")
    endif()

    if(DEFINED FFTW_linker)
    	message(" Using manually set FFTW linker: ${FFTW_linker}")
    else()
        set(FFTW_linker "-lfftw3") 
    	message(" Using default FFTW linker: ${FFTW_linker}")
    endif()

    try_compile(FOUND_FFTW "${CMAKE_BINARY_DIR}/temp" "${CMAKE_SOURCE_DIR}/cmake/tests/fftw.f90"
       CMAKE_FLAGS
                 "-DINCLUDE_DIRECTORIES=${FFTW_include_path}"
                 "-DLINK_DIRECTORIES=${FFTW_library_path}"
        LINK_LIBRARIES "${FFTW_linker}"
        OUTPUT_VARIABLE FFTW_test_output
        )
    if(FOUND_FFTW)
        message(" Success testing FFTW.\n Compiling with CPP_FFTW3")
        add_compile_definitions(CPP_FFTW3)
        set(USE_FFTW TRUE)
    elseif(USE_FFTW AND NOT FOUND_FFTW)
        message("Error information testing FFTW:")
        message("${FFTW_test_output}")
        message("\n")
        message( FATAL_ERROR "Unsuccessfull FFTW test.\n USE_FFTW has been set to TRUE, thus aborting.")
    else()
        message("Unsuccessfull FFTW test.\n Compiling without FFTW.")
        set(USE_FFTW FALSE)
    endif()

    try_compile(FFTW_threaded "${CMAKE_BINARY_DIR}/temp" "${CMAKE_SOURCE_DIR}/cmake/tests/fftw_thread.f90"
        LINK_LIBRARIES "${FFTW_linker}"
        OUTPUT_VARIABLE FFTW_test_output
        )
    if(FFTW_threaded AND USE_FFTW)
        message(" Success testing FFTW with threading .\n Compiling with CPP_FFTW3_THREAD")
        add_compile_definitions(CPP_FFTW3_THREAD)
    elseif(NOT FFTW_threaded AND USE_FFTW)
        message(" Failure to compile FFTW with threading, but normal fftw should work")
    endif()
elseif(FOUND_FFTW)
    message("FFTW implementation already found, skipping search of explicit FFTW implementation.")
else(NOT FOUND_FFTW AND NOT USE_FFTW)
    message("Skipping FFTW as USE_FFTW=FALSE.")
endif()


#check for FFTW3_MPI
message("\n\n Searching for FFTW3_MPI")
if(NOT DEFINED USE_FFTWMPI OR (DEFINED USE_FFTWMPI AND USE_FFTWMPI))

    if(DEFINED FFTWMPI_library_path)
            message(" Using manually set FFTW_MPI library path: ${FFTWMPI_library_path}")
    else()
        set(FFTWMPI_library_path "/usr/lib")
             message(" Using default FFTW library path: ${FFTWMPI_library_path}")
    endif()

    if(DEFINED FFTWMPI_include_path)
            message(" Using manually set FFTW_MPI include path: ${FFTWMPI_include_path}")
    else()
        set(FFTWMPI_include_path "/usr/include")
            message(" Using default FFTW include path: ${FFTWMPI_include_path}")
    endif()

    if(DEFINED FFTWMPI_linker)
        message(" Using manually set FFTW_MPI linker: ${FFTWMPI_linker}")
    else()
        set(FFTWMPI_linker "-lfftw3 -lfftw3_mpi -lm")
        message(" Using default FFTW_MPI linker: ${FFTWMPI_linker}")
    endif()

    try_compile(FOUND_FFTWMPI "${CMAKE_BINARY_DIR}/temp" "${CMAKE_SOURCE_DIR}/cmake/tests/fftw_mpi.f90"
       CMAKE_FLAGS
                 "-DINCLUDE_DIRECTORIES=${FFTW_include_path}"
                 "-DLINK_DIRECTORIES=${FFTW_library_path}"
        LINK_LIBRARIES "${FFTW_linker}"
        OUTPUT_VARIABLE FFTW_test_output
        )
    if(FOUND_FFTWMPI)
        message(" Success testing FFTW_MPI.\n Compiling with CPP_FFTWMPI")
        add_compile_definitions(CPP_FFTWMPI)
        set(USE_FFTWMPI TRUE)
    elseif(USE_FFTWMPI AND NOT FOUND_FFTWMPI)
        message("Error information testing FFTW_MPI:")
        message("${FFTW_test_output}")
        message("\n")
        message( FATAL_ERROR "Unsuccessfull FFTW_MPI test.\n USE_FFTW_MPI has been set to TRUE, thus aborting.")
    else()
        message("Unsuccessfull FFTW_MPI test.\n Compiling without FFTWMPI.")
        set(USE_FFTWMPI FALSE)
    endif()

elseif(FOUND_FFTWMPI)
    message("FFTW_MPI implementation already found, skipping search of explicit FFTW_MPI implementation.")
else(NOT FOUND_FFTWMPI AND NOT USE_FFTWMPI)
    message("Skipping FFTW_MPI as USE_FFTWMPI=FALSE.")
endif()




#check for BLAS
message("\n\n Search for BLAS")
find_package(BLAS)
if((USE_BLAS OR NOT DEFINED USE_BLAS) AND BLAS_FOUND)
    if(USE_MKL)
        message("FOUND BLAS library, but since MKL is already found that is used.\n\n" )
        set(USE_BLAS FALSE)
    else()
        message("FOUND BLAS library and using that implementation\n\n" )
        set(USE_BLAS TRUE)
    endif()
    add_compile_definitions(CPP_BLAS)
elseif(USE_BLAS AND NOT BLAS_FOUND AND NOT USE_MKL)
    message( FATAL_ERROR "Enforcing to use BLAS, but library can not be found.\n\n" )
else()
    message( "Skipping BLAS-part of compilation\n\n" )
endif()


#check for lapack
message("\n\n Search for LAPACK")
find_package(LAPACK)
if((USE_LAPACK OR NOT DEFINED USE_LAPACK) AND LAPACK_FOUND)
    if(USE_MKL)
        message("FOUND LAPACK library, but since MKL is already found that is used.\n\n" )
        set(USE_LAPACK FALSE)
    else()
        message("FOUND LAPACK library and using that implementation\n\n" )
        set(USE_LAPACK TRUE)
    endif()
    add_compile_definitions(CPP_LAPACK)
elseif(USE_LAPACK AND NOT LAPACK_FOUND AND NOT USE_MKL)
    message( FATAL_ERROR "Enforcing to use LAPACK, but library can not be found.\n\n" )
else()
    message( "Skipping LAPACK-part of compilation\n\n" )
endif()


#check for netcdf
message("\n\n Search for NETCDF")
if(NOT DEFINED USE_netCDF OR (DEFINED USE_netCDF AND USE_netCDF))
    if(DEFINED netCDF_library_path)
        message(" Using manually set netCDF library path: ${netCDF_library_path}")
    else()
        set(netCDF_library_path "/usr/lib")
         message(" Using default netCDF library path: ${netCDF_library_path}")
    endif()
    
    if(DEFINED netCDF_include_path)
        message(" Using manually set netCDF include path: ${netCDF_include_path}")
    else()
        set(netCDF_include_path "/usr/include")
        message(" Using default netCDF include path: ${netCDF_include_path}")
    endif()
    
    if(DEFINED netCDF_linker)
    	message(" Using manually set netCDF linker: ${netCDF_linker}")
    else()
        set(netCDF_linker "-lnetcdff") 
    	message(" Using default netCDF linker: ${netCDF_linker}")
    endif()
    
    try_compile(FOUND_netCDF "${CMAKE_BINARY_DIR}/temp" "${CMAKE_SOURCE_DIR}/cmake/tests/netcdf.f90"
        CMAKE_FLAGS
                  "-DINCLUDE_DIRECTORIES=${netCDF_include_path}"
                  "-DLINK_DIRECTORIES=${netCDF_library_path}"
        LINK_LIBRARIES "${netCDF_linker}"
        OUTPUT_VARIABLE netCDF_test_output)
    if((USE_netCDF OR NOT DEFINED USE_netCDF) AND FOUND_netCDF)
        message(" Success testing netCDF.\n")
        add_compile_definitions(CPP_NETCDF)
        set(USE_netCDF TRUE)
    elseif(USE_netCDF AND NOT FOUND_netCDF)
        message("Error output testing netcdf:")
        message("${netCDF_test_output}")
        message( FATAL_ERROR "Enforcing to use netCDF, but library can not be found.\nSet netCDF_library_path netCDF_include_path, netCDF_linker cmake variable in config.cmake with respecitive directory and library linker.\n\n")
    else()
        message( "Failed netcdf test.\nSkipping netCDF-part of compilation\n\n" )
    endif()
else()
    message("Skipping netCDF as USE_netCDF=FALSE.")
endif()

#check for gsl library
message("\n\n Search for GSL")
if(USE_GSL)
    message("using external GSL library\n\n" )
    add_compile_definitions(CPP_GSL)
    if(DEFINED GSL_library_path)
	message(" Using manually set GSL library path: ${GSL_library_path}")
    else()
	set(GSL_library_path "/usr/local/lib")
	message(" Using default GSL library path: ${GSL_library_path}")
    endif()

    if(DEFINED GSL_include_path)
        message(" Using manually set GSL include path: ${GSL_include_path}")
    else()
	set(GSL_include_path "/usr/local/include")
	message(" Using default GSL include path: ${GSL_include_path}")
    endif()

    if(DEFINED GSL_linker)
        message(" Using manually set GSL linker command: ${GSL_linker}")
    else()
        set(GSL_linker "-lgsl -lgslcblas")
        message(" Using default GSL linker command: ${GSL_linker}")
    endif()
else()
    message( "Skipping GSL-part of compilation\n\n" )
endif()
