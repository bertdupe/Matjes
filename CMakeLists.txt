cmake_minimum_required(VERSION 3.12) #find out actual minimum version

#Some stuff might be defined in config file
if (EXISTS "${CMAKE_BINARY_DIR}/config.cmake")
    message("\n\nReading config.cmake\n\n")
    include("${CMAKE_BINARY_DIR}/config.cmake")
else()
    message("\n\nMISSING config.cmake file, continue with default parameters.\nThis is probably very restrictive as many compilerflag options are not set by default\n\n")
endif()


project(Matjes NONE)

enable_language (Fortran CXX)

if(${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU" AND ${CMAKE_Fortran_COMPILER_VERSION} VERSION_LESS "9.0.0")
    message("\n")
    message(FATAL_ERROR " gfortran compiler with version ${CMAKE_Fortran_COMPILER_VERSION} detected.\n A gfortran version 9.0 or higher is necessary.")
endif()

#check which external libraries can be used and set some environment variables accordingly
include("cmake/external_lib_check.cmake")

#apply program branch options 
include("cmake/programoptions.cmake")

#set compiler options
include("cmake/compilerflags.cmake")

#set libraries and compile
include("cmake/Files_and_Targets.cmake")

install(TARGETS Matjes 
        CONFIGURATIONS Release
		RENAME Matjes)
#        DESTINATION bin)

