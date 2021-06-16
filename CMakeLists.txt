cmake_minimum_required(VERSION 3.0)

#Some stuff might be defined in config file
if (EXISTS "${CMAKE_BINARY_DIR}/config.cmake")
    message("\n\nReading config.cmake\n\n")
    include("${CMAKE_BINARY_DIR}/config.cmake")
else()
    message("\n\nMISSING config.cmake file, continue with default parameters.\nThis is probably very restrictive as many compilerflag options are not set by default\n\n")
endif()


project(Matjes NONE)

enable_language (Fortran CXX)

#check which external libraries can be used 
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

