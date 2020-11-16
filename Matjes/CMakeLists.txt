cmake_minimum_required(VERSION 3.0)


#Some stuff might be defined in config file
if (EXISTS "${CMAKE_BINARY_DIR}/config.cmake")
   include("${CMAKE_BINARY_DIR}/config.cmake")
endif()


project(Matjes NONE)

enable_language (Fortran CXX)

#if(NOT DEFINED custom_MPICCOMP_FORTRAN)
# find_package(MPI REQUIRED)
# if (MPI_FOUND)
# 	include_directories(SYSTEM ${MPI_INCLUDE_PATH})
# else (MPI_FOUND)
# 	 message(SEND_ERROR "This application cannot compile without MPI")
# endif (MPI_FOUND)
#endif()

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

