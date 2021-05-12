set(main_f90 "")
set(eigen_cpp "")
include(src/CMakeLists.txt)

set_source_files_properties(${main_f90} PROPERTIES Fortran_FORMAT FREE)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
add_library(eigen_cpp ${eigen_cpp})
add_library(eigen_H ${eigen_H})
add_library(cuda_fft ${cuda_fft})
add_library(cuda_H ${cuda_H})
target_link_libraries (eigen_cpp Eigen3::Eigen)
target_link_libraries (eigen_H Eigen3::Eigen)
#set_target_properties(eigen_cpp PROPERTIES LINKER_LANGUAGE CXX)
#set_target_properties(eigen_H PROPERTIES LINKER_LANGUAGE CXX)

#Serial executables
add_executable(Matjes ${main_f90})
target_link_libraries(Matjes eigen_cpp eigen_H cuda_fft cuda_H)
#set_target_properties(Matjes PROPERTIES LINKER_LANGUAGE Fortran)

target_compile_options(Matjes BEFORE PRIVATE "${COMPADD_only_own}")
if(DEFINED add_lib)
    target_link_libraries(Matjes ${add_lib})
endif()
if(DEFINED add_inc)
    message(${add_inc})
    target_include_directories(Matjes PUBLIC ${add_inc})
endif()
if(DEFINED add_cflags)
    message(${add_cflags})
    target_link_libraries(Matjes ${add_cflags})
endif()
