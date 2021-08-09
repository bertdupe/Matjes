include(src/CMakeLists.txt)   #read all the source-files

set_source_files_properties(${main_f90} PROPERTIES Fortran_FORMAT FREE)


add_executable(Matjes ${main_f90})
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

if(USE_EIGEN)
    add_library(eigen_cpp ${eigen_cpp})
    add_library(eigen_H ${eigen_H})
    target_link_libraries (eigen_cpp Eigen3::Eigen)
    target_link_libraries (eigen_H Eigen3::Eigen)
    target_link_libraries(Matjes eigen_cpp eigen_H)
endif()

if(USE_CUDA)
    add_library(cuda_fft ${cuda_fft})
    add_library(cuda_H ${cuda_H})
    target_link_libraries(cuda_fft cufft cudart)
    target_link_libraries(cuda_H cusparse cudart)
    target_link_libraries(Matjes cuda_fft cuda_H)
endif()

if(USE_BLAS)
    target_link_libraries(Matjes BLAS::BLAS)
endif()

if(USE_LAPACK)
    target_link_libraries(Matjes LAPACK::LAPACK)
endif()

if(DEFINED manual_linker_language)
    set_target_properties(Matjes PROPERTIES LINKER_LANGUAGE ${manual_linker_language})
endif()
