cmake_minimum_required(VERSION 3.0)
project(simpleFEM)

set(CMAKE_CXX_COMPILER /usr/bin/mpicxx)
set(CMAKE_CXX_STANDARD 17)
SET(CMAKE_CXX_FLAGS_DEBUG "-g")
#find_package (Eigen3 3.4 REQUIRED NO_MODULE)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})


#set(PETSC_DIR /home/red/gitclones/petsc/arch-python-linux-x86_64)
#set(PETSC_ARCH arch-python-linux-x86_64)
find_package(PkgConfig REQUIRED)
if (PKG_CONFIG_FOUND)
    pkg_check_modules(PETSC PETSc)
endif()
message(STATUS "Found PETSc ${PETSC_VERSION}")
message(STATUS "PETSC Dirs ${PETSC_INCLUDE_DIRS}")
if (PETSC_FOUND)
    list(APPEND COMPILE_OPTIONS ${PETSC_CFLAGS})

    include_directories(${PETSC_INCLUDE_DIRS})
    
   set(LINK_FLAGS "${LINK_FLAGS} ${PETSC_LDFLAGS}")   
    
    list(APPEND LIBRARIES ${PETSC_LINK_LIBRARIES})   
    
    set(CMAKE_REQUIRED_FLAGS ${PETSC_CFLAGS})   
    
    set(CMAKE_REQUIRED_INCLUDES "${PETSC_INCLUDE_DIRS}")
endif()

#include_directories(${PETSC_DIR}/include)
#include_directories(${PETSC_DIR}/arch-python-linux-x86_64/include)


add_subdirectory(src)



