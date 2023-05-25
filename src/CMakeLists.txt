cmake_minimum_required(VERSION 3.0)
project(simpleFEM)

set(CMAKE_CXX_COMPILER g++)
set(CMAKE_CXX_STANDARD 17)
SET(CMAKE_CXX_FLAGS_DEBUG "-g")
#find_package (Eigen3 3.4 REQUIRED NO_MODULE)

add_executable(${PROJECT_NAME} main.cpp)
add_subdirectory(core)

target_include_directories(${PROJECT_NAME} PUBLIC core)
target_link_directories(${PROJECT_NAME} PRIVATE core)
target_link_libraries(${PROJECT_NAME} CORE)
