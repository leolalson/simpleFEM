
add_executable(${PROJECT_NAME} main.cpp)

add_subdirectory(core)

target_link_libraries(${PROJECT_NAME} CORE)
target_link_libraries(${PROJECT_NAME} petsc)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})



