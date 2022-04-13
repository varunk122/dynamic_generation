cmake_minimum_required (VERSION 3.0)
project (dynamic_generator)
 
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

file(GLOB SOURCES
    ./*.hpp
    ./*.cpp
)

include_directories(${PROJECT_SOURCE_DIR}) 
add_executable (${PROJECT_NAME} ${SOURCES})
target_link_libraries (${PROJECT_NAME} Eigen3::Eigen )
