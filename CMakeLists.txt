cmake_minimum_required(VERSION 2.8)
project(ratyrace)

# set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/builds)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/libs)
set(CMAKE_CXX_STANDARD 14)
# set(CMAKE_CXX_F)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules")
set(CMAKE_CXX_FLAGS "-Ofast -march=native -ffast-math")
# set(CMAKE_CXX_FLAGS "-pg")

#Headers
include_directories(include)

#Sources
# AUX_SOURCE_DIRECTORY(src/ ALL_SOURCES)
set(SOURCES src/main.cu src/utils.cpp src/flags.cpp src/domainparser.cpp src/raytrace.cpp src/tree.cpp)
set(SFML_DIR "${CMAKE_SOURCE_DIR}/SFML/lib/cmake/SFML")

set(EXEC_NAME raytrace)

# set(SFML_DIR "F:/libs/SFML-2.5.1/lib/cmake/SFML")

# find_package(OpenCL REQUIRED)
# if(OpenCL_FOUND)
#   message(STATUS "OpenCL includes: ${OpenCL_INCLUDE_DIRS}")
#   message(STATUS "OpenCL libraries: ${OpenCL_LIBRARIES}")
#   include_directories(${OpenCL_INCLUDE_DIRS})
# endif()


find_package(SFML 2.5 COMPONENTS graphics window system)
if(SFML_FOUND)
  include_directories(${SFML_INCLUDE_DIR})
endif()

find_package(CUDA REQUIRED)
message("CUDA at ${CUDA_LIBRARIES}")

# if(SFML_FOUND AND OpenCL_FOUND)
#   add_executable(${EXEC_NAME} ${ALL_SOURCES})
#   target_link_libraries(${EXEC_NAME} sfml-window sfml-graphics sfml-system ${OpenCL_LIBRARIES})
# endif()
if(SFML_FOUND)
    cuda_add_executable(${EXEC_NAME} ${SOURCES} OPTIONS -maxrregcount=32)
    set_target_properties( ${EXEC_NAME} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    target_link_libraries(${EXEC_NAME} sfml-window sfml-graphics sfml-system ${CUDA_LIBRARIES})
endif()

# cuda_add_executable(${EXEC_NAME} ${SOURCES} OPTIONS -maxrregcount=32)
# target_link_libraries(${EXEC_NAME} )