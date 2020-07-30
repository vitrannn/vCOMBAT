cmake_minimum_required (VERSION 3.2)
project (TuberculosisSimulation)

# Version numbers
set (TuberculosisSimulation_VERSION_MAJOR 0)
set (TuberculosisSimulation_VERSION_MINOR 1)

# Local CMAKE modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
set(CMAKE_SYSTEM_PREFIX_PATH ${CMAKE_SYSTEM_PREFIX_PATH} /usr)
# configure a header file to pass some of the CMake settings to the source code
configure_file (
  "src/config.h.in"
  "${PROJECT_BINARY_DIR}/tuberculosis_simulation_config.h"
) 

#list of sources
set(sources src/main.c src/base_simulation.c src/full_model.c ${CMAKE_CURRENT_LIST_DIR}/arg_parser/carg_parser.c)

#Use release-level optimization
set(CMAKE_BUILD_TYPE Release)

find_package(Doxygen)
if(DOXYGEN_FOUND)
	configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
	add_custom_target(doc
		${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
		WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
		COMMENT "Generating API documentation with Doxygen" VERBATIM)
endif(DOXYGEN_FOUND)

# Include the GSL package
find_package(GSL REQUIRED)
set(LIBS ${LIBS} ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})
include_directories(${GSL_INCLUDE_DIRS} ${GSLCBLAS_INCLUDE_DIRS})
include_directories(${PROJECT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_LIST_DIR}/arg_parser)

# Include the math library
find_library(M_LIB m)
set(LIBS ${LIBS} ${M_LIB})

# Include the libYAML package
find_package(LIBYAML REQUIRED)
set(LIBS ${LIBS} ${LIBYAML_LIBRARIES})
include_directories(${LIBYAML_INCLUDE_DIRS})

# Output the executable file
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
add_executable(tuberculosis_simulation ${sources})
target_link_libraries(tuberculosis_simulation ${LIBS})
