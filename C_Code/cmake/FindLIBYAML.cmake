# CMake module to search for the libyaml library
# (library for parsing YAML files)
# If it's found it sets LIBYAML_FOUND to TRUE
# and following variables are set:
#    LIBYAML_INCLUDE_DIR
#    LIBYAML_LIBRARY

set(LIBYAML_DIR $ENV{LIBYAML_DIR})

if ($LIBYAML_DIR)
    set(LIBYAML "${LIBYAML_DIR}" CACHE PATH "Path to search for LIBYAML include and library files")
endif()

if ($LIBYAML_DIR)
    message(" Searching Yaml in ${LIBYAML_DIR}")
    find_path(LIBYAML_INCLUDE_DIRS 
                NAMES yaml.h
                HINTS "${LIBYAML_DIR}/include/yaml-cpp/")
                
    find_library(LIBYAML_LIBRARIES 
                NAMES yaml libyaml
                HINTS "${LIBYAML_DIR}/lib/")
else()
    FIND_PATH(LIBYAML_INCLUDE_DIRS NAMES yaml.h)
    FIND_LIBRARY(LIBYAML_LIBRARIES NAMES yaml libyaml)
endif()

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Yaml DEFAULT_MSG LIBYAML_LIBRARIES LIBYAML_INCLUDE_DIRS)
MARK_AS_ADVANCED(LIBYAML_INCLUDE_DIRS LIBYAML_LIBRARIES)
