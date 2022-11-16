#!For finding Boost specify only root directory-DBOOST_ROOT=...Software/boost_1_67_0/

#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16)
project(DataTableExtractor)

find_package(Boost 1.33 REQUIRED)

#----------------------------------------------------------------------------
# Setup include directory for this project
#
include_directories(${Boost_INCLUDE_DIR} ${PROJECT_SOURCE_DIR}/include)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
file(GLOB_RECURSE sources ${PROJECT_SOURCE_DIR}/source/*.cc)
file (GLOB_RECURSE sources1 ${PROJECT_SOURCE_DIR}/source/*.cpp)
list (APPEND sources ${sources1})
file(GLOB_RECURSE headers ${PROJECT_SOURCE_DIR}/include/*.hh)
file(GLOB_RECURSE headers1 ${PROJECT_SOURCE_DIR}/include/*.h)
list (APPEND headers ${headers1})

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

#----------------------------------------------------------------------------
# Set -rpath for binary. Allows it to find libraries in non-stardard locations
# (e.g. not in /usr/) without using LD_LIBRARY_PATH
#
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif("${isSystemDir}" STREQUAL "-1")

#----------------------------------------------------------------------------
# Add the executable, and link it to the libraries
#
add_executable(${PROJECT_NAME} ${sources} ${headers})
target_link_libraries(${PROJECT_NAME} PUBLIC ${Boost_LIBRARIES})
#message(STATUS "Geant4 libraries: ${Geant4_LIBRARIES}")

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build project. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(PROJECT_SCRIPTS
  "")

foreach(_script ${PROJECT_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${CMAKE_BUILD_TYPE}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to Debug/Release directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS ${PROJECT_NAME} DESTINATION ${PROJECT_BINARY_DIR}/${CMAKE_BUILD_TYPE})
