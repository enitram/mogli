########################################################################################################################
#     mogli - molecular graph library                                                                                  #
#                                                                                                                      #
#     Copyright (C) 2016-2019  Martin S. Engler                                                                        #
#                                                                                                                      #
#     This program is free software: you can redistribute it and/or modify                                             #
#     it under the terms of the GNU Lesser General Public License as published                                         #
#     by the Free Software Foundation, either version 3 of the License, or                                             #
#     (at your option) any later version.                                                                              #
#                                                                                                                      #
#     This program is distributed in the hope that it will be useful,                                                  #
#     but WITHOUT ANY WARRANTY; without even the implied warranty of                                                   #
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                                                     #
#     GNU General Public License for more details.                                                                     #
#                                                                                                                      #
#     You should have received a copy of the GNU Lesser General Public License                                         #
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.                                           #
########################################################################################################################

cmake_minimum_required(VERSION 3.12)
project(mogli)

set(CMAKE_CXX_STANDARD 17)

include(FetchContent)
include(ExternalProject)

###############
# Build options
###############

if(NOT DEFINED CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

OPTION(BUILD_PYTHON "Build python module." ON)
OPTION(BUILD_OFFLINE "Build without updating the external projects." OFF)
OPTION(BUILD_TESTS "Build unit tests" OFF)

if (BUILD_OFFLINE)
    set(OFFLINE 1)
else()
    set(OFFLINE 0)
endif()

###############
# Compile flags
###############

set(CMAKE_CXX_FLAGS_DEBUG "-g -ggdb -O0 -DDEBUG -std=c++17")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -std=c++17")

set(CMAKE_C_FLAGS_DEBUG "-g -ggdb -O0 -DDEBUG -std=c99")
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -std=c99")

if( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++ " )
endif()

#######################
# Dependency: msgpack-c
#######################

FetchContent_Declare(
        msgpack
        GIT_REPOSITORY https://github.com/msgpack/msgpack-c.git
        GIT_TAG        cpp-3.2.0
)

FetchContent_GetProperties(msgpack)
if(NOT msgpack_POPULATED)
    FetchContent_Populate(msgpack)
    add_subdirectory(${msgpack_SOURCE_DIR} ${msgpack_BINARY_DIR})
endif()

###################
# Dependency: lemon
###################

ExternalProject_Add(
        lemon
        PREFIX     external
        URL        http://lemon.cs.elte.hu/pub/sources/lemon-1.3.1.tar.gz
        URL_HASH   MD5=e89f887559113b68657eca67cf3329b5
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=lemon
        UPDATE_DISCONNECTED ${OFFLINE}
)

set(LIBLEMON_ROOT ${CMAKE_CURRENT_BINARY_DIR}/external/src/lemon-build/lemon/)

#######
# Files
#######

set( CommonLibs
        m
        emon
)

set(HEADER_FILES
        include/bronkerbosch.h
        include/canonization.h
        include/fcanonization.h
        include/fragment.h
        include/match.h
        include/mcf.h
        include/molecule.h
        include/orbits.h
        include/periodictable.h
        include/product.h
        include/subgraph_isomorphism.h
        include/types.h
        include/util/packing.h
        )

set(SOURCE_FILES
        src/bronkerbosch.cpp
        src/canonization.cpp
        src/fcanonization.cpp
        src/fragment.cpp
        src/match.cpp
        src/mcf.cpp
        src/molecule.cpp
        src/orbits.cpp
        src/product.cpp
        src/subgraph_isomorphism.cpp
        )

set(TEST_FILES
        tests/test_isomorphism.cpp
        tests/test_matching.cpp
        tests/test_packing.cpp
        tests/test_fixtures.h
        tests/test_io.cpp)

set(NAUTY_SRC
        lib/nauty/naugraph.c
        lib/nauty/naugroup.c
        lib/nauty/naurng.c
        lib/nauty/nautil.c
        lib/nauty/nauty.c
        lib/nauty/schreier.c
        )

set(LAD_SRC
        lib/lad/allDiff.c
        lib/lad/compatibility.c
        lib/lad/domains.c
        lib/lad/graph.c
        lib/lad/lad.c
        lib/lad/sublad.c
        lib/lad/sublad.h
        )

#############
# C++ library
#############

include_directories(
        "include"
        "lib/dynamic_bitset"
        "lib/lad"
        "lib/nauty"
        "${LIBLEMON_ROOT}/include"
        "${msgpack_SOURCE_DIR}/include"
)
link_directories( "${LIBLEMON_ROOT}/lib")

add_library(mogli SHARED ${HEADER_FILES} ${SOURCE_FILES} ${NAUTY_SRC} ${LAD_SRC})
target_link_libraries(mogli ${CommonLibs})
target_compile_definitions(mogli PRIVATE TR2="${USE_TR2}")
add_dependencies(mogli lemon)

add_executable(report tests/report.cpp ${SOURCE_FILES} ${NAUTY_SRC} ${LAD_SRC})
target_link_libraries(report ${CommonLibs})
add_dependencies(report lemon)

#add_executable(bbpmcs_compare tests/bbpmcs_compare.cpp ${SOURCE_FILES} ${NAUTY_SRC} ${LAD_SRC})
#target_link_libraries(bbpmcs_compare ${Boost_LIBRARIES} ${CommonLibs})
#
#add_executable(test_product tests/test_product.cpp ${SOURCE_FILES} ${NAUTY_SRC} ${LAD_SRC})
#target_link_libraries(test_product ${Boost_LIBRARIES} ${CommonLibs})
#
#add_executable(test_all_combinations tests/test_all_combinations.cpp ${SOURCE_FILES} ${NAUTY_SRC} ${LAD_SRC})
#target_link_libraries(test_all_combinations ${Boost_LIBRARIES} ${CommonLibs})
#
#add_executable(mcis tests/mcis.cpp ${SOURCE_FILES} ${NAUTY_SRC} ${LAD_SRC})
#target_link_libraries(mcis ${Boost_LIBRARIES} ${CommonLibs})
#
#add_executable(timeout tests/timeout.cpp ${SOURCE_FILES} ${NAUTY_SRC} ${LAD_SRC})
#target_link_libraries(timeout ${Boost_LIBRARIES} ${CommonLibs})

#################
# Python bindings
#################

if (BUILD_PYTHON)

    FetchContent_Declare(
            pybind11
            GIT_REPOSITORY https://github.com/pybind/pybind11
            GIT_TAG        v2.3.0
    )

    FetchContent_GetProperties(pybind11)
    if(NOT pybind11_POPULATED)
        FetchContent_Populate(pybind11)
        add_subdirectory(${pybind11_SOURCE_DIR} ${pybind11_BINARY_DIR})
    endif()

    pybind11_add_module(pymogli src/util/pymogli.cpp)
    target_link_libraries(pymogli PRIVATE mogli)
    add_dependencies(pymogli lemon)

endif()

############
# Unit tests
############

if (BUILD_TESTS)

    include(FetchContent)
    FetchContent_Declare(
            catch2
            GIT_REPOSITORY https://github.com/catchorg/Catch2.git
            GIT_TAG        v2.8.0
    )

    FetchContent_GetProperties(catch2)
    if(NOT catch2_POPULATED)
        FetchContent_Populate(catch2)
        add_subdirectory(${catch2_SOURCE_DIR} ${catch2_BINARY_DIR})
    endif()

    include(CTest)
    enable_testing()

    add_executable(unit_tests tests/test_main.cpp ${TEST_FILES})
    target_link_libraries(unit_tests PRIVATE mogli Catch2::Catch2)
    add_dependencies(unit_tests lemon)

    add_test(NAME test_c COMMAND unit_tests)

#    if (BUILD_PYTHON)
#
#        find_package(Python3 COMPONENTS Interpreter)
#
#        add_dependencies(unit_tests pymogli)
#        add_test(NAME test_py COMMAND Python3::Interpreter -m unittest ${PROJECT_SOURCE_DIR}/tests/test_python.py)
#
#
#    endif()

endif()

# TODO python unit tests
# TODO documentation
# TODO check if all methods (except methods with stream parameters) are exposed to python
