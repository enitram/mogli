cmake_minimum_required(VERSION 3.6)
project(mogli)

set( LIBLEMON_ROOT "$ENV{HOME}/lemon" CACHE PATH "Additional search directory for lemon library" )

IF(NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: Debug Release" FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

set(CMAKE_CXX_FLAGS_DEBUG "-g -ggdb -O0 -DDEBUG -std=c++11")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -std=c++11")

set(CMAKE_C_FLAGS_DEBUG "-g -ggdb -O0 -DDEBUG -std=c99")
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -std=c99")

find_package( Boost 1.60.0 REQUIRED COMPONENTS python)
set( Boost_USE_STATIC_LIBS OFF )
set( Boost_USE_MULTITHREADED ON )
set( Boost_USE_STATIC_RUNTIME OFF )
include_directories("${Boost_INCLUDE_DIRS}" "/usr/include/python2.7")

if( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++ " )
endif()

set( CommonLibs
        m
        emon
)

set(HEADER_FILES
        include/canonization.h
        include/molecule.h
        include/periodictable.h
        include/product.h
        include/bronkerbosch.h
        include/orbits.h
        include/fragment.h
        include/subgraph_isomorphism.h
        include/types.h
        include/util/packing.h
        include/mcf.h
        include/match.h
        include/fcanonization.h)

set(SOURCE_FILES
        src/canonization.cpp
        src/molecule.cpp
        src/product.cpp
        src/bronkerbosch.cpp
        src/orbits.cpp
        src/subgraph_isomorphism.cpp
        src/mcf.cpp)

set(NAUTY_SRC
        lib/nauty/nauty.c
        lib/nauty/nautil.c
        lib/nauty/naugraph.c
        lib/nauty/schreier.c
        lib/nauty/naurng.c
        lib/nauty/naugroup.c)

set(LAD_SRC
        lib/lad/allDiff.c
        lib/lad/domains.c
        lib/lad/graph.c
        lib/lad/lad.c
        lib/lad/sublad.c
        lib/lad/sublad.h
        lib/lad/compatibility.c)

include_directories(
        "include"
        "${LIBLEMON_ROOT}/include"
        "lib/nauty"
        "lib/lad"
        "lib/msgpack-c/include"
)
link_directories( "${LIBLEMON_ROOT}/lib")

add_library(mogli SHARED ${HEADER_FILES} ${SOURCE_FILES} src/util/boosting.cpp ${NAUTY_SRC} ${LAD_SRC})
target_link_libraries(mogli ${Boost_LIBRARIES} ${CommonLibs})

add_executable(report test/report.cpp ${SOURCE_FILES} ${NAUTY_SRC} ${LAD_SRC})
target_link_libraries(report ${Boost_LIBRARIES} ${CommonLibs})
