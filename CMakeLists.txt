cmake_minimum_required(VERSION 3.6)
project(mogli)

set( LIBLEMON_ROOT "$ENV{HOME}/lemon" CACHE PATH "Additional search directory for lemon library" )

IF(NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build, options are: Debug Release" FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

set(CMAKE_CXX_FLAGS_DEBUG "-g -ggdb -O0 -DDEBUG -std=c++11")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -std=c++11")

if( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++ " )
endif()

set( CommonLibs
        m
        emon
)

set(SOURCE_FILES
        src/canonization.h
        src/canonization.cpp
        src/isomorphism.h
        src/isomorphism.cpp
        src/molecule.h
        src/iacm.h
        src/product.h)

set(NAUTY_SRC
        include/nauty/nauty.c
        include/nauty/nautil.c
        include/nauty/naugraph.c
        include/nauty/schreier.c
        include/nauty/naurng.c)

include_directories( "src" "${LIBLEMON_ROOT}/include" "include/nauty" )
link_directories( "${LIBLEMON_ROOT}/lib" )

add_library(mogli ${SOURCE_FILES} ${NAUTY_SRC})
target_link_libraries(mogli ${CommonLibs} )
