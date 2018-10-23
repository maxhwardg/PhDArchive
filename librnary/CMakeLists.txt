cmake_minimum_required(VERSION 3.5)
project(librnary)

find_package(Threads)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra -pthread")


set(RNASTRUCTURE_SOURCES
        lib/RNAstructure/rna_library.cpp lib/RNAstructure/rna_library.h lib/RNAstructure/structure.cpp
        lib/RNAstructure/structure.h lib/RNAstructure/stackclass.cpp lib/RNAstructure/stackclass.h
        lib/RNAstructure/stackstruct.cpp lib/RNAstructure/stackstruct.h lib/RNAstructure/algorithm.cpp
        lib/RNAstructure/algorithm.h lib/RNAstructure/arrayclass.cpp lib/RNAstructure/arrayclass.h
        lib/RNAstructure/dotarray.cpp lib/RNAstructure/dotarray.h lib/RNAstructure/forceclass.cpp
        lib/RNAstructure/forceclass.h lib/RNAstructure/TProgressDialog.cpp lib/RNAstructure/TProgressDialog.h
        lib/RNAstructure/defines.h lib/RNAstructure/platform.h)

file(GLOB LIBRNARY_SOURCES
        "include/*.hpp"
        "include/folders/*.hpp"
        "include/models/*.hpp"
        "include/scorers/*.hpp"
        "include/training/*.hpp"
        "src/*.cpp"
        "src/folders/*.cpp"
        "src/models/*.cpp"
        "src/scorers/*.cpp"
        "src/training/*.cpp"
        )

set(SOURCE_FILES "${RNASTRUCTURE_SOURCES}" "${LIBRNARY_SOURCES}")

set(PROJECT_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/include" "${PROJECT_SOURCE_DIR}/lib")

add_library(${PROJECT_NAME} ${SOURCE_FILES})

target_include_directories(${PROJECT_NAME} PUBLIC "${PROJECT_INCLUDE_DIRS}")

target_link_libraries(${PROJECT_NAME} ${CMAKE_THREAD_LIBS_INIT})