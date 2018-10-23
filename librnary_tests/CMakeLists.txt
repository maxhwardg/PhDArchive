cmake_minimum_required(VERSION 3.5)
project(librnary_tests)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall")

add_subdirectory(lib/googletest)
enable_testing()


file(GLOB SOURCE_FILES
        "include/*.hpp"
        "src/*.cpp"
        )

set(PROJECT_INCLUDE_DIRS
        "${PROJECT_SOURCE_DIR}/lib/googletest/include" "${PROJECT_SOURCE_DIR}/lib/googletest"
        "${PROJECT_SOURCE_DIR}/include")


add_executable(${PROJECT_NAME} ${SOURCE_FILES})

target_link_libraries(${PROJECT_NAME} gtest gtest_main)
target_link_libraries(${PROJECT_NAME} librnary)
target_include_directories(${PROJECT_NAME} PRIVATE "${PROJECT_INCLUDE_DIRS}")