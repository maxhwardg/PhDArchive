cmake_minimum_required(VERSION 3.5)
project(programs)


set(INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/include" "${PROJECT_SOURCE_DIR}/lib")

set(LIB_SRC lib/cxxopts.hpp)

SET(PROGRAMS fold_linear fold_logarithmic fold_aalberts fold_avg_asym read_cts energy_linear energy_logarithmic
        energy_aalberts energy_avg_asym)

foreach (program ${PROGRAMS})
    add_executable(${program} src/${program}.cpp ${LIB_SRC})
    target_include_directories(${program} PRIVATE "${INCLUDE_DIRS}")
    target_link_libraries(${program} librnary)
endforeach (program ${PROGRAMS})
