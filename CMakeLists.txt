cmake_minimum_required(VERSION 3.20)
project(tsp)

set(CMAKE_CXX_STANDARD 14)

set(GUROBI_ROOT /Library/gurobi1000/macos_universal2)

message(STATUS "GUROBI_ROOT : ${GUROBI_ROOT}")

set(GUROBI_INCLUDE_DIR "${GUROBI_ROOT}/include")

message(STATUS "GUROBI_INCLUDE_DIR : ${GUROBI_INCLUDE_DIR}")

set(GUROBI_LIBRARY_1 "${GUROBI_ROOT}/lib/libgurobi_c++.a")

set(GUROBI_LIBRARY_2 "${GUROBI_ROOT}/lib/libgurobi100.dylib")

set(GUROBI_LIBRARIES ${GUROBI_LIBRARY_1} ${GUROBI_LIBRARY_2})

message(STATUS "GUROBI_LIBRARIES : ${GUROBI_LIBRARIES}")

include_directories(${GUROBI_INCLUDE_DIR})

include_directories(include)

file(GLOB_RECURSE SRCS src/*.cpp)

message(STATUS "SRCS : ${SRCS}")

add_executable(tsp ${SRCS})

target_link_libraries(tsp ${GUROBI_LIBRARIES})