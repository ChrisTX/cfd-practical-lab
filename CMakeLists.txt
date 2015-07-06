cmake_minimum_required(VERSION 3.0)
project(CFD_LAB)

set(CMAKE_CXX_FLAGS "-std=c++14")
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
	set(CMAKE_CXX_FLAGS "-std=c++14 -stdlib=libc++")
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -ffast-math")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")

find_library(CXXABI_LIB "c++abi")

add_executable(main-bfs main-bfs.cpp)
add_executable(main-ldc main-ldc.cpp)
add_executable(main-pipe main-pipe.cpp)
add_executable(main-karman main-karman.cpp)

find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})

add_executable(main-mpi main-mpi.cpp)
target_link_libraries(main-mpi ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(main-mpi PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(main-mpi PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
	target_link_libraries(main-bfs ${CXXABI_LIB})
	target_link_libraries(main-ldc ${CXXABI_LIB})
	target_link_libraries(main-pipe ${CXXABI_LIB})
	target_link_libraries(main-karman ${CXXABI_LIB})
	target_link_libraries(main-mpi ${CXXABI_LIB})
endif()