cmake_minimum_required(VERSION 3.10)

project(newtiger LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(HPX REQUIRED)
find_package(SILO REQUIRED)


#set(FLAGS "-Wall -march=native -ffast-math -fno-tree-vectorize -fno-tree-slp-vectorize")
set(FLAGS "-Wall -march=native -ffast-math -fopt-info-vec-all=vector.txt -fopt-info-vec-optimized  -fopt-info-vec-missed -fopt-info-vec-note ")
#set(FLAGS "-Wall -march=native -fsanitize=address -fno-omit-frame-pointer -g")
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS_RELEASE "${FLAGS} -DNDEBUG -Ofast")
set(CMAKE_CXX_FLAGS_DEBUG "${FLAGS} -g3 -O0")

add_hpx_executable(
  ${PROJECT_NAME}
  DEPENDENCIES
    SILO::silo
  SOURCES
    src/main.cpp
  HEADERS
  	include/Allocator.hpp
  	include/AutoDiff.hpp
  	include/Box.hpp
  	include/CoupledState.hpp
 	include/Concepts.hpp
 	include/Constants.hpp
 	include/Debug.hpp
  	include/GasState.hpp
  	include/Grid.hpp
  	include/Integer.hpp
  	include/Limiters.hpp
  	include/Material.hpp
  	include/Math.hpp
  	include/Matrix.hpp
  	include/RadiationState.hpp
  	include/Real.hpp
  	include/Select.hpp
  	include/Silo.hpp
  	include/Simd.hpp
  	include/Vector.hpp
)
target_link_options(${PROJECT_NAME} PRIVATE -rdynamic)

target_include_directories(${PROJECT_NAME} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
