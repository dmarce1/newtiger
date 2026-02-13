cmake_minimum_required(VERSION 3.10)

project(newtiger LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(HPX REQUIRED)
find_package(SILO REQUIRED)

set(FLAGS "-Wall -march=native")
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
    include/config.hpp
    include/Basis.hpp
    include/Constants.hpp
    include/Gas.hpp
    include/Interval.hpp
    include/IO.hpp
    include/Legendre.hpp
    include/Math.hpp
    include/Matrix.hpp
    include/Quadrature.hpp
    include/SubArray.hpp
    include/StateFactory.hpp
   	include/Units.hpp
    include/Vector.hpp
)
#target_link_libraries(yourTarget PRIVATE SILO::silo)

target_include_directories(${PROJECT_NAME} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
