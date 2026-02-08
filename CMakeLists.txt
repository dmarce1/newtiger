cmake_minimum_required(VERSION 3.10)

project(newtiger LANGUAGES CXX)

set(FLAGS "-Wall -march=native -fsanitize=address -fno-omit-frame-pointer -g")
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS_RELEASE "${FLAGS} -DNDEBUG -Ofast")
set(CMAKE_CXX_FLAGS_DEBUG "${FLAGS} -g3 -O0")

add_executable(${PROJECT_NAME})

target_sources(
		newtiger 
	PRIVATE
    	src/main.cpp
	PUBLIC
    	${CMAKE_CURRENT_BINARY_DIR}/config.hpp
    	include/Constants.hpp
    	include/IO.hpp
    	include/FourierLegendre.hpp
    	include/Legendre.hpp
    	include/Math.hpp
    	include/Matrix.hpp
    	include/Quadrature.hpp
   		include/Units.hpp
    	include/Vector.hpp
)

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/config.hpp.in
    ${CMAKE_CURRENT_BINARY_DIR}/config.hpp
)

target_include_directories(${PROJECT_NAME} PUBLIC "${CMAKE_CURRENT_BINARY_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/include")
