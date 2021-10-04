project(tri-mesh-hole-closer)

cmake_minimum_required(VERSION 3.13)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

add_library(${PROJECT_NAME} INTERFACE)

# libigl
option(LIBIGL_WITH_TRIANGLE "Use Triangle" ON)
target_include_directories(${PROJECT_NAME} INTERFACE external/libigl/include)
find_package(LIBIGL REQUIRED QUIET)

# eigen (header only library)
target_include_directories(${PROJECT_NAME} INTERFACE external/eigen)

target_link_libraries(${PROJECT_NAME} INTERFACE  igl::triangle)

# set project include dir
target_include_directories(${PROJECT_NAME} INTERFACE
                           "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
                           "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")