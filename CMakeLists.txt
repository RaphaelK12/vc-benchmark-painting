cmake_minimum_required(VERSION 2.8.12)

project(vc-benchmark-painting)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)

find_package(Qt5Core)
find_package(Vc 1.3.3)

include_directories(${Vc_INCLUDE_DIR})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Vc_DEFINITIONS}")

add_executable(${PROJECT_NAME} "main.cpp")

set_property(TARGET ${PROJECT_NAME} APPEND PROPERTY COMPILE_OPTIONS "${Vc_ARCHITECTURE_FLAGS}")

target_link_libraries(${PROJECT_NAME} Qt5::Core ${Vc_LIBRARIES})
