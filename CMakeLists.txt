cmake_minimum_required(VERSION 3.12)

project(mpm)

#add_definitions(-DDEBUG)

link_directories(SYSTEM)

set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

set(CMAKE_CXX_STANDARD 17)

include_directories(SYSTEM ${CMAKE_CURRENT_SOURCE_DIR} include)


add_executable(mpm_rigid_apic "mpm_rigid_apic.cpp")

add_executable(mpm_rigid_pic "mpm_rigid_pic.cpp")

add_executable(mpm_elastic_pic "mpm_elastic_pic.cpp")

add_executable(mpm_elastic_apic "mpm_elastic_apic.cpp")

add_executable(mpm_snow "mpm_snow.cpp")

add_executable(mpm_snow_elastic "mpm_snow_elastic.cpp")