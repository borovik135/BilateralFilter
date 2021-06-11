cmake_minimum_required(VERSION 2.8)

project(bilateralfilter)


find_package(OpenMP)
if(OPENMP_FOUND)
    message("OPENMP FOUND")
endif()
if(UNIX)
    set(CMAKE_CXX_FLAGS
        "-Wall -Wextra -Werror -std=c++98 -Wfatal-errors ${OpenMP_CXX_FLAGS}")
else()
    set(CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS}")
endif()


message(${CMAKE_CXX_COMPILER})
message(${CMAKE_CXX_FLAGS})


add_executable(bilateralfilter src/main.cpp
      src/Point.cpp
      src/Sample.cpp
      src/FileIO.cpp
      )

