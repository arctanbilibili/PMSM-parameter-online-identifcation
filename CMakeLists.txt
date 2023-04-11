cmake_minimum_required(VERSION 3.0.0)
project(bbb VERSION 0.1.0)

include(CTest)
enable_testing()

# Find the QtWidgets library

set(CMAKE_C_STANDARD 90)
set(CMAKE_C_STANDARD_REQUIRED ON)

add_executable(bbb main.c matrix.c RLS.c MARS.c EKF.c pmsm.c)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

#message( STATUS "1234" ${OpenCV_INCLUDE_DIRS})