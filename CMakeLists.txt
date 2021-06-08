cmake_minimum_required(VERSION 3.19)

project(POBR_Projekt)
find_package(OpenCV REQUIRED)
set(APP_TARGET pobr_starbucks)

set(CMAKE_CXX_STANDARD 20)

set(LIBRARIES ${OpenCV_LIBS})

add_executable(${APP_TARGET} src/main.cpp src/color.h src/feature.h src/coefficients.h src/processing.h src/structures.h)
target_link_libraries(${APP_TARGET} ${LIBRARIES})
