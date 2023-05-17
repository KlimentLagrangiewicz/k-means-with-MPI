cmake_minimum_required(VERSION 2.25)

project(k-means-with-MPI LANGUAGES C)

set(CMAKE_C_STANDARD 90)
find_package(MPI REQUIRED)
if (MPI_FOUND)
    include_directories(SYSTEM ${MPI_INCLUDE_PATH})
else (MPI_FOUND)
    message(SEND_ERROR "This application cannot compile without MPI")
endif (MPI_FOUND)

add_executable(k-means-with-MPI src/main.c src/help.h src/help.c src/kmeans.h src/kmeans.c)

target_link_libraries(k-means-with-MPI ${MPI_C_LIBRARIES} m)
