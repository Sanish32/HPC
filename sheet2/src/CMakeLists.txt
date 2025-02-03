#Run files
add_executable(quicksort-serial main.cpp quicksort.cpp)
add_executable(quicksort-distributed main.cpp quicksort_distributed.cpp)

# MPI
find_package(MPI REQUIRED)

target_link_libraries(quicksort-distributed MPI::MPI_C)


target_compile_options(quicksort-serial PRIVATE -DGPU=0 -DDISTRIBUTED=0)
target_compile_options(quicksort-distributed PRIVATE -DGPU=0 -DDISTRIBUTED=1)

