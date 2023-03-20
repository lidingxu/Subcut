# mkdir build
# cd build
# cmake ..
# make
include(ProcessorCount)
ProcessorCount(N)
if(NOT N EQUAL 0)
  set(CMAKE_MAKE_PROGRAM"${CMAKE_MAKE_PROGRAM} -j ${N}")      
endif()
cmake_minimum_required(VERSION 3.3)
set(CMAKE_BUILD_TYPE RELEASE)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g ")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}  -O3 -march=native -flto")
project(multilinear)


include_directories(${SCIP_DIR}/src  ${SCIP_DIR}/build) 
link_directories(${SCIP_DIR}/build/lib)

add_executable(multi
  src/main.cpp
  src/sepa_intersub.cpp
  src/sepa_interlattice.cpp
  src/utility.cpp
)

# link to math library if it is available
find_library(LIBM m)
if(NOT LIBM)
  set(LIBM "")
endif()


target_link_libraries(multi -ldl -lscip ${LIBM})