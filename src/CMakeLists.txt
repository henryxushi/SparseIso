cmake_minimum_required(VERSION 2.6)
project(SparseIso)

set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

set(SPARSEISO_DIR /home/cbil/Henry/RNA-seq/Assembly/SparseIso_package/SparseIso-github-V1.1)
set(BOOST_LIB_DIR /usr/local/lib)
set(BOOST_INCLUDE_DIR /usr/local/include)
set(BAMTOOLS_LIB_DIR /home/cbil/Henry/bamtools/lib)
set(BAMTOOLS_INCLUDE_DIR /home/cbil/Henry/bamtools/include)
set(EIGEN_INCLUDE_DIR /home/cbil/Downloads/eigen-eigen-bdd17ee3b1b3)

include_directories(. ${BOOST_INCLUDE_DIR} ${BAMTOOLS_INCLUDE_DIR} ${EIGEN_INCLUDE_DIR} ${SPARSEISO_DIR})
link_directories(${BOOST_LIB_DIR} ${BAMTOOLS_LIB_DIR})

add_executable(${PROJECT_NAME} SparseIso.cpp filter_bam.cpp getinsertsize.cpp Info.cpp options.cpp readinstance.cpp utility.cpp)

target_link_libraries(${PROJECT_NAME} pthread boost_math_c99 bamtools boost_program_options boost_regex boost_graph boost_thread boost_system boost_filesystem)
