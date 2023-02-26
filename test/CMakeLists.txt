cmake_minimum_required (VERSION 2.4)

if (COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

add_executable (testtype testtype.cpp)
target_link_libraries (testtype ellipse_detection )

add_executable (testlsd testlsd.cpp)
target_link_libraries (testlsd ellipse_detection )

add_executable (testgroup testgroup.cpp)
target_link_libraries (testgroup ellipse_detection )

add_executable (testinital testinital.cpp)
target_link_libraries (testinital ellipse_detection )

add_executable (testdetect testdetect.cpp)
target_link_libraries (testdetect ellipse_detection )