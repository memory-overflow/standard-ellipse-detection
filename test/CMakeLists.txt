cmake_minimum_required (VERSION 2.4)

if (COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

set (OpenCV_FOUND 1)
find_package (OpenCV REQUIRED)

project (tests)

add_executable (testtype testtype.cpp)
target_link_libraries (testtype ellipse_detection ${OpenCV_LIBS})

add_executable (testlsd testlsd.cpp)
target_link_libraries (testlsd ellipse_detection ${OpenCV_LIBS})

add_executable (testgroup testgroup.cpp)
target_link_libraries (testgroup ellipse_detection ${OpenCV_LIBS})

add_executable (testinital testinital.cpp)
target_link_libraries (testinital ellipse_detection ${OpenCV_LIBS})

add_executable (testdetect testdetect.cpp)
target_link_libraries (testdetect ellipse_detection ${OpenCV_LIBS})