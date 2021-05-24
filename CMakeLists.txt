## CMake configuration file of Ellipsedetcet project

cmake_minimum_required (VERSION 2.4)

if (POLICY CMP0048)
    cmake_policy (SET CMP0048 NEW)
endif ()


if (COMMAND cmake_policy)
  cmake_policy (SET CMP0003 NEW)
endif (COMMAND cmake_policy)


set (PACKAGE_NAME    "ellipse_detection")
set (PACKAGE_VERSION "1.0")


project (${PACKAGE_NAME} VERSION ${PACKAGE_VERSION} LANGUAGES CXX)

set (OpenCV_FOUND 1)
find_package (OpenCV REQUIRED)

# Manual set open library
#include_directories(/xxx/opencv_tag_v3.2.0)
#include_directories (/xxx/opencv_tag_v3.2.0/include)
#link_directories(/xxx/opencv_tag_v3.2.0/lib)
set (OTHER_LIBS -llapack -lblas -lgfortran)

set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRE ON)

set (EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin")
set (LIBRARY_OUTPUT_PATH "${PROJECT_BINARY_DIR}/lib")

include_directories ("${PROJECT_SOURCE_DIR}/include")

add_definitions ("-Wall")

set (DEPEND_FILE
  "src/unitily.cpp"
  "src/detect.cpp"
  "src/compute.cpp"
  "src/lsd.cpp"
  "src/cvcannyapi.cpp"
)

set (Opencv_LIBS
    -lopencv_core
    -lopencv_imgproc
    -lopencv_highgui
)

add_library (${PACKAGE_NAME} ${DEPEND_FILE})
target_link_libraries(${PACKAGE_NAME} ${Opencv_LIBS} ${OTHER_LIBS})

set (PUBLIC_HDRS
  include/types.hpp
  include/unitily.h
  include/defines.h
  include/detect.h
  include/compute.h
)

set (INCLUDE_INSTALL_DIR "include/${PACKAGE_NAME}")
set (RUNTIME_INSTALL_DIR "bin")
set (LIBRARY_INSTALL_DIR "lib")

install (FILES ${PUBLIC_HDRS} DESTINATION ${INCLUDE_INSTALL_DIR}/${GFLAGS_INCLUDE_DIR})

install (TARGETS ${PACKAGE_NAME}
  RUNTIME DESTINATION ${RUNTIME_INSTALL_DIR}
  LIBRARY DESTINATION ${LIBRARY_INSTALL_DIR}
  ARCHIVE DESTINATION ${LIBRARY_INSTALL_DIR}
)

if (BUILD_TESTING)
  add_subdirectory (test)
endif ()
