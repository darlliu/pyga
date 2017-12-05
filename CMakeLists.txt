cmake_minimum_required(VERSION 2.6)
project(mafslice)

#set(CMAKE_CXX_COMPILER "clang++")
#set(CMAKE_CXX_FLAGS "-O -std=c++1y")
set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_FLAGS "-O4 -std=gnu++14")
#set(_GLIBCXX_USE_CXX11_ABI 0)
set(SOURCE_DIR "../src")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "../bin")

#configure_file("${SOURCE_DIR}/base.hpp")

set(BOOST_ROOT)
message("Configuring core C++ binaries")
FIND_PACKAGE(Boost 1.60.0)
FIND_PACKAGE(PythonLibs 2.7.6)
add_library(ga SHARED ${SOURCE_DIR}/ga.cpp ${SOURCE_DIR}/ga.hpp)
target_link_libraries(ga m pthread stdc++)

#add_executable(test ${SOURCE_DIR}/main.cpp)
#target_link_libraries(test ga)
INCLUDE_DIRECTORIES(
    ${SOURCE_DIR}
    ${BOOST_ROOT}
)

LINK_DIRECTORIES(
    ${BOOST_LIB_DIR}
    )



IF(Boost_FOUND AND PYTHONLIBS_FOUND)
    message("BOOST AND PYTHON2 FOUND, configuring python interface")
    INCLUDE_DIRECTORIES ("${Boost_INCLUDE_DIRS}")
    INCLUDE_DIRECTORIES ( "${PYTHON_INCLUDE_DIRS}" )
    SET(Boost_USE_STATIC_LIBS OFF)
    SET(Boost_USE_MULTITHREADED ON)
    SET(Boost_USE_STATIC_RUNTIME OFF)
    FIND_PACKAGE(Boost COMPONENTS python)
    add_library(pyga SHARED ${SOURCE_DIR}/pyga.cpp)
    target_link_libraries(pyga ga boost_python)
    set_target_properties(pyga PROPERTIES PREFIX "")
ELSE()
    message("Error finding boost...")
ENDIF()
