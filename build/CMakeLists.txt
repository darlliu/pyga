cmake_minimum_required(VERSION 2.6)
project(mafslice)

#set(CMAKE_CXX_COMPILER "clang++")
#set(CMAKE_CXX_FLAGS "-O -std=c++1y")
set(CMAKE_CXX_COMPILER "g++")
#set(CMAKE_CXX_FLAGS "-O3 -fno-omit-frame-pointer -std=c++1y -lprofiler")
set(CMAKE_CXX_FLAGS "-O4 -std=c++1y")
#set(CMAKE_CXX_FLAGS "-O3 -std=c++1y  -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free ")
set(_GLIBCXX_USE_CXX11_ABI 0)
set(SOURCE_DIR "../src")
#set(BOOST_ROOT "/auto/igb-libs/linux/centos/6.x/x86_64/pkgs/boost/alt-1.59.0/boost_1_59_0")
set(BOOST_ROOT "/baldig/biotools/libraries/anaconda2")
#set(BOOST_LIB_DIR "${BOOST_ROOT}/stage/lib")
set(PYTHON_ROOT "/baldig/biotools/libraries/anaconda2")
set(PYTHON_LIBRARY "${PYTHON_ROOT}/lib")
set(PYTHON_INCLUDE_DIR "${PYTHON_ROOT}/include/python2.7")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "../bin")

#configure_file("${SOURCE_DIR}/base.hpp")

INCLUDE_DIRECTORIES(
    ${SOURCE_DIR}
    ${BOOST_ROOT}
)

LINK_DIRECTORIES(
    ${BOOST_LIB_DIR}
    )


message("Configuring core C++ binaries")
FIND_PACKAGE(Boost 1.60.0)
FIND_PACKAGE(PythonLibs 2.7.11)
add_library(ga SHARED ${SOURCE_DIR}/ga.cpp ${SOURCE_DIR}/ga.hpp)
target_link_libraries(ga m pthread stdc++)
#add_executable(test ${SOURCE_DIR}/main.cpp)
#target_link_libraries(test ga)

IF(Boost_FOUND AND PYTHONLIBS_FOUND)
    #message("BOOST 1.59 AND PYTHON 2.7.10 FOUND, configuring python interface")
    INCLUDE_DIRECTORIES ("${Boost_INCLUDE_DIRS}")
    INCLUDE_DIRECTORIES ( "${PYTHON_INCLUDE_DIRS}" )
    SET(Boost_USE_STATIC_LIBS OFF)
    SET(Boost_USE_MULTITHREADED ON)
    SET(Boost_USE_STATIC_RUNTIME OFF)
    FIND_PACKAGE(Boost 1.60.0 COMPONENTS python)
    add_library(pyga SHARED ${SOURCE_DIR}/pyga.cpp)
    target_link_libraries(pyga ga boost_python)
    set_target_properties(pyga PROPERTIES PREFIX "")
ELSE()
    message("Error finding boost...")
ENDIF()
