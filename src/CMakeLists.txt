# CMake
# Fernando de Goes (fdegoes@caltech.edu)
# Copyright @ 2012
project(smoother)

cmake_minimum_required(VERSION 2.6)
cmake_policy(VERSION 2.6)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake )

if(WIN32)
    add_definitions( "/W3 /D_CRT_SECURE_NO_WARNINGS /wd4503 /nologo" )
endif(WIN32)



###########
# OPTIONS #
###########

option(USE_DEBUG "Turn on the debug mode" OFF)
option(USE_64BIT_ARCH "Compile the 64bit executable" ON)
option(USE_EIGEN "Use Eigen-library" ON)
option(USE_IPOPT "Use Ipopt" ON)

if(USE_DEBUG)
    add_definitions(-DDEBUG)
    set(CMAKE_BUILD_TYPE Debug)
else(USE_DEBUG)
    set(CMAKE_BUILD_TYPE Release)
endif(USE_DEBUG)

if(USE_64BIT_ARCH)
    add_definitions(-m64)
    set(LINK_FLAGS -m64)
endif(USE_64BIT_ARCH)

############
# Packages #
############

set(QT_USE_QTXML TRUE)
set(QT_USE_QTOPENGL TRUE)

find_package(OpenGL)
find_package(Qt4)
find_package(QGLViewer)
find_package(CGAL COMPONENTS Qt4)

#########
# Files #
#########

set(
HDRS
third/array.h
third/energy.h
third/lbfgs.h
third/line_search.h
third/newton.h
third/sparse_array.h
third/sparse_matrix.h
third/suite_sparse_qr.h
console_color.h
render.h
scene.h
timer.h
viewer.h
window.h
)

set( 
SRCS
third/sparse_array.cpp
third/sparse_matrix.cpp
main.cpp
render.cpp
scene.cpp
window.cpp
)

set( 
MOCS 
moc_viewer.cxx 
moc_window.cxx
)

set(
UIS
window.ui
)

#########
# Build #
#########

# Includes
include_directories(BEFORE . ./build)
include(${QT_USE_FILE})
include_directories(${GLUT_INCLUDE_DIR})
include_directories(${OPENGL_INCLUDE_DIR})
include_directories(${QGLViewer_INCLUDE_DIR})
include(${CGAL_USE_FILE})

# UI files
qt4_wrap_ui( DT_UI_FILES ${UIS} )

# MOC files
qt4_generate_moc( viewer.h moc_viewer.cxx )
qt4_generate_moc( window.h moc_window.cxx )

# Executable
add_executable( ${PROJECT_NAME} ${SRCS} ${MOCS} ${DT_UI_FILES} ${DT_RESOURCE_FILES} )

# Link with Qt
target_link_libraries( ${PROJECT_NAME} ${QT_LIBRARIES} )

# Link with Glut and OpenGL
target_link_libraries( ${PROJECT_NAME} ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES}  )

# Link with QGLViewer
target_link_libraries( ${PROJECT_NAME} ${QGLViewer_LIBRARY} )

# Link with CGAL
target_link_libraries( ${PROJECT_NAME} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
