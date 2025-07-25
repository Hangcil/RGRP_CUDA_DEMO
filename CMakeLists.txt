cmake_minimum_required(VERSION 3.16)

project(RGRP_QML VERSION 0.1 LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(MSVC_CXX_FLAGS /Zc:__cplusplus /std:c++17 /EHsc)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CUDA_ARCHITECTURES 89)
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -O3 -use_fast_math --ptxas-options=-v -rdc=true")

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)


set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -rdc=true")

set(CMAKE_PREFIX_PATH "C:/libs\\Qt\\6.9.1\\msvc2022_64"
    )
find_package(Qt6 REQUIRED COMPONENTS Quick)
qt_standard_project_setup(REQUIRES 6.5)

set_source_files_properties(
  ${MATH_KERNELS}
  PROPERTIES LANGUAGE CUDA
)

file(GLOB_RECURSE MATH_KERNELS
    "${CMAKE_CURRENT_SOURCE_DIR}/mathSolvers/*.cu"
)
file(GLOB_RECURSE CXX_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/mathSolvers/*.cpp"
)
if(WIN32)
    set(APP_ICON_RESOURCE "${CMAKE_CURRENT_SOURCE_DIR}/images/icon.rc")
endif()
qt_add_executable(RGRP WIN32 main.cpp res.qrc ${APP_ICON_RESOURCE} ${CXX_SOURCES} ${MATH_KERNELS})
set_target_properties(RGRP PROPERTIES
  CUDA_PROPAGATE_HOST_FLAGS OFF
)
target_compile_options(RGRP PRIVATE
  # C++ files get the MSVC flag directly
  $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CXX_COMPILER_ID:MSVC>>:/Zc:__cplusplus>

  # CUDA files get it wrapped so nvcc passes it through
  $<$<AND:$<COMPILE_LANGUAGE:CUDA>,$<CXX_COMPILER_ID:MSVC>>:-Xcompiler="/Zc:__cplusplus">
)
qt_add_qml_module(RGRP
    URI RGRP_QML
    VERSION 1.0
    QML_FILES
        UI/Main.qml
)

set_target_properties(RGRP PROPERTIES
  CUDA_SEPARABLE_COMPILATION ON
  QT_QML_GENERATE_QMLLS_INI ON
)

target_include_directories(RGRP PRIVATE
    ${QT_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}/mathSolvers
)

target_link_libraries(RGRP PRIVATE
    Qt6::Quick
)