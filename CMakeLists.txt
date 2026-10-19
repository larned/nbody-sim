cmake_minimum_required(VERSION 3.16)
project(nbody-sim VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

if(MSVC)
  add_compile_options(/W4)
else()
  add_compile_options(-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wno-sign-conversion)
endif()

add_library(nbody
  src/system.cpp
  src/forces.cpp
  src/integrators.cpp
  src/scenarios.cpp
)
target_include_directories(nbody PUBLIC include)

add_executable(nbody_cli src/main.cpp)
target_link_libraries(nbody_cli PRIVATE nbody)
set_target_properties(nbody_cli PROPERTIES OUTPUT_NAME nbody)

add_executable(nbody_study tools/study.cpp)
target_link_libraries(nbody_study PRIVATE nbody)

include(CTest)
if(BUILD_TESTING)
  add_executable(nbody_tests tests/test_main.cpp)
  target_link_libraries(nbody_tests PRIVATE nbody)
  add_test(NAME nbody_tests COMMAND nbody_tests)
endif()
