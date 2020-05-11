cmake_minimum_required(VERSION 2.8)
project(optimization)

include(Sources.cmake)

# Target lib ===============================================
add_library(optimization STATIC
    ${OPTIMIZATION_HEADERS}
    ${OPTIMIZATION_SOURCES})

add_library(fieldopt::optimization ALIAS ${PROJECT_NAME})

target_link_libraries(optimization
    PUBLIC fieldopt::model
    PUBLIC fieldopt::simulation
    PUBLIC fieldopt::constraintmath
    PUBLIC fieldopt::reservoir
    PUBLIC fieldopt::runner
    PUBLIC ${gp}
    PUBLIC ${snopt}
    Qt5::Core
    PUBLIC ${Boost_LIBRARIES})

target_include_directories(optimization PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/constraints>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/objective>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/optimizers>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/solvers>)

add_compile_options(-std=c++11)

# Test exec ================================================
if (BUILD_TESTING)
  # Unit tests
  add_executable(test_optimization ${OPTIMIZATION_TESTS})

  target_link_libraries(test_optimization
      PUBLIC fieldopt::optimization
      PUBLIC fieldopt::runner
      ${gtest}
      ${gtest_main}
      ${CMAKE_THREAD_LIBS_INIT})

  target_include_directories(test_optimization PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/tests>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/constraints>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/objective>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/optimizers>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/solvers>)

  add_test(NAME test_optimization
      COMMAND $<TARGET_FILE:test_optimization>)
endif()

# Install ==================================================
install(TARGETS optimization
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib/static)
