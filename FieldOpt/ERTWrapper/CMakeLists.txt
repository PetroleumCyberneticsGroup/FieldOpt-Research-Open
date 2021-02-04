cmake_minimum_required(VERSION 2.8)
project(ertwrapper)

include(Sources.cmake)

# Target lib ===============================================
add_library(ertwrapper STATIC
  ${ERTWRAPPER_HEADERS}
  ${ERTWRAPPER_SOURCES})

add_library(fieldopt::ertwrapper ALIAS ${PROJECT_NAME})

include_directories(
  ${EIGEN3_INCLUDE_DIR}
  ${ERT_INCLUDE_DIR})

target_link_libraries(ertwrapper
  PUBLIC ri:ert_ecl)

target_include_directories(ertwrapper PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

add_compile_options(-std=c++11)

# Enable position independent code when compiling lib, so that
# shared library libwellindexcalculator.so may be linked with
# static reservoir and ertwrapper libraries
if (BUILD_WIC_ADGPRS)
    set_target_properties(ertwrapper
      PROPERTIES COMPILE_FLAGS "-fPIC")
endif()

# Test exec ================================================
# Unit tests
if (BUILD_TESTING)
    add_executable(test_ertwrapper ${ERTWRAPPER_TESTS})

    target_link_libraries(test_ertwrapper
      PUBLIC fieldopt::ertwrapper
      PUBLIC fieldopt::runner
      Qt5::Core
      ${Boost_LIBRARIES}
      ${gtest}
      ${gtest_main}
      ${CMAKE_THREAD_LIBS_INIT})

    set_target_properties(test_ertwrapper
      PROPERTIES COMPILE_FLAGS "-fPIC")

    target_include_directories(test_ertwrapper PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/tests>)

    add_test(NAME test_ertwrapper
      COMMAND $<TARGET_FILE:test_ertwrapper>)

endif()

# Install ==================================================
install(TARGETS ertwrapper
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib/static)
