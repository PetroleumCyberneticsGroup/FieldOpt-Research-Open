cmake_minimum_required(VERSION 3.2)
project(FieldOpt VERSION 0.3 LANGUAGES CXX)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Qt libraries
set(CMAKE_AUTOMOC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
find_package(Qt5Core REQUIRED)
find_package(Qt5Gui REQUIRED)
include_directories(${Qt5Core_INCLUDE_DIRS} ${Qt5Gui_INCLUDE_DIRS})

# Boost libraries
find_package(Boost)
if (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIR})
endif ()
set(Boost_USE_STATIC_LIBS OFF) # enable dynamic linking
set(Boost_USE_MULTITHREAD ON) # enable multithreading
find_package(Boost COMPONENTS REQUIRED program_options)

find_package(Armadillo)
find_package(Threads)

# Set this directory as the root dir
include_directories(${CMAKE_SOURCE_DIR})

include(CTest)
enable_testing()

add_subdirectory(ERTWrapper)
add_subdirectory(Model)
add_subdirectory(Utilities)
add_subdirectory(Optimization)
add_subdirectory(Simulation)
add_subdirectory(Runner)
add_subdirectory(WellIndexCalculator)
add_subdirectory(AdgprsResultsReader)

## Copy ADGPRS summary convertion python scripts
file(GLOB CONVERSION_SCRIPTS
        "AdgprsSummaryConverter/*.py")
file(COPY ${CONVERSION_SCRIPTS} DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/AdgprsSummaryConverter)

# Copy execution scripts
file(GLOB EXECUTION_SCRIPTS
        "execution_scripts/*.sh")
file(COPY ${EXECUTION_SCRIPTS} DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/execution_scripts)

## Copy the examples into the build dir
file(COPY ${CMAKE_SOURCE_DIR}/../examples DESTINATION ${CMAKE_BINARY_DIR}/)
#
## Copy the test driver file into the build dir
file(COPY ${CMAKE_SOURCE_DIR}/Utilities/tests/driver/driver.json DESTINATION ${CMAKE_BINARY_DIR}/examples)
#
## Create an empty directory for test output
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/fieldopt_output)