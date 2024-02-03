# CMakeLists.txt

cmake_minimum_required(VERSION 3.12)
project(libGA VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add header files to the library target
add_library(libGA INTERFACE)
target_include_directories(libGA INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Set compiler optimization flags
target_compile_options(libGA INTERFACE -O2)

# Add executable for testing
add_executable(testExecutable test/test.cpp)
target_link_libraries(testExecutable PRIVATE libGA)

# Set the output directory for the test executable
set_target_properties(testExecutable PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin
)

# You can add additional compiler options or flags here if needed

# Optional: Install the library
install(TARGETS libGA
    EXPORT libGATargets
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include
)

install(DIRECTORY include/ DESTINATION include)

# Optional: Export the library targets
export(TARGETS libGA FILE libGATargets.cmake)

# Optional: Create a package configuration file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    libGAConfigVersion.cmake
    VERSION 1.0.0
    COMPATIBILITY AnyNewerVersion
)

install(FILES libGAConfig.cmake ${CMAKE_CURRENT_BINARY_DIR}/libGAConfigVersion.cmake
    DESTINATION lib/cmake/libGA
)
