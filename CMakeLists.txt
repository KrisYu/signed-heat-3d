cmake_minimum_required(VERSION 3.12)

# ------------------------------------------------------------
# Project
# ------------------------------------------------------------
project(signed-heat-3d LANGUAGES CXX)

# ------------------------------------------------------------
# Build type
# ------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build: Debug or Release" FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# ------------------------------------------------------------
# C++ standard
# ------------------------------------------------------------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ------------------------------------------------------------
# IMPORTANT: Disable installation of third-party dependencies
#
# signed-heat-3d is a *tool*, not a system library.
# We only want to install the final executable, not its deps.
# ------------------------------------------------------------
set(CMAKE_SKIP_INSTALL_RULES ON CACHE BOOL "" FORCE)
set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY ON CACHE BOOL "" FORCE)
set(POLYSCOPE_INSTALL OFF CACHE BOOL "" FORCE)
set(GEOMETRY_CENTRAL_INSTALL OFF CACHE BOOL "" FORCE)
set(TETGEN_INSTALL OFF CACHE BOOL "" FORCE)
set(AMGCL_INSTALL OFF CACHE BOOL "" FORCE)
set(LIBIGL_INSTALL OFF CACHE BOOL "" FORCE)


# ------------------------------------------------------------
# Compiler flags (keep your original spirit, but simplified)
# ------------------------------------------------------------
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
  message(STATUS "Using clang/gcc compiler flags")
  add_compile_options(
    -Wall -Wextra
    -Wno-unused-parameter
    -Wno-unused-variable
    -Wno-unused-function
    -Wno-deprecated-declarations
    -Wno-narrowing 
  )

  if(CMAKE_BUILD_TYPE STREQUAL "Release")
    add_compile_options(-O3 -DNDEBUG)
  else()
    add_compile_options(-g)
  endif()

elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  message(STATUS "Using MSVC compiler flags")
  add_compile_options(/W4)
  add_definitions(-DNOMINMAX -D_USE_MATH_DEFINES)
endif()

# ------------------------------------------------------------
# Global compile definitions
# ------------------------------------------------------------
add_compile_definitions(NGC_SAFETY_CHECKS)

# ------------------------------------------------------------
# Third-party dependencies (build-only)
# ------------------------------------------------------------
add_subdirectory(deps/geometry-central)
add_subdirectory(deps/polyscope)
add_subdirectory(deps/tetgen)
add_subdirectory(deps/amgcl)

# ------------------------------------------------------------
# Sources
# ------------------------------------------------------------
set(SRCS
  src/signed_heat_tet_solver.cpp
  src/signed_heat_grid_solver.cpp
  src/signed_heat_3d.cpp
  src/main.cpp
)


# ------------------------------------------------------------
# Executable (IMPORTANT: do NOT call it "main")
# ------------------------------------------------------------
add_executable(signed_heat_3d
  ${SRCS}
)

# ------------------------------------------------------------
# Include directories
# ------------------------------------------------------------
target_include_directories(signed_heat_3d PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_CURRENT_SOURCE_DIR}/deps/tetgen
  ${CMAKE_CURRENT_SOURCE_DIR}/deps/polyscope/deps/args
)

# libigl is header-only here
target_include_directories(signed_heat_3d PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/deps/libigl/include
)

# ------------------------------------------------------------
# Link libraries
# ------------------------------------------------------------
target_link_libraries(signed_heat_3d
  geometry-central
  polyscope
  tetgen
  amgcl::amgcl
)

# ------------------------------------------------------------
# Install rule
# Install directly to <prefix>/signed_heat_3d
# ------------------------------------------------------------
install(TARGETS signed_heat_3d
  RUNTIME DESTINATION .
)
