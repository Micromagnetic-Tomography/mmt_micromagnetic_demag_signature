cmake_minimum_required(VERSION 3.15...3.26)
project(${SKBUILD_PROJECT_NAME} LANGUAGES C)

find_package(
  Python
  COMPONENTS Interpreter Development.Module NumPy
  REQUIRED)

# Parallelization
find_package(OpenMP REQUIRED COMPONENTS C)
message(STATUS "OpenMP found: C ${OpenMP_C_FLAGS}")

set(CLIB_DIR "${CMAKE_CURRENT_SOURCE_DIR}/mmt_micromagnetic_demag_signature/clib")

message(STATUS "Python: ${Python_VERSION} (${Python_EXECUTABLE})")
message(STATUS "NumPy include: ${Python_NumPy_INCLUDE_DIRS}")

add_custom_command(
  OUTPUT mds_clib.c
  COMMENT
  "Making mds_clib.c from ${CLIB_DIR}/mds_clib.pyx"
  COMMAND Python::Interpreter -m cython "${CLIB_DIR}/mds_clib.pyx" --output-file mds_clib.c
  DEPENDS "${CLIB_DIR}/mds_clib.pyx"
  VERBATIM)

# Note that we add clib.c to compile the C module, in addition to the pyx module
# The C module contains the functions being called in the pyx file
python_add_library(mds_clib MODULE mds_clib.c ${CLIB_DIR}/clib.c WITH_SOABI)
target_include_directories(mds_clib PRIVATE ${Python_NumPy_INCLUDE_DIRS})  # numpy header files dir
target_include_directories(mds_clib PRIVATE ${CLIB_DIR})  # Header files in the clib dir
# Link libraries
target_link_libraries(mds_clib PRIVATE OpenMP::OpenMP_C)
# Compiler options
target_compile_options(mds_clib PRIVATE -O3 -Wall -Wno-cpp -Wno-unused-function)

# Use relative path for the built modules
install(TARGETS mds_clib DESTINATION mmt_micromagnetic_demag_signature/clib)
