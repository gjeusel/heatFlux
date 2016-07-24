cmake_minimum_required(VERSION 3.3)
project(heatFlux)
enable_language(Fortran)

# Sanitize the build directory
set(heatFlux_BUILD_NAME     "" CACHE STRING "heatFlux configuration name")
set(heatFlux_BUILD_NAME_DIR .  CACHE STRING "heatFlux configuration directory")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib-heatFlux CACHE PATH "Single directory for all static libraries.")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib-heatFlux CACHE PATH "Single directory for all libraries.")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin CACHE PATH "Single directory for all executables.")

set(CMAKE_Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/modules" CACHE PATH "Single directory for all fortran modules.")

#-----------------------------------------------------------------------------
# - Compiling types availables

# CMAKE_BUILD_TYPE_STRINGS_EXTRA can be altered by the toolchain to add compile types
set(CMAKE_BUILD_TYPE_STRINGS "Debug" "Release" ${CMAKE_BUILD_TYPE_STRINGS_EXTRA})
set(CMAKE_BUILD_TYPE_PREFIX_STRINGS "")
foreach(type ${CMAKE_BUILD_TYPE_STRINGS})
  list (APPEND CMAKE_BUILD_TYPE_PREFIX_STRINGS "_${type}" )
endforeach()
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "" ${CMAKE_BUILD_TYPE_STRINGS})

foreach(type "" ${CMAKE_BUILD_TYPE_PREFIX_STRINGS})
  string(TOUPPER "${type}" type_UP)
  set(CMAKE_Fortran_FLAGS${type_UP} "${CMAKE_Fortran_FLAGS${type_UP}} ${CMAKE_Fortran_FLAGS${type_UP}_${CMAKE_Fortran_COMPILER_ID}}")
  message(STATUS "CMAKE_Fortran_FLAGS${type_UP}=${CMAKE_Fortran_FLAGS${type_UP}}")
endforeach()

#-----------------------------------------------------------------------------
# Compilers details

if(CMAKE_Fortran_COMPILER MATCHES "gfortran" OR CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
endif()

#-----------------------------------------------------------------------------
# Pré-compile options handling :

# macro to simplify options_list appending
macro(add_option opt default desc)
  list(APPEND options_list ${opt})
  set(${opt}_DEFAULT ${default} CACHE INTERNAL "")
  set(${opt}_DESC    ${desc} CACHE INTERNAL "")
  set(${opt}_VALUES  ${ARGN} CACHE INTERNAL "")
endmacro()

# pré-compile options :
add_option(heatFlux_CALC 3D "Dimensionality" 2D;3D)

# name's formalism about pre-compile options
include(FeatureSummary)
file(WRITE "${heatFlux_BINARY_DIR}/heatFlux_options.cmake" "")
foreach(opt ${options_list})
  if("${${opt}_VALUES}" MATCHES "ON;OFF")
    option(${opt} ${${opt}_DESC} ${${opt}_DEFAULT})
    add_feature_info(${opt} ${${opt}} "")
    file(APPEND "${heatFlux_BINARY_DIR}/heatFlux_options.cmake"
    "set(${opt} ${${opt}} CACHE BOOL \"${${opt}_DESC}\" FORCE)
")
  else() # combo box
    set(${opt} ${${opt}_DEFAULT} CACHE STRING "${${opt}_DESC} (${${opt}_VALUES})")
    set_property(CACHE ${opt} PROPERTY STRINGS ${${opt}_VALUES})
    add_feature_info(${opt}_${${opt}} ON "")
    file(APPEND "${heatFlux_BINARY_DIR}/heatFlux_options.cmake"
        "set(${opt} ${${opt}} CACHE STRING \"${${opt}_DESC} (${${opt}_VALUES})\" FORCE)
")

    # check wether ${opt} as a valid value
    set(${opt}_valid OFF)
    foreach(val ${${opt}_VALUES})
      if("${${opt}}" MATCHES ${val})
          set(${opt}_valid ON)
      endif()
    endforeach()

    if(NOT ${${opt}_valid})
      message(FATAL_ERROR "  ${opt}=${${opt}} should be one of ${${opt}_VALUES}  ")
    endif()
  endif()
endforeach()

FEATURE_SUMMARY(WHAT ENABLED_FEATURES  DESCRIPTION "Enable Features:")
FEATURE_SUMMARY(WHAT DISABLED_FEATURES DESCRIPTION "Disabled features:")

get_property(ENABLED_KEYWORDS GLOBAL PROPERTY ENABLED_FEATURES)
set(BOOLEAN_KEYWORDS_ON "")
foreach(keyword ${ENABLED_KEYWORDS})
  if("${${keyword}_VALUES}" MATCHES "ON;OFF")
    list(APPEND BOOLEAN_KEYWORDS_${${keyword}} ${keyword})
    list(REMOVE_ITEM ENABLED_KEYWORDS ${keyword})
  endif()
  list(APPEND definitions -D${keyword})
endforeach()

get_property(DISABLED_KEYWORDS GLOBAL PROPERTY DISABLED_FEATURES)
set(BOOLEAN_KEYWORDS_OFF "")
foreach(keyword ${DISABLED_KEYWORDS})
  if("${${keyword}_VALUES}" MATCHES "ON;OFF")
    list(APPEND BOOLEAN_KEYWORDS_${${keyword}} ${keyword})
    list(REMOVE_ITEM ENABLED_KEYWORDS ${keyword})
  endif()
endforeach()

# Adding -D[list of options ON]
add_definitions(${definitions})



#---------------------------------------------------------------------
# Source files handling

configure_file(${heatFlux_SOURCE_DIR}/config.h.in ${heatFlux_BINARY_DIR}/config.h  @ONLY)

# set up include-directories
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})
include_directories(${heatFlux_BINARY_DIR})
include_directories("${heatFlux_SOURCE_DIR}/src")


add_executable(heatFlux
  ./src/heatFlux.F90
  ./src/Parameters.F90
  ./src/Common_Variables.F90
  ./src/IO.F90
  ./src/Initialization.F90
  ./src/TemporalVariationRate.F90
  ./src/update_Temp.F90
  )
