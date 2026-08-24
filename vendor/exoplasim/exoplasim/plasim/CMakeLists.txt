# The PlaSim model build.
#
# Worldbuilding frame: this builds the Vesper project's climate model, a hard
# fork of ExoPlaSim. Nothing here concerns the simulated planet.
#
# Drive it through exoplasim/scripts/build_model.py, which is where the flags
# are read from config/planet.yaml and where the argument checking lives. This
# file is deliberately unfriendly to invoke by hand: every setting below is
# REQUIRED and has no default, because the build system this replaced defaulted
# silently and that is what it cost -- notes/audits/model-build-driver.md.
#
# WHY CMAKE AND NINJA. The rule set this replaces declared Fortran module
# dependencies BY HAND, one edge per line, and a missing edge does not fail: it
# works until the object order stops hiding it. CONS-13 is one such edge that
# survived because ${MPIMOD}.o happens to be first in the link list. CMake with
# Ninja scans the `use` statements and generates the edges, so that whole class
# of defect stops being possible rather than being fixed once.

cmake_minimum_required(VERSION 3.20)

# Nothing here is portable to another machine's compiler by accident: the flag
# line arrives whole from config/planet.yaml and CMake must not add to it.
set(CMAKE_BUILD_TYPE "" CACHE STRING "" FORCE)

project(plasim LANGUAGES Fortran C)

# ---------------------------------------------------------------------------
# Required settings. Each one is checked; none has a default.
# ---------------------------------------------------------------------------

set(PLASIM_NLAT      "" CACHE STRING "latitudes (32 48 64 96 128 160 192 256)")
set(PLASIM_NLEV      "" CACHE STRING "vertical levels")
set(PLASIM_NPRO      "" CACHE STRING "ranks (mpi) or threads (omp); 1 for serial")
set(PLASIM_PRECISION "" CACHE STRING "bytes per real: 4 or 8")
set(PLASIM_PARMODE   "" CACHE STRING "serial | mpi | omp")
set(PLASIM_FFLAGS    "" CACHE STRING "the whole Fortran flag line, from config/planet.yaml")
set(PLASIM_FFT       "fftmod"  CACHE STRING "fftmod | fft991mod")
set(PLASIM_PLANET    "p_earth" CACHE STRING "p_earth | p_mars")
set(PLASIM_SHTNS_PREFIX "" CACHE PATH "prefix holding include/shtns.f03 and lib/libshtns*.a")

function(require_one_of name value)
    if(NOT value IN_LIST ARGN)
        message(FATAL_ERROR
            "${name} is '${value}', which is not one of: ${ARGN}. "
            "It has no default: an unrecognised value is an error here, never a "
            "quiet fallback. notes/audits/model-build-driver.md.")
    endif()
endfunction()

foreach(v PLASIM_NLAT PLASIM_NLEV PLASIM_NPRO PLASIM_PRECISION PLASIM_PARMODE PLASIM_FFLAGS)
    if("${${v}}" STREQUAL "")
        message(FATAL_ERROR "${v} is not set, and has no default. Build through "
                            "exoplasim/scripts/build_model.py.")
    endif()
endforeach()

# THE LADDER, and it is `lib/rungs.py` that declares it. 48 is T31, which was
# missing here: `plasim/dat/T31` carries a complete surface set, `plasim.f90`
# has an mpstep branch for it and `radmod.f90` a tuning branch, and
# build_model.py resolved --res T31 cleanly and then died at this line.
# world-ajx.
require_one_of(PLASIM_NLAT      "${PLASIM_NLAT}"      32 48 64 96 128 160 192 256)
require_one_of(PLASIM_PRECISION "${PLASIM_PRECISION}" 4 8)
require_one_of(PLASIM_PARMODE   "${PLASIM_PARMODE}"   serial mpi omp)
require_one_of(PLASIM_FFT       "${PLASIM_FFT}"       fftmod fft991mod)
require_one_of(PLASIM_PLANET    "${PLASIM_PLANET}"    p_earth p_mars)

if(NOT PLASIM_NLEV MATCHES "^[0-9]+$" OR PLASIM_NLEV LESS 1)
    message(FATAL_ERROR "PLASIM_NLEV is '${PLASIM_NLEV}'; expected a positive integer.")
endif()
if(NOT PLASIM_NPRO MATCHES "^[0-9]+$" OR PLASIM_NPRO LESS 1)
    message(FATAL_ERROR "PLASIM_NPRO is '${PLASIM_NPRO}'; expected a positive integer.")
endif()

# NLPP = NLAT/NPRO is a `parameter` in plasimmod, so a non-dividing count builds
# a model whose bands do not tile the globe. Integer division would hide it.
math(EXPR _nlpp_rem "${PLASIM_NLAT} % ${PLASIM_NPRO}")
if(NOT _nlpp_rem EQUAL 0)
    message(FATAL_ERROR
        "PLASIM_NPRO=${PLASIM_NPRO} does not divide PLASIM_NLAT=${PLASIM_NLAT}. "
        "NLPP is NLAT/NPRO and is a parameter, so the bands would not tile the globe.")
endif()

math(EXPR PLASIM_NLON "${PLASIM_NLAT} * 2")

# ---------------------------------------------------------------------------
# Flags. The line arrives whole and CMake adds nothing to it.
# ---------------------------------------------------------------------------

separate_arguments(_fflags NATIVE_COMMAND "${PLASIM_FFLAGS}")
if(PLASIM_PRECISION EQUAL 8)
    # -fdefault-real-8 last, exactly where the rule set it had it: it is
    # MOST_PREC, appended after MOST_F90_OPTS on every compile line.
    list(APPEND _fflags -fdefault-real-8)
endif()
# The source path gfortran embeds in its runtime diagnostics -- "In file
# '<path>/radmod.f90', around line N" -- is the path it was handed. Building out
# of tree makes that absolute, which put 100 characters into every one of 1640
# such strings and grew .rodata by 35,776 bytes. Mapping it away restores the
# bare filename the old in-tree build produced, which is both a better message
# and what lets this build be compared byte for byte against that one.
list(APPEND _fflags "-ffile-prefix-map=${CMAKE_CURRENT_SOURCE_DIR}/src/=")
list(APPEND _fflags "-ffile-prefix-map=${CMAKE_CURRENT_BINARY_DIR}/=")

string(JOIN " " CMAKE_Fortran_FLAGS ${_fflags})
set(CMAKE_C_FLAGS "-O3")

# CMake puts -rdynamic on every executable link so that a program can resolve
# its own symbols at runtime. The model does not, and it costs an exported
# dynamic symbol table -- 36 KB here, and a binary that cannot be compared byte
# for byte against one the old rule set produced.
set(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "")
set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "")

# CMake appends per-configuration flags to the ones above unless they are empty.
foreach(cfg DEBUG RELEASE RELWITHDEBINFO MINSIZEREL)
    set(CMAKE_Fortran_FLAGS_${cfg} "")
    set(CMAKE_C_FLAGS_${cfg} "")
endforeach()

# ---------------------------------------------------------------------------
# The configuration-dependent source variants.
# ---------------------------------------------------------------------------

if(PLASIM_PARMODE STREQUAL "omp")
    set(_mpimod  mpimod_omp)
    set(_utilmod utilities_omp)
elseif(PLASIM_PARMODE STREQUAL "mpi")
    set(_mpimod  mpimod)
    set(_utilmod utilities)
else()
    set(_mpimod  mpimod_stub)
    set(_utilmod utilities_stub)
endif()

# resmod.f90 is generated, so it lives in the build directory and never in the
# source tree. The build system this replaced wrote it next to the sources,
# which is why a build for one resolution could be read as another's.
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/resmod.f90"
"      module resmod ! generated by CMakeLists.txt
      parameter(NLAT_ATM = ${PLASIM_NLAT})
      parameter(NLEV_ATM = ${PLASIM_NLEV})
      parameter(NPRO_ATM = ${PLASIM_NPRO})
      end module resmod

")

# ORDER IS THE LINK ORDER and is kept as the rule set had it, so the executable
# this produces can be compared byte for byte against one the old build made.
# The dependency EDGES are not listed anywhere here: Ninja scans the `use`
# statements and derives them, which is the point of the change.
set(_sources
    shtnsmod.f90
    ${_mpimod}.f90
    ${PLASIM_FFT}.f90
    guimod_stub.f90
    rainmod.f90
    simba.f90
    ${PLASIM_PLANET}.f90
    carbonmod.f90
    hurricanemod.f90
    ${_utilmod}.f90
    "${CMAKE_CURRENT_BINARY_DIR}/resmod.f90"
    plasim.f90
    plasimmod.f90
    specblock.f90
    calmod.f90
    gaussmod.f90
    legmod.f90
    outmod.f90
    miscmod.f90
    fluxmod.f90
    radmod.f90
    surfmod.f90
    landmod.f90
    glaciermod.f90
    seamod.f90
    icemod.f90
    oceanmod.f90
    restartmod.f90
    tracermod.f90
    tpcore.f90
    trc_routines.f90
    aeromod.f90
    aerocore.f90
    pumax_stub.c
    lsgmod.f90
    cpl_stub.f90
)
list(TRANSFORM _sources PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/src/"
     REGEX "^[^/]+$")

add_executable(plasim.x ${_sources})
set_target_properties(plasim.x PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules"
    OUTPUT_NAME plasim.x
    SUFFIX "")

# ---------------------------------------------------------------------------
# SHTns. shtnsmod includes shtns.f03, and the model links it in every
# configuration -- the MPI build cannot USE it, having no single address space
# holding every latitude, but it still has to resolve the symbols.
# ---------------------------------------------------------------------------

if(PLASIM_SHTNS_PREFIX STREQUAL "")
    message(FATAL_ERROR
        "PLASIM_SHTNS_PREFIX is not set. shtnsmod.f90 includes shtns.f03 and the "
        "model links SHTns in every configuration. exoplasim/scripts/build_shtns.sh "
        "installs it.")
endif()
if(NOT EXISTS "${PLASIM_SHTNS_PREFIX}/include/shtns.f03")
    message(FATAL_ERROR
        "no shtns.f03 under ${PLASIM_SHTNS_PREFIX}/include. Run "
        "exoplasim/scripts/build_shtns.sh. Failing here is deliberate: the build "
        "this replaced carried on without it and the model then failed at nshtns=1, "
        "which is midnight rather than now.")
endif()
file(GLOB _shtns_lib "${PLASIM_SHTNS_PREFIX}/lib/libshtns*.a")
if(NOT _shtns_lib)
    message(FATAL_ERROR "no libshtns*.a under ${PLASIM_SHTNS_PREFIX}/lib.")
endif()
list(GET _shtns_lib 0 _shtns_lib)

target_include_directories(plasim.x PRIVATE "${PLASIM_SHTNS_PREFIX}/include")
# -lgomp explicitly: libshtns_omp pulls fftw3_omp, which needs it, and the MPI
# build does not compile with -fopenmp so it would not otherwise be on the link
# line.
target_link_libraries(plasim.x PRIVATE "${_shtns_lib}" fftw3_omp fftw3 gomp m)

message(STATUS "plasim: NLAT=${PLASIM_NLAT} NLON=${PLASIM_NLON} NLEV=${PLASIM_NLEV} "
               "NPRO=${PLASIM_NPRO} ${PLASIM_PARMODE} real*${PLASIM_PRECISION}")
message(STATUS "plasim: ${CMAKE_Fortran_FLAGS}")
