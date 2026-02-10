# cmake/FindSILO.cmake
#
# Finds the Silo library (LLNL).
#
# Result variables:
#   SILO_FOUND
#   SILO_VERSION
#   SILO_INCLUDE_DIRS
#   SILO_LIBRARIES
#
# Imported targets:
#   SILO::silo
# Optionally provides:
#   hpx_SILO::silo  (alias to SILO::silo) if SILO_CREATE_HPX_ALIAS is ON
#
# User hints:
#   SILO_ROOT        - installation prefix (preferred)
#   SILO_DIR         - installation prefix
#   SILO_INCLUDE_DIR - explicit include directory
#   SILO_LIBRARY     - explicit library path

cmake_minimum_required(VERSION 3.16)

include(FindPackageHandleStandardArgs)

set(_silo_hints "")
if (DEFINED SILO_ROOT)
  list(APPEND _silo_hints "${SILO_ROOT}")
endif()
if (DEFINED SILO_DIR)
  list(APPEND _silo_hints "${SILO_DIR}")
endif()

# --- First try pkg-config (best for apt installs) ---
find_package(PkgConfig QUIET)
if (PkgConfig_FOUND)
  # Ubuntu provides /usr/lib/x86_64-linux-gnu/pkgconfig/silo.pc
  pkg_check_modules(PC_SILO QUIET silo)

  if (PC_SILO_FOUND)
    set(SILO_VERSION "${PC_SILO_VERSION}")

    # pkg-config may return only -lsiloh5; keep both the flags and the dirs.
    set(SILO_INCLUDE_DIRS "${PC_SILO_INCLUDE_DIRS}")
    set(SILO_LIBRARY_DIRS "${PC_SILO_LIBRARY_DIRS}")

    # Prefer explicit library discovery so we can make an imported target cleanly.
    # Try common names. On Ubuntu the library is typically "siloh5".
    find_library(SILO_LIBRARY
      NAMES siloh5 silo
      HINTS ${SILO_LIBRARY_DIRS} ${_silo_hints}
      PATH_SUFFIXES lib lib64 lib/x86_64-linux-gnu
    )

    # If pkg-config didn't provide include dirs, still try to locate silo.h.
    if (NOT SILO_INCLUDE_DIRS)
      find_path(SILO_INCLUDE_DIR
        NAMES silo.h
        HINTS ${PC_SILO_INCLUDEDIR} ${_silo_hints}
        PATH_SUFFIXES include
      )
      if (SILO_INCLUDE_DIR)
        set(SILO_INCLUDE_DIRS "${SILO_INCLUDE_DIR}")
      endif()
    endif()
  endif()
endif()

# --- Fallback: direct search (manual builds, non-pkg-config systems) ---
if (NOT SILO_INCLUDE_DIRS)
  if (DEFINED SILO_INCLUDE_DIR)
    set(SILO_INCLUDE_DIRS "${SILO_INCLUDE_DIR}")
  else()
    find_path(SILO_INCLUDE_DIR
      NAMES silo.h
      HINTS ${_silo_hints}
      PATH_SUFFIXES include
    )
    if (SILO_INCLUDE_DIR)
      set(SILO_INCLUDE_DIRS "${SILO_INCLUDE_DIR}")
    endif()
  endif()
endif()

if (NOT SILO_LIBRARY)
  if (DEFINED SILO_LIBRARY)
    # user provided an explicit path; keep it
  else()
    find_library(SILO_LIBRARY
      NAMES siloh5 silo
      HINTS ${_silo_hints}
      PATH_SUFFIXES lib lib64 lib/x86_64-linux-gnu
    )
  endif()
endif()

set(SILO_LIBRARIES "${SILO_LIBRARY}")

# Handle standard args and mark advanced variables
find_package_handle_standard_args(
  SILO
  REQUIRED_VARS SILO_LIBRARY SILO_INCLUDE_DIRS
  VERSION_VAR SILO_VERSION
)

mark_as_advanced(SILO_LIBRARY SILO_INCLUDE_DIR SILO_ROOT SILO_DIR)

# --- Create imported target SILO::silo ---
if (SILO_FOUND AND NOT TARGET SILO::silo)
  add_library(SILO::silo UNKNOWN IMPORTED)

  set_target_properties(SILO::silo PROPERTIES
    IMPORTED_LOCATION "${SILO_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${SILO_INCLUDE_DIRS}"
  )

  # If pkg-config reported extra link flags (rarely needed), propagate them.
  # Note: Avoid dumping "-L..." and "-l..." into INTERFACE_LINK_LIBRARIES.
  if (PC_SILO_FOUND AND PC_SILO_LDFLAGS_OTHER)
    set_property(TARGET SILO::silo APPEND PROPERTY
      INTERFACE_LINK_OPTIONS "${PC_SILO_LDFLAGS_OTHER}"
    )
  endif()

  # Some distros list private deps that may be needed for static linking.
  if (PC_SILO_FOUND AND PC_SILO_LIBS_PRIVATE)
    # PC_SILO_LIBS_PRIVATE may contain -l/-L flags; best-effort:
    separate_arguments(_silo_private_libs NATIVE_COMMAND "${PC_SILO_LIBS_PRIVATE}")
    set_property(TARGET SILO::silo APPEND PROPERTY
      INTERFACE_LINK_LIBRARIES "${_silo_private_libs}"
    )
  endif()
endif()

# Optional HPX alias expected by some HPX helper macros
option(SILO_CREATE_HPX_ALIAS "Create hpx_SILO::silo alias target" ON)
if (SILO_FOUND AND SILO_CREATE_HPX_ALIAS AND NOT TARGET hpx_SILO::silo)
  add_library(hpx_SILO::silo ALIAS SILO::silo)
endif()

