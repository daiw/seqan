# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /core/apps/index_saBwt_performance/CMakeLists.txt
#
# CMakeLists.txt file for index_saBwt_performance.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_core_apps_index_saBwt_performance)
message (STATUS "Configuring core/apps/index_saBwt_performance")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
#set (SEQAN_FIND_DEPENDENCIES NONE)
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Update the list of file names below if you add source files to your application.
add_executable (index_saBwt_performance index_saBwt_performance.cpp index_saBwt_performance.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (index_saBwt_performance ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (index_saBwt_performance)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install index_saBwt_performance in ${PREFIX}/bin directory
install (TARGETS index_saBwt_performance
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/index_saBwt_performance for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES example/fasta1.fa
#               example/fasta2.fa
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

#seqan_add_app_test (index_saBwt_performance)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:index_saBwt_performance")
  set (CPACK_PACKAGE_NAME "index_saBwt_performance")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "index_saBwt_performance")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "David Iwanowitsch <iwanowit@inf.fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "David Iwanowitsch <iwanowit@inf.fu-berlin.de>")

  seqan_configure_cpack_app (index_saBwt_performance "index_saBwt_performance")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:index_saBwt_performance")

