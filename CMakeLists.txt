cmake_minimum_required(VERSION 2.6)

# Example usage ############################################################
#
# To build a .deb (or .rpm) package
#  mkdir build; cd build; cmake ..; make; cpack -G {DEB|RPM}
#
# To do a debug build
#  mkdir build; cd build; cmake .. -DCMAKE_BUILD_TYPE=debug; make
#
# Note: in source builds are not allowed (so we don't clobber bwa's existing
# makefile).

project(bwa)

# Versioning ###############################################################
set(UPSTREAM_VERSION 0.5.9)
set(INTERNAL_VERSION_FULL  0.3.0)
set(INTERNAL_VERSION_SHORT 0.3)

# use these to created a versioned copy of the bwa binary, e.g., bwa-0.5.9
set(EXE_VERSION_SUFFIX "-${UPSTREAM_VERSION}-${INTERNAL_VERSION_SHORT}")
set(BWA_EXECUTABLE_NAME ibwa${EXE_VERSION_SUFFIX})

# Compiler options and libraries ###########################################
if (NOT CMAKE_BUILD_TYPE)
    message("CMAKE_BUILD_TYPE not specified, defaulting to RelWithDebInfo")
    set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif ()

if (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR
        "\nIn source builds are not allowed. Please run cmake from another directory, e.g.,\n"
        "mkdir build; cd build; cmake ..; make\n\n"
        "Before proceeding you should remove CMakeCache.txt and CMakeFiles/ from the source directory.")
endif ()

if (CMAKE_COMPILER_IS_GNUCC)
    set(CMAKE_C_FLAGS "-Wall")
endif ()

find_package(ZLIB REQUIRED)
include_directories(${ZLIB_INCLUDE_DIRS})
find_package(Threads)
if (CMAKE_USE_PTHREADS_INIT)
    add_definitions(-DHAVE_PTHREAD)
endif ()

if(UNIX)
    set(MATH_LIBRARY m)
endif()

# Build ####################################################################
add_subdirectory(bwt_gen)
set(LIB_SOURCES
    bamlite.c bamlite.h bntseq.c bntseq.h bwape.c bwase.c bwase.h bwaseqio.c
    bwt.c bwt.h bwt_lite.c bwt_lite.h bwtaln.c bwtaln.h bwtcache.c bwtcache.h
    bwtgap.c bwtgap.h bwtindex.c bwtio.c bwtmisc.c bwtsw2.h bwtsw2_aux.c
    bwtsw2_chain.c bwtsw2_core.c bwtsw2_main.c cs2nt.c is.c
    khash.h kseq.h ksort.h kstring.c kstring.h kvec.h
    simple_dp.c stdaln.c stdaln.h threadblock.c threadblock.h utils.c utils.h
    dbset.c dbset.h saiset.c saiset.h 
    byteorder.c byteorder.h
    bwapair.c bwapair.h
    bwasw.c bwasw.h
    bwaremap.c bwaremap.h
)
set(MAIN_SOURCES main.c main.h)

add_library(bwa ${LIB_SOURCES})
add_executable(${BWA_EXECUTABLE_NAME} ${MAIN_SOURCES})
target_link_libraries(${BWA_EXECUTABLE_NAME} bwa
    bwtgen ${MATH_LIBRARY} ${ZLIB_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
install(TARGETS ${BWA_EXECUTABLE_NAME} DESTINATION bin/)
install(FILES bwa.1 DESTINATION share/man/man1/
	RENAME ${BWA_EXECUTABLE_NAME}.1)


# Packaging ################################################################
set(ARCH "i686")
if(${CMAKE_SIZEOF_VOID_P} MATCHES 8)
    set(ARCH "x86_64")
endif ()

set(CPACK_PACKAGE_VERSION ${INTERNAL_VERSION_FULL})
set(CPACK_PACKAGE_NAME "ibwa-${UPSTREAM_VERSION}-${INTERNAL_VERSION_SHORT}")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "A sequence aligner based on the BWT.")
set(CPACK_PACKAGE_VENDOR "Travis Abbott")
#set(CPACK_PACKAGE_VERSION ${FULL_VERSION})
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Travis Abbott <tabbott@genome.wustl.edu>>")
set(CPACK_SYSTEM_NAME "Linux-${ARCH}")
set(CPACK_TOPLEVEL_TAG "Linux-${ARCH}")
set(CPACK_DEBIAN_PACKAGE_SECTION science)
set(CPACK_DEBIAN_PACKAGE_PRIORITY optional)
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libc6 (>= 2.4), libgcc1 (>= 1:4.1.1-21)")

include(CPack)
