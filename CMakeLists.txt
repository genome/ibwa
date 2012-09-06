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
set(INTERNAL_VERSION_FULL 0.5.0)
set(INTERNAL_VERSION_SHORT 0.5)
set(FULL_VERSION "${INTERNAL_VERSION_FULL}")
configure_file(version.h.in version.h)
include_directories(${PROJECT_BINARY_DIR})

# use these to created a versioned copy of the bwa binary, e.g., bwa-0.5.9
set(EXE_VERSION_SUFFIX "${INTERNAL_VERSION_SHORT}")
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

if (CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "-Wall -std=c++0x")
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
    bwaremap.cpp bwaremap.h
    translate_cigar.cpp translate_cigar.h
    filter_alignments.cpp filter_alignments.h
)
set(MAIN_SOURCES main.c main.h)

add_library(bwa ${LIB_SOURCES})
target_link_libraries(bwa ${MATH_LIBRARY} ${ZLIB_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
add_executable(${BWA_EXECUTABLE_NAME} ${MAIN_SOURCES})
target_link_libraries(${BWA_EXECUTABLE_NAME} bwa
    bwtgen ${MATH_LIBRARY} ${ZLIB_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
install(TARGETS ${BWA_EXECUTABLE_NAME} DESTINATION bin/)
install(FILES bwa.1 DESTINATION share/man/man1/
	RENAME ${BWA_EXECUTABLE_NAME}.1)

# unit tests
#find_package(GTest)
#if(GTEST_FOUND)
    #message("Google Test framework found, building unit tests")
    #enable_testing(true)
    #add_subdirectory(test)
#else()
    #message("Google Test framework not found, no tests will be built")
    #message("GTest is available at http://code.google.com/p/googletest/")
    #message("Ubuntu users can likely sudo apt-get install libgtest-dev")
#endif()



# Packaging ################################################################
set(ARCH "i686")
set(DEB_ARCH "i686")
if(${CMAKE_SIZEOF_VOID_P} MATCHES 8)
    set(ARCH "x86_64")
    set(DEB_ARCH "amd64")
endif ()

set(CPACK_PACKAGE_VERSION "${INTERNAL_VERSION_FULL}")
set(CPACK_PACKAGE_NAME "ibwa${INTERNAL_VERSION_SHORT}")
set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}_${INTERNAL_VERSION_FULL}_${DEB_ARCH}")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "A sequence aligner based on the BWT.")
set(CPACK_PACKAGE_VENDOR "TGI")
#set(CPACK_PACKAGE_VERSION ${FULL_VERSION})
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "genome-dev <genome-dev@genome.wustl.edu>>")
set(CPACK_SYSTEM_NAME "Linux-${ARCH}")
set(CPACK_TOPLEVEL_TAG "Linux-${ARCH}")
set(CPACK_DEBIAN_PACKAGE_SECTION science)
set(CPACK_DEBIAN_PACKAGE_PRIORITY optional)
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libc6 (>= 2.4), libgcc1 (>= 1:4.1.1-21)")

include(CPack)
