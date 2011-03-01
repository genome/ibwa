cmake_minimum_required(VERSION 2.6)

project(bwa)

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

add_subdirectory(bwt_gen)

set(SOURCES
    bamlite.c
    bamlite.h
    bntseq.c
    bntseq.h
    bwape.c
    bwase.c
    bwase.h
    bwaseqio.c
    bwt.c
    bwt.h
    bwt_lite.c
    bwt_lite.h
    bwtaln.c
    bwtaln.h
    bwtcache.c
    bwtcache.h
    bwtgap.c
    bwtgap.h
    bwtindex.c
    bwtio.c
    bwtmisc.c
    bwtsw2.h
    bwtsw2_aux.c
    bwtsw2_chain.c
    bwtsw2_core.c
    bwtsw2_main.c
    cs2nt.c
    is.c
    khash.h
    kseq.h
    ksort.h
    kstring.c
    kstring.h
    kvec.h
    main.c
    main.h
    simple_dp.c
    stdaln.c
    stdaln.h
    threadblock.c
    threadblock.h
    utils.c
    utils.h
)

if (CMAKE_COMPILER_IS_GNUCC)
    set(CMAKE_C_FLAGS "-Wall")
endif ()

find_package(Threads)
if (CMAKE_USE_PTHREADS_INIT)
    add_definitions(-DHAVE_PTHREAD)
endif ()


add_executable(bwa ${SOURCES})
target_link_libraries(bwa m z ${CMAKE_THREAD_LIBS_INIT} bwtgen)
