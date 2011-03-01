project(libbwtgen)

set(SOURCES
    QSufSort.c
    QSufSort.h
    bwt_gen.c
    bwt_gen.h
)

add_library(bwtgen ${SOURCES})
