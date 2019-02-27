project(cephes)

add_library(cephes STATIC
  bdtr.c
  chdtr.c
  const.c
  fdtr.c
  gamma.c
  igam.c
  igami.c
  incbet.c
  incbi.c
  mtherr.c
  ndtr.c
  ndtri.c
  polevl.c
  stdtr.c
)
