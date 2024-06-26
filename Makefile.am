bin_PROGRAMS = conf_check \
               debug_rng \
               debug_son \
               debug_su2 \
               debug_su2_adj \
               debug_su2sun \
               debug_sun \
               debug_u1 \
               debug_vecs \
               yang_mills_higgs \
               yang_mills_local \
               yang_mills_local_fundadj \
               yang_mills_polycorr \
               yang_mills_polycorr_long \
               yang_mills_polycorr_higgs \
               yang_mills_polycorr_higgs_long \
               yang_mills_polycorradj \
               yang_mills_polycorradj_long \
               yang_mills_t0 \
               yang_mills_tracedef \
			   yang_mills_tracedef_polyakov \
               yang_mills_tracedef_plaquettes \
               yang_mills_tube_conn \
               yang_mills_tube_conn_long \
               yang_mills_tube_disc \
               yang_mills_tube_disc_long \
               yang_mills_tubeadj_disc \
               yang_mills_tubeadj_disc_long \
	       	   test_lexeo_orth 

align = ./include/aligncheck.h \
        ./lib/aligncheck.c

dSFMT = ./include/dSFMT-common.h \
        ./include/dSFMT.h \
        ./include/dSFMT-params11213.h \
        ./include/dSFMT-params1279.h \
        ./include/dSFMT-params132049.h \
        ./include/dSFMT-params19937.h \
        ./include/dSFMT-params216091.h \
        ./include/dSFMT-params2203.h \
        ./include/dSFMT-params4253.h \
        ./include/dSFMT-params44497.h \
        ./include/dSFMT-params521.h \
        ./include/dSFMT-params86243.h \
        ./include/dSFMT-params.h \
        ./lib/dSFMT.c

endian = ./include/endianness.h  \
         ./lib/endianness.c 

flavourmatrix = ./include/flavour_matrix.h \
                ./lib/flavour_matrix.c

funcpoint = ./include/function_pointers.h \
            ./lib/function_pointers.c

gauge_conf = ./include/gauge_conf.h \
             ./lib/gauge_conf_def.c \
             ./lib/gauge_conf_meas.c \
             ./lib/gauge_conf_meas_multilevel.c \
             ./lib/gauge_conf_multilevel.c \
             ./lib/gauge_conf_upd.c 

geometry = ./include/geometry.h \
           ./lib/geometry.c

gparam = ./include/gparam.h \
         ./lib/gparam.c

son = ./include/son.h \
      ./include/son_monopoles.h \
      ./include/son_upd.h \
      ./lib/son.c \
      ./lib/son_monopoles.c \
      ./lib/son_upd.c

random = ./include/myrngparam.h \
         ./include/random.h \
         ./lib/random.c 

su2 = ./include/su2.h \
      ./include/su2_monopoles.h \
      ./include/su2_upd.h \
      ./lib/su2.c \
      ./lib/su2_monopoles.c \
      ./lib/su2_upd.c

sun = ./include/sun.h \
      ./include/sun_monopoles.h \
      ./include/sun_upd.h \
      ./lib/sun.c \
      ./lib/sun_monopoles.c \
      ./lib/sun_upd.c 

tens_prod = ./include/tens_prod.h \
            ./include/tens_prod_adj.h \
            ./lib/tens_prod.c \
            ./lib/tens_prod_adj.c

u1 = ./include/u1.h \
     ./include/u1_upd.h \
     ./include/u1_monopoles.h \
     ./lib/u1.c \
     ./lib/u1_upd.c \
     ./lib/u1_monopoles.c

everything =  $(align) \
              $(dSFMT)  \
              $(endian) \
              $(flavourmatrix) \
              $(funcpoint) \
              $(gauge_conf) \
              $(geometry) \
              $(gparam) \
              ./include/macro.h \
              $(son) \
              $(random) \
              $(su2) \
              $(sun) \
              $(tens_prod) \
              $(u1) 

conf_check_SOURCES=./src/conf_check.c \
                   $(align) \
                   $(dSFMT) \
                   $(endian) \
                   $(flavourmatrix) \
                   $(funcpoint) \
                   $(gparam) \
                   $(random) \
                   $(son) \
                   $(su2) \
                   $(sun) \
                   $(tens_prod) \
                   $(u1) 

debug_rng_SOURCES=./src/debug_rng.c \
                  $(align) \
                  $(dSFMT) \
                  $(random)

debug_son_SOURCES=./src/debug_son.c \
                  $(align) \
                  $(dSFMT) \
                  $(endian) \
                  $(flavourmatrix) \
                  ./include/macro.h \
                  $(random) \
                  $(son) \
                  $(tens_prod) 

debug_su2_SOURCES=./src/debug_su2.c \
                  $(align) \
                  $(dSFMT) \
                  $(endian) \
                  $(flavourmatrix) \
                  $(gparam) \
                  ./include/macro.h \
                  $(random) \
                  $(su2) \
                  $(tens_prod) 

debug_su2_adj_SOURCES=./src/debug_su2_adj.c \
                      $(align) \
                      $(dSFMT) \
                      $(endian) \
                      $(flavourmatrix) \
                      ./include/macro.h \
                      $(random) \
                      $(su2) \
                      $(tens_prod) 

debug_su2sun_SOURCES=./src/debug_su2sun.c \
                      $(align) \
                      $(dSFMT) \
                      $(endian) \
                      $(flavourmatrix) \
                      ./include/macro.h \
                      $(random) \
                      $(su2) \
                      $(sun) \
                      $(tens_prod) 

debug_sun_SOURCES=./src/debug_sun.c \
                  $(align) \
                  $(dSFMT) \
                  $(endian) \
                  $(flavourmatrix) \
                  ./include/macro.h \
                  $(random) \
                  $(su2) \
                  $(sun) \
                  $(tens_prod) 

debug_u1_SOURCES=./src/debug_u1.c \
                  $(align) \
                  $(dSFMT) \
                  $(endian) \
                  $(gparam) \
                  ./include/macro.h \
                  $(random) \
                  $(tens_prod) \
                  $(u1)

debug_vecs_SOURCES=./src/debug_vecs.c \
                   $(align) \
                   $(dSFMT) \
                   $(endian) \
                   $(funcpoint) \
                   $(gparam) \
                   ./include/macro.h \
                   $(random) \
                   $(son) \
                   $(su2) \
                   $(sun) \
                   $(tens_prod) \
                   $(u1)

yang_mills_higgs_SOURCES=./src/yang_mills_higgs.c \
                          $(everything)

yang_mills_local_SOURCES=./src/yang_mills_local.c \
                          $(everything)

yang_mills_local_fundadj_SOURCES=./src/yang_mills_local_fundadj.c \
                                 $(everything)

yang_mills_polycorr_SOURCES=./src/yang_mills_polycorr.c \
                            $(everything)

yang_mills_polycorr_long_SOURCES=./src/yang_mills_polycorr_long.c \
                                 $(everything)

yang_mills_polycorr_higgs_SOURCES=./src/yang_mills_polycorr_higgs.c \
                                  $(everything)

yang_mills_polycorr_higgs_long_SOURCES=./src/yang_mills_polycorr_higgs_long.c \
                                  $(everything)

yang_mills_polycorradj_SOURCES=./src/yang_mills_polycorradj.c \
                               $(everything)

yang_mills_polycorradj_long_SOURCES=./src/yang_mills_polycorradj_long.c \
                                    $(everything)

yang_mills_t0_SOURCES=./src/yang_mills_t0.c \
                      $(everything)

yang_mills_tracedef_SOURCES=./src/yang_mills_tracedef.c \
                            $(everything)

yang_mills_tracedef_polyakov_SOURCES=./src/yang_mills_tracedef_polyakov.c \
							$(everything)

yang_mills_tracedef_plaquettes_SOURCES=./src/yang_mills_tracedef_plaquettes.c \
                                          $(everything)
                                          
yang_mills_tube_conn_SOURCES=./src/yang_mills_tube_conn.c \
                             $(everything)

yang_mills_tube_conn_long_SOURCES=./src/yang_mills_tube_conn_long.c \
                                  $(everything)

yang_mills_tube_disc_SOURCES=./src/yang_mills_tube_disc.c \
                             $(everything)
 
yang_mills_tube_disc_long_SOURCES=./src/yang_mills_tube_disc_long.c \
                                  $(everything)

yang_mills_tubeadj_disc_SOURCES=./src/yang_mills_tubeadj_disc.c \
                                $(everything)

yang_mills_tubeadj_disc_long_SOURCES=./src/yang_mills_tubeadj_disc_long.c \
                                     $(everything)

test_lexeo_orth_SOURCES=./tests/test_lexeo_orth.c \
				$(everything)
   
AM_CFLAGS = -std=c99 -O2 -Wall -Wextra -Werror -pedantic -Wconversion @OPENMP_CFLAGS@

CLEANFILES = dati.dat log.dat dati.log conf.dat* ml.dat template_input.in higgs_conf.dat* mon.dat
