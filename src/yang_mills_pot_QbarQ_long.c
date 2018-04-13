#ifndef YM_POLYCORR_LONG_C
#define YM_POLYCORR_LONG_C

#include"../include/macro.h"

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

#ifdef OPENMP_MODE
  #include<omp.h>
#endif

#include"../include/function_pointers.h"
#include"../include/gauge_conf.h"
#include"../include/geometry.h"
#include"../include/gparam.h"
#include"../include/random.h"

void real_main(char *in_file)
    {
    Gauge_Conf GC;
    Geometry geo;
    GParam param;

    int count;
    FILE *datafilep;
    time_t time1, time2;

    // to disable nested parallelism
    #ifdef OPENMP_MODE
      omp_set_nested(0);
    #endif

    // read input file
    readinput(in_file, &param);

    int tmp=param.d_size[1];
    for(count=2; count<STDIM; count++)
       {
       if(tmp!= param.d_size[count])
         {
         fprintf(stderr, "When using yang_mills_pot_QbarQ_long all the spatial sizes have to be of equal length.\n");
         exit(EXIT_FAILURE);
         }
       }

    // initialize random generator
    initrand(param.d_randseed);

    // open data_file
    init_data_file(&datafilep, &param);

    // initialize function_pointers
    init_function_pointers();

    // initialize geometry
    init_indexing_lexeo();
    init_geometry(&geo, &param);

    // initialize gauge configuration
    init_gauge_conf(&GC, &param);

    // initialize ml_polycorr arrays
    init_gauge_conf_polycorr(&GC, &param);

    // montecarlo starts
    time(&time1);
    if(param.d_start != 2) // NEW SIMULATION
      {
      for(count=1; count<param.d_measevery; count++)
         {
         update(&GC, &geo, &param);
         }

      // save configuration
      save_on_file(&GC, &param);
      // backup copy
      save_on_file_back(&GC, &param);

      // save ml polycorr arrays
      save_polycorr_on_file(&GC, &param, 0, 0);
      }
    else // CONTINUATION OF PREVIOUS SIMULATION
      {
      int count, tstart, iteration;

      // read multilevel stuff
      read_polycorr_from_file(&GC, &param, &tstart, &iteration);

      if(tstart<0) // update the conf, no multilevel
        {
        for(count=1; count<param.d_measevery; count++)
           {
           update(&GC, &geo, &param);
           }

        // save configuration
        save_on_file(&GC, &param);
        // backup copy
        save_on_file_back(&GC, &param);

        // save multilevel stuff
        save_polycorr_on_file(&GC, &param, 0, 0);
        }
      else // tstart>=0, perform multilevel
        {
        multilevel_pot_QbarQ_long(&GC,
                                  &geo,
                                  &param,
                                  tstart,
                                  param.d_ml_step[0],
                                  iteration);

        iteration+=1;
        if(iteration==param.d_ml_level0_repeat)
          {
          iteration=0;
          tstart+=param.d_ml_step[0];
          }

        if(tstart==param.d_size[0])
          {
          // print the measure
          perform_measures_pot_QbarQ_long(&GC, &param, datafilep);

          tstart=-1; // next time the conf will be updated, no multilevel
          }

        // save multilevel stuff
        save_polycorr_on_file(&GC, &param, tstart, iteration);
        }
      }
    time(&time2);
    // montecarlo end

    // close data file
    fclose(datafilep);

    // save configuration
    if(param.d_saveconf_back_every!=0)
      {
      save_on_file(&GC, &param);
      }

    // print simulation details
    print_parameters_polycorr(&param, time1, time2);

    // free gauge configuration
    end_gauge_conf(&GC, &param);

    // free ml_polycorr
    end_gauge_conf_polycorr(&GC);

    // free geometry
    free_geometry(&geo, &param);

    exit(EXIT_SUCCESS);
    }


int main (int argc, char **argv)
    {
    char in_file[50];

    if(argc != 2)
      {
      printf("\nPackage %s version %s\n", PACKAGE_NAME, PACKAGE_VERSION);
      printf("Claudio Bonati %s\n", PACKAGE_BUGREPORT);
      printf("Usage: %s input_file\n\n", argv[0]);

      printf("Compilation details:\n");
      printf("\tN_c (number of colors): %d\n", NCOLOR);
      printf("\tST_dim (space-time dimensionality): %d\n", STDIM);
      printf("\tNum_levels (number of levels): %d\n", NLEVELS);
      printf("\n");
      printf("\tINT_ALIGN: %s\n", QUOTEME(INT_ALIGN));
      printf("\tDOUBLE_ALIGN: %s\n", QUOTEME(DOUBLE_ALIGN));

      #ifdef DEBUG
        printf("\n\tDEBUG mode\n");
      #endif

      #ifdef OPENMP_MODE
        printf("\n\tusing OpenMP with %d threads\n", NTHREADS);
      #endif

      #ifdef OPT_MULTIHIT
        printf("\tcompiled for multihit optimization\n");
      #endif

      #ifdef OPT_MULTILEVEL
        printf("\tcompiled for multilevel optimization\n");
      #endif


      printf("\n");

      #ifdef __INTEL_COMPILER
        printf("\tcompiled with icc\n");
      #elif defined(__clang__)
        printf("\tcompiled with clang\n");
      #elif defined( __GNUC__ )
        printf("\tcompiled with gcc version: %d.%d.%d\n",
                __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
      #endif

      return EXIT_SUCCESS;
      }
    else
      {
      if(strlen(argv[1]) >= STD_STRING_LENGTH)
        {
        fprintf(stderr, "File name too long. Increse STD_STRING_LENGTH in include/macro.h\n");
        }
      else
        {
        strcpy(in_file, argv[1]);
        }
      }

    real_main(in_file);

    return EXIT_SUCCESS;
    }

#endif
