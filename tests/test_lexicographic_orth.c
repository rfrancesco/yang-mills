#ifndef YM_TRACEDEF_C
#define YM_TRACEDEF_C

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


    // to disable nested parallelism
    #ifdef OPENMP_MODE
      omp_set_nested(0);
    #endif

    // read input file
    readinput(in_file, &param);

    // initialize random generator
    initrand(param.d_randseed);

    // initialize geometry
    init_indexing_lexeo();
    init_geometry(&geo, &param);

    // initialize gauge configuration
    init_gauge_conf(&GC, &param);

    printf("Testing orthogonal volumes...");
    long prod=1;
    long prod_test=1;
    for(int axis=0; axis<STDIM; axis++)
      {
      prod*=param.d_orth_vol[axis];
      }
    for(int i=0; i<STDIM-1; i++)
      {
      prod_test*=param.d_volume;
      }
    if(prod == prod_test)
      {
      printf("OK!\n");
      }
    else
      {
      printf("ERROR! \n");
      exit(EXIT_FAILURE);
      }

    printf("Testing consistency of the (rsp, t, axis) <-> r <-> (x_1...x_n, t) mapping. \n");
    printf("This should loop over every combination of (axis, rsp, t). This should be checked visually. \n");
    printf("VALID/ERROR: consistency check, by transforming coordinates both ways. ERROR should never be displayed, by design. \n");
    printf("OK_SISP: consistency check for axis=0, by checking that rsp is indeed the sisp used in the original code. \n");

    for(int axis=0; axis<STDIM; axis++)
      {
      printf("BEGINNING TEST: (rsp, t, axis=%d)\n", axis);
      for(long rsp=0; rsp<param.d_orth_vol[axis]; rsp++)
        {
        printf("Testing rsp = %ld, relative to axis %d. Varying t...\n", rsp, axis);
        for(int t=0; t<param.d_size[axis]; t++)
          {
          int cartcoord[STDIM];
          long r=siorth_and_par_to_si(&geo, rsp, t, axis);
          si_to_cart(cartcoord, r, &param);
          printf("(%ld, %d, %d) - Cart.:", rsp, t, axis);
          for(int i=0; i<STDIM; i++)
            {
            printf("%d ",cartcoord[i]);
            }
          printf(" | ");
          long rsp_test;
          int t_test;
          si_to_siorth_and_par(&rsp_test, &t_test, axis, r, &geo);
          if (rsp == rsp_test && t == t_test)
            {
            printf("VALID ");
            }
          else
            {
            printf("ERROR\n");
            exit(EXIT_FAILURE);
            }
	  if (axis == 0)
	    {
            long rsp_test_bonati;
	    int t_test_bonati;
	    si_to_sisp_and_t(&rsp_test_bonati, &t_test_bonati, &geo, r);
	    if (rsp == rsp_test_bonati && t == t_test_bonati)
               {
	       printf("OK_SISP: (sisp, t) == (%ld, %d) ", rsp_test_bonati, t_test_bonati);
	       }
	    else
	       {
	       printf("ERROR_SISP\n");
	       exit(EXIT_FAILURE);
	       }
            }
	  printf("\n");
          }
        }
      printf("COMPLETED test for axis=%d\n", axis);
      }


    // free gauge configuration
    free_gauge_conf(&GC, &param);

    // free geometry
    free_geometry(&geo, &param);
    }


void print_template_input(void)
  {
  FILE *fp;

  fp=fopen("template_input.in", "w");

  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file template_input.in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    fprintf(fp, "size 4 4 4 4\n");
    fprintf(fp,"\n");
    fprintf(fp, "beta 5.705\n");
    fprintf(fp, "htracedef  1.1\n");
    fprintf(fp, "theta 1.5\n");
    fprintf(fp,"\n");
    fprintf(fp, "sample    10\n");
    fprintf(fp, "thermal   0\n");
    fprintf(fp, "overrelax 5\n");
    fprintf(fp, "measevery 1\n");
    fprintf(fp, "monomeas 0   # 1=monopoles measures are performed\n");
    fprintf(fp,"\n");
    fprintf(fp, "start                   0  # 0=ordered  1=random  2=from saved configuration\n");
    fprintf(fp, "saveconf_back_every     5  # if 0 does not save, else save backup configurations every ... updates\n");
    fprintf(fp, "saveconf_analysis_every 5  # if 0 does not save, else save configurations for analysis every ... updates\n");
    fprintf(fp, "\n");
    fprintf(fp, "epsilon_metro    0.25  #distance from the identity of the random matrix for metropolis\n");
    fprintf(fp,"\n");
    fprintf(fp, "coolsteps  3     # number of cooling steps to be used\n");
    fprintf(fp, "coolrepeat 5     # number of times 'coolsteps' are repeated\n");
    fprintf(fp,"\n");
    fprintf(fp, "#output files\n");
    fprintf(fp, "conf_file  conf.dat\n");
    fprintf(fp, "data_file  dati.dat\n");
    fprintf(fp, "mon_file   mon.dat\n");
    fprintf(fp, "log_file   log.dat\n");
    fprintf(fp, "\n");
    fprintf(fp, "randseed 0    #(0=time)\n");
    fclose(fp);
    }
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
      printf("\tGGROUP: %s\n", QUOTEME(GGROUP));
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

      #ifdef THETA_MODE
        printf("\n\tusing imaginary theta\n");
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

      print_template_input();


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
