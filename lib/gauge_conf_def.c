#ifndef GAUGE_CONF_DEF_C
#define GAUGE_CONF_DEF_C

#include"../include/macro.h"

#ifdef HASH_MODE
  #include<openssl/md5.h>
#endif
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include"../include/function_pointers.h"
#include"../include/gparam.h"
#include"../include/geometry.h"
#include"../include/gauge_conf.h"
#include"../include/tens_prod.h"

void init_gauge_conf(Gauge_Conf *GC, GParam const * const param)
  {
  long r, j;
  int err;

  // allocate the local lattice
  err=posix_memalign((void**) &(GC->lattice), (size_t) DOUBLE_ALIGN, (size_t) param->d_volume * sizeof(GAUGE_GROUP *));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating the lattice! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  for(r=0; r<(param->d_volume); r++)
     {
     err=posix_memalign((void**)&(GC->lattice[r]), (size_t) DOUBLE_ALIGN, (size_t )param->d_stdim * sizeof(GAUGE_GROUP));
     if(err!=0)
       {
       fprintf(stderr, "Problems in allocating the lattice! (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     }

  #ifdef THETA_MODE
    init_clover_array(GC, param);
  #endif

  // initialize lattice
  if(param->d_start==0) // ordered start
    {
    GAUGE_GROUP aux1, aux2;
    one(&aux1);

    GC->update_index=0;

    for(r=0; r<(param->d_volume); r++)
       {
       for(j=0; j<param->d_stdim; j++)
          {
          rand_matrix(&aux2);
          times_equal_real(&aux2, 0.001);
          plus_equal(&aux2, &aux1);
          unitarize(&aux2);
          equal_dag(&(GC->lattice[r][j]), &aux2);
          }
       }
    }
  if(param->d_start==1)  // random start
    {
    GAUGE_GROUP aux1;

    GC->update_index=0;

    for(r=0; r<(param->d_volume); r++)
       {
       for(j=0; j<param->d_stdim; j++)
          {
          rand_matrix(&aux1);
          equal(&(GC->lattice[r][j]), &aux1);
          }
       }
    }

  if(param->d_start==2) // initialize from stored conf
    {
    read_gauge_conf(GC, param);
    }
  }


void read_gauge_conf(Gauge_Conf *GC, GParam const * const param)
  {
  FILE *fp;
  int i, dimension, tmp_i;
  int err, mu;
  long lex, si;
  GAUGE_GROUP matrix;
  #ifdef HASH_MODE
    char md5sum_new[2*MD5_DIGEST_LENGTH+1];
    char md5sum_old[2*MD5_DIGEST_LENGTH+1];
  #else
    char md5sum_old[2*STD_STRING_LENGTH+1]={0};
  #endif

  fp=fopen(param->d_conf_file, "r"); // open the configuration file
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_conf_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else // read the txt header of the configuration
    {
    err=fscanf(fp, "%d", &dimension);
    if(err!=1)
      {
      fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_conf_file, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
    if(dimension != param->d_stdim)
      {
      fprintf(stderr, "The space time dimension of the configuration (%d) does not coincide with the one of the global parameter (%d)\n",
              dimension, param->d_stdim);
      exit(EXIT_FAILURE);
      }

    for(i=0; i<param->d_stdim; i++)
       {
       err=fscanf(fp, "%d", &tmp_i);
       if(err!=1)
         {
         fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_conf_file, __FILE__, __LINE__);
         exit(EXIT_FAILURE);
         }
       if(tmp_i != param->d_size[i])
         {
         fprintf(stderr, "The size of the configuration lattice does not coincide with the one of the global parameter\n");
         exit(EXIT_FAILURE);
         }
       }

    err=fscanf(fp, "%ld %s\n", &(GC->update_index), md5sum_old);
    if(err!=2)
      {
      fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_conf_file, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }

    fclose(fp);
    }

  fp=fopen(param->d_conf_file, "rb"); // open the configuration file in binary
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_conf_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    // read again the header
    err=0;
    while(err!='\n')
         {
         err=fgetc(fp);
         }

    for(lex=0; lex<param->d_volume; lex++)
       {
       si=lex_to_si(lex, param);
       for(mu=0; mu<param->d_stdim; mu++)
          {
          read_from_binary_file_bigen(fp, &matrix);

          equal(&(GC->lattice[si][mu]), &matrix);
          }
       }
    fclose(fp);

    #ifdef HASH_MODE
      // compute the new md5sum and check for consistency
      compute_md5sum(md5sum_new, GC, param);
      if(strncmp(md5sum_old, md5sum_new, 2*MD5_DIGEST_LENGTH+1)!=0)
        {
        fprintf(stderr, "The computed md5sum %s does not match the stored %s (%s, %d)\n", md5sum_new, md5sum_old, __FILE__, __LINE__);
        exit(EXIT_FAILURE);
        }
    #endif
    }
  }


void end_gauge_conf(Gauge_Conf *GC, GParam const * const param)
  {
  long i;

  for(i=0; i<(param->d_volume); i++)
     {
     free(GC->lattice[i]);
     }
  free(GC->lattice);

  #ifdef THETA_MODE
    end_clover_array(GC, param);
  #endif
  }


// save a configuration in ILDG-like format
void write_conf_on_file_with_name(Gauge_Conf const * const GC,
                            GParam const * const param,
                            char const * const namefile)
  {
  long si, lex;
  int i, mu;
  #ifdef HASH_MODE
    char md5sum[2*MD5_DIGEST_LENGTH+1];
  #else
     char md5sum[2*STD_STRING_LENGTH+1]={0};
  #endif
  FILE *fp;

  #ifdef HASH_MODE
    compute_md5sum(md5sum, GC, param);
  #endif

  fp=fopen(namefile, "w"); // open the configuration file
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", namefile, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    fprintf(fp, "%d ", param->d_stdim);
    for(i=0; i<param->d_stdim; i++)
       {
       fprintf(fp, "%d ", param->d_size[i]);
       }
    fprintf(fp, "%ld %s\n", GC->update_index, md5sum);
    }
  fclose(fp);

  fp=fopen(namefile, "ab"); // open the configuration file in binary mode
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", namefile, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    for(lex=0; lex<param->d_volume; lex++)
       {
       si=lex_to_si(lex, param);
       for(mu=0; mu<param->d_stdim; mu++)
          {
          print_on_binary_file_bigen(fp, &(GC->lattice[si][mu]) );
          }
       }
    fclose(fp);
    }
  }


void write_conf_on_file(Gauge_Conf const * const GC, GParam const * const param)
  {
  write_conf_on_file_with_name(GC, param, param->d_conf_file);
  }


void write_conf_on_file_back(Gauge_Conf const * const GC, GParam const * const param)
  {
  char name[STD_STRING_LENGTH], aux[STD_STRING_LENGTH];
  static int counter=0;

  strcpy(name, param->d_conf_file);
  if(counter==0)
    {
    sprintf(aux, "_back0");
    }
  else
    {
    sprintf(aux, "_back1");
    }
  strcat(name, aux);

  write_conf_on_file_with_name(GC, param, name);

  counter=1-counter;
  }


// allocate GC and initialize with GC2
void init_gauge_conf_from_gauge_conf(Gauge_Conf *GC, Gauge_Conf const * const GC2, GParam const * const param) 
  {
  long r;
  int mu, err;

  // allocate the lattice
  err=posix_memalign((void**)&(GC->lattice), (size_t) DOUBLE_ALIGN, (size_t) param->d_volume * sizeof(GAUGE_GROUP *));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating the lattice! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  for(r=0; r<(param->d_volume); r++)
     {
     err=posix_memalign((void**)&(GC->lattice[r]), (size_t) DOUBLE_ALIGN, (size_t) param->d_stdim * sizeof(GAUGE_GROUP));
     if(err!=0)
       {
       fprintf(stderr, "Problems in allocating the lattice! (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     }

  #ifdef THETA_MODE
    init_clover_array(GC, param);
  #endif

  // initialize GC
  for(r=0; r<(param->d_volume); r++)
     {
     for(mu=0; mu<param->d_stdim; mu++)
        {
        equal(&(GC->lattice[r][mu]), &(GC2->lattice[r][mu]) );
        }
     }

  GC->update_index=GC2->update_index;
  }


// compute the md5sum of the configuration and save it in res, that is a char[2*MD5_DIGEST_LENGTH]
void compute_md5sum(char *res, Gauge_Conf const * const GC, GParam const * const param)
  {
  #ifdef HASH_MODE
    MD5_CTX mdContext;
    unsigned char c[MD5_DIGEST_LENGTH];
    long si, lex;
    GAUGE_GROUP matrix;
    int mu, k;

    MD5_Init(&mdContext);
    for(lex=0; lex<param->d_volume; lex++)
       {
       si=lex_to_si(lex, param);
       for(mu=0; mu<param->d_stdim; mu++)
          {
          equal(&matrix, &(GC->lattice[si][mu]));

          #if NCOLOR==1
            MD5_Update(&mdContext, &(matrix.comp), sizeof(double complex));
          #elif NCOLOR==2
            for(k=0; k<4; k++)
               {
               MD5_Update(&mdContext, &(matrix.comp[k]), sizeof(double));
               }
          #else
            for(k=0; k<NCOLOR*NCOLOR; k++)
               {
               MD5_Update(&mdContext, &(matrix.comp[k]), sizeof(double complex));
               }
          #endif
          }
       }
    MD5_Final(c, &mdContext);

    for(k = 0; k < MD5_DIGEST_LENGTH; k++)
       {
       sprintf(&(res[2*k]), "%02x", c[k]);
       }
  #else
    // just to avoid warning at compile time
    (void) res;
    (void) GC;
    (void) param;
  #endif
  }


// allocate the ml_polycorr arrays
void init_polycorr(Gauge_Conf *GC,
                              GParam const * const param)
  {
  int i, err;

  err=posix_memalign((void**)&(GC->ml_polycorr_ris), (size_t) DOUBLE_ALIGN, (size_t) NLEVELS *sizeof(TensProd *));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating ml_polycorr_ris (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    for(i=0; i<NLEVELS; i++)
       {
       err=posix_memalign((void**)&(GC->ml_polycorr_ris[i]), (size_t) DOUBLE_ALIGN, (size_t) param->d_space_vol *sizeof(TensProd));
       if(err!=0)
         {
         fprintf(stderr, "Problems in allocating ml_polycorr_ris[%d] (%s, %d)\n", i, __FILE__, __LINE__);
         exit(EXIT_FAILURE);
         }
       }
    }

  err=posix_memalign((void**)&(GC->ml_polycorr_tmp), (size_t) DOUBLE_ALIGN, (size_t) NLEVELS *sizeof(TensProd *));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating ml_polycorr_tmp (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    for(i=0; i<NLEVELS; i++)
       {
       err=posix_memalign((void**)&(GC->ml_polycorr_tmp[i]), (size_t) DOUBLE_ALIGN, (size_t) param->d_space_vol *sizeof(TensProd));
       if(err!=0)
         {
         fprintf(stderr, "Problems in allocating ml_polycorr_tmp[%d] (%s, %d)\n", i, __FILE__, __LINE__);
         exit(EXIT_FAILURE);
         }
       }
    }
  }


// free the ml_polycorr arrays
void end_polycorr(Gauge_Conf *GC)
  {
  int i;

  for(i=0; i<NLEVELS; i++)
     {
     free(GC->ml_polycorr_ris[i]);
     free(GC->ml_polycorr_tmp[i]);
     }
  free(GC->ml_polycorr_ris);
  free(GC->ml_polycorr_tmp);
  }


// save ml_polycorr[0] arrays on file
void write_polycorr_on_file(Gauge_Conf const * const GC,
                           GParam const * const param,
                           int tstart,
                           int iteration)
  {
  const int ml=0;
  long i;
  #ifdef HASH_MODE
    char md5sum[2*MD5_DIGEST_LENGTH+1];
  #else
    char md5sum[2*STD_STRING_LENGTH+1]={0};
  #endif
  FILE *fp;

  #ifdef HASH_MODE
    compute_md5sum_polycorr(md5sum, GC, param);
  #endif

  fp=fopen(param->d_ml_file, "w"); // open the configuration file
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    fprintf(fp, "%ld %d %d %s\n",
                param->d_space_vol,
                tstart,
                iteration,
                md5sum);
    }
  fclose(fp);

  fp=fopen(param->d_ml_file, "ab"); // open the configuration file in binary mode
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    for(i=0; i<(param->d_space_vol); i++)
       {
       print_on_binary_file_bigen_TensProd(fp, &(GC->ml_polycorr_ris[ml][i]));
       }
    for(i=0; i<(param->d_space_vol); i++)
       {
       print_on_binary_file_bigen_TensProd(fp, &(GC->ml_polycorr_tmp[ml][i]));
       }

    fclose(fp);
    }
  }


// read ml_polycorr[0] arrays from file
void read_polycorr_from_file(Gauge_Conf const * const GC,
                             GParam const * const param,
                             int *tstart,
                             int *iteration)
  {
  const int ml=0;
  long i, loc_space_vol;
  FILE *fp;
  #ifdef HASH_MODE
    char md5sum_new[2*MD5_DIGEST_LENGTH+1];
    char md5sum_old[2*MD5_DIGEST_LENGTH+1];
  #else
    char md5sum_old[2*STD_STRING_LENGTH+1]={0};
  #endif

  fp=fopen(param->d_ml_file, "r"); // open the multilevel file
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    i=fscanf(fp, "%ld %d %d %s\n", &loc_space_vol, tstart, iteration, md5sum_old);
    if(i!=4)
      {
      fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
    if(loc_space_vol != param->d_space_vol)
      {
      fprintf(stderr, "Error: space_vol in the multilevel file %s is different from the one in the input (%s, %d)\n",
              param->d_ml_file, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
    }
  fclose(fp);

  fp=fopen(param->d_ml_file, "rb"); // open the multilevel file in binary mode
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    // read again the header: loc_space_vol, tstart, iteration, hash
    i=0;
    while(i!='\n')
         {
         i=fgetc(fp);
         }

    for(i=0; i<(param->d_space_vol); i++)
       {
       read_from_binary_file_bigen_TensProd(fp, &(GC->ml_polycorr_ris[ml][i]));
       }
    for(i=0; i<(param->d_space_vol); i++)
       {
       read_from_binary_file_bigen_TensProd(fp, &(GC->ml_polycorr_tmp[ml][i]));
       }

    fclose(fp);
    }

  #ifdef HASH_MODE
    // compute the new md5sum and check for consistency
    compute_md5sum_polycorr(md5sum_new, GC, param);
    if(strncmp(md5sum_old, md5sum_new, 2*MD5_DIGEST_LENGTH+1)!=0)
      {
      fprintf(stderr, "The computed md5sum %s of the multilevel file does not match the stored %s\n", md5sum_new, md5sum_old);
      exit(EXIT_FAILURE);
      }
  #endif
  }


// compute the md5sum of the ml_polycorr[0] arrays and save it in res, that is a char[2*MD5_DIGEST_LENGTH]
void compute_md5sum_polycorr(char *res, Gauge_Conf const * const GC, GParam const * const param)
  {
  #ifdef HASH_MODE
    MD5_CTX mdContext;
    unsigned char c[MD5_DIGEST_LENGTH];
    long i;
    int n1, n2, n3, n4;
    const int ml=0;

    MD5_Init(&mdContext);

    for(i=0; i<(param->d_space_vol); i++)
       {
       for(n1=0; n1<NCOLOR; n1++)
          {
          for(n2=0; n2<NCOLOR; n2++)
             {
             for(n3=0; n3<NCOLOR; n3++)
                {
                for(n4=0; n4<NCOLOR; n4++)
                   {
                   MD5_Update(&mdContext, &((GC->ml_polycorr_ris[ml][i]).comp[n1][n2][n3][n4]), sizeof(double complex));
                   }
                }
             }
          }
       }

    for(i=0; i<(param->d_space_vol); i++)
       {
       for(n1=0; n1<NCOLOR; n1++)
          {
          for(n2=0; n2<NCOLOR; n2++)
             {
             for(n3=0; n3<NCOLOR; n3++)
                {
                for(n4=0; n4<NCOLOR; n4++)
                   {
                   MD5_Update(&mdContext, &((GC->ml_polycorr_tmp[ml][i]).comp[n1][n2][n3][n4]), sizeof(double complex));
                   }
                }
             }
          }
       }

    MD5_Final(c, &mdContext);

    for(i = 0; i < MD5_DIGEST_LENGTH; i++)
       {
       sprintf(&(res[2*i]), "%02x", c[i]);
       }
  #else
    // just to avoid warning at compile time
    (void) res;
    (void) GC;
    (void) param;
  #endif
  }


// allocate the ml_polycorr arrays
void init_polycorr_and_polyplaq(Gauge_Conf *GC,
                                GParam const * const param)
  {
  const unsigned long numplaq=(unsigned long) (param->d_stdim*(param->d_stdim-1))/2;
  int i, err;
  long r;

  init_polycorr(GC, param);

  err=posix_memalign((void**)&(GC->ml_polyplaq_ris), (size_t)DOUBLE_ALIGN, (size_t) NLEVELS *sizeof(TensProd **));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating ml_polyplaq_ris (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    for(i=0; i<NLEVELS; i++)
       {
       err=posix_memalign((void**)&(GC->ml_polyplaq_ris[i]), (size_t)DOUBLE_ALIGN, (size_t) param->d_space_vol *sizeof(TensProd *));
       if(err!=0)
         {
         fprintf(stderr, "Problems in allocating ml_polyplaq_ris[%d] (%s, %d)\n", i, __FILE__, __LINE__);
         exit(EXIT_FAILURE);
         }
       for(r=0; r<param->d_space_vol; r++)
          {
          err=posix_memalign((void**)&(GC->ml_polyplaq_ris[i][r]), (size_t)DOUBLE_ALIGN, (size_t) numplaq *sizeof(TensProd));
          if(err!=0)
            {
            fprintf(stderr, "Problems in allocating ml_polyplaq_ris[%d][%ld] (%s, %d)\n", i, r, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
            }
          }
       }
    }

  err=posix_memalign((void**)&(GC->ml_polyplaq_tmp), (size_t)DOUBLE_ALIGN, (size_t) NLEVELS *sizeof(TensProd **));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating ml_polyplaq_tmp (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    for(i=0; i<NLEVELS; i++)
       {
       err=posix_memalign((void**)&(GC->ml_polyplaq_tmp[i]), (size_t)DOUBLE_ALIGN, (size_t) param->d_space_vol *sizeof(TensProd *));
       if(err!=0)
         {
         fprintf(stderr, "Problems in allocating ml_polyplaq_tmp[%d] (%s, %d)\n", i, __FILE__, __LINE__);
         exit(EXIT_FAILURE);
         }
       for(r=0; r<param->d_space_vol; r++)
          {
          err=posix_memalign((void**)&(GC->ml_polyplaq_tmp[i][r]), (size_t)DOUBLE_ALIGN, (size_t) numplaq *sizeof(TensProd));
          if(err!=0)
            {
            fprintf(stderr, "Problems in allocating ml_polyplaq_tmp[%d][%ld] (%s, %d)\n", i, r, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
            }
          }
       }
    }
  }


// free the ml_polycorr arrays
void end_polycorr_and_polyplaq(Gauge_Conf *GC,
                               GParam const * const param)
  {
  int i;
  long r;

  end_polycorr(GC);

  for(i=0; i<NLEVELS; i++)
     {
     for(r=0; r<param->d_space_vol; r++)
        {
        free(GC->ml_polyplaq_ris[i][r]);
        free(GC->ml_polyplaq_tmp[i][r]);
        }
     free(GC->ml_polyplaq_ris[i]);
     free(GC->ml_polyplaq_tmp[i]);
     }
  free(GC->ml_polyplaq_ris);
  free(GC->ml_polyplaq_tmp);
  }


// allocate the clovers arrays
void init_clover_array(Gauge_Conf *GC,
                       GParam const * const param)
  {
  int i, err;
  long r;

  err=posix_memalign((void**)&(GC->clover_array), (size_t)DOUBLE_ALIGN, (size_t) param->d_volume *sizeof(GAUGE_GROUP **));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating clovers (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    for(r=0; r<param->d_volume; r++)
       {
       err=posix_memalign((void**)&(GC->clover_array[r]), (size_t)DOUBLE_ALIGN, 4*sizeof(GAUGE_GROUP *));
       if(err!=0)
         {
         fprintf(stderr, "Problems in allocating clovers[%ld] (%s, %d)\n", r, __FILE__, __LINE__);
         exit(EXIT_FAILURE);
         }
       for(i=0; i<4; i++)
          {
          err=posix_memalign((void**)&(GC->clover_array[r][i]), (size_t)DOUBLE_ALIGN, 4*sizeof(GAUGE_GROUP));
          if(err!=0)
            {
            fprintf(stderr, "Problems in allocating clovers[%ld][%d] (%s, %d)\n", r, i, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
            }
          }
       }
    }
  }


// free the clovers arrays
void end_clover_array(Gauge_Conf *GC,
                      GParam const * const param)
  {
  int i;
  long r;

  for(r=0; r<param->d_space_vol; r++)
     {
     for(i=0; i<4; i++)
        {
        free(GC->clover_array[r][i]);
        }
     free(GC->clover_array[r]);
     }
  free(GC->clover_array);
  }


#endif
