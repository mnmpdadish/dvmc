/*
*/
#include <complex.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <signal.h>

//
void interrupt_handler(int sig) {exit(0);}
//

int mergeOutput(int in_ncomp, int n_file, const char ** string_list, int NExcitation, int L, 
           int W, double *phys_CA_averaged, double *phys_AC_averaged, double *phys_CHA_averaged, double *phys_AHC_averaged) {
           
  // to be able to use ctrl-c
  signal(SIGINT, interrupt_handler); 
  //
    
  int NExcitation2 = NExcitation*NExcitation;
  
  int ii,dr,Nsite,nn,mm;
  FILE *fp;
  char *cerr;
  char ctmp[1024];
  char ctmp2[1024];
  //char fileName[] = "output/zvo_nCHAm_nAHCm_001.dat";
  int dummy1,dummy2,dummy3;
  double factor = 1./n_file;
  for(ii = 0; ii<n_file; ii++){
    printf("%s\n",string_list[ii]);
  }
  
  Nsite = W*L;
  
  fprintf(stdout, "\n\n%d sites to read per file\n\n",Nsite);
  for(ii = 0; ii<n_file; ii++){
    fprintf(stdout, "reading file '%s'\n",string_list[ii]);
    fp = fopen(string_list[ii], "r");
    if (fp == NULL) {
      fprintf(stdout, "error: no '%s' found.\n",string_list[ii]); 
      fclose(fp);
      exit(-1);
    }
    cerr = fgets(ctmp, sizeof(ctmp) / sizeof(char), fp); //skip header
    if (cerr != NULL) {
      for (dr = 0; dr < Nsite; dr++) {
       printf("%d  ", dr+1); fflush(stdout);
       for (nn = 0; nn < NExcitation; nn++) {
        for (mm = 0; mm < NExcitation; mm++) {
         cerr = fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp);
         double val1,val2,val3,val0;
         sscanf(ctmp2, "%d %d %d %le %le %le %le\n", &dummy1,&dummy2,&dummy3, &val0,&val1,&val2,&val3);
         if((dummy1!=dr)||(dummy2!=nn)||(dummy3!=mm)){
          fprintf(stdout, "\n\nerror: file wrong.\n"); 
          fclose(fp);
          exit(-1);           
         }
         phys_CA_averaged [nn+NExcitation*mm + dr*NExcitation2] += factor*val0; 
         phys_AC_averaged [nn+NExcitation*mm + dr*NExcitation2] += factor*val1; 
         phys_CHA_averaged[nn+NExcitation*mm + dr*NExcitation2] += factor*val2; 
         phys_AHC_averaged[nn+NExcitation*mm + dr*NExcitation2] += factor*val3; 
         //fprintf(stdout,"%d %d %d %f %f %f\n", dummy1,dummy2,dummy3, val1,val2,val3,val4);
        }
       }
      }
      printf("\n\n");
    }
    fclose(fp);
  }
  
  return 0;

}





//simpler, faster than function above (require binary output (*.bin)
int mergeOutputBin(int n_file, const char ** string_list, int *n_exc, int *exc_L, int *exc_W, 
                   float *phys_nCAm_averaged, float *phys_nACm_averaged, float *phys_nCHAm_averaged, float *phys_nAHCm_averaged) {
           
  // to be able to use ctrl-c
  signal(SIGINT, interrupt_handler); 
  //
  
  double factor = 1./n_file;
  int ii;
  
  for(ii = 0; ii<n_file; ii++){
    printf("%s\n",string_list[ii]);
  }

  FILE *fp;
  fp = fopen(string_list[0], "rb");
  if (fp == NULL) {
    fprintf(stdout, "error: no '%s' found.\n",string_list[ii]); 
    fclose(fp);
    exit(-1);
  }
  
  int NExcitation, Excitation_L, Excitation_W;
  if(fread(&NExcitation,  sizeof(int), 1, fp) ==0){
    fprintf(stdout, "error: reading NExcitation in file '%s'.\n",string_list[ii]); 
    exit(-1);
  }
  if(fread(&Excitation_L,  sizeof(int), 1, fp) ==0){
    fprintf(stdout, "error: reading Excitation_L in file '%s'.\n",string_list[ii]); 
    exit(-1);
  }
  if(fread(&Excitation_W,  sizeof(int), 1, fp) ==0){
    fprintf(stdout, "error: reading Excitation_W in file '%s'.\n",string_list[ii]); 
    exit(-1);
  }
  
  fclose(fp);
  
  long int jj,size = NExcitation*NExcitation*Excitation_L*Excitation_W;
  int Nsite = Excitation_L * Excitation_W;
  float * data_read  = (float *)calloc((4*size),sizeof(float));
  
  //fread(data_read, sizeof(double), 4*size, fp);
  
  fprintf(stdout, "\n\n%d sites to read per file\n\n",Nsite);
  for(ii = 0; ii<n_file; ii++){
    fprintf(stdout, "reading file '%s'\n",string_list[ii]);
    fp = fopen(string_list[ii], "rb");
    if (fp == NULL) {
      fprintf(stdout, "error: no '%s' found.\n",string_list[ii]); 
      fclose(fp);
      exit(-1);
    }
    
    int NExcitation_tmp, Excitation_L_tmp, Excitation_W_tmp;
    if(fread(&NExcitation_tmp,  sizeof(int), 1, fp) ==0){
      fprintf(stdout, "error: reading NExcitation in file '%s'.\n",string_list[ii]); 
      exit(-1);
    }
    if(fread(&Excitation_L_tmp,  sizeof(int), 1, fp) ==0){
      fprintf(stdout, "error: reading Excitation_L in file '%s'.\n",string_list[ii]); 
      exit(-1);
    }
    if(fread(&Excitation_W_tmp,  sizeof(int), 1, fp) ==0){
      fprintf(stdout, "error: reading Excitation_W in file '%s'.\n",string_list[ii]); 
      exit(-1);
    }
    
    if(NExcitation_tmp!=NExcitation || Excitation_L_tmp!=Excitation_L || Excitation_W_tmp!=Excitation_W){
      fprintf(stdout, "error: incoherent header in file '%s'.\n",string_list[ii]); 
      exit(-1);
    }

    if(fread(data_read, sizeof(float), 4*size, fp) ==0){
      fprintf(stdout, "error: reading Excitation_W in file '%s'.\n",string_list[ii]); 
      exit(-1);
    }
    
    for(jj=0;jj<size;jj++){
      phys_nCAm_averaged [jj] += factor*data_read[jj+size*0]; 
      phys_nACm_averaged [jj] += factor*data_read[jj+size*1];
      phys_nCHAm_averaged[jj] += factor*data_read[jj+size*2]; 
      phys_nAHCm_averaged[jj] += factor*data_read[jj+size*3]; 
    }
    fclose(fp);
  }
  
  
  free(data_read);

  printf("\nData read:\n");
  for(ii=0;ii<10;ii++) printf("% 4.5f % 4.5f % 4.5f % 4.5f \n",phys_nCAm_averaged[ii], phys_nACm_averaged[ii], phys_nCHAm_averaged[ii], phys_nAHCm_averaged[ii]);
  printf("...\n");
  
  return 0;
}


float sumEigenEquation(int n_exc, float *e_re, float *e_im,  
                                  float *u_re, float *u_im,  
                                  float *us_re,float *us_im, 
                                  float eta, float omega)
{
  int nn;
  double complex g=0.0; 
  for(nn=0;nn<n_exc;nn++){
    
    double complex e = e_re[nn] + I*e_im[nn];
    double complex u = u_re[nn] + I*u_im[nn];
    double complex us = us_re[nn] + I*us_im[nn];
    double complex z = omega+I*eta;
    g += u * us / (z + e);
  }    
  //printf("% 4.5f \n",cimag(g));
  return cimag(g);
}




