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


//simpler, faster than function above (require binary output (*.bin)
int mergeOutputBin(int n_file, const char ** string_list, int *n_exc, int *exc_L, int *exc_W, 
                   float *S_CA, float *S_AC, float *H_CA, float *H_AC, int verbose) {
           
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
  
  int Nsite = Excitation_L*Excitation_W;
  long int jj,size = NExcitation*NExcitation*Nsite*Nsite;
  float * data_read  = (float *)calloc((4*size),sizeof(float));
  
  float *phys_nCAm_averaged   = (float *)malloc((Nsite*Nsite*NExcitation*NExcitation)*sizeof(float));
  float *phys_nACm_averaged   = (float *)malloc((Nsite*Nsite*NExcitation*NExcitation)*sizeof(float));
  float *phys_nAHCm_averaged  = (float *)malloc((Nsite*Nsite*NExcitation*NExcitation)*sizeof(float));
  float *phys_nCHAm_averaged  = (float *)malloc((Nsite*Nsite*NExcitation*NExcitation)*sizeof(float));
  
  //fread(data_read, sizeof(double), 4*size, fp);
  
  if(verbose) fprintf(stdout, "\n\n%d sites to read per file\n\n",Nsite);
  for(ii = 0; ii<n_file; ii++){
    if(verbose) fprintf(stdout, "reading file '%s'\n",string_list[ii]);
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
  
  // untangling
  int nn,mm;
  for(ii=0;ii<Nsite;ii++) {
   for(jj=0;jj<Nsite;jj++) {
    for(nn=0;nn<NExcitation;nn++) {
     for(mm=0;mm<NExcitation;mm++) {
       S_CA[ii + Nsite*(mm + NExcitation*(jj + Nsite*nn))] += phys_nCAm_averaged[ii + Nsite*(jj + Nsite*(mm + NExcitation*nn))]; 
       S_AC[ii + Nsite*(mm + NExcitation*(jj + Nsite*nn))] += phys_nACm_averaged[ii + Nsite*(jj + Nsite*(mm + NExcitation*nn))]; 
       H_CA[ii + Nsite*(mm + NExcitation*(jj + Nsite*nn))] += phys_nCHAm_averaged[ii + Nsite*(jj + Nsite*(mm + NExcitation*nn))]; 
       H_AC[ii + Nsite*(mm + NExcitation*(jj + Nsite*nn))] += phys_nAHCm_averaged[ii + Nsite*(jj + Nsite*(mm + NExcitation*nn))]; 
     }
    }
   }
  }
  
  

  
  free(phys_nCAm_averaged);
  free(phys_nACm_averaged);
  free(phys_nAHCm_averaged);
  free(phys_nCHAm_averaged); 
  free(data_read);

  if(verbose) {
    printf("\nData read:\n");
    for(ii=0;ii<10;ii++) printf("% 4.5f % 4.5f % 4.5f % 4.5f \n",S_CA[ii], S_AC[ii], H_CA[ii], H_AC[ii]);
    printf("...\n");
  }
  
  return 0;
}



int greenFrom_e_U_Uinv_S(int sign, int Nw, double *w, int dim, double complex *e, double complex *U, double complex *Uinv_S, double Omega, double mu, double eta, double complex *g) //)H, double *S,)
{
  int ii,nn;
   
  for(ii=0;ii<Nw;ii++) {
    g[ii]=0.0;
    double complex z = w[ii] + I*eta + sign*Omega + mu;
    for(nn=0;nn<dim;nn++){
      g[ii] += U[nn]*Uinv_S[dim*nn] /(z-sign*e[nn]);
    }
  }

  return 0;
}



