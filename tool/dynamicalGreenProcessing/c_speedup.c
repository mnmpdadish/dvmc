/*
*/
#include <complex.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <signal.h>

//#include "matrix/mat.c"





extern void zgetrf_(int*, int* , double complex*, int*, int*, int*);
extern void zgetri_(int*, double complex*, int*, int*, double complex*, int*, int*);
extern void zgeev_(char*, char*, int *, double complex*, int*, double complex*, double complex*, int*, double complex*,int*, double complex*, int*, double*, int*);
extern void zgemm_(char*, char*, int*, int*, int*, double*, double complex*, int*, double complex*, int*, double*, double complex*, int*);


//compute eigenvector and right eigenvectors:
int eigen(int DIM, double complex *H, double complex *e, double complex *ur)
{
  int INFO, LWORK = 10*DIM;
  double complex *WORK = (double complex *)calloc((LWORK),sizeof(double complex));
  double *RWORK = (double  *)calloc((2*DIM),sizeof(double));
  char JOBVL='N', JOBVR='V';
  
  zgeev_( &JOBVL, &JOBVR, &DIM, H, &DIM, e, NULL,&DIM, ur,&DIM, WORK, &LWORK, RWORK, &INFO);
  if( !(INFO == 0) ) {
      printf( "The algorithm failed to compute eigenvalues.\n" );
      exit( 1 );
  }
  free(RWORK);
  free(WORK);
  return 0;
}



// A=>A^-1
void invert_matrix(int N, double complex *A, double complex *Ainv) {
  int INFO1=0;
  int INFO2=0;
  int nEntry=N*N;
  int * IPIV;
  double complex * WORK;
  WORK = (double complex *) malloc(nEntry*sizeof(double complex));
  IPIV = (int *) malloc(N*sizeof(int));
  
  int ii;
  for(ii=0;ii<N*N;ii++) Ainv[ii] = A[ii];
  
  //for(ii=0;ii<N*N;ii++) printf("% 4.5f  % 4.5f \n",creal(Ainv[ii]),cimag(Ainv[ii]));
  
  zgetrf_(&N, &N, Ainv, &N, IPIV, &INFO1);
  zgetri_(&N, Ainv, &N, IPIV, WORK, &nEntry, &INFO2);
  if( !(INFO1 == 0) || !(INFO2 == 0) ) {
    printf( "The algorithm failed to invert the matrix. %d %d\n", INFO1, INFO2);
    exit( 1 );
  }
  free(WORK);
  free(IPIV);
}


//C=A*B
int MatMatMul(int N, double complex * A, double complex * B, double complex * C) {
  double one=1.0;
  double zero=0.0;
  char no = 'n';
  zgemm_(&no,&no,&N,&N,&N, &one, A, &N, B, &N, &zero, C, &N); 
  return 0;
}


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


float sumEigenEquation0(int n_exc, double complex *e, double complex *u, double complex *us, double eta, double omega)
{
  int nn;
  double complex g=0.0; 
  double complex z = omega+I*eta;
  for(nn=0;nn<n_exc;nn++){
    g += u[nn] * us[nn] / (z + e[nn]);
  }    
  //printf("% 4.5f \n",cimag(g));
  return cimag(g);
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



int smallTest(int DIM, float *w, int n, float complex *H) //)H, float *S,)
{
  int ii;
  printf("\n");
  //for(ii=0;ii<DIM;ii++){
  //  printf(" % 4.5f \n",w[ii]);
  //}
  
  for(ii=0;ii<n*n;ii++){
    printf(" % 4.5f  % 4.5f \n",creal(H[ii]),cimag(H[ii]));
  }
  H[1] = 9.42341+I*132.3;
  return 0;
}

int smallTest0(int DIM, float *w, int n, float *H) //)H, float *S,)
{
  int ii;
  printf("\n");
  //for(ii=0;ii<DIM;ii++){
  //  printf(" % 4.5f \n",w[ii]);
  //}
  
  for(ii=0;ii<n*n;ii++){
    printf(" % 4.5f \n",H[ii]);
  }
  
  return 0;
}



int greenFrom_H_and_S(int sign, int Nw, double *w, int dim, double complex *H, double complex *S, double Omega, double mu, double eta, double complex *g) //)H, double *S,)
{
  int ii,nn;
  //printf("salut1\n"); fflush(stdout);
  
  double complex *Sinv   = (double complex *)calloc((dim*dim),sizeof(double complex));
  double complex *H_Sinv = (double complex *)calloc((dim*dim),sizeof(double complex));
  double complex *Uinv_S = (double complex *)calloc((dim*dim),sizeof(double complex));

  double complex *e    = (double complex *)calloc((dim),sizeof(double complex));
  double complex *U    = (double complex *)calloc((dim*dim),sizeof(double complex));
  double complex *Uinv = (double complex *)calloc((dim*dim),sizeof(double complex));
  
  //for(nn=0;nn<dim*dim;nn++){
  //  printf("% 4.5f ",creal(S[nn]));
  //}
  invert_matrix(dim, S, Sinv);
  MatMatMul(dim, H, Sinv, H_Sinv);
  eigen(dim,H_Sinv,e,U);
  invert_matrix(dim, U, Uinv);
  MatMatMul(dim, Uinv, S, Uinv_S);  
  
  for(ii=0;ii<Nw;ii++) {
    g[ii]=0.0;
    double complex z = w[ii] + I*eta + sign*Omega + mu;
    for(nn=0;nn<dim;nn++){
      g[ii] += U[dim*nn]*Uinv_S[nn] /(z-sign*e[nn]);
    }
  }
  free(Sinv);
  free(H_Sinv);
  free(Uinv_S);
  
  free(e);
  free(U);
  free(Uinv);
  return 0;
}

  //from python, this is equal to:
  //e,u = la.eig(np.dot(H,la.inv(S)))
  //uinv = la.inv(u)
  //us = np.dot(uinv,S)
  
  //  z = w + 1.0j*eta + Omega + mu
  //  for nn in range(n_exc_choice):
  //   g += u[0,nn]*us[nn,0] / (ze - e[nn])



