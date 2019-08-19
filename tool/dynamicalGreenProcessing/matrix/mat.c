
#include <complex.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>


extern void sgetrf_(int*, int* , float*, int*, int*, int*);
extern void sgetri_(int*, float*, int*, int*, float*, int*, int*);
extern void sgeev_(char*, char*, int *, float*, int*, float*, float*, float*,int*, float*,int*, float*, int*, int*);
extern void sgemm_(char*, char*, int*, int*, int*, float*, float*, int*, float*, int*, float*, float*, int*);


//compute eigenvector and right eigenvectors:
int eigen(int DIM, float *H, float complex *e, float complex *u)
{
  int INFO, LWORK = 10*DIM;
  float *WORK = (float *)calloc((LWORK),sizeof(float));
  char JOBVL='N', JOBVR='V';
  
  float *e_re = (float *)calloc((DIM),sizeof(float));
  float *e_im = (float *)calloc((DIM),sizeof(float));
  
  float *vl = (float *)calloc((DIM*DIM),sizeof(float));
  float *vr = (float *)calloc((DIM*DIM),sizeof(float));
  
  sgeev_( &JOBVL, &JOBVR, &DIM, H, &DIM, e_re, e_im, vl,&DIM, vr,&DIM, WORK, &LWORK, &INFO);
  if( !(INFO == 0) ) {
      printf( "The algorithm failed to compute eigenvalues.\n" );
      exit( 1 );
  }
  
  
  int ii,jj;
  for(ii=0;ii<DIM;ii++){
    e[ii] = e_re[ii] + I * e_im[ii];
    if(e_im[ii]==0.0){
      //printf("%d % 10.5f  % 10.5f\n",ii,e_re[ii],e_im[ii]);
      for(jj=0;jj<DIM;jj++) u[jj+DIM*ii] = vr[jj+DIM*ii];
    }
    else{
      //printf("%d % 12.5f  % 12.5f\n",ii,e_re[ii],e_im[ii]);
      for(jj=0;jj<DIM;jj++) {
        u[jj+DIM*ii]     = vr[jj+DIM*ii] + I*vr[jj+DIM*(ii+1)];
        u[jj+DIM*(ii+1)] = vr[jj+DIM*ii] - I*vr[jj+DIM*(ii+1)];
      }
      ii++;
      //printf("%d % 12.5f  % 12.5f\n",ii,e_re[ii],e_im[ii]);
      e[ii] = e_re[ii] + I*e_im[ii];
    }
  }
  
  free(e_re);
  free(e_im);
  free(vl);//useless
  free(vr); 
  free(WORK);
  return 0;
}



// A=>A^-1
void invert_matrix(int N, float *A) {
  int INFO1=0;
  int INFO2=0;
  int nEntry=N*N;
  int * IPIV;
  float * WORK;
  WORK = (float *) malloc(nEntry*sizeof(float));
  IPIV = (int *) malloc(N*sizeof(int));
  
  sgetrf_(&N, &N, A, &N, IPIV, &INFO1);
  sgetri_(&N, A, &N, IPIV, WORK, &nEntry, &INFO2);
  if( !(INFO1 == 0) || !(INFO2 == 0) ) {
    printf( "The algorithm failed to invert the fMatrix. %d %d\n", INFO1, INFO2);
    exit( 1 );
  }
  free(WORK);
  free(IPIV);
}

//C=A*B
int MatMatMul(int N, float * A, float * B, float * C) {
  float one=1.0;
  float zero=0.0;
  char no = 'n';
  sgemm_(&no,&no,&N,&N,&N, &one, A, &N, B, &N, &zero, C, &N); 
  return 0;
}

void main()
{
  
  int ii,jj,dim=3;
  float *H = (float *)calloc(dim*dim,sizeof(float));
  float *H2 = (float *)calloc(dim*dim,sizeof(float));
  float *H3 = (float *)calloc(dim*dim,sizeof(float));
  
  float complex *U = (float complex *)calloc(dim*dim,sizeof(float complex));
  float complex *U2 = (float complex *)calloc(dim*dim,sizeof(float complex));

  float complex *e = (float complex *)calloc(dim,sizeof(float complex));
  float complex *e2 = (float complex *)calloc(dim,sizeof(float complex));
  
  
  /*  
  for(ii=0;ii<dim;ii++){
    for(jj=0;jj<dim;jj++){
      H[ii+jj*dim] = (ii+jj*dim);
      H2[ii+jj*dim] = (ii-jj*dim)*(ii+jj*dim);
    }
  }*/

  H[0] = 1.0; H[0+dim] = 2.0; H[0+dim*2] = 6.0;
  H[1] = 3.0; H[1+dim] = 4.0; H[1+dim*2] = 0.0;
  H[2] = 1.0; H[2+dim] = 2.0; H[2+dim*2] = 1.0;

  printf("\nH\n");
  for(ii=0;ii<dim;ii++){
    for(jj=0;jj<dim;jj++){
      printf(" % 4.5f",H[ii+jj*dim]);
    }
    printf("\n");
  }

  H2[0] = 4.0; H2[0+dim] = 1.0; H2[0+dim*2] = 3.0;
  H2[1] = 4.0; H2[1+dim] = 3.0; H2[1+dim*2] = 1.0;
  H2[2] = 4.0; H2[2+dim] = 1.0; H2[2+dim*2] = 2.0;
  
  printf("\nmultiplication\n");
  MatMatMul(dim,H,H2,H3);
  
  for(ii=0;ii<dim;ii++){
    for(jj=0;jj<dim;jj++){
      printf(" % 4.5f",H3[ii+jj*dim]);
    }
    printf("\n");
  }

  
  printf("\ninv:\n");
  H[0] = 1.0; H[0+dim] = 2.0; H[0+dim*2] = 6.0;
  H[1] = 3.0; H[1+dim] = 4.0; H[1+dim*2] = 0.0;
  H[2] = 1.0; H[2+dim] = 2.0; H[2+dim*2] = 1.0;
  
  invert_matrix(dim, H);
  for(ii=0;ii<dim;ii++){
    for(jj=0;jj<dim;jj++){
      printf(" % 4.5f",H[ii+jj*dim]);
    }
    printf("\n");
  }
  

  printf("\neigen:\n");
  H[0] = 1.0; H[0+dim] = 2.0; H[0+dim*2] = 6.0;
  H[1] = 3.0; H[1+dim] = 4.0; H[1+dim*2] = 0.0;
  H[2] = 1.0; H[2+dim] = 2.0; H[2+dim*2] = 1.0;
  
  eigen(dim,  H,  e,  U);
  //eigen(dim, H2, e2, U2);
  printf("\nre ");
  for(ii=0;ii<dim;ii++){
    printf(" % 4.5f",creal(e[ii]));
  }
  printf("\nim ");
  for(ii=0;ii<dim;ii++){
    printf(" % 4.5f",cimag(e[ii]));
  }

  printf("\n\nre");
  for(ii=0;ii<dim;ii++){
    for(jj=0;jj<dim;jj++){
      printf(" % 4.5f",creal(U[ii+jj*dim]));
    }
    printf("\n  ");
  }

  printf("\nim");
  for(ii=0;ii<dim;ii++){
    for(jj=0;jj<dim;jj++){
      printf(" % 4.5f",cimag(U[ii+jj*dim]));
    }
    printf("\n  ");
  }

  printf("\n% 4.5f\n",creal(U[0+dim]));
  printf("\ncompare these values with mat.py in the same directory.\nthis is a very basic benchmark.\n"); 
  return;
}
