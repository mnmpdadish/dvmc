
#include <complex.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>


extern void cgetrf_(int*, int* , float complex*, int*, int*, int*);
extern void cgetri_(int*, float complex*, int*, int*, float complex*, int*, int*);
extern void cgeev_(char*, char*, int *, float complex*, int*, float complex*, float complex*, int*, float complex*,int*, float complex*, int*, float*, int*);
extern void cgemm_(char*, char*, int*, int*, int*, float*, float complex*, int*, float complex*, int*, float*, float complex*, int*);


//compute eigenvector and right eigenvectors:
int eigen(int DIM, float complex *H, float complex *e, float complex *ur)
{
  int INFO, LWORK = 10*DIM;
  float complex *WORK = (float complex *)calloc((LWORK),sizeof(float complex));
  float *RWORK = (float  *)calloc((2*DIM),sizeof(float));
  char JOBVL='N', JOBVR='V';
  
  cgeev_( &JOBVL, &JOBVR, &DIM, H, &DIM, e, NULL,&DIM, ur,&DIM, WORK, &LWORK, RWORK, &INFO);
  if( !(INFO == 0) ) {
      printf( "The algorithm failed to compute eigenvalues.\n" );
      exit( 1 );
  }
  free(RWORK);
  free(WORK);
  return 0;
}



// A=>A^-1
void invert_matrix(int N, float complex *A) {
  int INFO1=0;
  int INFO2=0;
  int nEntry=N*N;
  int * IPIV;
  float complex * WORK;
  WORK = (float complex *) malloc(nEntry*sizeof(float));
  IPIV = (int *) malloc(N*sizeof(int));
  
  cgetrf_(&N, &N, A, &N, IPIV, &INFO1);
  cgetri_(&N, A, &N, IPIV, WORK, &nEntry, &INFO2);
  if( !(INFO1 == 0) || !(INFO2 == 0) ) {
    printf( "The algorithm failed to invert the fMatrix. %d %d\n", INFO1, INFO2);
    exit( 1 );
  }
  free(WORK);
  free(IPIV);
}

//C=A*B
int MatMatMul(int N, float complex * A, float complex * B, float complex * C) {
  float one=1.0;
  float zero=0.0;
  char no = 'n';
  cgemm_(&no,&no,&N,&N,&N, &one, A, &N, B, &N, &zero, C, &N); 
  return 0;
}

void main()
{
  
  int ii,jj,dim=3;
  float complex *H = (float complex *)calloc(dim*dim,sizeof(float complex));
  float complex *H2 = (float complex *)calloc(dim*dim,sizeof(float complex));
  float complex *H3 = (float complex *)calloc(dim*dim,sizeof(float complex));
  
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

  H[0] = 1.0; H[0+dim] = 2.0; H[0+dim*2] = 6.0+I*5.0;
  H[1] = 3.0; H[1+dim] = 4.0; H[1+dim*2] = 0.0;
  H[2] = 1.0; H[2+dim] = 2.0; H[2+dim*2] = 1.0;

  printf("\nH\n");
  for(ii=0;ii<dim;ii++){
    for(jj=0;jj<dim;jj++){
      printf(" % 7.3f +i% 7.3f  ",creal(H[ii+jj*dim]),cimag(H[ii+jj*dim]));
    }
    printf("\n");
  }

  H2[0] = 4.0; H2[0+dim] = 1.0; H2[0+dim*2] = 3.0;
  H2[1] = 4.0; H2[1+dim] = 3.0; H2[1+dim*2] = 1.0;
  H2[2] = 4.0; H2[2+dim] = 1.0; H2[2+dim*2] = 2.0+I*1.0;
  
  printf("\nmultiplication\n");
  MatMatMul(dim,H,H2,H3);
  
  for(ii=0;ii<dim;ii++){
    for(jj=0;jj<dim;jj++){
      printf(" % 7.3f +i% 7.3f  ",creal(H3[ii+jj*dim]),cimag(H3[ii+jj*dim]));
    }
    printf("\n");
  }

  
  printf("\ninv:\n");
  H[0] = 1.0; H[0+dim] = 2.0; H[0+dim*2] = 6.0+I*5.0;
  H[1] = 3.0; H[1+dim] = 4.0; H[1+dim*2] = 0.0;
  H[2] = 1.0; H[2+dim] = 2.0; H[2+dim*2] = 1.0;
  
  invert_matrix(dim, H);
  for(ii=0;ii<dim;ii++){
    for(jj=0;jj<dim;jj++){
      printf(" % 7.3f +i% 7.3f  ",creal(H[ii+jj*dim]),cimag(H[ii+jj*dim]));
    }
    printf("\n");
  }
  

  printf("\neigen:\n");
  H[0] = 1.0; H[0+dim] = 2.0; H[0+dim*2] = 6.0+I*5.0;
  H[1] = 3.0; H[1+dim] = 4.0; H[1+dim*2] = 0.0;
  H[2] = 1.0; H[2+dim] = 2.0; H[2+dim*2] = 1.0;
  
  eigen(dim,  H,  e,  U);
  //eigen(dim, H2, e2, U2);
  printf("\nre ");
  for(ii=0;ii<dim;ii++){
    printf(" % 7.3f +i% 7.3f  ",creal(e[ii]),cimag(e[ii]));
  }

  printf("\n\n");
  for(ii=0;ii<dim;ii++){
    for(jj=0;jj<dim;jj++){
      printf(" % 7.3f +i% 7.3f  ",creal(U[ii+jj*dim]),cimag(U[ii+jj*dim]));
    }
    printf("\n");
  }

  printf("\n% 4.5f\n",creal(U[0+dim]));
  printf("\ncompare these values with mat.py in the same directory.\nthis is a very basic benchmark.\n"); 
  return;
}
