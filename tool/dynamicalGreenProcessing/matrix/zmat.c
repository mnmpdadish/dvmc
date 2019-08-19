
#include <complex.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>



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




void main()
{
  
  int ii,jj,dim=3;
  double complex *H = (double complex *)calloc(dim*dim,sizeof(double complex));
  double complex *H2 = (double complex *)calloc(dim*dim,sizeof(double complex));
  double complex *H3 = (double complex *)calloc(dim*dim,sizeof(double complex));
  
  double complex *U = (double complex *)calloc(dim*dim,sizeof(double complex));
  double complex *U2 = (double complex *)calloc(dim*dim,sizeof(double complex));

  double complex *e = (double complex *)calloc(dim,sizeof(double complex));
  double complex *e2 = (double complex *)calloc(dim,sizeof(double complex));
  
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
  
  invert_matrix(dim, H,H2);
  for(ii=0;ii<dim;ii++){
    for(jj=0;jj<dim;jj++){
      printf(" % 7.3f +i% 7.3f  ",creal(H2[ii+jj*dim]),cimag(H2[ii+jj*dim]));
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
