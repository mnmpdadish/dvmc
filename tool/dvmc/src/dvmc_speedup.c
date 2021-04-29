/*
*/
#include <complex.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <signal.h>
#include <math.h>

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
  
  float *phys_nCAm_averaged   = (float *)calloc((Nsite*Nsite*NExcitation*NExcitation),sizeof(float));
  float *phys_nACm_averaged   = (float *)calloc((Nsite*Nsite*NExcitation*NExcitation),sizeof(float));
  float *phys_nAHCm_averaged  = (float *)calloc((Nsite*Nsite*NExcitation*NExcitation),sizeof(float));
  float *phys_nCHAm_averaged  = (float *)calloc((Nsite*Nsite*NExcitation*NExcitation),sizeof(float));
  
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

  int nn,mm;
  
  for(ii=0;ii<Nsite;ii++) 
   for(jj=0;jj<Nsite;jj++) 
    for(mm=0;mm<NExcitation;mm++) {
     for(nn=0;nn<NExcitation;nn++) {
      printf("%d %d %d %d  %d %d % 4.2e % 4.2e % 4.2e % 4.2e \n",ii,0,jj,0,nn,mm,
                                    phys_nACm_averaged[nn + NExcitation*(mm + NExcitation*(ii + Nsite*jj))],
                                    phys_nCAm_averaged[nn + NExcitation*(mm + NExcitation*(ii + Nsite*jj))],
                                    phys_nAHCm_averaged[nn + NExcitation*(mm + NExcitation*(ii + Nsite*jj))],
                                    phys_nCHAm_averaged[nn + NExcitation*(mm + NExcitation*(ii + Nsite*jj))]);
      //S_CA[ii + Nsite*(mm + NExcitation*(jj + Nsite*nn))] += phys_nCAm_averaged[nn + NExcitation*(mm + NExcitation*(ii + Nsite*jj))]; 
    }
  }
  
  // untangling
  for(ii=0;ii<Nsite;ii++) {
   for(jj=0;jj<Nsite;jj++) {
    for(nn=0;nn<NExcitation;nn++) {
     for(mm=0;mm<NExcitation;mm++) {
       S_CA[ii + Nsite*(mm + NExcitation*(jj + Nsite*nn))] += phys_nCAm_averaged[nn + NExcitation*(mm + NExcitation*(ii + Nsite*jj))]; 
       S_AC[ii + Nsite*(mm + NExcitation*(jj + Nsite*nn))] += phys_nACm_averaged[nn + NExcitation*(mm + NExcitation*(ii + Nsite*jj))]; 
       H_CA[ii + Nsite*(mm + NExcitation*(jj + Nsite*nn))] += phys_nCHAm_averaged[nn + NExcitation*(mm + NExcitation*(ii + Nsite*jj))]; 
       H_AC[ii + Nsite*(mm + NExcitation*(jj + Nsite*nn))] += phys_nAHCm_averaged[nn + NExcitation*(mm + NExcitation*(ii + Nsite*jj))]; 
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

int greenFrom_e_U_Uinv_S_general(int sign, int ll, int pp, int Ns, int Nw, double *w, int dim, double complex *e, double complex *U, double complex *Uinv_S, double Omega, double mu, double eta, double complex *g) //)H, double *S,)
{
  int ii,nn,ll2,pp2,Nex;
  int ind1, ind2;

  Nex = dim;

  for(ii=0;ii<Ns*Ns*Nw;ii++) g[ii]=0.0;

  for(ii=0;ii<Nw;ii++) {
    double complex z = w[ii] + I*eta + sign*Omega + mu;

    // G_lp(w) matrix structure
    // [ G_00(w_0)    G_10(w_0)    ... G_01(w_0)    G_11(w_0)    ... G_NsNs(w_0)    ]
    // [ G_00(w_1)    G_10(w_1)    ... G_01(w_1)    G_11(w_1)    ... G_NsNs(w_1)    ]
    // [   .                                                              .         ]
    // [   .                                                              .         ]
    // [   .                                                              .         ]
    // [ G_00(w_Nw-1) G_10(w_Nw-1) ... G_01(w_Nw-1) G_11(w_Nw-1) ... G_NsNs(w_Nw-1) ]
    // [ G_00(w_Nw)   G_10(w_Nw)   ... G_01(w_Nw)   G_11(w_Nw)   ... G_NsNs(w_Nw)   ]
    //for(ll=0;ll<Ns;ll++){
    //for(pp=0;pp<Ns;pp++){
    //for(ll2=0;ll2<Ns;ll2++){
    for(pp2=0;pp2<Ns;pp2++){
      for(nn=0;nn<Nex;nn++){
	  //g[ll+Ns*pp+Ns*Ns*ii] += U[nn+Ns*Nex*ll]*Uinv_S[Ns*Nex*pp+Ns*Nex*nn]/(z-sign*e[nn]);
	  //g[ll+Ns*pp+Ns*Ns*ii] += U[nn+Ns*Nex*ll]*Uinv_S[pp+Ns*Nex*nn]/(z-sign*e[nn]);
	      //ind1 = ll + Ns*(0+Nex*(pp2+Ns*mm));
	ind1 = pp2 + Ns*(nn+Nex*(ll+Ns*0));
	ind2 = pp + Ns*(0+Nex*(pp2+Ns*nn));
	      //g[ll+Ns*pp+Ns*Ns*ii] += U[nn+Ns*Nex*ll]*Uinv_S[pp+Ns*Nex*nn]/(z-sign*e[nn]);
	g[ii] += U[ind1]*Uinv_S[ind2]/(z-sign*e[nn]);
	  //if(ii==0 && ll==0 && pp==0) printf("%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",nn,creal(g[ll+Ns*pp+Ns*Ns*ii]),cimag(g[ll+Ns*pp+Ns*Ns*ii]),creal(U[nn+Ns*Nex*ll]), cimag(U[nn+Ns*Nex*ll]), creal(Uinv_S[pp+Ns*Nex*nn]), cimag(Uinv_S[pp+Ns*Nex*nn]),creal(z),cimag(z),creal(e[nn]),cimag(e[nn]));
	  //printf("%d\t%d\t%d\t%f\t%f\t%f\t%f\n",ll,pp,nn,creal(U[nn+Ns*Nex*ll]), cimag(U[nn+Ns*Nex*ll]), creal(Uinv_S[pp+Ns*Nex*nn]), cimag(Uinv_S[pp+Ns*Nex*nn]) );
	  //printf("%d\t%d\t%d\t%f\t%f\n",ll,pp,nn,creal(U[nn+Ns*Nex*ll]), cimag(U[nn+Ns*Nex*ll]));
	}
	//if(ll<5 && pp<5 && ii<5) printf("%d\t%d\t%d\t%f\t%f\n",ll,pp,ii,creal(g[ll+Ns*pp+Ns*Ns*ii]),cimag(g[ll+Ns*pp+Ns*Ns*ii]));
      }
    }
    //}

  return 0;
}

int fourier_Green(int Lx, int Ly, double complex *gr, double complex *gk)
{

  int ii, jj, mm, nn, ll, pp;
  int Ns = Lx*Ly;
  int N[2] ={Lx,Ly} ;

  double kx[Ns], ky[Ns], x[Ns], y[Ns];
  double pi = acos(-1.0);

  int ntemp, den;

  // get coordinates for lattice sites
  for(ii=0;ii<Ns;ii++){
    ntemp=ii;
    for(jj=1;jj>=0;jj--){
      den=1;
      for(ll=0;ll<jj;ll++){
	den=den*N[ll];
      }
      if(jj==0) {
	x[ii]=ntemp/den;
	ntemp=ntemp-x[ii]*den;
	x[ii]=x[ii];//+1;
      }
      if(jj==1) {
	y[ii]=ntemp/den;
	ntemp=ntemp-y[ii]*den;
	y[ii]=y[ii];//+1;
      }
    }
    //printf("%d\t%f\t%f\n",ii,x[ii],y[ii]);
  }

  for(ii=0;ii<Ns;ii++){
    kx[ii] = (2*pi*x[ii])/Lx;
    ky[ii] = (2*pi*y[ii])/Ly;
    //printf("%d\t%f\t%f\t%f\t%f\n",ii,(2*pi*x[ii]),(2*pi*y[ii]),kx[ii],ky[ii]);
  }

  for(ii=0;ii<Ns*Ns;ii++) gk[ii] = 0.0;

  /*
  for(ii=00;ii<Ns;ii++){
    for(jj=00;jj<Ns;jj++){
      printf("%d\t%d\t%f\t%f\n",ii,jj,creal(gr[Ns*ii+jj]),cimag(gr[Ns*ii+jj]));
    }
  }
  */

  // fourier transform G_{k,kp} = 1/Ns * sum_{ij} e^ik.r_i*e^ikp.r_j * G_{ij}
  for(mm=0;mm<Ns;mm++){
    for(nn=0;nn<Ns;nn++){
      for(ii=0;ii<Ns;ii++){
	for(jj=0;jj<Ns;jj++){
	  double complex phase = cexp(I*(kx[mm]*x[ii]+ky[mm]*y[ii]))*cexp(-1.0*I*(kx[nn]*x[jj]+ky[nn]*y[jj]));
	  //double complex phase = I*(kx[mm]*x[ii]+ky[mm]*y[ii]);
	  //if(mm==nn) printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",kx[mm],ky[mm],x[ii],y[ii],creal(phase),cimag(phase),creal(gr[Ns*ii+jj]),cimag(gr[Ns*ii+jj]));
	  gk[Ns*mm+nn] += (1.0/Ns)*phase*gr[Ns*ii+jj];
	}
      }
      //if(mm==nn) 
      //printf("%d\t%d\t%d\t%f\t%f\n",mm,nn,Ns*mm+nn,creal(gk[Ns*mm+nn]),cimag(gk[Ns*mm+nn]));
    }
  }

  return 0;

}

