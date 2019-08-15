/*
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>



int ReadExcitationDef(int *iNExcitation, int *iExcitation_L, int *iExcitation_W) {

  FILE *fp;
  char *cerr;
  char ctmp[1024];
  char ctmp2[1024];

  fprintf(stdout, "  Read File 'excitation.def'.\n");
  fp = fopen("excitation.def", "r");
  if (fp == NULL) {
    fprintf(stdout, "error: no 'excitation.def' found.\n"); 
    fclose(fp);
    exit(-1);
  }
  cerr = fgets(ctmp, sizeof(ctmp) / sizeof(char), fp);
  if (cerr != NULL) {
    cerr = fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp);
    sscanf(ctmp2, "%s %d\n", ctmp, iNExcitation); //2
    cerr = fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp);
    sscanf(ctmp2, "%s %d\n", ctmp, iExcitation_L); //3
    cerr = fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp);
    sscanf(ctmp2, "%s %d\n", ctmp, iExcitation_W); //4
  }
  fclose(fp);
  return 0;
}





int main(int argc, char *argv[])
{
  
  int ii;
  if(argc>1) 
  {
      for(ii = 0; ii<argc-1; ii++){
        fprintf(stdout,"%s\n",argv[1+ii]);
      }
    
  }
  else{
    fprintf(stdout, "error: need input files\n");
  }
  
  int n_file = argc-1;
  double factor = 1.0/n_file;
  int NExcitation, Excitation_L, Excitation_W;
  ReadExcitationDef(&NExcitation, &Excitation_L, &Excitation_W);
  fprintf(stdout, "NExcitation, Excitation_L, Excitation_W\n%11d %13d %13d\n", NExcitation, Excitation_L, Excitation_W);
  int Nsite = Excitation_L * Excitation_W;
  int dr,nn,mm;
  double *Phys_averaged;  
  Phys_averaged  = (double *)calloc((Nsite*NExcitation*NExcitation*4),sizeof(double));
  int NExcitation2 = NExcitation*NExcitation;
  
  FILE *fp;
  char *cerr;
  char ctmp[1024];
  char ctmp2[1024];
  //char fileName[] = "output/zvo_nCHAm_nAHCm_001.dat";
  int dummy1,dummy2,dummy3;
  
  for(ii = 0; ii<argc-1; ii++){
    fprintf(stdout, "  Read File '%s'.\n",argv[1+ii]);
    fp = fopen(argv[1+ii], "r");
    if (fp == NULL) {
      fprintf(stdout, "error: no '%s' found.\n",argv[1+ii]); 
      fclose(fp);
      exit(-1);
    }
    cerr = fgets(ctmp, sizeof(ctmp) / sizeof(char), fp); //skip header
    if (cerr != NULL) {
      for (dr = 0; dr < Nsite; dr++) {
       printf("\n%d  ", dr);
       for (nn = 0; nn < NExcitation; nn++) {
        for (mm = 0; mm < NExcitation; mm++) {
         cerr = fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp);
         /*if (cerr != NULL) {
           fprintf(stdout, "error in file '%s' found. \n",fileName); 
           fclose(fp);
           exit(-1);
         }*/
         double val1,val2,val3,val0;
         sscanf(ctmp2, "%d %d %d %le %le %le %le\n", &dummy1,&dummy2,&dummy3, &val0,&val1,&val2,&val3);
         Phys_averaged[4*(nn+NExcitation*mm + dr*NExcitation2)+0] += factor*val0; 
         Phys_averaged[4*(nn+NExcitation*mm + dr*NExcitation2)+1] += factor*val1; 
         Phys_averaged[4*(nn+NExcitation*mm + dr*NExcitation2)+2] += factor*val2; 
         Phys_averaged[4*(nn+NExcitation*mm + dr*NExcitation2)+3] += factor*val3; 
         //fprintf(stdout,"%d %d %d %f %f %f\n", dummy1,dummy2,dummy3, val1,val2,val3,val4);
        }
       }
      }
    }
    fclose(fp);
  }
  return 0;
}


