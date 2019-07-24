/*
mVMC - A numerical solver package for a wide range of quantum lattice models based on many-variable Variational Monte Carlo method
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is developed based on the mVMC-mini program
(https://github.com/fiber-miniapp/mVMC-mini)
which follows "The BSD 3-Clause License".

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details. 

You should have received a copy of the GNU General Public License 
along with this program. If not, see http://www.gnu.org/licenses/. 
*/
/*-------------------------------------------------------------
 * Variational Monte Carlo
 * Cauculate Green Functions
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/
#include "calgrn.h"
#ifndef _CALGRN_SRC
#define _CALGRN_SRC

int read_StdFace_L_W(){
  FILE * tmpfile;
  tmpfile = fopen("StdFace.def",  "r");
  //rewind(file);
    
  char tmpbuff[200];
  while(!feof(tmpfile)) 
  {
    if (fgets(tmpbuff,200,tmpfile))
    {
      if(sscanf (tmpbuff,"L = %d",&StdFace_L)) { printf("found L=%d\n", StdFace_L);} 
      if(sscanf (tmpbuff,"W = %d",&StdFace_W)) { printf("found W=%d\n", StdFace_W);} 
    }
  }
  assert(StdFace_L*StdFace_W==Nsite);
  
  fclose(tmpfile);
  return 0;
}




void CalculateGreenFuncMoments(const double w, const double complex ip, 
                               int *eleIdx, int *eleCfg,
                               int *eleNum, int *eleProjCnt) {

  int idx,idx0,idx1;
  int ri,rj,s,rk,rl,t;
  double complex tmp;
  int *myEleIdx, *myEleNum, *myProjCntNew;
  double complex *myBuffer;
  
  RequestWorkSpaceThreadInt(Nsize+Nsite2+NProj);
  //RequestWorkSpaceThreadComplex(NQPFull+2*Nsize);
  RequestWorkSpaceThreadComplex(NQPFull + 2*Nsize);
  // GreenFunc1: NQPFull, GreenFunc2: NQPFull+2*Nsize 

  //#pragma omp parallel default(shared)                \
  private(myEleIdx,myEleNum,myProjCntNew,myBuffer,idx)
  {
    myEleIdx = GetWorkSpaceThreadInt(Nsize);
    myEleNum = GetWorkSpaceThreadInt(Nsite2);
    myProjCntNew = GetWorkSpaceThreadInt(NProj);
    myBuffer = GetWorkSpaceThreadComplex(NQPFull+2*Nsize);

    //#pragma loop noalias
    for(idx=0;idx<Nsize;idx++) myEleIdx[idx] = eleIdx[idx];
    //#pragma loop noalias
    for(idx=0;idx<Nsite2;idx++) myEleNum[idx] = eleNum[idx];

    //#pragma omp master
    {
      //printf("start:\n");
      StartTimer(50);
    }

    //#pragma omp for private(idx,ri,rj,s,tmp) schedule(dynamic) nowait
    for(idx=0;idx<NCisAjs;idx++) {
      ri = CisAjsIdx[idx][0];
      rj = CisAjsIdx[idx][2];
      s  = CisAjsIdx[idx][3];


      
      //Doublon-Holon TJS
      PhysN2[idx+NCisAjs*0] += w*myEleNum[ri+s*Nsite]*myEleNum[ri+(1-s)*Nsite]
                                   *(1.0-myEleNum[rj+s*Nsite])*(1.0-myEleNum[rj+(1-s)*Nsite]);
                                   
      //Doublon-Doublon TJS
      PhysN2[idx+NCisAjs*1] += w*myEleNum[ri+s*Nsite]*myEleNum[ri+(1-s)*Nsite]
                                      *myEleNum[rj+s*Nsite]*myEleNum[rj+(1-s)*Nsite];

      //Charge-Doublon TJS
      PhysN2[idx+NCisAjs*2] += w*myEleNum[ri+s*Nsite] *myEleNum[rj+s*Nsite]*myEleNum[rj+(1-s)*Nsite];
      
      //n_sigma (1-n_sigma) MC
      PhysN2[idx+NCisAjs*3] += w*myEleNum[ri+s*Nsite] *(1.0-myEleNum[rj+s*Nsite]);
    }

    //#pragma omp for private(ri,s) schedule(dynamic) nowait
    for(ri=0;ri<Nsite;ri++) {
      s=0;
      //Doublon MC
      PhysN1[ri+Nsite*0] += w*myEleNum[ri]*myEleNum[ri+Nsite];

      //Holon MC
      PhysN1[ri+Nsite*1] += w*(1.0-myEleNum[ri])*(1.0-myEleNum[ri+Nsite]);
            
      //Density_up (s==0) MC
      PhysN1[ri+Nsite*2] += w*myEleNum[ri];

      //Density_down (s==1) MC
      PhysN1[ri+Nsite*3] += w*myEleNum[ri+Nsite];           
    }
    
    //#pragma omp master
    {StopTimer(50);}
  }
  ReleaseWorkSpaceThreadInt();
  ReleaseWorkSpaceThreadComplex();
  return;
}




int del(int i,int j){
  return ((i==j)? 1:0);
}

int Commute_Nat_(commuting_with commuting, int ra, int rb, int t, int ri, int rj, int s, int rm, int rn, int u, int *eleNum) {
  int sign;
  if((commuting==with_CisCmuAnuAjs) || (commuting==with_AisCmuAnuCjs)) {
    if(commuting==with_CisCmuAnuAjs) {sign = 1;}
    else if(commuting==with_AisCmuAnuCjs) {sign = -1;}
    
    if(t==0){ 
      return 1;
    }
    else if(t==1){ 
      return eleNum[ra+(1-s)*Nsite]  + del((1-s),u) * (del(ra,rm)-del(ra,rn));
    }
    else if(t==2){
      if(ra==rb){
        return (eleNum[ra +s *Nsite] + (del(ra,rm)-del(ra,rn)) * del(u,s)    + sign * (del(ra,ri)-del(ra,rj)));      
      }
      return (eleNum[ra+   s *Nsite] + (del(ra,rm)-del(ra,rn)) * del(u,s)    + sign * (del(ra,ri)-del(ra,rj)) )
            *(eleNum[rb+   s *Nsite] + (del(rb,rm)-del(rb,rn)) * del(u,s)    + sign * (del(rb,ri)-del(rb,rj)) );
    }
    else if(t==3){
      return (eleNum[ra+(1-s)*Nsite] + (del(ra,rm)-del(ra,rn)) * del(u,1-s))
            *(eleNum[rb+(1-s)*Nsite] + (del(rb,rm)-del(rb,rn)) * del(u,1-s));
    }
    else if(t==4){
      return (eleNum[ra+(1-s)*Nsite] + (del(ra,rm)-del(ra,rn)) * del(u,1-s))
            *(eleNum[rb+   s *Nsite] + (del(rb,rm)-del(rb,rn)) * del(u,s)    + sign * (del(rb,ri)-del(rb,rj)) )   ;
    }
    else{
      printf("oups, error %d\n", t);
      exit(1);
      return -1;
    }
  
  }
  else if((commuting==with_CisAjs) || (commuting==with_AisCjs)) {
  
    if(commuting==with_CisAjs) {sign = 1;}
    else if(commuting==with_AisCjs) {sign = -1;}
    
    if(t==0){ //no charge
      return 1;
    }
    else if(t==1){ // electron reverse-spin
      return eleNum[ra+(1-s)*Nsite] ;
    }
    else if(t==2){ 
      if(ra==rb) {
        return (eleNum[ra+    s*Nsite] + sign * (del(ra,ri)-del(ra,rj) ) );
      }
        return (eleNum[ra+    s*Nsite] + sign * (del(ra,ri)-del(ra,rj) ) )
              *(eleNum[rb+    s*Nsite] + sign * (del(rb,ri)-del(rb,rj) ) );      
    }
    else if(t==3){ 
      return (eleNum[ra+(1-s)*Nsite] )
            *(eleNum[rb+(1-s)*Nsite] );
    }
    else if(t==4){ 
      return (eleNum[ra+(1-s)*Nsite] )
            *(eleNum[rb+    s*Nsite] + sign * (del(rb,ri)-del(rb,rj) ) );
    }
    else{
      printf("oups, error %d\n", t);
      exit(1);
      return -1;
    }
  }
  else if(commuting==with_nothing) {
    if(t==0){ //no charge
      return 1;
    }    
    else if(t==1){ //1 electron reverse-spin
      return eleNum[ra+(1-s)*Nsite];
    }
    else if(t==2){ 
      return (eleNum[ra+s*Nsite]*eleNum[rb+s*Nsite]);
    }
    else if(t==3){ 
      return (eleNum[ra+(1-s)*Nsite]*eleNum[rb+(1-s)*Nsite]);
    }
    else if(t==4){ 
      return (eleNum[ra+(1-s)*Nsite]*eleNum[rb+s*Nsite]);
    }
    else{
      printf("oups, error %d\n", t);
      exit(1);
      return -1;
    }
  }
  
}




// retunr i%N, in the same way that python treat the negative numbers
int moduloPython(int i,int N){
  return ((i % N) + N) % N;
}

int find_neighbor_site(int r,int dx,int dy){
  assert(Dimension_1*Dimension_2==Nsite);
  int r_x   = moduloPython(r + dx            , Dimension_1);  
  int r_y   = moduloPython(r + Dimension_1*dy, Nsite);  
  int r_out = moduloPython(r_x + Dimension_1*r_y,Nsite);    
  return r_out;
}

//C=C+weight*A*B
unsigned int C_ADD_AxB(double * C, double const * A, double const * B, int N, double weight, int sampleSize) {
  char transA= 'N', transB= 'C';
  double beta=1.0;
  //int ONE = 1;
  //M_DGEMM(&transA,&transB,&N,&N,&sampleSize, &weight, &A[0], &N, &B[0], &N, &beta, &C[0], &N); 
  //*
  int ii, jj, kk;
  for(ii=0; ii<N; ii++){
    for(jj=0; jj<N; jj++){
      for(kk=0; kk<sampleSize; kk++){
        C[ii+jj*N] += weight * A[ii+kk*N] * B[jj+kk*N];
      }
    }
  }//*/
  return 0;
}


void CalculateGreenFuncMoments2_real(const double w, const double ip, 
                                      int *eleIdx, int *eleCfg,
                                      int *eleNum, int *eleProjCnt, int sample) {
  int idx,idx0,idx1;
  int ri,rj,s,rk,rl,t;
  double tmp;
  int *myEleIdx, *myEleNum, *myEleCfg, *myProjCntNew;
  double *myBuffer_real;
  
  double AC_tmp, CA_tmp;
  
  //#pragma omp parallel default(shared)                \
  //private(myEleIdx,myEleNum,myProjCntNew,myBuffer,idx, O_AC_vec, O_CA_vec, O0_vec, H_vec)//, O_vec)
  {
    myEleIdx = (int*)malloc(sizeof(int) * Nsize );
    myEleNum = (int*)malloc(sizeof(int) * Nsite2 );
    //myEleCfg = (int*)malloc(sizeof(int) * Nsite2 );
    myProjCntNew = (int*)malloc(sizeof(int) * NProj );
    myBuffer_real = (double *)malloc(sizeof(double) * (NQPFull+2*Nsize) );

    for(idx=0;idx<Nsize;idx++)  myEleIdx[idx] = eleIdx[idx];
    for(idx=0;idx<Nsite2;idx++) myEleNum[idx] = eleNum[idx];
    //for(idx=0;idx<Nsite2;idx++) myEleCfg[idx] = eleCfg[idx];
    
    double factor;
    int rm, rn, u;
    int idx_int, idx_trans;

    double multiplicity = (double) (Nsite*Nsite) / (double) NDynamicalGIdx;
    double f0 = 1.0/multiplicity;
    //printf("multiplicity: %f,  %f\n", multiplicity, f0);
    //exit(0);
    
        

    int ri, rj;
//#pragma omp for private(idx,ri,rj,s,tmp) schedule(dynamic) nowait
    for(ri=0;ri<Nsite;ri++) {
     for(rj=0;rj<Nsite;rj++) {
      for(s=0;s<2;s++) {
       Local_CA[ri+Nsite*rj+Nsite*Nsite*s] = GreenFunc1_real(ri,rj,s,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,
                                                              myProjCntNew,myBuffer_real); 

       //Local_AC[ri+Nsite*rj] = GreenFunc1_real_AisCjs(ri,rj,spin,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,
       //                                                   myProjCntNew,myBuffer_real); 
      }
      //printf("%d %d   %d  / %d \n", ri,rj,DynamicalGIdx[ri][rj],NDynamicalGIdx); fflush(stdout);
     }
    }
//#pragma omp for private(idx) nowait
    for(ri=0;ri<Nsite;ri++) {
     for(rj=0;rj<Nsite;rj++) {
      for(s=0;s<2;s++) {
       Phys_CA[DynamicalGIdx[ri][rj] + s*NDynamicalGIdx] += w*f0*Local_CA[ri+Nsite*rj+Nsite*Nsite*s];
      }
     }
    }

//#pragma omp for private(idx,idx_trans,rk,rl,t,rm,rn,u) schedule(dynamic) nowait
    for(ri=0;ri<Nsite;ri++) {
     for(rj=0;rj<Nsite;rj++) {
      s=0;

      for(idx_trans=0;idx_trans<NTransfer;idx_trans++) {
      
        rm = Transfer[idx_trans][0];
        rn = Transfer[idx_trans][2];
        u  = Transfer[idx_trans][3];
        
        Local_CisAjsCmuAnu[idx_trans][ri+Nsite*rj] = GreenFunc2_real(ri,rj,rm,rn,s,u,ip,
                        myEleIdx,eleCfg,myEleNum,eleProjCnt,myProjCntNew,myBuffer_real);

        //Local_AisCjsCmuAnu[idx_trans][idx] = GreenFunc2_real_AisCjsCktAlt(ri,rj,rm,rn,s,u,ip,
        //                      myEleIdx,eleCfg,myEleNum,eleProjCnt,myProjCntNew,myBuffer_real);
      }
     }
    }
    
///*
//////////////////////////////////////////////////////////////////////////////////////////////////////
//#pragma omp for private(idx,ri,rj,s,rk,rl,t,tmp) schedule(dynamic) nowait    
    for(ri=0;ri<Nsite;ri++) {
     for(rj=0;rj<Nsite;rj++) {
      s=0;
      //int idx_green1 = OrbitalIdx[ri][rj];
      //int idx_green2 = OrbitalIdx[rj][ri];
      int idx1 = ri+Nsite*rj;
      int idx2 = rj+Nsite*ri;
      
      int idx_exc;
      for(idx_exc=0;idx_exc<NExcitation;idx_exc++){
        //int idx_reverse = ijst_to_idx[rj+s*Nsite][ri+s*Nsite];
        //int idx_reverse = ijst_to_idx[ri+s*Nsite][(2*ri-rj+Nsite)%Nsite+s*Nsite];  // BEWARE: this line impose invariance under translation and assume periodic boundary condition
        
        int t     = ChargeExcitationIdx[idx_exc][0];  // type
        int dr1_x = ChargeExcitationIdx[idx_exc][1];  // position 1
        int dr1_y = ChargeExcitationIdx[idx_exc][2];  // position 1
        int dr2_x = ChargeExcitationIdx[idx_exc][3];  // position 2
        int dr2_y = ChargeExcitationIdx[idx_exc][4];  // position 2
        
        int ra1 = find_neighbor_site(ri,dr1_x,dr1_y);
        int ra2 = find_neighbor_site(ri,dr2_x,dr2_y);
        int rb1 = find_neighbor_site(rj,dr1_x,dr1_y);
        int rb2 = find_neighbor_site(rj,dr2_x,dr2_y);
        
        int idx_vector1 = idx_exc + NExcitation*(sample + sampleChunk*idx1);
        int idx_vector2 = idx_exc + NExcitation*(sample + sampleChunk*idx2);
        
        //printf("%d %d %d %d %d\n", ChargeExcitationIdx[idx_exc][0],ChargeExcitationIdx[idx_exc][1],ChargeExcitationIdx[idx_exc][2],ChargeExcitationIdx[idx_exc][3],ChargeExcitationIdx[idx_exc][4]);
        //printf("%d %d \n", Dimension_1, Dimension_2); fflush(stdout);
        //printf("%d %d %d %d \n", ra1,ra1,rb1,rb2);
        
        //printf("%d %d %d \n", ChargeExcitationIdx[idx_exc][0],ChargeExcitationIdx[idx_exc][1],ChargeExcitationIdx[idx_exc][2]);
        
        // <phi|ca|x> / <phi|x>
        CA_tmp = Local_CA[idx1] * Commute_Nat_(with_CisAjs,  ra1, ra2, t, ri, rj, s, 0,0,0, myEleNum);        
        O_CA_vec1[idx_vector1] = CA_tmp;        
        O_CA_vec2[idx_vector2] = CA_tmp;//conj(CA_tmp);        

        // <phi|ac|x> / <phi|x> = delta_{ri,rj} * <phi|x> / <phi|x> - <phi|ca|x> / <phi|x>           <-- need to reverse indices
        AC_tmp  = del(ri,rj);
        AC_tmp -= Local_CA[idx2];
        AC_tmp *= Commute_Nat_(with_AisCjs, ra1, ra2, t, ri, rj, s, 0,0,0, myEleNum);
        //AC_tmp = Local_AC[idx] * Commute_Nat_(with_AisCjs,  ra1, ra2, t, ri, rj, s, 0,0,0, myEleNum);        
        O_AC_vec1[idx_vector1] = AC_tmp;
        O_AC_vec2[idx_vector2] = AC_tmp;//conj(AC_tmp);
        
        // <phi|x>
        O0_vec1[idx_vector1]  = w * ((double) (Commute_Nat_(with_nothing, rb1, rb2, t,  0,  0, s, 0,0,0, myEleNum)));
        O0_vec2[idx_vector2]  = w * ((double) (Commute_Nat_(with_nothing, rb1, rb2, t,  0,  0, s, 0,0,0, myEleNum)));
        //
        // <phi|H_U|x> 
        double tmp_int_AHC=0.0;
        double tmp_int_CHA=0.0;
        
        ///*
        for(idx_int=0;idx_int<NCoulombIntra;idx_int++) {
          rm = CoulombIntra[idx_int];
          factor = ParaCoulombIntra[idx_int] *
                   ( 1.*del(rj,rm) + myEleNum[rm+s*Nsite])*myEleNum[rm+(1-s)*Nsite];
          tmp_int_AHC += factor;        

          factor = ParaCoulombIntra[idx_int] *
                   ( -1.*del(rj,rm) + myEleNum[rm+s*Nsite])*myEleNum[rm+(1-s)*Nsite];
          tmp_int_CHA += factor;
        }
        
        H_AC_vec1[idx_vector1] = tmp_int_AHC * AC_tmp ;
        H_CA_vec1[idx_vector1] = tmp_int_CHA * CA_tmp ;
        H_AC_vec2[idx_vector2] = tmp_int_AHC * AC_tmp ;
        H_CA_vec2[idx_vector2] = tmp_int_CHA * CA_tmp ;
        
        // <phi|H_T|x> / <phi|x>
        for(idx_trans=0;idx_trans<NTransfer;idx_trans++) {
          //*
          
          rm = Transfer[idx_trans][0];
          rn = Transfer[idx_trans][2];
          u  = Transfer[idx_trans][3];
          
          int idx_green0 = ri+Nsite*rn;
          tmp = -1.0 * ParaTransfer[idx_trans] 
                     * (Local_CisAjsCmuAnu[idx_trans][idx1] 
                        - del(rm,rj) * del(s,u) * Local_CA[idx_green0] ) 
                     * Commute_Nat_(with_CisCmuAnuAjs, ra1, ra2, t, ri, rj, s, rm, rn, u, myEleNum) ;
          H_CA_vec1[idx_vector1] += tmp;
          H_CA_vec2[idx_vector2] += tmp;//conj(tmp);
          
          int idx_green1 = rj+Nsite*ri;
          int idx_green2 = rm+Nsite*rn+Nsite*Nsite*u;
          int idx_green3 = rm+Nsite*ri;
          //int idx_green1 = ijst_to_idx[ri+s*Nsite][rm+s*Nsite];
          tmp = -1.0 * ParaTransfer[idx_trans]
                     * ( - Local_CisAjsCmuAnu[idx_trans][idx_green1] 
                         + del(ri,rj) * Local_CA[idx_green2] 
                         + del(rn,rj) * del(s,u) * ( del(rm,ri) - Local_CA[idx_green3] ))
                     * Commute_Nat_(with_AisCmuAnuCjs, ra1, ra2, t, ri, rj, s, rm, rn, u, myEleNum) ;
          //tmp = -1.0 * ParaTransfer[idx_trans] 
          //           * (Local_AisCjsCmuAnu[idx_trans][idx] 
          //              + del(rn,rj) * del(s,u) * Local_AC[idx_green1] ) 
          //           * Commute_Nat_(with_AisCmuAnuCjs, ra1, ra2, t, ri, rj, s, rm, rn, u, myEleNum) ;
          H_AC_vec1[idx_vector1] += tmp;
          H_AC_vec2[idx_vector2] += tmp;//conj(tmp);
          //
        }
        ///
      }
     }
    }//*/
    free(myEleIdx);
    free(myEleNum);
    //free(myEleCfg);
    free(myProjCntNew);
    free(myBuffer_real);
  } //
  
  return;
}







void CalculateGreenFunc(const double w, const double complex ip, int *eleIdx, int *eleCfg,
                        int *eleNum, int *eleProjCnt) {

  int idx,idx0,idx1;
  int ri,rj,s,rk,rl,t;
  double complex tmp;
  int *myEleIdx, *myEleNum, *myProjCntNew;
  double complex *myBuffer;
  
  RequestWorkSpaceThreadInt(Nsize+Nsite2+NProj);
  RequestWorkSpaceThreadComplex(NQPFull+2*Nsize);
  /* GreenFunc1: NQPFull, GreenFunc2: NQPFull+2*Nsize */

  #pragma omp parallel default(shared)                \
  private(myEleIdx,myEleNum,myProjCntNew,myBuffer,idx)
  {
    myEleIdx = GetWorkSpaceThreadInt(Nsize);
    myEleNum = GetWorkSpaceThreadInt(Nsite2);
    myProjCntNew = GetWorkSpaceThreadInt(NProj);
    myBuffer = GetWorkSpaceThreadComplex(NQPFull+2*Nsize);

    #pragma loop noalias
    for(idx=0;idx<Nsize;idx++) myEleIdx[idx] = eleIdx[idx];
    #pragma loop noalias
    for(idx=0;idx<Nsite2;idx++) myEleNum[idx] = eleNum[idx];

    #pragma omp master
    {StartTimer(50);}

    #pragma omp for private(idx,ri,rj,s,tmp) schedule(dynamic) nowait
    for(idx=0;idx<NCisAjs;idx++) {
      ri = CisAjsIdx[idx][0];
      rj = CisAjsIdx[idx][2];
      s  = CisAjsIdx[idx][3];
      tmp = GreenFunc1(ri,rj,s,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,
                       myProjCntNew,myBuffer);
      LocalCisAjs[idx] = tmp;
    }
    #pragma omp master
    {StopTimer(50);StartTimer(51);}
    
    #pragma omp for private(idx,ri,rj,s,rk,rl,t,tmp) schedule(dynamic)
    for(idx=0;idx<NCisAjsCktAltDC;idx++) {
      ri = CisAjsCktAltDCIdx[idx][0];
      rj = CisAjsCktAltDCIdx[idx][2];
      s  = CisAjsCktAltDCIdx[idx][1];
      rk = CisAjsCktAltDCIdx[idx][4];
      rl = CisAjsCktAltDCIdx[idx][6];
      t  = CisAjsCktAltDCIdx[idx][5];

      tmp = GreenFunc2(ri,rj,rk,rl,s,t,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,
                       myProjCntNew,myBuffer);
      PhysCisAjsCktAltDC[idx] += w*tmp;
    }
    
    #pragma omp master
    {StopTimer(51);StartTimer(52);}

    #pragma omp for private(idx) nowait
    for(idx=0;idx<NCisAjs;idx++) {
      PhysCisAjs[idx] += w*LocalCisAjs[idx];
    }

    #pragma omp master
    {StopTimer(52);StartTimer(53);}

    #pragma omp for private(idx,idx0,idx1) nowait
    for(idx=0;idx<NCisAjsCktAlt;idx++) {
      idx0 = CisAjsCktAltIdx[idx][0];
      idx1 = CisAjsCktAltIdx[idx][1];
      PhysCisAjsCktAlt[idx] += w*LocalCisAjs[idx0]*conj(LocalCisAjs[idx1]);// TBC conj ok?
    }

    #pragma omp master
    {StopTimer(53);}
  }

  ReleaseWorkSpaceThreadInt();
  ReleaseWorkSpaceThreadComplex();
  return;
}


void CalculateGreenFuncBF(const double w, const double ip, int *eleIdx, int *eleCfg,
                          int *eleNum, int *eleProjCnt, const int *eleProjBFCnt) {

  int idx,idx0,idx1;
  int ri,rj,s,rk,rl,t;
  double complex tmp;
  int *myEleIdx, *myEleNum, *myEleCfg, *myProjCntNew, *myProjBFCntNew;
  double complex* mySltBFTmp;
  double complex* myBuffer;

  RequestWorkSpaceThreadInt(Nsize+2*Nsite2+NProj+16*Nsite*Nrange);
  RequestWorkSpaceThreadComplex(NQPFull+2*Nsize+NQPFull*Nsite2*Nsite2);
  /* GreenFunc1: NQPFull, GreenFunc2: NQPFull+2*Nsize */

#pragma omp parallel default(shared)\
  private(myEleIdx,myEleNum,myEleCfg,myProjCntNew,myProjBFCntNew,myBuffer,mySltBFTmp)
  {
    myEleIdx = GetWorkSpaceThreadInt(Nsize);
    myEleNum = GetWorkSpaceThreadInt(Nsite2);
    myEleCfg = GetWorkSpaceThreadInt(Nsite2);
    myProjCntNew   = GetWorkSpaceThreadInt(NProj);
    myProjBFCntNew = GetWorkSpaceThreadInt(16*Nsite*Nrange);
    myBuffer   = GetWorkSpaceThreadComplex(NQPFull+2*Nsize);
    mySltBFTmp = GetWorkSpaceThreadComplex(NQPFull*Nsite2*Nsite2);

#pragma loop noalias
    for(idx=0;idx<Nsize;idx++) myEleIdx[idx] = eleIdx[idx];
#pragma loop noalias
    for(idx=0;idx<Nsite2;idx++) myEleNum[idx] = eleNum[idx];
#pragma loop noalias
    for(idx=0;idx<Nsite2;idx++) myEleCfg[idx] = eleCfg[idx];

    StoreSlaterElmBF_fcmp(mySltBFTmp);
#pragma omp master
    {StartTimer(50);}

#pragma omp for private(idx,ri,rj,s,tmp) schedule(dynamic) nowait
    for(idx=0;idx<NCisAjs;idx++) {
      ri = CisAjsIdx[idx][0];
      rj = CisAjsIdx[idx][2];
      s  = CisAjsIdx[idx][3];
      tmp = GreenFunc1BF(ri,rj,s,ip,mySltBFTmp,myEleIdx,myEleCfg,myEleNum,eleProjCnt,myProjCntNew,eleProjBFCnt,myProjBFCntNew,myBuffer);
      LocalCisAjs[idx] = tmp;
    }

#pragma omp master
    {StopTimer(50);StartTimer(51);}

#pragma omp for private(idx,ri,rj,s,rk,rl,t,tmp) schedule(dynamic)
    for(idx=0;idx<NCisAjsCktAltDC;idx++) {
      ri = CisAjsCktAltDCIdx[idx][0];
      rj = CisAjsCktAltDCIdx[idx][2];
      s  = CisAjsCktAltDCIdx[idx][1];
      rk = CisAjsCktAltDCIdx[idx][4];
      rl = CisAjsCktAltDCIdx[idx][6];
      t  = CisAjsCktAltDCIdx[idx][5];
      tmp = GreenFunc2(ri,rj,rk,rl,s,t,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,
                       myProjCntNew,myBuffer);
      PhysCisAjsCktAltDC[idx] += w*tmp;
    }

#pragma omp master
    {StopTimer(51);StartTimer(52);}

#pragma omp for private(idx) nowait
    for(idx=0;idx<NCisAjs;idx++) {
      PhysCisAjs[idx] += w*LocalCisAjs[idx];
    }

#pragma omp master
    {StopTimer(52);StartTimer(53);}

#pragma omp for private(idx,idx0,idx1) nowait
    for(idx=0;idx<NCisAjsCktAlt;idx++) {
      idx0 = CisAjsCktAltIdx[idx][0];
      idx1 = CisAjsCktAltIdx[idx][1];
      PhysCisAjsCktAlt[idx] += w*LocalCisAjs[idx0]*LocalCisAjs[idx1];
    }

#pragma omp master
    {StopTimer(53);}
  }

  ReleaseWorkSpaceThreadInt();
  ReleaseWorkSpaceThreadComplex();
  return;
}
#endif
