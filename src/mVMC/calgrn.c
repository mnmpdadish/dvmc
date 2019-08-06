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



void CalculateStaticQuantities_real(const double w, 
                               int *eleIdx, int *eleCfg,
                               int *eleNum, int *eleProjCnt) {

  int idx;
  int ri,rj,s;
  int *myEleNum;
  
  RequestWorkSpaceThreadInt(Nsize+Nsite2+NProj);
  {
    myEleNum = GetWorkSpaceThreadInt(Nsite2);
    for(idx=0;idx<Nsite2;idx++) myEleNum[idx] = eleNum[idx];

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
  }
  ReleaseWorkSpaceThreadInt();
  return;
}




int del(int i,int j){
  return ((i==j)? 1:0);
}


// This function defines every type of excitation (0 to 7) and its 
// anticommutation with different combination of creation operator (C)
// and annihilation operator (A)

int Commute_Nat_(commuting_with commuting, int ra, int rb, int t, int ri, int rj, int s, int rm, int rn, int u, int *eleNum) {
  int sign;
  if((commuting==with_CisCmuAnuAjs) || (commuting==with_AisCmuAnuCjs)) {
    if(commuting==with_CisCmuAnuAjs) {sign = 1;}
    else if(commuting==with_AisCmuAnuCjs) {sign = -1;}
    
    if(t==0){ 
      return 1;
    }
    else if(t==1){ 
      return (eleNum[ra+(1-s)*Nsite] + del((1-s),u) * (del(ra,rm)-del(ra,rn)));
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
    else if(t==5){
      return (eleNum[ri+(1-s)*Nsite] + del((1-s),u) * (del(ri,rm)-del(ri,rn)))
            *(eleNum[ra+   s *Nsite] + (del(ra,rm)-del(ra,rn)) * del(u,s)    + sign * (del(ra,ri)-del(ra,rj)) )
            *(eleNum[rb+   s *Nsite] + (del(rb,rm)-del(rb,rn)) * del(u,s)    + sign * (del(rb,ri)-del(rb,rj)) );
    }
    else if(t==6){
      return (eleNum[ri+(1-s)*Nsite] + del((1-s),u) * (del(ri,rm)-del(ri,rn)))
            *(eleNum[ra+(1-s)*Nsite] + (del(ra,rm)-del(ra,rn)) * del(u,1-s))
            *(eleNum[rb+(1-s)*Nsite] + (del(rb,rm)-del(rb,rn)) * del(u,1-s));
    }
    else if(t==7){
      return (eleNum[ri+(1-s)*Nsite] + del((1-s),u) * (del(ri,rm)-del(ri,rn)))
            *(eleNum[ra+(1-s)*Nsite] + (del(ra,rm)-del(ra,rn)) * del(u,1-s))
            *(eleNum[rb+   s *Nsite] + (del(rb,rm)-del(rb,rn)) * del(u,s)    + sign * (del(rb,ri)-del(rb,rj)) )   ;
    }
//    else if(t==5){
//      return (eleNum[ra+(1-s)*Nsite] + (del(ra,rm)-del(ra,rn)) * del(u,1-s)) * (eleNum[ra+s*Nsite] + (del(ra,rm)-del(ra,rn)) * del(u,s) + sign * (del(ra,ri)-del(ra,rj)))
//            *(eleNum[rb+(1-s)*Nsite] + (del(rb,rm)-del(rb,rn)) * del(u,1-s)) * (eleNum[rb+s*Nsite] + (del(rb,rm)-del(rb,rn)) * del(u,s) + sign * (del(rb,ri)-del(rb,rj)));
//    }
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
      return (eleNum[ra+(1-s)*Nsite]) ;
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
    else if(t==5){ 
      return (eleNum[ri+(1-s)*Nsite])
            *(eleNum[ra+    s*Nsite] + sign * (del(ra,ri)-del(ra,rj) ) )
            *(eleNum[rb+    s*Nsite] + sign * (del(rb,ri)-del(rb,rj) ) );      
    }
    else if(t==6){ 
      return (eleNum[ri+(1-s)*Nsite])
            *(eleNum[ra+(1-s)*Nsite] )
            *(eleNum[rb+(1-s)*Nsite] );
    }
    else if(t==7){ 
      return (eleNum[ri+(1-s)*Nsite])
            *(eleNum[ra+(1-s)*Nsite] )
            *(eleNum[rb+    s*Nsite] + sign * (del(rb,ri)-del(rb,rj) ) );
    }
//    else if(t==5){ 
//      return (eleNum[ra+s*Nsite] + sign * (del(ra,ri)-del(ra,rj) ) ) * (eleNum[ra+(1-s)*Nsite])
//            *(eleNum[rb+s*Nsite] + sign * (del(rb,ri)-del(rb,rj) ) ) * (eleNum[rb+(1-s)*Nsite]);      
//    }
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
      return (eleNum[ra+(1-s)*Nsite]);
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
    else if(t==5){ 
      return (eleNum[ri+(1-s)*Nsite])
            *(eleNum[ra+s*Nsite]*eleNum[rb+s*Nsite]);
    }
    else if(t==6){ 
      return (eleNum[ri+(1-s)*Nsite])
            *(eleNum[ra+(1-s)*Nsite]*eleNum[rb+(1-s)*Nsite]);
    }
    else if(t==7){ 
      return (eleNum[ri+(1-s)*Nsite])
            *(eleNum[ra+(1-s)*Nsite]*eleNum[rb+s*Nsite]);
    }
//    else if(t==5){ 
//      return (eleNum[ra+(1-s)*Nsite]*eleNum[ra+s*Nsite]) * (eleNum[rb+(1-s)*Nsite]*eleNum[rb+s*Nsite]);
//    }
    else{
      printf("oups, error %d\n", t);
      exit(1);
      return -1;
    }
  }
  exit(1);
  return -2; 
}

// retunr i%N, in the same way that python treat the negative numbers
int moduloPython(int i,int N){
  return ((i % N) + N) % N;
}

int find_neighbor_difference(int ri,int rj){
  int ri_x = moduloPython(ri,              Excitation_L);
  int ri_y = moduloPython(ri/Excitation_L, Excitation_W);
  int rj_x = moduloPython(rj,              Excitation_L);
  int rj_y = moduloPython(rj/Excitation_L, Excitation_W);

  int dr_out = moduloPython(rj_x-ri_x,Excitation_L) + Excitation_L*moduloPython(rj_y-ri_y,Excitation_W);
  //int dr_out = moduloPython(moduloPython(rj_x-ri_x,Nsite) + Excitation_L*moduloPython(rj_y-ri_y,Nsite),Nsite);    
  return dr_out;
}

int find_neighbor_site(int r,int dx,int dy){
  int r_x   = moduloPython(r + dx             , Excitation_L);  
  int r_y   = moduloPython(r/Excitation_L+dy  , Excitation_W);  
  int r_out = moduloPython(r_x + Excitation_L*r_y,Nsite);    
  return r_out;
}

int find_neighbor_site2(int r,int dr){
  int dx = moduloPython(dr,Excitation_W);
  int dy = moduloPython(dr/Excitation_W,Excitation_L);
  return find_neighbor_site(r,dx,dy);
}


//C=C+weight*A*B
unsigned int C_ADD_AxB(double * C, double const * A, double const * B, int N, double weight, int sampleSize) {
  char transA= 'N', transB= 'C';
  double beta=1.0;
  //int ONE = 1;
  M_DGEMM(&transA,&transB,&N,&N,&sampleSize, &weight, &A[0], &N, &B[0], &N, &beta, &C[0], &N); 
  // The function M_DGEMM calculate the same as the following 8 lines
  // but faster, providing that the user set 
  /*
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


void CalculateDynamicalGreenFunc_real(const double w, const double ip, 
                                      int *eleIdx, int *eleCfg,
                                      int *eleNum, int *eleProjCnt, int sample) {
  int idx;
  int ri,rj,s;
  int *myEleIdx, *myEleNum, *myProjCntNew; //, *myEleCfg
  double *myBuffer_real;
  
  RequestWorkSpaceThreadInt(Nsize+Nsite2+NProj);
  RequestWorkSpaceThreadDouble(NQPFull+2*Nsize);
  
#pragma omp parallel default(shared)\
  private(myEleIdx,myEleNum,myProjCntNew,myBuffer_real)
  {
    myEleIdx = GetWorkSpaceThreadInt(Nsize);
    myEleNum = GetWorkSpaceThreadInt(Nsite2);
    myProjCntNew   = GetWorkSpaceThreadInt(NProj);
    myBuffer_real  = GetWorkSpaceThreadDouble(NQPFull+2*Nsize);
    
    #pragma loop noalias
    for(idx=0;idx<Nsize; idx++) myEleIdx[idx] = eleIdx[idx];
    #pragma loop noalias
    for(idx=0;idx<Nsite2;idx++) myEleNum[idx] = eleNum[idx];
    
    int rm, rn, u;
    int idx_int, idx_trans;


    #pragma omp for private(idx,ri,rj,s) schedule(dynamic) 
    for(ri=0;ri<Nsite;ri++) {
     for(rj=0;rj<Nsite;rj++) {
      for(s=0;s<2;s++) {
       Local_CA[ri+Nsite*rj+Nsite*Nsite*s] = GreenFunc1_real(ri,rj,s,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,
                                                              myProjCntNew,myBuffer_real); 
      }
     }
    }
    #pragma omp for private(idx,ri,rj,s) schedule(dynamic) nowait
    for(ri=0;ri<Nsite;ri++) {
     for(rj=0;rj<Nsite;rj++) {
      for(s=0;s<2;s++) {
       Phys_CA[ri+Nsite*rj+Nsite*Nsite*s] += w*Local_CA[ri+Nsite*rj+Nsite*Nsite*s];
      }
     }
    }

    #pragma omp for private(idx,ri,rj,s,idx_trans,rm,rn,u) schedule(dynamic) 
    for(ri=0;ri<Nsite;ri++) {
     for(rj=0;rj<Nsite;rj++) {
      s=0;

      for(idx_trans=0;idx_trans<NTransfer;idx_trans++) {
      
        rm = Transfer[idx_trans][0];
        rn = Transfer[idx_trans][2];
        u  = Transfer[idx_trans][3];
        
        Local_CisAjsCmuAnu[idx_trans][ri+Nsite*rj] = GreenFunc2_real(ri,rj,rm,rn,s,u,ip,
                        myEleIdx,eleCfg,myEleNum,eleProjCnt,myProjCntNew,myBuffer_real);
      }
     }
    }
    
    #pragma omp for private(idx,ri,rj,s,idx_trans,rm,rn,u) schedule(dynamic) nowait
    for(ri=0;ri<Nsite;ri++) {
     for(rj=0;rj<Nsite;rj++) {
      s=0;
      int idx1 = ri+Nsite*rj;
      int idx2 = rj+Nsite*ri;
      int idx_exc;
          
      for(idx_exc=0;idx_exc<NExcitation;idx_exc++){
        
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
        //printf("%d %d \n", Excitation_L, Excitation_W); fflush(stdout);
        //printf("%d %d %d %d \n", ra1,ra1,rb1,rb2);
        
        //printf("%d %d %d \n", ChargeExcitationIdx[idx_exc][0],ChargeExcitationIdx[idx_exc][1],ChargeExcitationIdx[idx_exc][2]);
        
        // <phi|ca|x> / <phi|x>
        double CA_tmp = Local_CA[idx1] * Commute_Nat_(with_CisAjs,  ra1, ra2, t, ri, rj, s, 0,0,0, myEleNum);        
        O_CA_vec1[idx_vector1] = CA_tmp;        
        O_CA_vec2[idx_vector2] = CA_tmp;//conj(CA_tmp);        

        // <phi|ac|x> / <phi|x> = delta_{ri,rj} * <phi|x> / <phi|x> - <phi|ca|x> / <phi|x>           <-- need to reverse indices
        double AC_tmp  = del(ri,rj);
        AC_tmp -= Local_CA[idx2];
        AC_tmp *= Commute_Nat_(with_AisCjs, ra1, ra2, t, ri, rj, s, 0,0,0, myEleNum);
        //AC_tmp = Local_AC[idx] * Commute_Nat_(with_AisCjs,  ra1, ra2, t, ri, rj, s, 0,0,0, myEleNum);        
        O_AC_vec1[idx_vector1] = AC_tmp;
        O_AC_vec2[idx_vector2] = AC_tmp;//conj(AC_tmp);
        
        // <phi|x>
        O0_vec1[idx_vector1]  = w * ((double) (Commute_Nat_(with_nothing, rb1, rb2, t,  rj, ri, s, 0,0,0, myEleNum)));
        O0_vec2[idx_vector2]  = w * ((double) (Commute_Nat_(with_nothing, rb1, rb2, t,  rj, ri, s, 0,0,0, myEleNum)));
        //
        // <phi|H_U|x> 
        double tmp_int_AHC=0.0;
        double tmp_int_CHA=0.0;
        
        ///*
        for(idx_int=0;idx_int<NCoulombIntra;idx_int++) {
          rm = CoulombIntra[idx_int];
          double factor = ParaCoulombIntra[idx_int] *
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
          double tmp = -1.0 * ParaTransfer[idx_trans] 
                     * (Local_CisAjsCmuAnu[idx_trans][idx1] 
                        - del(rm,rj) * del(s,u) * Local_CA[idx_green0] ) 
                     * Commute_Nat_(with_CisCmuAnuAjs, ra1, ra2, t, ri, rj, s, rm, rn, u, myEleNum) ;
          H_CA_vec1[idx_vector1] += tmp;
          H_CA_vec2[idx_vector2] += tmp;//conj(tmp);
          
          int idx_green1 = rj+Nsite*ri;
          int idx_green2 = rm+Nsite*rn+Nsite*Nsite*u;
          int idx_green3 = rm+Nsite*ri;
          tmp = -1.0 * ParaTransfer[idx_trans]
                     * ( - Local_CisAjsCmuAnu[idx_trans][idx_green1] 
                         + del(ri,rj) * Local_CA[idx_green2] 
                         + del(rn,rj) * del(s,u) * ( del(rm,ri) - Local_CA[idx_green3] ))
                     * Commute_Nat_(with_AisCmuAnuCjs, ra1, ra2, t, ri, rj, s, rm, rn, u, myEleNum) ;
          H_AC_vec1[idx_vector1] += tmp;
          H_AC_vec2[idx_vector2] += tmp;//conj(tmp);
        }
      }
     }
    }//*/
    //free(myEleIdx);
    //free(myEleNum);
    ////free(myEleCfg);
    //free(myProjCntNew);
    //free(myBuffer_real);
  } //

  ReleaseWorkSpaceThreadInt();
  ReleaseWorkSpaceThreadComplex();  
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
