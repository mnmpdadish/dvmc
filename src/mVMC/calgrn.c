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



/*Function to Multiply charge or doublon Factors ("N", "D") from left or right TJS
double complex MultiplyFactor(double complex original,
			      int n, int *rsk, int *rsl,
			      int LeftFactor, int ri, int p,
			      int RightFactor, int rj, int q, int *eleNum) {
    
  int Factor = 1;
  // Multiply Factor from Right
  if(RightFactor==1){
    Factor = Factor * eleNum[rj+q*Nsite];
  }
  
  int i;
  int rpi = ri+p*Nsite;

  // Multiply Factor from Left
  if(LeftFactor == 1){
    int modification = 0;
    for(i=0;i<n;i++){
      if(rsk[i]==rpi){modification += 1;}
      if(rsl[i]==rpi){modification -= 1;}
    }
    Factor = Factor * (eleNum[rpi]+modification);
  }
      
  return ((double)Factor *original);//TBC
  
}

/* Function to calculate GFs with charge and doublon factors 
void MultiplyFactor2GFs(int idx, int NumElem, int n, int *rsk, int *rsl,int *eleNum,
			int flip, int PH,
			double complex original,
			double complex ACM,  //M is like N, but for opposite spin of C
			double complex MAC,  
			double complex MACM) {

  int rk,rl,s;
  int ri,p;
  int rj,q;
  
  rk = CisAjsIdx[idx][0];
  rl = CisAjsIdx[idx][2];
  s  = CisAjsIdx[idx][3];
  
  int idx_ex = ijst_to_idx[rl+s*Nsite][rk+s*Nsite]; 

  if (PH==1){
    rk = CisAjsIdx[idx][2];
    rl = CisAjsIdx[idx][0];
  }
  
  double complex tmp;

  tmp = MultiplyFactor(original,n,rsk,rsl,0,0,0,1,rk,1-s,eleNum);
  if(flip==0){
    AU[idx] += tmp;
    UC[idx_ex] += conj(tmp);  
    //  UC[idx] += MultiplyFactor(original,n,rsk,rsl,1,rl,1-s,0,0,0,eleNum);
    UU[idx] +=    MultiplyFactor(original,n,rsk,rsl,1,rl,1-s,1,rk,1-s,eleNum);
  }
  else{
    UC[idx_ex] += tmp;
    AU[idx] += conj(tmp);  
    UU[idx_ex] += MultiplyFactor(original,n,rsk,rsl,1,rl,1-s,1,rk,1-s,eleNum);
  }  
} // end of MultiplyFactor2GFs

/**/


void CalculateGreenFuncMoments(const double w, const double complex ip, int *eleIdx, int *eleCfg,
                        int *eleNum, int *eleProjCnt) {

  int idx,idx0,idx1;
  int ri,rj,s,rk,rl,t;
  double complex tmp;
  int *myEleIdx, *myEleNum, *myProjCntNew;
  double complex *myBuffer;

  RequestWorkSpaceThreadInt(Nsize+Nsite2+NProj);
  RequestWorkSpaceThreadComplex(NQPFull+2*Nsize);
  /* GreenFunc1: NQPFull, GreenFunc2: NQPFull+2*Nsize */

  #pragma omp parallel default(shared)		\
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
    {
      //printf("start:\n");
      StartTimer(50);
    }

    #pragma omp for private(idx,ri,rj,s,tmp) schedule(dynamic) nowait
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
      
      //n_sigma (1-n_sigma) Maxime
      PhysN2[idx+NCisAjs*3] += w*myEleNum[ri+s*Nsite] *(1.0-myEleNum[rj+s*Nsite]);
    }

    #pragma omp for private(ri,s) schedule(dynamic) nowait
    for(ri=0;ri<Nsite;ri++) {
      s=0;
      //Doublon Maxime
      PhysN1[ri+Nsite*0] += w*myEleNum[ri]*myEleNum[ri+Nsite];

	    //Holon Maxime
      PhysN1[ri+Nsite*1] += w*(1.0-myEleNum[ri])*(1.0-myEleNum[ri+Nsite]);
	    
      //Density_up (s==0) Maxime
      PhysN1[ri+Nsite*2] += w*myEleNum[ri];

      //Density_down (s==1) Maxime
      PhysN1[ri+Nsite*3] += w*myEleNum[ri+Nsite];	   
    }
    
    #pragma omp master
    {StopTimer(50);}

  }


  /*Local two-body Green's fuction LocalCktAltCmuAnu TJS*/
  int idx_int, idx_trans;
  int rm, rn, u;

/*
#pragma omp for private(idx,idx_trans,rk,rl,t,rm,rn,u) schedule(dynamic) nowait
  for(idx=0;idx<NCisAjs;idx++) {
    rk = CisAjsIdx[idx][0];
    rl = CisAjsIdx[idx][2];
    t  = CisAjsIdx[idx][3];

    for(idx_trans=0;idx_trans<NTransfer;idx_trans++) {
    
      rm = Transfer[idx_trans][0];
      rn = Transfer[idx_trans][2];
      u  = Transfer[idx_trans][3];
	
      LocalCktAltCmuAnu[idx_trans][idx] = GreenFunc2(rk,rl,rm,rn,t,u,ip,
             myEleIdx,eleCfg,myEleNum,eleProjCnt,myProjCntNew,myBuffer);
    }
  }

  /*Composite Fermion Energies TJS + MC
  double factor=0.0;
  double complex tmp_int = 0.0;
  double complex tmp_trans = 0.0;
  
#pragma omp for private(idx,idx_int,idx_trans,s,rk,rl,t,rm,rn,u,tmp,tmp_int,tmp_trans,factor) schedule(dynamic) nowait
  for(idx=0;idx<NCisAjs;idx++) {
      
    rk = CisAjsIdx[idx][0];
    rl = CisAjsIdx[idx][2];
    s  = CisAjsIdx[idx][3];
    
    //assuming the same spin for Ck and Al
    t = s;

    // <phi|Al|H_U|Ck|phi>    
    tmp_int=0.0;
    for(idx_int=0;idx_int<NCoulombIntra;idx_int++) {	

      rm = CoulombIntra[idx_int];
      factor = ParaCoulombIntra[idx_int]
              *(((rk==rm)? 1.:0.) + myEleNum[rm+s*Nsite])*myEleNum[rm+(1-s)*Nsite];
      tmp_int += factor*(((rk==rl)? 1.:0.) - LocalCisAjs[idx]);
    }
      
    // <phi|Al|H_T|Ck|phi>    
    tmp_trans = 0.0;
    for(idx_trans=0;idx_trans<NTransfer;idx_trans++) {
    
      rm = Transfer[idx_trans][0];
      rn = Transfer[idx_trans][2];
      u  = Transfer[idx_trans][3];

      factor = ParaTransfer[idx_trans];
      tmp = factor*LocalCktAltCmuAnu[idx_trans][idx];

      if(rk==rl){
        tmp += -factor*LocalCisAjs[ijst_to_idx[rm+u*Nsite][rn+u*Nsite]];
      }
      if((rk==rn)&&(s==u)){
        tmp += factor*LocalCisAjs[ijst_to_idx[rm+u*Nsite][rl+s*Nsite]];
        if((rl==rm)&&(s==u)){
          tmp += -factor;
        }
      }
      tmp_trans += tmp;
      LocalAHTCijsklm[idx_trans][idx] = tmp;
      MultiplyFactor2GFs(idx,NCisAjs,1,rsk,rsl,myEleNum,0,0,w*tmp,PhysACM,PhysMAC,PhysMACM);
    }
    Phys_G_Moments[idx] += w*(tmp_trans + tmp_int);
    
    
  }

*/
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

  #pragma omp parallel default(shared)		\
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
