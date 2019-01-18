/* mcparameters.h -- This file is part of Archimedes release 1.1.0.
   Archimedes is a simulator for Submicron 2D III-V semiconductor
   Devices. It implements the Monte Carlo method
   for the simulation of the semiclassical Boltzmann equation for both
   electrons and holes. It also includes the quantum effects by means 
   of effective potential method. It is now able to simulate applied
   magnetic fields along with self consistent Faraday equation.

   Copyright (C) 2004-2011 Jean Michel Sellier <jeanmichel.sellier@gmail.com>
 
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


// ######################################################
// Created on 05 sep.2004, Siracusa, J.M.Sellier
// Last modif. : 18 Aug.2011, Carry le Rouet, France, J.M.Sellier
// ######################################################

// Definition of the variables needed for taking into
// account scatterings from acoustic and optical,non-polar
// phonons (which are the most relevant scatterings in common semiconductors like Silicon)

void
MCparameters(int Material)
{
 real wo,no,w[8],n[8],dos,aco[5],oge[7],oga[7];
 real dos1,dos2,dos3,am1,am2,am3,cl,deq,dij,dde;
 // added conduction mass by Suzuki 18/9/23
 real am1_dos,am2_dos,am3_dos;
 real hwe,hwij,wij,we,ne,nij;
 real poe,poa,ope,opa,eqe,eqa,ode,oda,qmin,qmax;
 real initialenergy,sei,finalenergy,sef;
 real eps,epf,ep,bimp,cimp,cimpn,qd;
 real ak,qq,wk;
 int ie,i,j;
 real z2=4.;

 ISEED = 38467.;  //  initial value for random number generator

// These definitions are valid for every material
 BKTQ=KB*TL/Q; // in eV
 QH=Q/HBAR;

// Material with 2 valleys
// #######################
 if(NOVALLEY[Material]==2){
// Effective mass for the GAMMA and L-Valley, respectively
  am1=MSTAR[Material][1][0]*M;
  am2=MSTAR[Material][2][0]*M;
  eps=EPSR[Material]*EPS0;
  epf=EPF[Material]*EPS0;
  ep=1./(1./epf-1./eps);

// Parameters for Phonon Scattering
  cl=RHO[Material]*pow(UL[Material],2.);
  deq=dij=DTK[Material][0]*Q;
  hwe=hwij=HWO[Material][0];

  SMH[Material][1]=sqrt(2.*am1*Q)/HBAR;
  SMH[Material][2]=sqrt(2.*am2*Q)/HBAR;
  HHM[Material][1]=HBAR*HBAR/(2.*am1*Q);
  HHM[Material][2]=HBAR*HBAR/(2.*am2*Q);
  HM[Material][1]=HBAR/am1;
  HM[Material][2]=HBAR/am2;

  wo=HWO[Material][0]*Q/HBAR;
  wij=hwij*Q/HBAR;
  we=hwe*Q/HBAR;

  no=1./(exp(HWO[Material][0]/BKTQ)-1.);
  nij=1./(exp(hwij/BKTQ)-1.);
  ne=1./(exp(hwe/BKTQ)-1.);

  dos1=pow(sqrt(2.*am1*Q)/HBAR,3.)/pow(2.*PI,2.);
  dos2=pow(sqrt(2.*am2*Q)/HBAR,3.)/pow(2.*PI,2.);

  poe=Q/8./PI/ep*Q*wo*(no+1.);
  poa=poe*no/(1.+no);
  aco[0]=2.*PI*DA[Material][0]/Q*DA[Material][0]*BKTQ/HBAR*Q/cl;
  ope=PI*dij/wij*dij/RHO[Material]/Q*(nij+1.);
  opa=ope*nij/(1.+nij);
  eqe=PI*deq/we*deq/RHO[Material]/Q*(ne+1.);
  eqa=eqe*ne/(1.+ne);

// Parameters for impurity scatterings
  cimp=CIMP; // <--- impurity concentration
  qd=sqrt(Q*cimp/BKTQ/eps);
  QD2=qd*qd;
  bimp=2.*PI*cimp*Q*Q/HBAR*Q/eps/eps;

// Calculation of scattering rates
  for(ie=1;ie<=DIME;ie++){
   initialenergy=DE*((real)(ie));
   sei=sqrt(initialenergy);
// GAMMA-valley
// ============
   if(OPTICALPHONONS==ON){
// Polar optical phonon
// Emission - Parabolic
    finalenergy=initialenergy-HWO[Material][0];
    if(finalenergy>0.){
     sef=sqrt(finalenergy);
     qmax=sef+sei;
     qmin=sei-sef;
     SWK[Material][1][1][ie]=poe*SMH[Material][1]*sei/initialenergy/Q*log(qmax/qmin);
    }
    else SWK[Material][1][1][ie]=0.;
// Absorption - Parabolic
    finalenergy=initialenergy+HWO[Material][0];
    sef=sqrt(finalenergy);
    qmax=sef+sei;
    qmin=sef-sei;
    SWK[Material][1][2][ie]=SWK[Material][1][1][ie]
                 +poa*SMH[Material][1]*sei/initialenergy/Q*log(qmax/qmin);
// Emission
    finalenergy=initialenergy-hwij+EMIN[Material][1]-EMIN[Material][2];
    if(finalenergy>0.){
     sef=sqrt(finalenergy*(1.+alphaK[Material][2]*finalenergy));
     SWK[Material][1][3][ie]=SWK[Material][1][2][ie]
                            +z2*ope*sef*dos2*(1.+2.*alphaK[Material][2]*finalenergy);
    }
    else SWK[Material][1][3][ie]=SWK[Material][1][2][ie];
// Absorption
    finalenergy=initialenergy+hwij+EMIN[Material][1]-EMIN[Material][2];
    if(finalenergy>0.){
     sef=sqrt(finalenergy*(1.+alphaK[Material][2]*finalenergy));
     SWK[Material][1][4][ie]=SWK[Material][1][3][ie]
                            +z2*opa*sef*dos2*(1.+2.*alphaK[Material][2]*finalenergy);
    }
    else SWK[Material][1][4][ie]=SWK[Material][1][3][ie];
   }
   else{
    // NO OPTICAL SCATTERING
    SWK[Material][1][1][ie]=0.0;
    SWK[Material][1][2][ie]=0.0;
    SWK[Material][1][3][ie]=0.0;
    SWK[Material][1][4][ie]=0.0;
   }
   if(ACOUSTICPHONONS==ON){
// Acoustic Phonon
    finalenergy=initialenergy;
    sef=sqrt(finalenergy*(1.+alphaK[Material][1]*finalenergy));
    SWK[Material][1][5][ie]=SWK[Material][1][4][ie]
                           +aco[0]*sef*dos1*(1.+2.*alphaK[Material][1]*finalenergy);
   }
   else {
    // NO ACOUSTIC PHONON
    SWK[Material][1][5][ie]=SWK[Material][1][4][ie]+0.0;
   }
// Impurity scattering
   if(IMPURITYPHONONS==ON){
    finalenergy=initialenergy;
    sef=sqrt(finalenergy*(1.+alphaK[Material][1]*finalenergy));
    ak=SMH[Material][1]*sef;
    qq=QD2*(4.*ak*ak+QD2);
    wk=bimp/qq*sef*dos1*(1.+2.*alphaK[Material][1]*finalenergy);
    if(wk>1.e14) wk=1.e14;
    SWK[Material][1][6][ie]=SWK[Material][1][5][ie]+wk;
   }
   else {
    // NO IMPURITY SCATTERING
    SWK[Material][1][6][ie]=SWK[Material][1][5][ie]+0.0;
   }
// L-valley
// ========
   if(OPTICALPHONONS==ON){
// Polar optical phonon
// Emission
    finalenergy=initialenergy-HWO[Material][0];
    if(finalenergy>0.){
      sef=sqrt(finalenergy);
      qmax=sef+sei;
      qmin=sei-sef;
      SWK[Material][2][1][ie]=poe*SMH[Material][2]*sei/initialenergy/Q*log(qmax/qmin);
    }
    else SWK[Material][2][1][ie]=0.;
// Absorption
    finalenergy=initialenergy+HWO[Material][0];
    sef=sqrt(finalenergy);
    qmax=sef+sei;
    qmin=sef-sei;
    SWK[Material][2][2][ie]=SWK[Material][2][1][ie]
                 +poa*SMH[Material][2]*sei/initialenergy/Q*log(qmax/qmin);
// Non-polar optical phonon
// Emission
    finalenergy=initialenergy-hwe;
    if(finalenergy>0.){
      sef=sqrt(finalenergy*(1.+alphaK[Material][2]*finalenergy));
      SWK[Material][2][3][ie]=SWK[Material][2][2][ie]
                   +(z2-1.)*eqe*sef*dos2*(1.+2.*alphaK[Material][2]*finalenergy);
    }
    else SWK[Material][2][3][ie]=SWK[Material][2][2][ie];
// Absorption
    finalenergy=initialenergy+hwe;
    sef=sqrt(finalenergy*(1.+alphaK[Material][2]*finalenergy));
    SWK[Material][2][4][ie]=SWK[Material][2][3][ie]
                 +(z2-1.)*eqa*sef*dos2*(1.+2.*alphaK[Material][2]*finalenergy);

// Emission
    finalenergy=initialenergy-hwij+EMIN[Material][2]-EMIN[Material][1];
    if(finalenergy>0.){
     sef=sqrt(finalenergy*(1.+alphaK[Material][2]*finalenergy));
     SWK[Material][2][5][ie]=SWK[Material][2][4][ie]
                            +ope*sef*dos1*(1.+2.*alphaK[Material][2]*finalenergy);
    }
    else SWK[Material][2][5][ie]=SWK[Material][2][4][ie];
// Absorption
    finalenergy=initialenergy+hwij+EMIN[Material][2]-EMIN[Material][1];
    if(finalenergy>0.){
     sef=sqrt(finalenergy*(1.+alphaK[Material][2]*finalenergy));
     SWK[Material][2][6][ie]=SWK[Material][2][5][ie]
                            +opa*sef*dos1*(1.+2.*alphaK[Material][2]*finalenergy);
    }
    else SWK[Material][2][6][ie]=SWK[Material][2][5][ie];
   }
   else {
      // NO OPTICAL SCATTERING
      SWK[Material][2][1][ie]=0.0;
      SWK[Material][2][2][ie]=0.0;
      SWK[Material][2][3][ie]=0.0;
      SWK[Material][2][4][ie]=0.0;
      SWK[Material][2][5][ie]=0.0;
      SWK[Material][2][6][ie]=0.0;
   }
   if(ACOUSTICPHONONS==ON){
// Acoustic phonon
    finalenergy=initialenergy;
    sef=sqrt(finalenergy*(1.+alphaK[Material][2]*finalenergy));
    SWK[Material][2][7][ie]=SWK[Material][2][6][ie]
                           +aco[0]*sef*dos2*(1.+2.*alphaK[Material][2]*finalenergy);
   }
   else {
      // NO ACOUSTIC SCATTERING
      SWK[Material][2][7][ie]=SWK[Material][2][6][ie]+0.0;
   }
// ==================
// --->  } <--- This is a bug! Removed by J.M.Sellier on 24 dec.2006 (SR)
// (Thanks to Kun-Yuan Xu who helped me to find it!)
// ==================
   if(IMPURITYPHONONS==ON){
// Impurity scattering
    finalenergy=initialenergy;
    sef=sqrt(finalenergy*(1.+alphaK[Material][2]*finalenergy));
    ak=SMH[Material][2]*sef;
    qq=QD2*(4.*ak*ak+QD2);
    wk=bimp/qq*sef*dos2*(1.+2.*alphaK[Material][2]*finalenergy);
    if(wk>1.e14) wk=1.e14;
    SWK[Material][2][8][ie]=SWK[Material][2][7][ie]+wk;
   }
   else {
    // NO IMPURITY SCATTERING
    SWK[Material][2][8][ie]=SWK[Material][2][7][ie]+0.0;
   }
  } // <--- this is the correct place for the symbol } (J.M.Sellier,24/12/2006)
// Evalutation of gamma
  GM[Material]=SWK[Material][1][6][1];
  for(ie=1;ie<=DIME;ie++){
    if(SWK[Material][1][6][ie]>GM[Material]) GM[Material]=SWK[Material][1][6][ie];
    if(SWK[Material][2][8][ie]>GM[Material]) GM[Material]=SWK[Material][2][8][ie];
  }
  printf("GAMMA[%d] = %g \n",Material, GM[Material]);
  for(i=1;i<=6;i++)
    for(ie=1;ie<=DIME;ie++)
      SWK[Material][1][i][ie]/=GM[Material];
  for(i=1;i<=8;i++)
    for(ie=1;ie<=DIME;ie++)
      SWK[Material][2][i][ie]/=GM[Material];
 }
// End of two-valleys materials
// ############################

// Material = one valley
// #####################
 if(NOVALLEY[Material]==1){
  SMH[Material][0]=sqrt(2.*MSTAR[Material][1][0]*M*Q)/HBAR;
  HHM[Material][0]=HBAR*HBAR/(2.*MSTAR[Material][1][0]*M*Q);
  HM[Material][0]=HBAR/(MSTAR[Material][1][0]*M);
// Density of states
  dos=pow((sqrt(2.*MSTAR[Material][1][0]*M)*sqrt(Q)/HBAR),3.)/(4.*PI*PI);
// constant for the acoustic phonon
  aco[0]=2.*PI*(DA[Material][0]/Q)*DA[Material][0]*(BKTQ/HBAR)
     *(Q/(RHO[Material]*UL[Material]*UL[Material]));
// Constants for the 6 Silicon-like non-polar optical phonons
   for(i=1;i<=6;i++){
// i-th Optical Phonon
    oge[i]=0.;
    oga[i]=0.;
    if(ZF[Material][i-1]!=0.){
      wo=HWO[Material][i-1]*Q/HBAR; // frequency of phonon
      no=1./(exp(HWO[Material][i-1]/BKTQ) - 1.); // population of phonons
      oge[i]=ZF[Material][i-1]*PI*(DTK[Material][i-1]*Q/wo)
            *((DTK[Material][i-1]*Q/RHO[Material])/Q)*(no+1.);
      oga[i]=oge[i]*no/(1.+no);
    }
   }
// Calculation of scattering rates
   for(ie=1; ie<=DIME; ++ie) SWK[Material][0][0][ie]=0.;
   for(ie=1; ie<=DIME; ++ie){
    initialenergy=DE*((real) ie);
    sei=sqrt(initialenergy);
    if(OPTICALPHONONS==ON){
// non polar optical phonons
     for(i=1;i<=6;i++){
      finalenergy=initialenergy-HWO[Material][i-1];
      SWK[Material][0][i*2-1][ie]=SWK[Material][0][i*2-2][ie];
      if(finalenergy>0.){
       sef=sqrt(finalenergy*(1.+alphaK[Material][1]*finalenergy));
       SWK[Material][0][i*2-1][ie]=SWK[Material][0][i*2-2][ie]
                     +oge[i]*sef*dos*(1.+2.*alphaK[Material][1]*finalenergy);
      }
      finalenergy=initialenergy+HWO[Material][i-1];
      SWK[Material][0][i*2][ie]=SWK[Material][0][i*2-1][ie];
      if(finalenergy>0.){
       sef=sqrt(finalenergy*(1.+alphaK[Material][1]*finalenergy));
       SWK[Material][0][i*2][ie]=SWK[Material][0][i*2-1][ie]
                   +oga[i]*sef*dos*(1.+2.*alphaK[Material][1]*finalenergy);
      }
     }
    }
    else {
     // NO OPTICAL SCATTERING
     SWK[Material][0][1][ie]=0.0;
     SWK[Material][0][2][ie]=0.0;
     SWK[Material][0][3][ie]=0.0;
     SWK[Material][0][4][ie]=0.0;
     SWK[Material][0][5][ie]=0.0;
     SWK[Material][0][6][ie]=0.0;
     SWK[Material][0][7][ie]=0.0;
     SWK[Material][0][8][ie]=0.0;
     SWK[Material][0][9][ie]=0.0;
     SWK[Material][0][10][ie]=0.0;
     SWK[Material][0][11][ie]=0.0;
     SWK[Material][0][12][ie]=0.0;
    }
    if(ACOUSTICPHONONS==ON){
// Acoustic phonon
    finalenergy=initialenergy;
    sef=sqrt(finalenergy*(1.+alphaK[Material][1]*finalenergy));
    SWK[Material][0][13][ie]=SWK[Material][0][12][ie]
                  +aco[0]*sef*dos*(1.+2.*alphaK[Material][1]*finalenergy);
    }
    else {
    // NO ACOUSTIC PHONONS 
    SWK[Material][0][13][ie]=SWK[Material][0][12][ie]+0.0;
    }
   }
// Evaluation of gamma
  GM[Material]=SWK[Material][0][13][1];
  for(ie=1;ie<=DIME;++ie)
    if(SWK[Material][0][13][ie]>GM[Material]) GM[Material]=SWK[Material][0][13][ie];
  printf("GAMMA[%d] = %g \n",Material, GM[Material]);
  for(ie=1;ie<=DIME;ie++)
    for(i=1;i<=13;i++)
      SWK[Material][0][i][ie]/=GM[Material];
 }
// End of one-valley material
// Germanium-simulation
// #####################
 if(NOVALLEY[Material]==3){
    // add tensor of Herring - Vogt Transformation by Suzuki
    //########################################################
    am1_dos=MSTAR[Material][1][1]*M;
    am2_dos=MSTAR[Material][2][1]*M;
    am3_dos=MSTAR[Material][3][1]*M;
    dos1=pow(sqrt(2.*am1_dos*Q)/HBAR,3.)/pow(2.*PI,2.);
    dos2=pow(sqrt(2.*am2_dos*Q)/HBAR,3.)/pow(2.*PI,2.);
    dos3=pow(sqrt(2.*am3_dos*Q)/HBAR,3.)/pow(2.*PI,2.);
    eps=EPSR[Material]*EPS0;
    epf=EPF[Material]*EPS0;
    ep=1./(1./epf-1./eps);

    // Parameters for Phonon Scattering
    cl=RHO[Material]*pow(UL[Material],2.);
    deq=dij=DTK[Material][0]*Q;
    hwe=hwij=HWO[Material][0];
    SMH[Material][1]=sqrt(2.*am1_dos*Q)/HBAR;
    SMH[Material][2]=sqrt(2.*am2_dos*Q)/HBAR;
    SMH[Material][3]=sqrt(2.*am3_dos*Q)/HBAR;
    HHM[Material][1]=HBAR*HBAR/(2.*am1_dos*Q);
    HHM[Material][2]=HBAR*HBAR/(2.*am2_dos*Q);
    HHM[Material][3]=HBAR*HBAR/(2.*am3_dos*Q);
    HM[Material][1]=HBAR/am1_dos;
    HM[Material][2]=HBAR/am2_dos;
    HM[Material][3]=HBAR/am3_dos;

    // Parameters for Optical Deformation Potential Scatterng
    dde=5.5e10*Q; //
    wo=HWO[Material][0]*Q/HBAR;
    no=1./(exp(HWO[Material][0]/BKTQ)-1.);
    ode=PI*dde/wo*dde/RHO[Material]/Q*(no+1.);
    oda=ode*no/(1.+no);

    // Parameters for Inter Valley Scattering
    real zl = 4.0;
    real zd = 2.0;
    for(j=0; j<=6; j++){
        w[j]=HIV[j]*Q/HBAR;
        n[j]=1./(exp(HIV[j]/BKTQ)-1.);
    }
    // Parameters for Acoustic Phonon Scattering
    aco[1]=2.*PI*DA[Material][0]/Q*DA[Material][0]*BKTQ/HBAR*Q/cl;
    aco[2]=2.*PI*DA[Material][1]/Q*DA[Material][1]*BKTQ/HBAR*Q/cl;
    aco[3]=2.*PI*DA[Material][2]/Q*DA[Material][2]*BKTQ/HBAR*Q/cl;

    // Parameters for impurity scatterings
    cimp=undope; // <--- impurity concentration (undope)
    qd=sqrt(Q*cimp/BKTQ/eps);
    QD2=qd*qd;  
    bimp=2.*PI*cimp*Q*Q/HBAR*Q/eps/eps;

    real temp=0.;
    // ###########################################

    for(ie=1; ie<=DIME; ie++){
        initialenergy=DE*((real)(ie));
        sei=sqrt(initialenergy);
        // L-Valley
        // =======================================
        if(OPTICALPHONONS==ON){
            // Optical Deformation Potential Scattering 

            // Emission
            finalenergy=initialenergy-HWO[Material][0];
            if(finalenergy>0.){
                sef=sqrt(finalenergy*(1.+alphaK[Material][1]*finalenergy));
                SWK[Material][1][1][ie]=ode*sef*dos1*(1.+2.*alphaK[Material][1]*finalenergy);
            }
            else SWK[Material][1][1][ie]=0.;

            // Absorption
            finalenergy=initialenergy+HWO[Material][0];
            sef=sqrt(finalenergy*(1.+alphaK[Material][1]*finalenergy));
            SWK[Material][1][2][ie]=SWK[Material][1][1][ie]+oda*sef*dos1*(1.+2.*alphaK[Material][1]*finalenergy);

            // Intervalley phonon Scattering (iv)

            // L -> L (1)
            ope=PI*DIV[0]/w[0]*DIV[0]/RHO[Material]/Q*(n[0]+1.);
            opa=ope*n[0]/(1.+n[0]);
            // Emission
            finalenergy=initialenergy-HIV[0]+EMIN[Material][1]-EMIN[Material][1];
            if(finalenergy>0.){
                sef=sqrt(finalenergy*(1.+alphaK[Material][1]*finalenergy));
                SWK[Material][1][3][ie]=SWK[Material][1][2][ie]+(zl-1.)*ope*sef*dos1*(1.+2.*alphaK[Material][1]*finalenergy);
            }
            else SWK[Material][1][3][ie]=SWK[Material][1][2][ie];

            // Absorption
            finalenergy=initialenergy+HIV[0]+EMIN[Material][1]-EMIN[Material][1];
            if(finalenergy>0.){
                sef=sqrt(finalenergy*(1.+alphaK[Material][1]*finalenergy));
                SWK[Material][1][4][ie]=SWK[Material][1][3][ie]+(zl-1.)*opa*sef*dos1*(1.+2.*alphaK[Material][1]*finalenergy);
            }
            else SWK[Material][1][4][ie]=SWK[Material][1][3][ie];

            // L -> L (2)
            ope=PI*DIV[1]/w[1]*DIV[1]/RHO[Material]/Q*(n[1]+1.);
            opa=ope*n[1]/(1.+n[1]);
            // Emission
            finalenergy=initialenergy-HIV[1]+EMIN[Material][1]-EMIN[Material][1];
            if(finalenergy>0.){
                sef=sqrt(finalenergy*(1.+alphaK[Material][1]*finalenergy));
                SWK[Material][1][5][ie]=SWK[Material][1][4][ie]+(zl-1.)*ope*sef*dos1*(1.+2.*alphaK[Material][1]*finalenergy);
            }
            else SWK[Material][1][5][ie]=SWK[Material][1][4][ie];

            // Absorption
            finalenergy=initialenergy+HIV[1]+EMIN[Material][1]-EMIN[Material][1];
            if(finalenergy>0.){
                sef=sqrt(finalenergy*(1.+alphaK[Material][1]*finalenergy));
                SWK[Material][1][6][ie]=SWK[Material][1][5][ie]+(zl-1.)*opa*sef*dos1*(1.+2.*alphaK[Material][1]*finalenergy);
            }
            else SWK[Material][1][6][ie]=SWK[Material][1][5][ie];

            // L -> G (3)
            ope=PI*DIV[2]/w[2]*DIV[2]/RHO[Material]/Q*(n[2]+1.);
            opa=ope*n[2]/(1.+n[2]);
            // Emission
            finalenergy=initialenergy-HIV[2]+EMIN[Material][1]-EMIN[Material][2];
            if(finalenergy>0.){
                sef=sqrt(finalenergy*(1.+alphaK[Material][2]*finalenergy));
                SWK[Material][1][7][ie]=SWK[Material][1][6][ie]+ope*sef*dos2*(1.+2.*alphaK[Material][2]*finalenergy);
            }
            else SWK[Material][1][7][ie]=SWK[Material][1][6][ie];

            // Absorption
            finalenergy=initialenergy+HIV[2]+EMIN[Material][1]-EMIN[Material][2];
            if(finalenergy>0.){
                sef=sqrt(finalenergy*(1.+alphaK[Material][2]*finalenergy));
                SWK[Material][1][8][ie]=SWK[Material][1][7][ie]+opa*sef*dos2*(1.+2.*alphaK[Material][2]*finalenergy);
            }
            else SWK[Material][1][8][ie]=SWK[Material][1][7][ie];

            // L -> D (4)
            ope=PI*DIV[3]/w[3]*DIV[3]/RHO[Material]/Q*(n[3]+1.);
            opa=ope*n[3]/(1.+n[3]);
            // Emission
            finalenergy=initialenergy-HIV[3]+EMIN[Material][1]-EMIN[Material][3];
            if(finalenergy>0.){
                sef=sqrt(finalenergy*(1.+alphaK[Material][3]*finalenergy));
                SWK[Material][1][9][ie]=SWK[Material][1][8][ie]+zd*ope*sef*dos3*(1.+2.*alphaK[Material][3]*finalenergy);
            }
            else SWK[Material][1][9][ie]=SWK[Material][1][8][ie];

            // Absorption
            finalenergy=initialenergy+HIV[3]+EMIN[Material][1]-EMIN[Material][3];
            if(finalenergy>0.){
                sef=sqrt(finalenergy*(1.+alphaK[Material][3]*finalenergy));
                SWK[Material][1][10][ie]=SWK[Material][1][9][ie]+zd*opa*sef*dos3*(1.+2.*alphaK[Material][3]*finalenergy);
            }
            else SWK[Material][1][10][ie]=SWK[Material][1][9][ie];
        }
        else {
            // NO OPTICAL PHONON
            SWK[Material][1][1][ie]=0.0;
            SWK[Material][1][2][ie]=0.0;
            SWK[Material][1][3][ie]=0.0;
            SWK[Material][1][4][ie]=0.0;
            SWK[Material][1][5][ie]=0.0;
            SWK[Material][1][6][ie]=0.0;
            SWK[Material][1][7][ie]=0.0;
            SWK[Material][1][8][ie]=0.0;
            SWK[Material][1][9][ie]=0.0;
            SWK[Material][1][10][ie]=0.0;
        }

        if(ACOUSTICPHONONS==ON){
            // Acoustic Phonon Scattering
            finalenergy=initialenergy;
            sef=sqrt(finalenergy*(1.+alphaK[Material][1]*finalenergy));
            SWK[Material][1][11][ie]=SWK[Material][1][10][ie]
                                    +aco[1]*sef*dos1*(1.+2.*alphaK[Material][1]*finalenergy);
            }
        else {
            // NO ACOUSTIC PHONON
            SWK[Material][1][11][ie]=SWK[Material][1][10][ie]+0.0;
        }

        // Impurity scattering
        if(IMPURITYPHONONS==ON){
            finalenergy=initialenergy;
            sef=sqrt(finalenergy*(1.+alphaK[Material][1]*finalenergy));
            ak=SMH[Material][1]*sef;
            qq=QD2*(4.*ak*ak+QD2);
            wk=bimp/qq*sef*dos1*(1.+2.*alphaK[Material][1]*finalenergy);
            if(wk>1.e14) wk=1.e14;
            SWK[Material][1][12][ie]=SWK[Material][1][11][ie]+wk;
        }
        else {
            // NO IMPURITY SCATTERING
            SWK[Material][1][12][ie]=SWK[Material][1][11][ie]+0.0;
        }
        // Scattering rate calculation before normalization
        Scat_Rate[Material][1][1][ie]=SWK[Material][1][1][ie];
        for(i=2;i<=12;i++){
            Scat_Rate[Material][1][i][ie]=SWK[Material][1][i][ie]-SWK[Material][1][i-1][ie];
        }

        // G-Valley
        // =======================================
        if(OPTICALPHONONS==ON){
        // Intervalley phonon Scattering (iv)

            // G -> L (3)
            ope=PI*DIV[2]/w[2]*DIV[2]/RHO[Material]/Q*(n[2]+1.);
            opa=ope*n[2]/(1.+n[2]);
            // Emission
            finalenergy=initialenergy-HIV[2]+EMIN[Material][2]-EMIN[Material][1];
            if(finalenergy>0.){
                sef=sqrt(finalenergy*(1.+alphaK[Material][1]*finalenergy));
                SWK[Material][2][1][ie]=zl*ope*sef*dos1*(1.+2.*alphaK[Material][1]*finalenergy);
            }
            else SWK[Material][2][1][ie]=0.;

            // Absorption
            finalenergy=initialenergy+HIV[2]+EMIN[Material][2]-EMIN[Material][1];
            if(finalenergy>0.){
                sef=sqrt(finalenergy*(1.+alphaK[Material][1]*finalenergy));
                SWK[Material][2][2][ie]=SWK[Material][2][1][ie]+zl*opa*sef*dos1*(1.+2.*alphaK[Material][1]*finalenergy);
            }
            else SWK[Material][2][2][ie]=SWK[Material][2][1][ie];

            // G -> D (7)
            ope=PI*DIV[6]/w[6]*DIV[6]/RHO[Material]/Q*(n[6]+1.);
            opa=ope*n[6]/(1.+n[6]);
            // Emission
            finalenergy=initialenergy-HIV[6]+EMIN[Material][2]-EMIN[Material][3];
            if(finalenergy>0.){
                sef=sqrt(finalenergy*(1.+alphaK[Material][3]*finalenergy));
                SWK[Material][2][3][ie]=SWK[Material][2][2][ie]+zd*ope*sef*dos3*(1.+2.*alphaK[Material][3]*finalenergy);
            }
            else SWK[Material][2][3][ie]=SWK[Material][2][2][ie];

            // Absorption
            finalenergy=initialenergy+HIV[6]+EMIN[Material][2]-EMIN[Material][3];
            if(finalenergy>0.){
                sef=sqrt(finalenergy*(1.+alphaK[Material][3]*finalenergy));
                SWK[Material][2][4][ie]=SWK[Material][2][3][ie]+zd*opa*sef*dos3*(1.+2.*alphaK[Material][3]*finalenergy);
            }
            else SWK[Material][2][4][ie]=SWK[Material][2][3][ie];
        }
        else {
            // NO OPTICAL PHONON
            SWK[Material][2][1][ie]=0.0;
            SWK[Material][2][2][ie]=0.0;
            SWK[Material][2][3][ie]=0.0;
            SWK[Material][2][4][ie]=0.0;
        }

        if(ACOUSTICPHONONS==ON){
            // Acoustic Phonon Scattering
            finalenergy=initialenergy;
            sef=sqrt(finalenergy*(1.+alphaK[Material][2]*finalenergy));
            SWK[Material][2][5][ie]=SWK[Material][2][4][ie]
                                    +aco[2]*sef*dos2*(1.+2.*alphaK[Material][2]*finalenergy);
            }
        else {
            // NO ACOUSTIC PHONON
            SWK[Material][2][5][ie]=SWK[Material][2][4][ie]+0.0;
        }

        if(IMPURITYPHONONS==ON){
            // Impurity scattering
            finalenergy=initialenergy;
            sef=sqrt(finalenergy*(1.+alphaK[Material][2]*finalenergy));
            ak=SMH[Material][2]*sef;
            qq=QD2*(4.*ak*ak+QD2);
            wk=bimp/qq*sef*dos2*(1.+2.*alphaK[Material][2]*finalenergy);
            if(wk>1.e14) wk=1.e14;
            SWK[Material][2][6][ie]=SWK[Material][2][5][ie]+wk;
        }
        else {
            // NO IMPURITY SCATTERING
            SWK[Material][2][6][ie]=SWK[Material][2][5][ie]+0.0;
        }
        // Scattering rate calculation before normalization
        Scat_Rate[Material][2][1][ie]=SWK[Material][2][1][ie];
        for(i=2;i<=6;i++){
            Scat_Rate[Material][2][i][ie]=SWK[Material][2][i][ie]-SWK[Material][2][i-1][ie];
        }

        // Delta-Valley
        // =======================================

        if(OPTICALPHONONS==ON){
        // Intervalley phonon Scattering (iv)

            // D -> L (4)
            ope=PI*DIV[3]/w[3]*DIV[3]/RHO[Material]/Q*(n[3]+1.);
            opa=ope*n[3]/(1.+n[3]);
            // Emission
            finalenergy=initialenergy-HIV[3]+EMIN[Material][3]-EMIN[Material][1];
            if(finalenergy>0.){
                sef=sqrt(finalenergy*(1.+alphaK[Material][1]*finalenergy));
                SWK[Material][3][1][ie]=zl*ope*sef*dos1*(1.+2.*alphaK[Material][1]*finalenergy);
            }
            else SWK[Material][3][1][ie]=0.;

            // Absorption
            finalenergy=initialenergy+HIV[3]+EMIN[Material][3]-EMIN[Material][1];
            if(finalenergy>0.){
                sef=sqrt(finalenergy*(1.+alphaK[Material][1]*finalenergy));
                SWK[Material][3][2][ie]=SWK[Material][3][1][ie]+zl*opa*sef*dos1*(1.+2.*alphaK[Material][1]*finalenergy);
            }
            else SWK[Material][3][2][ie]=SWK[Material][3][1][ie];

            // D -> G (7)
            ope=PI*DIV[6]/w[6]*DIV[6]/RHO[Material]/Q*(n[6]+1.);
            opa=ope*n[6]/(1.+n[6]);
            // Emission
            finalenergy=initialenergy-HIV[6]+EMIN[Material][3]-EMIN[Material][2];
            if(finalenergy>0.){
                sef=sqrt(finalenergy*(1.+alphaK[Material][2]*finalenergy));
                SWK[Material][3][3][ie]=SWK[Material][3][2][ie]+ope*sef*dos2*(1.+2.*alphaK[Material][2]*finalenergy);
            }
            else SWK[Material][3][3][ie]=SWK[Material][3][2][ie];

            // Absorption
            finalenergy=initialenergy+HIV[6]+EMIN[Material][3]-EMIN[Material][2];
            if(finalenergy>0.){
                sef=sqrt(finalenergy*(1.+alphaK[Material][2]*finalenergy));
                SWK[Material][3][4][ie]=SWK[Material][3][3][ie]+opa*sef*dos2*(1.+2.*alphaK[Material][2]*finalenergy);
            }
            else SWK[Material][3][4][ie]=SWK[Material][3][3][ie];

            // D -> D (5)
            ope=PI*DIV[4]/w[4]*DIV[4]/RHO[Material]/Q*(n[4]+1.);
            opa=ope*n[4]/(1.+n[4]);
            // Emission
            finalenergy=initialenergy-HIV[4]+EMIN[Material][3]-EMIN[Material][3];
            if(finalenergy>0.){
                sef=sqrt(finalenergy*(1.+alphaK[Material][3]*finalenergy));
                SWK[Material][3][5][ie]=SWK[Material][3][4][ie]+(zd-1.)*ope*sef*dos3*(1.+2.*alphaK[Material][3]*finalenergy);
            }
            else SWK[Material][3][5][ie]=SWK[Material][3][4][ie];

            // Absorption
            finalenergy=initialenergy+HIV[4]+EMIN[Material][3]-EMIN[Material][3];
            if(finalenergy>0.){
                sef=sqrt(finalenergy*(1.+alphaK[Material][3]*finalenergy));
                SWK[Material][3][6][ie]=SWK[Material][3][5][ie]+(zd-1.)*opa*sef*dos3*(1.+2.*alphaK[Material][3]*finalenergy);
            }
            else SWK[Material][3][6][ie]=SWK[Material][3][5][ie];

            // D -> D (6)
            ope=PI*DIV[5]/w[5]*DIV[5]/RHO[Material]/Q*(n[5]+1.);
            opa=ope*n[5]/(1.+n[5]);
            // Emission
            finalenergy=initialenergy-HIV[5]+EMIN[Material][3]-EMIN[Material][3];
            if(finalenergy>0.){
                sef=sqrt(finalenergy*(1.+alphaK[Material][3]*finalenergy));
                SWK[Material][3][7][ie]=SWK[Material][3][6][ie]+(zd-1.)*ope*sef*dos3*(1.+2.*alphaK[Material][3]*finalenergy);
            }
            else SWK[Material][3][7][ie]=SWK[Material][3][6][ie];

            // Absorption
            finalenergy=initialenergy+HIV[5]+EMIN[Material][3]-EMIN[Material][3];
            if(finalenergy>0.){
                sef=sqrt(finalenergy*(1.+alphaK[Material][3]*finalenergy));
                SWK[Material][3][8][ie]=SWK[Material][3][7][ie]+(zd-1.)*opa*sef*dos3*(1.+2.*alphaK[Material][3]*finalenergy);
            }
            else SWK[Material][3][8][ie]=SWK[Material][3][7][ie];

        }
        else {
            // NO OPTICAL PHONON
            SWK[Material][1][1][ie]=0.0;
            SWK[Material][1][2][ie]=0.0;
            SWK[Material][1][3][ie]=0.0;
            SWK[Material][1][4][ie]=0.0;
            SWK[Material][1][5][ie]=0.0;
            SWK[Material][1][6][ie]=0.0;
            SWK[Material][1][7][ie]=0.0;
            SWK[Material][1][8][ie]=0.0;
        }

        if(ACOUSTICPHONONS==ON){
            // Acoustic Phonon Scattering
            finalenergy=initialenergy;
            sef=sqrt(finalenergy*(1.+alphaK[Material][3]*finalenergy));
            SWK[Material][3][9][ie]=SWK[Material][3][8][ie]
                                +aco[3]*sef*dos3*(1.+2.*alphaK[Material][3]*finalenergy);
        }
        else {
            // NO ACOUSTIC PHONON
            SWK[Material][3][9][ie]=SWK[Material][3][8][ie]+0.0;
        }

        if(IMPURITYPHONONS==ON){
            // Impurity scattering
            finalenergy=initialenergy;
            sef=sqrt(finalenergy*(1.+alphaK[Material][3]*finalenergy));
            ak=SMH[Material][3]*sef;
            qq=QD2*(4.*ak*ak+QD2);
            wk=bimp/qq*sef*dos3*(1.+2.*alphaK[Material][3]*finalenergy);
            if(wk>1.e14) wk=1.e14;
            SWK[Material][3][10][ie]=SWK[Material][3][9][ie]+wk;
        }
        else {
            // NO IMPURITY SCATTERING
            SWK[Material][3][10][ie]=SWK[Material][3][9][ie]+0.0;
        }
        // Scattering rate calculation before normalization
        Scat_Rate[Material][3][1][ie]=SWK[Material][3][1][ie];
        for(i=2;i<=10;i++){
            Scat_Rate[Material][3][i][ie]=SWK[Material][3][i][ie]-SWK[Material][3][i-1][ie];
        }
    }
    // Evaluation scattering rate
    gm[Material][1][0]=SWK[Material][1][12][1];
    for(ie=1;ie<=DIME;ie++){
        if(SWK[Material][1][12][ie]>gm[Material][1][0]) {
            gm[Material][1][0]=SWK[Material][1][12][ie];
        }
    }   
    printf("L[%d]undope = %g \n",Material, gm[Material][1][0]);
    for(i=1;i<=12;i++){
        for(ie=1;ie<=DIME;ie++){
            SWK[Material][1][i][ie]/=gm[Material][1][0];
            if(SWK[Material][1][i][ie] > 1.0){
                SWK[Material][1][i][ie] = 1.0;
            }
        }
    }

    gm[Material][2][0]=SWK[Material][2][6][1];
    for(ie=1;ie<=DIME;ie++){
        if(SWK[Material][2][6][ie]>gm[Material][2][0]){
            gm[Material][2][0]=SWK[Material][2][6][ie];
        }
    }
    printf("GAMMA[%d]undope = %g \n",Material, gm[Material][2][0]);
    for(i=1;i<=6;i++){
        for(ie=1;ie<=DIME;ie++){
            SWK[Material][2][i][ie]/=gm[Material][2][0];
            if(SWK[Material][2][i][ie] > 1.0){
                SWK[Material][2][i][ie] = 1.0;
            }
        }
    }

    gm[Material][3][0]=SWK[Material][3][10][1];
    for(ie=1;ie<=DIME;ie++){
        if(SWK[Material][3][10][ie]>gm[Material][3][0]){
            gm[Material][3][0]=SWK[Material][3][10][ie];
        }
    }
    printf("DELTA[%d]undope = %g \n",Material, gm[Material][3][0]);
    for(i=1;i<=10;i++){
        for(ie=1;ie<=DIME;ie++){
            SWK[Material][3][i][ie]/=gm[Material][3][0];
            if(SWK[Material][3][i][ie] > 1.0){
                SWK[Material][3][i][ie] = 1.0;
            }
        }
    }
    GM[Material]=gm[Material][3][0];
    // End of undope layer
    // ###################
  
    // Start of N+ layer
    am1_dos=MSTAR[Material][1][1]*M;
    am2_dos=MSTAR[Material][2][1]*M;
    am3_dos=MSTAR[Material][3][1]*M;
    dos1=pow(sqrt(2.*am1_dos*Q)/HBAR,3.)/pow(2.*PI,2.);
    dos2=pow(sqrt(2.*am2_dos*Q)/HBAR,3.)/pow(2.*PI,2.);
    dos3=pow(sqrt(2.*am3_dos*Q)/HBAR,3.)/pow(2.*PI,2.);
    eps=EPSR[Material]*EPS0;
    epf=EPF[Material]*EPS0;
    ep=1./(1./epf-1./eps);

    // Parameters for Phonon Scattering
    cl=RHO[Material]*pow(UL[Material],2.);
    deq=dij=DTK[Material][0]*Q;
    hwe=hwij=HWO[Material][0];
    SMH[Material][1]=sqrt(2.*am1_dos*Q)/HBAR;
    SMH[Material][2]=sqrt(2.*am2_dos*Q)/HBAR;
    SMH[Material][3]=sqrt(2.*am3_dos*Q)/HBAR;
    HHM[Material][1]=HBAR*HBAR/(2.*am1_dos*Q);
    HHM[Material][2]=HBAR*HBAR/(2.*am2_dos*Q);
    HHM[Material][3]=HBAR*HBAR/(2.*am3_dos*Q);
    HM[Material][1]=HBAR/am1_dos;
    HM[Material][2]=HBAR/am2_dos;
    HM[Material][3]=HBAR/am3_dos;

    // Parameters for Optical Deformation Potential Scatterng
    dde=5.5e10*Q; //
    wo=HWO[Material][0]*Q/HBAR;
    no=1./(exp(HWO[Material][0]/BKTQ)-1.);
    ode=PI*dde/wo*dde/RHO[Material]/Q*(no+1.);
    oda=ode*no/(1.+no);

    // Parameters for Inter Valley Scattering
    for(j=0; j<=6; j++){
        w[j]=HIV[j]*Q/HBAR;
        n[j]=1./(exp(HIV[j]/BKTQ)-1.);
    }
    // Parameters for Acoustic Phonon Scattering
    aco[1]=2.*PI*DA[Material][0]/Q*DA[Material][0]*BKTQ/HBAR*Q/cl;
    aco[2]=2.*PI*DA[Material][1]/Q*DA[Material][1]*BKTQ/HBAR*Q/cl;
    aco[3]=2.*PI*DA[Material][2]/Q*DA[Material][2]*BKTQ/HBAR*Q/cl;

    // Parameters for impurity scatterings
    cimp=Nplus; // <--- impurity concentration (N+)
    qd=sqrt(Q*cimp/BKTQ/eps);
    QD2=qd*qd;  
    bimp=2.*PI*cimp*Q*Q/HBAR*Q/eps/eps;

    // ###########################################

    for(ie=1; ie<=DIME; ie++){
        initialenergy=DE*((real)(ie));
        sei=sqrt(initialenergy);
        // L-Valley
        // =======================================
        if(OPTICALPHONONS==ON){
            // Optical Deformation Potential Scattering 

            // Emission
            finalenergy=initialenergy-HWO[Material][0];
            if(finalenergy>0.){
                sef=sqrt(finalenergy*(1.+alphaK[Material][1]*finalenergy));
                SWK_N[Material][1][1][ie]=ode*sef*dos1*(1.+2.*alphaK[Material][1]*finalenergy);
            }
            else SWK_N[Material][1][1][ie]=0.;

            // Absorption
            finalenergy=initialenergy+HWO[Material][0];
            sef=sqrt(finalenergy*(1.+alphaK[Material][1]*finalenergy));
            SWK_N[Material][1][2][ie]=SWK_N[Material][1][1][ie]+oda*sef*dos1*(1.+2.*alphaK[Material][1]*finalenergy);

            // Intervalley phonon Scattering (iv)

            // L -> L (1)
            ope=PI*DIV[0]/w[0]*DIV[0]/RHO[Material]/Q*(n[0]+1.);
            opa=ope*n[0]/(1.+n[0]);
            // Emission
            finalenergy=initialenergy-HIV[0]+EMIN[Material][1]-EMIN[Material][1];
            if(finalenergy>0.){
                sef=sqrt(finalenergy*(1.+alphaK[Material][1]*finalenergy));
                SWK_N[Material][1][3][ie]=SWK_N[Material][1][2][ie]+(zl-1.)*ope*sef*dos1*(1.+2.*alphaK[Material][1]*finalenergy);
            }
            else SWK_N[Material][1][3][ie]=SWK_N[Material][1][2][ie];

            // Absorption
            finalenergy=initialenergy+HIV[0]+EMIN[Material][1]-EMIN[Material][1];
            if(finalenergy>0.){
                sef=sqrt(finalenergy*(1.+alphaK[Material][1]*finalenergy));
                SWK_N[Material][1][4][ie]=SWK_N[Material][1][3][ie]+(zl-1.)*opa*sef*dos1*(1.+2.*alphaK[Material][1]*finalenergy);
            }
            else SWK_N[Material][1][4][ie]=SWK_N[Material][1][3][ie];

            // L -> L (2)
            ope=PI*DIV[1]/w[1]*DIV[1]/RHO[Material]/Q*(n[1]+1.);
            opa=ope*n[1]/(1.+n[1]);
            // Emission
            finalenergy=initialenergy-HIV[1]+EMIN[Material][1]-EMIN[Material][1];
            if(finalenergy>0.){
                sef=sqrt(finalenergy*(1.+alphaK[Material][1]*finalenergy));
                SWK_N[Material][1][5][ie]=SWK_N[Material][1][4][ie]+(zl-1.)*ope*sef*dos1*(1.+2.*alphaK[Material][1]*finalenergy);
            }
            else SWK_N[Material][1][5][ie]=SWK_N[Material][1][4][ie];

            // Absorption
            finalenergy=initialenergy+HIV[1]+EMIN[Material][1]-EMIN[Material][1];
            if(finalenergy>0.){
                sef=sqrt(finalenergy*(1.+alphaK[Material][1]*finalenergy));
                SWK_N[Material][1][6][ie]=SWK_N[Material][1][5][ie]+(zl-1.)*opa*sef*dos1*(1.+2.*alphaK[Material][1]*finalenergy);
            }
            else SWK_N[Material][1][6][ie]=SWK_N[Material][1][5][ie];

            // L -> G (3)
            ope=PI*DIV[2]/w[2]*DIV[2]/RHO[Material]/Q*(n[2]+1.);
            opa=ope*n[2]/(1.+n[2]);
            // Emission
            finalenergy=initialenergy-HIV[2]+EMIN[Material][1]-EMIN[Material][2];
            if(finalenergy>0.){
                sef=sqrt(finalenergy*(1.+alphaK[Material][2]*finalenergy));
                SWK_N[Material][1][7][ie]=SWK_N[Material][1][6][ie]+ope*sef*dos2*(1.+2.*alphaK[Material][2]*finalenergy);
            }
            else SWK_N[Material][1][7][ie]=SWK_N[Material][1][6][ie];

            // Absorption
            finalenergy=initialenergy+HIV[2]+EMIN[Material][1]-EMIN[Material][2];
            if(finalenergy>0.){
                sef=sqrt(finalenergy*(1.+alphaK[Material][2]*finalenergy));
                SWK_N[Material][1][8][ie]=SWK_N[Material][1][7][ie]+opa*sef*dos2*(1.+2.*alphaK[Material][2]*finalenergy);
            }
            else SWK_N[Material][1][8][ie]=SWK_N[Material][1][7][ie];

            // L -> D (4)
            ope=PI*DIV[3]/w[3]*DIV[3]/RHO[Material]/Q*(n[3]+1.);
            opa=ope*n[3]/(1.+n[3]);
            // Emission
            finalenergy=initialenergy-HIV[3]+EMIN[Material][1]-EMIN[Material][3];
            if(finalenergy>0.){
                sef=sqrt(finalenergy*(1.+alphaK[Material][3]*finalenergy));
                SWK_N[Material][1][9][ie]=SWK_N[Material][1][8][ie]+zd*ope*sef*dos3*(1.+2.*alphaK[Material][3]*finalenergy);
            }
            else SWK_N[Material][1][9][ie]=SWK_N[Material][1][8][ie];

            // Absorption
            finalenergy=initialenergy+HIV[3]+EMIN[Material][1]-EMIN[Material][3];
            if(finalenergy>0.){
                sef=sqrt(finalenergy*(1.+alphaK[Material][3]*finalenergy));
                SWK_N[Material][1][10][ie]=SWK_N[Material][1][9][ie]+zd*opa*sef*dos3*(1.+2.*alphaK[Material][3]*finalenergy);
            }
            else SWK_N[Material][1][10][ie]=SWK_N[Material][1][9][ie];
        }
        else {
            // NO OPTICAL PHONON
            SWK_N[Material][1][1][ie]=0.0;
            SWK_N[Material][1][2][ie]=0.0;
            SWK_N[Material][1][3][ie]=0.0;
            SWK_N[Material][1][4][ie]=0.0;
            SWK_N[Material][1][5][ie]=0.0;
            SWK_N[Material][1][6][ie]=0.0;
            SWK_N[Material][1][7][ie]=0.0;
            SWK_N[Material][1][8][ie]=0.0;
            SWK_N[Material][1][9][ie]=0.0;
            SWK_N[Material][1][10][ie]=0.0;
        }

        if(ACOUSTICPHONONS==ON){
            // Acoustic Phonon Scattering
            finalenergy=initialenergy;
            sef=sqrt(finalenergy*(1.+alphaK[Material][1]*finalenergy));
            SWK_N[Material][1][11][ie]=SWK_N[Material][1][10][ie]
                                    +aco[1]*sef*dos1*(1.+2.*alphaK[Material][1]*finalenergy);
            }
        else {
            // NO ACOUSTIC PHONON
            SWK_N[Material][1][11][ie]=SWK_N[Material][1][10][ie]+0.0;
        }

        // Impurity scattering
        if(IMPURITYPHONONS==ON){
            finalenergy=initialenergy;
            sef=sqrt(finalenergy*(1.+alphaK[Material][1]*finalenergy));
            ak=SMH[Material][1]*sef;
            qq=QD2*(4.*ak*ak+QD2);
            wk=bimp/qq*sef*dos1*(1.+2.*alphaK[Material][1]*finalenergy);
            if(wk>1.e14) wk=1.e14;
            SWK_N[Material][1][12][ie]=SWK_N[Material][1][11][ie]+wk;
        }
        else {
            // NO IMPURITY SCATTERING
            SWK_N[Material][1][12][ie]=SWK_N[Material][1][11][ie]+0.0;
        }
        // Scattering rate calculation before normalization
        Scat_Rate_N[Material][1][1][ie]=SWK_N[Material][1][1][ie];
        for(i=2;i<=12;i++){
            Scat_Rate_N[Material][1][i][ie]=SWK_N[Material][1][i][ie]-SWK_N[Material][1][i-1][ie];
        }

        // G-Valley
        // =======================================
        if(OPTICALPHONONS==ON){
        // Intervalley phonon Scattering (iv)

            // G -> L (3)
            ope=PI*DIV[2]/w[2]*DIV[2]/RHO[Material]/Q*(n[2]+1.);
            opa=ope*n[2]/(1.+n[2]);
            // Emission
            finalenergy=initialenergy-HIV[2]+EMIN[Material][2]-EMIN[Material][1];
            if(finalenergy>0.){
                sef=sqrt(finalenergy*(1.+alphaK[Material][1]*finalenergy));
                SWK_N[Material][2][1][ie]=zl*ope*sef*dos1*(1.+2.*alphaK[Material][1]*finalenergy);
            }
            else SWK_N[Material][2][1][ie]=0.;

            // Absorption
            finalenergy=initialenergy+HIV[2]+EMIN[Material][2]-EMIN[Material][1];
            if(finalenergy>0.){
                sef=sqrt(finalenergy*(1.+alphaK[Material][1]*finalenergy));
                SWK_N[Material][2][2][ie]=SWK_N[Material][2][1][ie]+zl*opa*sef*dos1*(1.+2.*alphaK[Material][1]*finalenergy);
            }
            else SWK_N[Material][2][2][ie]=SWK_N[Material][2][1][ie];

            // G -> D (7)
            ope=PI*DIV[6]/w[6]*DIV[6]/RHO[Material]/Q*(n[6]+1.);
            opa=ope*n[6]/(1.+n[6]);
            // Emission
            finalenergy=initialenergy-HIV[6]+EMIN[Material][2]-EMIN[Material][3];
            if(finalenergy>0.){
                sef=sqrt(finalenergy*(1.+alphaK[Material][3]*finalenergy));
                SWK_N[Material][2][3][ie]=SWK_N[Material][2][2][ie]+zd*ope*sef*dos3*(1.+2.*alphaK[Material][3]*finalenergy);
            }
            else SWK_N[Material][2][3][ie]=SWK_N[Material][2][2][ie];

            // Absorption
            finalenergy=initialenergy+HIV[6]+EMIN[Material][2]-EMIN[Material][3];
            if(finalenergy>0.){
                sef=sqrt(finalenergy*(1.+alphaK[Material][3]*finalenergy));
                SWK_N[Material][2][4][ie]=SWK_N[Material][2][3][ie]+zd*opa*sef*dos3*(1.+2.*alphaK[Material][3]*finalenergy);
            }
            else SWK_N[Material][2][4][ie]=SWK_N[Material][2][3][ie];
        }
        else {
            // NO OPTICAL PHONON
            SWK_N[Material][2][1][ie]=0.0;
            SWK_N[Material][2][2][ie]=0.0;
            SWK_N[Material][2][3][ie]=0.0;
            SWK_N[Material][2][4][ie]=0.0;
        }

        if(ACOUSTICPHONONS==ON){
            // Acoustic Phonon Scattering
            finalenergy=initialenergy;
            sef=sqrt(finalenergy*(1.+alphaK[Material][2]*finalenergy));
            SWK_N[Material][2][5][ie]=SWK_N[Material][2][4][ie]
                                    +aco[2]*sef*dos2*(1.+2.*alphaK[Material][2]*finalenergy);
            }
        else {
            // NO ACOUSTIC PHONON
            SWK_N[Material][2][5][ie]=SWK_N[Material][2][4][ie]+0.0;
        }

        if(IMPURITYPHONONS==ON){
            // Impurity scattering
            finalenergy=initialenergy;
            sef=sqrt(finalenergy*(1.+alphaK[Material][2]*finalenergy));
            ak=SMH[Material][2]*sef;
            qq=QD2*(4.*ak*ak+QD2);
            wk=bimp/qq*sef*dos2*(1.+2.*alphaK[Material][2]*finalenergy);
            if(wk>1.e14) wk=1.e14;
            SWK_N[Material][2][6][ie]=SWK_N[Material][2][5][ie]+wk;
        }
        else {
            // NO IMPURITY SCATTERING
            SWK_N[Material][2][6][ie]=SWK_N[Material][2][5][ie]+0.0;
        }
        // Scattering rate calculation before normalization
        Scat_Rate_N[Material][2][1][ie]=SWK_N[Material][2][1][ie];
        for(i=2;i<=6;i++){
            Scat_Rate_N[Material][2][i][ie]=SWK_N[Material][2][i][ie]-SWK_N[Material][2][i-1][ie];
        }

        // Delta-Valley
        // =======================================

        if(OPTICALPHONONS==ON){
        // Intervalley phonon Scattering (iv)

            // D -> L (4)
            ope=PI*DIV[3]/w[3]*DIV[3]/RHO[Material]/Q*(n[3]+1.);
            opa=ope*n[3]/(1.+n[3]);
            // Emission
            finalenergy=initialenergy-HIV[3]+EMIN[Material][3]-EMIN[Material][1];
            if(finalenergy>0.){
                sef=sqrt(finalenergy*(1.+alphaK[Material][1]*finalenergy));
                SWK_N[Material][3][1][ie]=zl*ope*sef*dos1*(1.+2.*alphaK[Material][1]*finalenergy);
            }
            else SWK_N[Material][3][1][ie]=0.;

            // Absorption
            finalenergy=initialenergy+HIV[3]+EMIN[Material][3]-EMIN[Material][1];
            if(finalenergy>0.){
                sef=sqrt(finalenergy*(1.+alphaK[Material][1]*finalenergy));
                SWK_N[Material][3][2][ie]=SWK_N[Material][3][1][ie]+zl*opa*sef*dos1*(1.+2.*alphaK[Material][1]*finalenergy);
            }
            else SWK_N[Material][3][2][ie]=SWK_N[Material][3][1][ie];

            // D -> G (7)
            ope=PI*DIV[6]/w[6]*DIV[6]/RHO[Material]/Q*(n[6]+1.);
            opa=ope*n[6]/(1.+n[6]);
            // Emission
            finalenergy=initialenergy-HIV[6]+EMIN[Material][3]-EMIN[Material][2];
            if(finalenergy>0.){
                sef=sqrt(finalenergy*(1.+alphaK[Material][2]*finalenergy));
                SWK_N[Material][3][3][ie]=SWK_N[Material][3][2][ie]+ope*sef*dos2*(1.+2.*alphaK[Material][2]*finalenergy);
            }
            else SWK_N[Material][3][3][ie]=SWK_N[Material][3][2][ie];

            // Absorption
            finalenergy=initialenergy+HIV[6]+EMIN[Material][3]-EMIN[Material][2];
            if(finalenergy>0.){
                sef=sqrt(finalenergy*(1.+alphaK[Material][2]*finalenergy));
                SWK_N[Material][3][4][ie]=SWK_N[Material][3][3][ie]+opa*sef*dos2*(1.+2.*alphaK[Material][2]*finalenergy);
            }
            else SWK_N[Material][3][4][ie]=SWK_N[Material][3][3][ie];

            // D -> D (5)
            ope=PI*DIV[4]/w[4]*DIV[4]/RHO[Material]/Q*(n[4]+1.);
            opa=ope*n[4]/(1.+n[4]);
            // Emission
            finalenergy=initialenergy-HIV[4]+EMIN[Material][3]-EMIN[Material][3];
            if(finalenergy>0.){
                sef=sqrt(finalenergy*(1.+alphaK[Material][3]*finalenergy));
                SWK_N[Material][3][5][ie]=SWK_N[Material][3][4][ie]+(zd-1.)*ope*sef*dos3*(1.+2.*alphaK[Material][3]*finalenergy);
            }
            else SWK_N[Material][3][5][ie]=SWK_N[Material][3][4][ie];

            // Absorption
            finalenergy=initialenergy+HIV[4]+EMIN[Material][3]-EMIN[Material][3];
            if(finalenergy>0.){
                sef=sqrt(finalenergy*(1.+alphaK[Material][3]*finalenergy));
                SWK_N[Material][3][6][ie]=SWK_N[Material][3][5][ie]+(zd-1.)*opa*sef*dos3*(1.+2.*alphaK[Material][3]*finalenergy);
            }
            else SWK_N[Material][3][6][ie]=SWK_N[Material][3][5][ie];

            // D -> D (6)
            ope=PI*DIV[5]/w[5]*DIV[5]/RHO[Material]/Q*(n[5]+1.);
            opa=ope*n[5]/(1.+n[5]);
            // Emission
            finalenergy=initialenergy-HIV[5]+EMIN[Material][3]-EMIN[Material][3];
            if(finalenergy>0.){
                sef=sqrt(finalenergy*(1.+alphaK[Material][3]*finalenergy));
                SWK_N[Material][3][7][ie]=SWK_N[Material][3][6][ie]+(zd-1.)*ope*sef*dos3*(1.+2.*alphaK[Material][3]*finalenergy);
            }
            else SWK_N[Material][3][7][ie]=SWK_N[Material][3][6][ie];

            // Absorption
            finalenergy=initialenergy+HIV[5]+EMIN[Material][3]-EMIN[Material][3];
            if(finalenergy>0.){
                sef=sqrt(finalenergy*(1.+alphaK[Material][3]*finalenergy));
                SWK_N[Material][3][8][ie]=SWK_N[Material][3][7][ie]+(zd-1.)*opa*sef*dos3*(1.+2.*alphaK[Material][3]*finalenergy);
            }
            else SWK_N[Material][3][8][ie]=SWK_N[Material][3][7][ie];

        }
        else {
            // NO OPTICAL PHONON
            SWK_N[Material][1][1][ie]=0.0;
            SWK_N[Material][1][2][ie]=0.0;
            SWK_N[Material][1][3][ie]=0.0;
            SWK_N[Material][1][4][ie]=0.0;
            SWK_N[Material][1][5][ie]=0.0;
            SWK_N[Material][1][6][ie]=0.0;
            SWK_N[Material][1][7][ie]=0.0;
            SWK_N[Material][1][8][ie]=0.0;
        }

        if(ACOUSTICPHONONS==ON){
            // Acoustic Phonon Scattering
            finalenergy=initialenergy;
            sef=sqrt(finalenergy*(1.+alphaK[Material][3]*finalenergy));
            SWK_N[Material][3][9][ie]=SWK_N[Material][3][8][ie]
                                +aco[3]*sef*dos3*(1.+2.*alphaK[Material][3]*finalenergy);
        }
        else {
            // NO ACOUSTIC PHONON
            SWK_N[Material][3][9][ie]=SWK_N[Material][3][8][ie]+0.0;
        }

        if(IMPURITYPHONONS==ON){
            // Impurity scattering
            finalenergy=initialenergy;
            sef=sqrt(finalenergy*(1.+alphaK[Material][3]*finalenergy));
            ak=SMH[Material][3]*sef;
            qq=QD2*(4.*ak*ak+QD2);
            wk=bimp/qq*sef*dos3*(1.+2.*alphaK[Material][3]*finalenergy);
            if(wk>1.e14) wk=1.e14;
            SWK_N[Material][3][10][ie]=SWK_N[Material][3][9][ie]+wk;
        }
        else {
            // NO IMPURITY SCATTERING
            SWK_N[Material][3][10][ie]=SWK_N[Material][3][9][ie]+0.0;
        }
        // Scattering rate calculation before normalization
        Scat_Rate_N[Material][3][1][ie]=SWK_N[Material][3][1][ie];
        for(i=2;i<=10;i++){
            Scat_Rate_N[Material][3][i][ie]=SWK_N[Material][3][i][ie]-SWK_N[Material][3][i-1][ie];
        }
    }
    // Evaluation scattering rate
    gm[Material][1][1]=SWK_N[Material][1][12][1];
    for(ie=1;ie<=DIME;ie++){
        if(SWK_N[Material][1][12][ie]>gm[Material][1][1]) {
            gm[Material][1][1]=SWK_N[Material][1][12][ie];
        }
    }   
    printf("L[%d]N+ = %g \n",Material, gm[Material][1][1]);
    for(i=1;i<=12;i++){
        for(ie=1;ie<=DIME;ie++){
            SWK_N[Material][1][i][ie]/=gm[Material][1][1];
            if(SWK_N[Material][1][i][ie] > 1.0){
                SWK_N[Material][1][i][ie] = 1.0;
            }
        }
    }

    gm[Material][2][1]=SWK_N[Material][2][6][1];
    for(ie=1;ie<=DIME;ie++){
        if(SWK_N[Material][2][6][ie]>gm[Material][2][1]){
            gm[Material][2][1]=SWK_N[Material][2][6][ie];
        }
    }
    printf("GAMMA[%d]N+ = %g \n",Material, gm[Material][2][1]);
    for(i=1;i<=6;i++){
        for(ie=1;ie<=DIME;ie++){
            SWK_N[Material][2][i][ie]/=gm[Material][2][1];
            if(SWK_N[Material][2][i][ie] > 1.0){
                SWK_N[Material][2][i][ie] = 1.0;
            }
        }
    }

    gm[Material][3][1]=SWK_N[Material][3][10][1];
    for(ie=1;ie<=DIME;ie++){
        if(SWK_N[Material][3][10][ie]>gm[Material][3][1]){
            gm[Material][3][1]=SWK_N[Material][3][10][ie];
        }
    }
    printf("DELTA[%d]N+ = %g \n",Material, gm[Material][3][1]);
    for(i=1;i<=10;i++){
        for(ie=1;ie<=DIME;ie++){
            SWK_N[Material][3][i][ie]/=gm[Material][3][1];
            if(SWK_N[Material][3][i][ie] > 1.0){
                SWK_N[Material][3][i][ie] = 1.0;
            }
        }
    }
    for(i=0;i<=4;i++){
        for(j=0;j<=13;j++){
            SC[i][j]=0;
        }
    }
    // End of Three-valleys materials
    // ############################
  }
}

// =======================================================
