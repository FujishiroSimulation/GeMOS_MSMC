/* saveoutput2dgnuplot.h -- This file is part of Archimedes release 0.1.0.
   Archimedes is a simulator for Submicron 2D III-V semiconductor
   Devices. It implements the Monte Carlo method and Hybrid MEP model
   for the simulation of the semiclassical Boltzmann equation for both
   electrons and holes. It also includes the quantum effects by means 
   of effective potential method. It is now able to simulate applied
   magnetic fields along with self consistent Faraday equation.

   Copyright (C) 2004, 2005, 2006, 2007 Jean Michel Sellier <sellier@dmi.unict.it>
 
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
// Created on 10 Sep.2004, Siracusa (Italy), J.M.Sellier
// Last modif. : 10 Sep.2007, Siracusa (Italy), J.M.Sellier
// ######################################################

// Here we save all the macroscopic variables 
// related to the electrons.
// The format of saving is compatible with GNUPLOT
// and is as follows (constituted of 3 columns):
// x1 y1 u(x1,y1)
// x1 y2 u(x1,y2)
// ...
// xn yn u(xn,yn)
// Actually this file format is completly compatible
// with the "xyz" format of xd3d package. It is higly
// recomended that you use this program in order to
// obtain a very good scientifical visualisation of
// the solutions. In this way the post-processing become
// very easy and clear, so let use this program! It
// is free and you download it from :
// www.gnu.org ---> Free Directory ---> Science ---> Scientific visualisation

// INPUT : je (integer). It is the integer which will be the index of
// the solution. In this way you can easily create a succession of
// solutions like this:
// density000.xyz, density001.xyz, density002.xyz, ..., density010.xyz
// and so on, in order to create movies or animated gif.

void
SaveOutput2DGNUPLOT(int je)
{
 char s[150];
 FILE *fp;
 FILE *xv;
 FILE *xvt;
 FILE *up;
 FILE *vp;
 FILE *lp;
 FILE *xet;
 FILE *lxp;
 FILE *lyp;
 FILE *ep;
 FILE *qp;
 FILE *zp;
 FILE *fM;
 FILE *eg;
 FILE *oc;
 FILE *egl;
 FILE *egg;
 FILE *egx;
 FILE *scbl;
 FILE *scbg;
 FILE *scbx;
 FILE *scbln;
 FILE *scbgn;
 FILE *scbxn;
 FILE *scl;
 FILE *scg;
 FILE *scx;
 FILE *scln;
 FILE *scgn;
 FILE *scxn;
 FILE *sc;
 FILE *su;
 register int i,j;
 real x, y;

//if(Material==SILICON || Material==GERMANIUM){
 if(je<=9){
   sprintf(s,"density00%d.xyz",je);
   fp=fopen(s,"w");
   sprintf(s,"x_velocity_xcoordinate00%d.xyz",je);
   xv=fopen(s,"w");
   sprintf(s,"x_velocity_xtransit00%d.xyz",je);
   xvt=fopen(s,"w");
   sprintf(s,"x_velocity00%d.xyz",je);
   up=fopen(s,"w");
   sprintf(s,"y_velocity00%d.xyz",je);
   vp=fopen(s,"w");
   sprintf(s,"energy00%d.xyz",je);
   ep=fopen(s,"w");
   sprintf(s,"potential00%d.xyz",je);
   lp=fopen(s,"w");
   sprintf(s,"magnetic_field00%d.xyz",je);
   fM=fopen(s,"w");
   sprintf(s,"x_Efield_xtransit00%d.xyz",je);
   xet=fopen(s,"w");
   sprintf(s,"x_Efield00%d.xyz",je);
   lxp=fopen(s,"w");
   sprintf(s,"y_Efield00%d.xyz",je);
   lyp=fopen(s,"w");
   sprintf(s,"quantum_potential00%d.xyz",je);
   qp=fopen(s,"w");
   sprintf(s,"z_potential00%d.xyz",je);
   zp=fopen(s,"w");
   sprintf(s,"occupancy00%d.xyz",je);
   oc=fopen(s,"w");
   sprintf(s,"z_energy_distribution00%d.xyz",je);
   eg=fopen(s,"w");
   sprintf(s,"energy_distribution_L00%d.xyz",je);
   egl=fopen(s,"w");
   sprintf(s,"energy_distribution_G00%d.xyz",je);
   egg=fopen(s,"w");
   sprintf(s,"energy_distribution_X00%d.xyz",je);
   egx=fopen(s,"w");
   sprintf(s, "Scat_Rate_L00%d.xyz",je);
   scbl=fopen(s, "w");
   sprintf(s, "Scat_Rate_G00%d.xyz",je);
   scbg=fopen(s, "w");
   sprintf(s, "Scat_Rate_X00%d.xyz",je);
   scbx=fopen(s, "w");
   sprintf(s, "Scat_Rate_L_N+00%d.xyz",je);
   scbln=fopen(s, "w");
   sprintf(s, "Scat_Rate_G_N+00%d.xyz",je);
   scbgn=fopen(s, "w");
   sprintf(s, "Scat_Rate_X_N+00%d.xyz",je);
   scbxn=fopen(s, "w");
   sprintf(s,"scattering_rate_L00%d.xyz",je);
   scl=fopen(s,"w");
   sprintf(s,"scattering_rate_G00%d.xyz",je);
   scg=fopen(s,"w");
   sprintf(s,"scattering_rate_X00%d.xyz",je);
   scx=fopen(s,"w");
   sprintf(s,"scattering_rate_L_N+00%d.xyz",je);
   scln=fopen(s,"w");
   sprintf(s,"scattering_rate_G_N+00%d.xyz",je);
   scgn=fopen(s,"w");
   sprintf(s,"scattering_rate_X_N+00%d.xyz",je);
   scxn=fopen(s,"w");
   sprintf(s,"scatterign_count00%d.xyz",je);
   sc=fopen(s,"w");
   sprintf(s,"Summary00%d.xyz",je);
   su=fopen(s,"w");
 }
 else if(je>9 && je<=99){
   sprintf(s,"density0%d.xyz",je);
   fp=fopen(s,"w");
   sprintf(s,"x_velocity0%d.xyz",je);
   up=fopen(s,"w");
   sprintf(s,"y_velocity0%d.xyz",je);
   vp=fopen(s,"w");
   sprintf(s,"energy0%d.xyz",je);
   ep=fopen(s,"w");
   sprintf(s,"potential0%d.xyz",je);
   lp=fopen(s,"w");
   sprintf(s,"magnetic_field0%d.xyz",je);
   fM=fopen(s,"w");
   sprintf(s,"x_Efield0%d.xyz",je);
   lxp=fopen(s,"w");
   sprintf(s,"y_Efield0%d.xyz",je);
   lyp=fopen(s,"w");
   sprintf(s,"quantum_potential0%d.xyz",je);
   qp=fopen(s,"w");
 }
 else{
   sprintf(s,"density%d.xyz",je);
   fp=fopen(s,"w");
   sprintf(s,"x_velocity%d.xyz",je);
   up=fopen(s,"w");
   sprintf(s,"y_velocity%d.xyz",je);
   vp=fopen(s,"w");
   sprintf(s,"energy%d.xyz",je);
   ep=fopen(s,"w");
   sprintf(s,"potential%d.xyz",je);
   lp=fopen(s,"w");
   sprintf(s,"magnetic_field%d.xyz",je);
   fM=fopen(s,"w");
   sprintf(s,"x_Efield%d.xyz",je);
   lxp=fopen(s,"w");
   sprintf(s,"y_Efield%d.xyz",je);
   lyp=fopen(s,"w");
   sprintf(s,"quantum_potential%d.xyz",je);
   qp=fopen(s,"w");
 }
// Save the Monte Carlo results
 if(Model_Number==MCE || Model_Number==MCEH){
// Electron Density Output
// =======================
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(fp,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,u2d[i][j][1]);
      fprintf(fp,"\n");
  }
// X-velocity of X-coordinate output
// =================================
  float TV;
  for(i=1;i<=nx+1;i++){
    TV=0.0;
    for(j=1;j<=ny+1;j++){
      TV = TV + u2d[i][j][2]/MEDIA;
    }
    fprintf(xv,"%g %g\n",1.e6*(i-1.)*dx, TV/(ny+1));
  }
// X-velocity of X-transit output (34nm to 166nm portion in the case of 200nm)
// ==============================
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++){
      if(i >= 18 && i <= 84){
        fprintf(xvt,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,
              u2d[i][j][2]/MEDIA);
      }
    }  
  }
// X-component of electronic velocity output
// =========================================
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(up,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,
              u2d[i][j][2]/MEDIA);
      fprintf(up,"\n");
  }
// Y-component of electronic velocity output
// =========================================
  for (j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(vp,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,
              u2d[i][j][3]/MEDIA);
      fprintf(vp,"\n");
  }
// Electrostatic Potential
// =======================
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(lp,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,-PSI[i][j]);
      fprintf(lp,"\n");
  }
// Magnetic Field
// ==============
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(fM,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,B[i][j]);
      fprintf(fM,"\n");
  }
// X-field of X-transit output (34nm to 166nm portion in the case of 200nm)
// ==============================
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++){
      if(i >= 18 && i <= 84){
        fprintf(xet,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,E[i][j][0]);
      }
    }
  }
// X-component of electric field
// =============================
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(lxp,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,E[i][j][0]);
      fprintf(lxp,"\n");
  }
// Y-component of electric field
// =============================
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(lyp,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,E[i][j][1]);
      fprintf(lyp,"\n");
  }
// Electron Energy (in eV)
// ===============
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(ep,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,
              u2d[i][j][4]/MEDIA);
      fprintf(ep,"\n");
  }
// Quantum Effective Potential
// ===========================
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(qp,"%g %g %g\n",
          1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,-u2d[i][j][0]);
      fprintf(qp,"\n");
  }

// z-direction Effective Potential
// ===========================
  for(j=1;j<=ny+1;j++){
    fprintf(zp,"%g %g\n",
        1.e6*(j-1.)*dy,PSI[60][j]);
  }
  
// occupancy
// =========
  int CL, CG, CX, CT,SCL,SCG,SCX;
  float OL, OG, OX;
  SCL=SCG=SCX=0;
  for (j=0; j<=nx-1; j++){
    CL=0;
    CG=0;
    CX=0;
    CT=0;
    OL=0;
    OG=0;
    OX=0;
    for(i=1; i<=INUM; i++){
      if(int(P[i][5]/dx) == j){
        if(P[i][0] == 1){
          CL++;
          if(j>=int(undope_INI/dx) && j<=int(undope_FIN/dx)){
            SCL++;
          }
        }
        else if(P[i][0] == 2){
          CG++;
          if(j>=int(undope_INI/dx) && j<=int(undope_FIN/dx)){
            SCG++;
          }
        }
        else if(P[i][0] == 3){
          CX++;
          if(j>=int(undope_INI/dx) && j<=int(undope_FIN/dx)){
            SCX++;
          }
        }
      }
    }
    CT = CL + CG + CX;
    OL = float(CL) / float(CT);
    OG = float(CG) / float(CT);
    OX = float(CX) / float(CT);
    fprintf(oc,"%g %g %g %g\n",
        1.e6*j*dx, OL, OG, OX);
  }
// Energy distribution
// ===========================
  float energy_L, energy_G, energy_X;
  for(i=1; i<=INUM; i++){
    if(P[i][0] == 1){
      x = P[i][5]/dx;
      y = P[i][6]/dy;
      energy_L = pow(HBAR, 2)*(pow(P[i][1], 2) + pow(P[i][2], 2) + pow(P[i][3], 2)) / (2.*MSTAR[GERMANIUM][1][0] * M)/Q;
      fprintf(egl,"%g %g\n",
          P[i][5], energy_L-PSI[int(x)][int(y)]);
    }
    else if(P[i][0] == 2){
      x = P[i][5]/dx;
      y = P[i][6]/dy;
      energy_G = pow(HBAR, 2)*(pow(P[i][1], 2) + pow(P[i][2], 2) + pow(P[i][3], 2)) / (2.*MSTAR[GERMANIUM][2][0] * M)/Q;
      fprintf(egg,"%g %g\n",
          P[i][5], energy_G-PSI[int(x)][int(y)]+0.14194);
    }
    else{
      x = P[i][5]/dx;
      y = P[i][6]/dy;
      energy_X = pow(HBAR, 2)*(pow(P[i][1], 2) + pow(P[i][2], 2) + pow(P[i][3], 2)) / (2.*MSTAR[GERMANIUM][3][0] * M)/Q;
      fprintf(egx,"%g %g\n",
          P[i][5], energy_X-PSI[int(x)][int(y)]+0.14946);
    }
  }
// Energy distribution (no energy)
// =====================
  for(i=1; i<=INUM; i++){
      fprintf(eg,"%g %g\n",
          P[i][5],P[i][6]);
  }
// Scattering rate before normalization
// ====================================
// undope layer
// Valley L => 12 G => 6 X => 10
  int ie=1;
  for(ie=1;ie<DIME;ie++){
    fprintf(scbl, "%g %g %g %g %g %g %g %g %g %g %g %g %g\n",
      DE*((real)(ie)), Scat_Rate[2][1][1][ie], Scat_Rate[2][1][2][ie], Scat_Rate[2][1][3][ie], Scat_Rate[2][1][4][ie], Scat_Rate[2][1][5][ie], Scat_Rate[2][1][6][ie], Scat_Rate[2][1][7][ie], Scat_Rate[2][1][8][ie], Scat_Rate[2][1][9][ie], Scat_Rate[2][1][10][ie], Scat_Rate[2][1][11][ie], Scat_Rate[2][1][12][ie]);
  }
  for(ie=1;ie<DIME;ie++){
    fprintf(scbg, "%g %g %g %g %g %g %g\n",
      DE*((real)(ie)), Scat_Rate[2][2][1][ie], Scat_Rate[2][2][2][ie], Scat_Rate[2][2][3][ie], Scat_Rate[2][2][4][ie], Scat_Rate[2][2][5][ie], Scat_Rate[2][2][6][ie]);
  }
  for(ie=1;ie<DIME;ie++){
    fprintf(scbx, "%g %g %g %g %g %g %g %g %g %g %g\n",
      DE*((real)(ie)), Scat_Rate[2][3][1][ie], Scat_Rate[2][3][2][ie], Scat_Rate[2][3][3][ie], Scat_Rate[2][3][4][ie], Scat_Rate[2][3][5][ie], Scat_Rate[2][3][6][ie], Scat_Rate[2][3][7][ie], Scat_Rate[2][3][8][ie], Scat_Rate[2][3][9][ie], Scat_Rate[2][3][10][ie]);
  }
// n+ layer
// Valley L => 12 G => 6 X => 10
  for(ie=1;ie<DIME;ie++){
    fprintf(scbln, "%g %g %g %g %g %g %g %g %g %g %g %g %g\n",
      DE*((real)(ie)), Scat_Rate_N[2][1][1][ie], Scat_Rate_N[2][1][2][ie], Scat_Rate_N[2][1][3][ie], Scat_Rate_N[2][1][4][ie], Scat_Rate_N[2][1][5][ie], Scat_Rate_N[2][1][6][ie], Scat_Rate_N[2][1][7][ie], Scat_Rate_N[2][1][8][ie], Scat_Rate_N[2][1][9][ie], Scat_Rate_N[2][1][10][ie], Scat_Rate_N[2][1][11][ie], Scat_Rate_N[2][1][12][ie]);
  }
  for(ie=1;ie<DIME;ie++){
    fprintf(scbgn, "%g %g %g %g %g %g %g\n",
      DE*((real)(ie)), Scat_Rate_N[2][2][1][ie], Scat_Rate_N[2][2][2][ie], Scat_Rate_N[2][2][3][ie], Scat_Rate_N[2][2][4][ie], Scat_Rate_N[2][2][5][ie], Scat_Rate_N[2][2][6][ie]);
  }
  for(ie=1;ie<DIME;ie++){
    fprintf(scbxn, "%g %g %g %g %g %g %g %g %g %g %g\n",
      DE*((real)(ie)), Scat_Rate_N[2][3][1][ie], Scat_Rate_N[2][3][2][ie], Scat_Rate_N[2][3][3][ie], Scat_Rate_N[2][3][4][ie], Scat_Rate_N[2][3][5][ie], Scat_Rate_N[2][3][6][ie], Scat_Rate_N[2][3][7][ie], Scat_Rate_N[2][3][8][ie], Scat_Rate_N[2][3][9][ie], Scat_Rate_N[2][3][10][ie]);
  }

// Scattering rate
// ===============
// undope layer
// Valley  L => 12 G => 6 X => 10
  for(ie=1;ie<=DIME;ie++){
    fprintf(scl, "%g %g %g %g %g %g %g %g %g %g %g %g %g\n",
      DE*((real)(ie)),  SWK[2][1][1][ie], SWK[2][1][2][ie], SWK[2][1][3][ie], SWK[2][1][4][ie], SWK[2][1][5][ie], SWK[2][1][6][ie], SWK[2][1][7][ie], SWK[2][1][8][ie], SWK[2][1][9][ie], SWK[2][1][10][ie], SWK[2][1][11][ie], SWK[2][1][12][ie]);
  }
// Ganmma Valley
  for(ie=1;ie<=DIME;ie++){
    fprintf(scg, "%g %g %g %g %g %g %g\n",
      DE*((real)(ie)),  SWK[2][2][1][ie], SWK[2][2][2][ie], SWK[2][2][3][ie], SWK[2][2][4][ie], SWK[2][2][5][ie],  SWK[2][2][6][ie]);
  }

// X Valley
  for(ie=1;ie<=DIME;ie++){
    fprintf(scx, "%g %g %g %g %g %g %g %g %g %g %g\n",
      DE*((real)(ie)),  SWK[2][3][1][ie], SWK[2][3][2][ie], SWK[2][3][3][ie], SWK[2][3][4][ie], SWK[2][3][5][ie], SWK[2][3][6][ie], SWK[2][3][7][ie], SWK[2][3][8][ie], SWK[2][3][9][ie], SWK[2][3][10][ie]);
  }

// N+ layer
// Valley  L => 12 G => 6 X => 10
  for(ie=1;ie<=DIME;ie++){
    fprintf(scln, "%g %g %g %g %g %g %g %g %g %g %g %g %g\n",
      DE*((real)(ie)),  SWK_N[2][1][1][ie], SWK_N[2][1][2][ie], SWK_N[2][1][3][ie], SWK_N[2][1][4][ie], SWK_N[2][1][5][ie], SWK_N[2][1][6][ie], SWK_N[2][1][7][ie], SWK_N[2][1][8][ie], SWK_N[2][1][9][ie], SWK_N[2][1][10][ie], SWK_N[2][1][11][ie], SWK_N[2][1][12][ie]);
  }
// Ganmma Valley
  for(ie=1;ie<=DIME;ie++){
    fprintf(scgn, "%g %g %g %g %g %g %g\n",
      DE*((real)(ie)),  SWK_N[2][2][1][ie], SWK_N[2][2][2][ie], SWK_N[2][2][3][ie], SWK_N[2][2][4][ie], SWK_N[2][2][5][ie],  SWK_N[2][2][6][ie]);
  }

// X Valley
  for(ie=1;ie<=DIME;ie++){
    fprintf(scxn, "%g %g %g %g %g %g %g %g %g %g %g\n",
      DE*((real)(ie)),  SWK_N[2][3][1][ie], SWK_N[2][3][2][ie], SWK_N[2][3][3][ie], SWK_N[2][3][4][ie], SWK_N[2][3][5][ie], SWK_N[2][3][6][ie], SWK_N[2][3][7][ie], SWK_N[2][3][8][ie], SWK_N[2][3][9][ie], SWK_N[2][3][10][ie]);
  }

// Scattering count
  fprintf(sc, "ode E=%d, A=%d, L(1) E=%d, A=%d, L(2) E=%d, A=%d, G E=%d, A=%d, X E=%d, A=%d, Aco=%d, Imu=%d\n",
    SC[1][1], SC[1][2], SC[1][3], SC[1][4], SC[1][5], SC[1][6], SC[1][7], SC[1][8], SC[1][9], SC[1][10], SC[1][11], SC[1][12]);
  fprintf(sc, "L E=%d, A=%d, X E=%d, A=%d, Aco=%d, Imu=%d\n",
    SC[2][1], SC[2][2], SC[2][3], SC[2][4], SC[2][5], SC[2][6]);
  fprintf(sc, "L E=%d, A=%d, G E=%d, A=%d, X(1) E=%d, A=%d, X(2) E=%d, A=%d, Aco=%d, Imu=%d\n",
    SC[3][1], SC[3][2], SC[3][3], SC[3][4], SC[3][5], SC[3][6], SC[3][7], SC[3][8], SC[3][9], SC[3][10]);

// Summary
// =======
  int count;
  float SE, SV;
  count = 0;
  SE = 0.0;
  SV = 0.0;
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++){
      if(i >= 18 && i <= 84){
        SE = SE + E[i][j][0];
        SV = SV + u2d[i][j][2]/MEDIA;
        count++;
      }
    }
  }
  fprintf(su,"Efield:xvelocity(n+ region 1/6) %g %g\n",SE/float(count),SV/float(count));

  count = 0;
  SE = 0.0;
  SV = 0.0;
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++){
      if(i >= 31 && i <= 71){
        SE = SE + E[i][j][0];
        SV = SV + u2d[i][j][2]/MEDIA;
        count++;
      }
    }
  }
  fprintf(su,"Efield:xvelocity(n+ region 3/10) %g %g\n",SE/float(count),SV/float(count));

  count = 0;
  SE = 0.0;
  SV = 0.0;
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++){
      if(i >= 51 && i <= 91){
        SE = SE + E[i][j][0];
        SV = SV + u2d[i][j][2]/MEDIA;
        count++;
      }
    }
  }
  fprintf(su,"Efield:xvelocity(for140nm,N+:50nm,undope:40nm) %g %g\n",SE/float(count),SV/float(count));
  
  count = 0;
  SE = 0.0;
  SV = 0.0;
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++){
      if(i >= 46 && i <= 86){
        SE = SE + E[i][j][0];
        SV = SV + u2d[i][j][2]/MEDIA;
        count++;
      }
    }
  }
  fprintf(su,"Efield:xvelocity(for130nm,N+:45nm,undope:40nm) %g %g\n",SE/float(count),SV/float(count));
  
  count = 0;
  SE = 0.0;
  SV = 0.0;
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++){
      if(i >= 41 && i <= 81){
        SE = SE + E[i][j][0];
        SV = SV + u2d[i][j][2]/MEDIA;
        count++;
      }
    }
  }
  fprintf(su,"Efield:xvelocity(for120nm,N+:40nm,undope:40nm) %g %g\n",SE/float(count),SV/float(count));
  
  //ocupancy in chanel
  fprintf(su,"ocupancy in chanel L:G:X %d %d %d\n",SCL,SCG,SCX);

// Closure of output files
// =======================
   fclose(fp);
   fclose(xv);
   fclose(xvt);
   fclose(up);
   fclose(vp);
   fclose(lp);
   fclose(fM);
   fclose(xet);
   fclose(lxp);
   fclose(lyp);
   fclose(ep);
   fclose(qp);
   fclose(zp);
   fclose(eg);
   fclose(oc);
   fclose(egl);
   fclose(egg);
   fclose(egx);
   fclose(scbl);
   fclose(scbg);
   fclose(scbx);
   fclose(scbln);
   fclose(scbgn);
   fclose(scbxn);
   fclose(scl);
   fclose(scg);
   fclose(scx);
   fclose(scln);
   fclose(scgn);
   fclose(scxn);
   fclose(sc);
   fclose(su);
  }
// Save the Hybrid MEP results
 if(Model_Number==MEPE || Model_Number==MEPEH){
// Electron Density Output
// =======================
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(fp,"%g %g %g\n",1.e6*(i-1.)*dx,
              1.e6*(j-1.)*dy,u2d[i+2][j+2][1]);
      fprintf(fp,"\n");
  }
// X-component of electronic velocity output
// =========================================
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(up,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,
              u2d[i+2][j+2][2]/u2d[i+2][j+2][1]);
      fprintf(up,"\n");
  }
// Y-component of electronic velocity output
// =========================================
  for (j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(vp,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,
              u2d[i+2][j+2][3]/u2d[i+2][j+2][1]);
      fprintf(vp,"\n");
  }
// Electrostatic Potential
// =======================
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(lp,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,PSI[i][j]);
      fprintf(lp,"\n");
  }
// X-component of electric field
// =============================
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(lxp,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,E[i][j][0]);
      fprintf(lxp,"\n");
  }
// Y-component of electric field
// =============================
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(lyp,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,E[i][j][1]);
      fprintf(lyp,"\n");
  }
// Electron Energy (in eV)
// ===============
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(ep,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,
              u2d[i+2][j+2][4]/u2d[i+2][j+2][1]/Q);
      fprintf(ep,"\n");
  }
// Quantum Effective Potential
// ===========================
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(qp,"%g %g %g\n",
          1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,u2d[i][j][0]);
      fprintf(qp,"\n");
  }
// Closure of output files
// =======================
   fclose(fp);
   fclose(up);
   fclose(vp);
   fclose(lp);
   fclose(lxp);
   fclose(lyp);
   fclose(ep);
   fclose(qp);
  }
// }
// ****************************************************************
/*if(Material==GAAS){
 if(je<=9){
   sprintf(s,"density00%d.xyz",je);
   fp=fopen(s,"w");
   sprintf(s,"densityGAMMA_valley00%d.xyz",je);
   fG=fopen(s,"w");
   sprintf(s,"densityL_valley00%d.xyz",je);
   fL=fopen(s,"w");
   sprintf(s,"x_velocity00%d.xyz",je);
   up=fopen(s,"w");
   sprintf(s,"y_velocity00%d.xyz",je);
   vp=fopen(s,"w");
   sprintf(s,"energy00%d.xyz",je);
   ep=fopen(s,"w");
   sprintf(s,"potential00%d.xyz",je);
   lp=fopen(s,"w");
   sprintf(s,"x_Efield00%d.xyz",je);
   lxp=fopen(s,"w");
   sprintf(s,"y_Efield00%d.xyz",je);
   lyp=fopen(s,"w");
   sprintf(s,"quantum_potential00%d.xyz",je);
   qp=fopen(s,"w");
 }
 else if(je>9 && je<=99){
   sprintf(s,"density0%d.xyz",je);
   fp=fopen(s,"w");
   sprintf(s,"densityGAMMA_valley0%d.xyz",je);
   fG=fopen(s,"w");
   sprintf(s,"densityL_valley0%d.xyz",je);
   fL=fopen(s,"w");
   sprintf(s,"x_velocity0%d.xyz",je);
   up=fopen(s,"w");
   sprintf(s,"y_velocity0%d.xyz",je);
   vp=fopen(s,"w");
   sprintf(s,"energy0%d.xyz",je);
   ep=fopen(s,"w");
   sprintf(s,"potential0%d.xyz",je);
   lp=fopen(s,"w");
   sprintf(s,"x_Efield0%d.xyz",je);
   lxp=fopen(s,"w");
   sprintf(s,"y_Efield0%d.xyz",je);
   lyp=fopen(s,"w");
   sprintf(s,"quantum_potential0%d.xyz",je);
   qp=fopen(s,"w");
 }
 else{
   sprintf(s,"density%d.xyz",je);
   fp=fopen(s,"w");
   sprintf(s,"densityGAMMA_valley%d.xyz",je);
   fG=fopen(s,"w");
   sprintf(s,"densityL_valley%d.xyz",je);
   fL=fopen(s,"w");
   sprintf(s,"x_velocity%d.xyz",je);
   up=fopen(s,"w");
   sprintf(s,"y_velocity%d.xyz",je);
   vp=fopen(s,"w");
   sprintf(s,"energy%d.xyz",je);
   ep=fopen(s,"w");
   sprintf(s,"potential%d.xyz",je);
   lp=fopen(s,"w");
   sprintf(s,"x_Efield%d.xyz",je);
   lxp=fopen(s,"w");
   sprintf(s,"y_Efield%d.xyz",je);
   lyp=fopen(s,"w");
   sprintf(s,"quantum_potential%d.xyz",je);
   qp=fopen(s,"w");
 }

// Electron Density Output
// =======================
 for(j=1;j<=ny+1;j++){
   for(i=1;i<=nx+1;i++)
     fprintf(fp,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,u2d[i][j][1]);
     fprintf(fp,"\n");
 }
// Electron Density GAMMA-valley Output
// ====================================
 for(j=1;j<=ny+1;j++){
   for(i=1;i<=nx+1;i++)
     fprintf(fG,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,DG[i][j]);
     fprintf(fG,"\n");
 }
// Electron Density L-valley Output
// ================================
 for(j=1;j<=ny+1;j++){
   for(i=1;i<=nx+1;i++)
     fprintf(fL,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,DL[i][j]);
     fprintf(fL,"\n");
 }
// X-component of electronic velocity output
// =========================================
 for(j=1;j<=ny+1;j++){
   for(i=1;i<=nx+1;i++)
     fprintf(up,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,
             u2d[i][j][2]/MEDIA);
     fprintf(up,"\n");
 }
// Y-component of electronic velocity output
// =========================================
 for (j=1;j<=ny+1;j++){
   for(i=1;i<=nx+1;i++)
     fprintf(vp,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,
             u2d[i][j][3]/MEDIA);
     fprintf(vp,"\n");
 }
// Electrostatic Potential
// =======================
 for(j=1;j<=ny+1;j++){
   for(i=1;i<=nx+1;i++)
     fprintf(lp,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,PSI[i][j]);
     fprintf(lp,"\n");
 }
// X-component of electric field
// =============================
 for(j=1;j<=ny+1;j++){
   for(i=1;i<=nx+1;i++)
     fprintf(lxp,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,E[i][j][0]);
     fprintf(lxp,"\n");
 }
// Y-component of electric field
// =============================
 for(j=1;j<=ny+1;j++){
   for(i=1;i<=nx+1;i++)
     fprintf(lyp,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,E[i][j][1]);
     fprintf(lyp,"\n");
 }
// Electron Energy (in eV)
// ===============
 for(j=1;j<=ny+1;j++){
   for(i=1;i<=nx+1;i++)
     fprintf(ep,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,
             u2d[i][j][4]/MEDIA);
     fprintf(ep,"\n");
 }
// Quantum Effective Potential
// ===========================
 for(j=1;j<=ny+1;j++){
   for(i=1;i<=nx+1;i++)
     fprintf(qp,"%g %g %g\n",
         1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,u2d[i][j][0]);
     fprintf(qp,"\n");
 }

// Closure of output files
// =======================
 fclose(fp);
 fclose(up);
 fclose(vp);
 fclose(lp);
 fclose(lxp);
 fclose(lyp);
 fclose(ep);
 fclose(qp);
 fclose(fG);
 fclose(fL);
 }*/
}

// ============================================================
