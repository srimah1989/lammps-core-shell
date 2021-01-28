/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_core_shell.h"
#include "atom.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define SMALL 0.001

/*---------------------------------------------------------------- */

PairCoreShell::PairCoreShell(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;
  single_enable = 0;
  manybody_flag = 1;
}

/* ---------------------------------------------------------------------- */

PairCoreShell::~PairCoreShell()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut_lj);
    memory->destroy(cut_ljsq);
    memory->destroy(a);
    memory->destroy(rho);
    memory->destroy(c);
    memory->destroy(rhoinv);
    memory->destroy(buck1);
    memory->destroy(buck2);
    memory->destroy(offset);
    memory->destroy(kspring);	//Destroy the array for the kspring (eVA**-2)
    memory->destroy(kspring4);   //Higher Order Sprink Term 
    memory->destroy(cut_spring);
    memory->destroy(cut_springsq);
    memory->destroy(cut_3body);
    memory->destroy(cut_3bodysq);
  }
}

/* ---------------------------------------------------------------------- */

void PairCoreShell::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair;
  double rsq,r2inv,r6inv,forcecoul,forcebuck,factor_coul,factor_lj;
  double prefactor;
  double r,rexp,rinv;
  int *ilist,*jlist,*numneigh,**firstneigh;
  
  double onebysix, oneby24;
  
  onebysix = 1/6;
  oneby24 = 1/24;  
  
  //Reducing FULL to HALF N List 
  
  int itag,jtag;
  int *tag = atom->tag;
  
  // Pair Style Coul Wolf variables
  double erfcc,erfcd,v_sh,dvdrr,e_self,e_shift,f_shift,qisq;

  evdwl = ecoul = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

// Three Body term
  
  double cs, sn, theta, del_theta;
  double E1, fs, ft, fu;
  double fj[3], fk[3], delr1[3], delr2[3], delr3[3];
  
  double rinvsq1, rinvsq2, rinv12, rinv13, rinv23, r1, r2, r3, rsq1, rsq2, rsq3;
  
  int k, kk, ktype, jnumm1;

//-----------------------------//

double minus;
minus  = -1.0;

//-----------------------------------------------//

// self and shifted coulombic energy - Wolf Summation

  e_self = v_sh = 0.0;
  e_shift = erfc(alf*cut_coul)/cut_coul;
  f_shift = -(e_shift+ 2.0*alf/MY_PIS * exp(-alf*alf*cut_coul*cut_coul))/cut_coul;

//------------------------------//

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    itag = tag[i];

    qisq = q_ion[itype]*q_ion[itype];
    //Self Energy calculation	 
    e_self = -(e_shift/2.0 + alf/MY_PIS) * qisq*qqrd2e;
    if (evflag){
     ev_tally(i,i,nlocal,0,0.0,e_self,0.0,0.0,0.0,0.0);
            }

    for (jj = 0; jj < jnum; jj++) {   // Start of 2 body term 
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      
      j &= NEIGHMASK;
      
      jtag = tag[j];
      
      if (itag > jtag) {
        if ((itag+jtag) % 2 == 0) continue;
      } else if (itag < jtag) {
        if ((itag+jtag) % 2 == 1) continue;
      } else {
        if (x[j][2] < ztmp) continue;
        if (x[j][2] == ztmp && x[j][1] < ytmp) continue;
        if (x[j][2] == ztmp && x[j][1] == ytmp && x[j][0] < xtmp) continue;
      }

      delx =  xtmp - x[j][0];
      dely =  ytmp - x[j][1];
      delz =  ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;		//Compute distance between i and j
      jtype = type[j];
      
      if(rsq < cut_springsq[itype][jtype]){
       r = sqrt(rsq);
       rinv = 1/r;
       fpair = minus * kspring[itype][jtype] + minus * kspring4[itype][jtype] * rsq * onebysix ;	//Force field spring term 
       
        //Update spring term to force
    	f[i][0] += delx * fpair;
    	f[i][1] += dely * fpair;
    	f[i][2] += delz * fpair;

	f[j][0] -= delx * fpair;
    	f[j][1] -= dely * fpair;
    	f[j][2] -= delz * fpair;

    	
    	if(eflag){
    		 // Compute energy of Harmonic spring term
    		 evdwl =  0.5 * kspring[itype][jtype] * rsq + kspring4[itype][jtype] * rsq * rsq * oneby24 ;
    	}
    	
    	if (evflag){
//     	 ev_tally_full(i,evdwl,0.0,fpair,delx,dely,delz);
	 ev_tally(i,j,nlocal,newton_pair, evdwl,0.0,fpair,delx,dely,delz);
	}
      }		//End of 'IF' command after checking spring
      else{
      
      if (rsq < cutsq[itype][jtype]) {		// Wolf summation two body term
 	
        r2inv = 1.0/rsq;          
        if (rsq < cut_coulsq) {
          r = sqrt(rsq);
          prefactor = qqrd2e*qtmp*q[j]/r;
          erfcc = erfc(alf*r);
          erfcd = exp(-alf*alf*r*r);
          v_sh = (erfcc - e_shift*r) * prefactor;
          dvdrr = (erfcc/rsq + 2.0*alf/MY_PIS * erfcd/r) + f_shift;
          forcecoul = dvdrr*rsq*prefactor;
          if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;
        } else forcecoul = 0.0;
	
        if (rsq < cut_ljsq[itype][jtype]) {	// Buckingham term
          r = sqrt(rsq);
          r6inv = r2inv*r2inv*r2inv;
          rexp = exp(-r*rhoinv[itype][jtype]);
          forcebuck = buck1[itype][jtype]*r*rexp - buck2[itype][jtype]*r6inv;
          
        } else forcebuck = 0.0;
	// Update coulomb and buckingham term
        fpair = (forcecoul + factor_lj*forcebuck) * r2inv;
	
	
        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
	
        f[j][0] -= delx*fpair;
        f[j][1] -= dely*fpair;
        f[j][2] -= delz*fpair;

        if (eflag) {			// Update two body energy term
          if (rsq < cut_coulsq) {		
              ecoul =   v_sh;
            if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor;
          } else ecoul = 0.0;
          if (rsq < cut_ljsq[itype][jtype]) {
            evdwl = a[itype][jtype]*rexp - c[itype][jtype]*r6inv -
              offset[itype][jtype];
            evdwl *= factor_lj;
          } else evdwl = 0.0;
        }
        if (evflag) {
// 	      ev_tally_full(i,evdwl,ecoul,fpair,delx,dely,delz);
	      ev_tally(i,j,nlocal,newton_pair,evdwl,ecoul,fpair,delx,dely,delz);
	  
	}
          
            
      }
      } 	//End of 'ELSE' command after checking spring cutoff
    }
    
     //Three Body term 
    if(si_tag == 0 ) continue;	// To find the Si atom for 3 body calculation
    jnumm1 = jnum -1;
    
    for(jj = 0; jj < jnumm1; jj++){	
    	j = jlist[jj];
      	j &= NEIGHMASK;
      	jtype = type[j];
	      
	delr1[0] =  x[j][0] - xtmp;
      	delr1[1] =  x[j][1] - ytmp;
      	delr1[2] =  x[j][2] - ztmp;
      	rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];
      	
      	if (rsq1 > cut_3bodysq[itype][jtype]) continue; //Check distance between i and jj < cut off
      	
      	for (kk = jj+1; kk < jnum; kk++) {	
      	
        	k = jlist[kk];
       	 	k &= NEIGHMASK;
        	ktype = type[k];
      	        	
        	delr2[0] = x[k][0] - xtmp;
	        delr2[1] = x[k][1] - ytmp;
        	delr2[2] = x[k][2] - ztmp;
        	rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
    
		delr3[0] = x[j][0] - x[k][0];
	        delr3[1] = x[j][1] - x[k][1];
        	delr3[2] = x[j][2] - x[k][2];
        	rsq3 = delr3[0]*delr3[0] + delr3[1]*delr3[1] + delr3[2]*delr3[2];
        	
        	r1 = sqrt(rsq1);
  		rinvsq1 = 1.0/rsq1;
  
		r2 = sqrt(rsq2);
  		rinvsq2 = 1.0/rsq2;
  
		r3 = sqrt(rsq3);
  
  		rinv12 = 1.0/(r1*r2);
  
  		rinv13 = 1.0/(r1*r3);
  
  		rinv23 = 1.0/(r2*r3);
  		
  		if (rsq2 > cut_3bodysq[itype][ktype]) continue;	//To restrict the 3 body computation within a single silicate
  		
  		if (rsq3 > cut_3bodysq[jtype][ktype]) continue;
    		
    		if(itype == si_tag)
    		{
       		threebody(rsq1,rsq2,delr1,delr2,fj,fk,eflag,evdwl);	// Update force fields
       		// For debugging force terms
       		//fprintf(screen,"\ni = %d j = %d k = %d ",i,j,k);
       		
       		//fprintf(logfile,"\n Out side :: fj = %g fk = %g fi = %g\n",fj[2],fk[2],fj[2]+fk[2]);
       		//fprintf(screen,"\n %d \t %d \t %d  \t %lf \n ",itype,jtype,ktype,evdwl);
		   		  f[i][0] -= fj[0] + fk[0];
        			  f[i][1] -= fj[1] + fk[1];
        			  f[i][2] -= fj[2] + fk[2];
			          f[j][0] += fj[0];
				  f[j][1] += fj[1];
				  f[j][2] += fj[2];
				  f[k][0] += fk[0];
				  f[k][1] += fk[1];
				  f[k][2] += fk[2];
				  
        			if (evflag) ev_tally3(i,j,k,evdwl,0.0,fj,fk,delr1,delr2);
       		}
       		else continue;
    	}
    } //End of jj loop of the Three body term ;)
  }
//   for(i=0;i<nlocal;i++)  fprintf(logfile,"\n%d\t 1 %lf 2 %lf 3 %lf\n",i,f[i][0],f[i][1],f[i][2]);
  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairCoreShell::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut_lj,n+1,n+1,"pair:cut_lj");
  memory->create(cut_ljsq,n+1,n+1,"pair:cut_ljsq");
  memory->create(a,n+1,n+1,"pair:a");
  memory->create(rho,n+1,n+1,"pair:rho");
  memory->create(c,n+1,n+1,"pair:c");
  memory->create(rhoinv,n+1,n+1,"pair:rhoinv");
  memory->create(buck1,n+1,n+1,"pair:buck1");
  memory->create(buck2,n+1,n+1,"pair:buck2");
  memory->create(offset,n+1,n+1,"pair:offset");
  memory->create(kspring,n+1,n+1,"pair:kspring"); //Allocate size for the array for the kspring (eVA**-2)
  memory->create(kspring4,n+1,n+1,"pair:kspring4");//Allocate size for the array for the kspring (eVA**-2) - 21/05/14
  memory->create(cut_spring,n+1,n+1,"pair:cut_spring");
  memory->create(cut_springsq,n+1,n+1,"pair:cut_springsq");
  memory->create(cut_3body,n+1,n+1,"pair:cut_3body");
  memory->create(cut_3bodysq,n+1,n+1,"pair:cut_3bodysq");	
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairCoreShell::settings(int narg, char **arg)
{
  
  int crit = 0;
    
  int n_types = atom->ntypes;
  
  crit =  5 + n_types ;
  
  if (narg < crit || narg > crit+1) error->all(FLERR,"Illegal pair_style command");

  si_tag  = force->inumeric(FLERR,arg[0]);
  k3body = force->numeric(FLERR,arg[1]);
  th_theta = force->numeric(FLERR,arg[2]);
  alf = force->numeric(FLERR,arg[3]);
  
  memory->create(q_ion,n_types+1,"pair:q_ion");
  
  int ii, jj;
  jj = 1;
  for(int ii = 4;ii<crit-1;ii++){
    q_ion[jj] = force->numeric(FLERR,arg[ii]);
    jj++;
  }  
    
  cut_lj_global = force->numeric(FLERR,arg[crit-1]);
  if (narg == crit) cut_coul = cut_lj_global;
  else cut_coul = force->numeric(FLERR,arg[crit]);
  
  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut_lj[i][j] = cut_lj_global;
  }
  

}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairCoreShell::coeff(int narg, char **arg)
{
  if (narg < 9 || narg > 10) 
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  double a_one = force->numeric(FLERR,arg[2]);
  double rho_one = force->numeric(FLERR,arg[3]);
  if (rho_one <= 0) error->all(FLERR,"Incorrect args for pair coefficients");
  double c_one = force->numeric(FLERR,arg[4]);
  double kspring_one = force->numeric(FLERR,arg[5]);	//csm to get spring constant for shells
  double kspring4_one = force->numeric(FLERR,arg[6]);	//csm to get spring constant for shells 
  double cut_spring_one = force->numeric(FLERR,arg[7]);
  double cut_3body_one = force->numeric(FLERR,arg[8]);

  double cut_lj_one = cut_lj_global;
  if (narg == 10) cut_lj_one = force->numeric(FLERR,arg[9]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      a[i][j] = a_one;
      rho[i][j] = rho_one;
      c[i][j] = c_one;
      cut_lj[i][j] = cut_lj_one;
      kspring[i][j] = kspring_one;
      kspring4[i][j] = kspring4_one;
      cut_spring[i][j] = cut_spring_one;
      cut_3body[i][j] = cut_3body_one;
      setflag[i][j] = 1;
      count++;
    }
  }
  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairCoreShell::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");
  
  double cut = MAX(cut_lj[i][j],cut_coul);
  cut_ljsq[i][j] = cut_lj[i][j] * cut_lj[i][j];

  rhoinv[i][j] = 1.0/rho[i][j];
  buck1[i][j] = a[i][j]/rho[i][j];
  buck2[i][j] = 6.0*c[i][j];

  if (offset_flag) {
    double rexp = exp(-cut_lj[i][j]/rho[i][j]);
    offset[i][j] = a[i][j]*rexp - c[i][j]/pow(cut_lj[i][j],6.0);
  } else offset[i][j] = 0.0;

  cut_ljsq[j][i] = cut_ljsq[i][j];
  a[j][i] = a[i][j];
  c[j][i] = c[i][j];
  rhoinv[j][i] = rhoinv[i][j];
  buck1[j][i] = buck1[i][j];
  buck2[j][i] = buck2[i][j];
  offset[j][i] = offset[i][j];
  kspring[j][i] = kspring[i][j];
  kspring4[j][i] = kspring4[i][j];
  
  cut_springsq[i][j] = cut_spring[i][j] * cut_spring[i][j];
  cut_springsq[j][i] = cut_springsq[i][j];
  
  cut_3bodysq[i][j] = cut_3body[i][j] * cut_3body[i][j];
  cut_3bodysq[j][i] = cut_3bodysq[i][j];
	
  // compute I,J contribution to long-range tail correction
  // count total # of atoms of type I and J via Allreduce

  if (tail_flag) {
    int *type = atom->type;
    int nlocal = atom->nlocal;

    double count[2],all[2];
    count[0] = count[1] = 0.0;
    for (int k = 0; k < nlocal; k++) {
      if (type[k] == i) count[0] += 1.0;
      if (type[k] == j) count[1] += 1.0;
    }
    MPI_Allreduce(count,all,2,MPI_DOUBLE,MPI_SUM,world);

    double rho1 = rho[i][j];
    double rho2 = rho1*rho1;
    double rho3 = rho2*rho1;
    double rc = cut_lj[i][j];
    double rc2 = rc*rc;
    double rc3 = rc2*rc;
    etail_ij = 2.0*MY_PI*all[0]*all[1]*
      (a[i][j]*exp(-rc/rho1)*rho1*(rc2 + 2.0*rho1*rc + 2.0*rho2) -
       c[i][j]/(3.0*rc3));
    ptail_ij = (-1/3.0)*2.0*MY_PI*all[0]*all[1]*
      (-a[i][j]*exp(-rc/rho1)*
       (rc3 + 3.0*rho1*rc2 + 6.0*rho2*rc + 6.0*rho3) + 2.0*c[i][j]/rc3);
  }

  return cut;
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairCoreShell::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR,"Pair style core/shell requires atom attribute q");
    
   if (force->newton_pair == 0)
    error->all(FLERR,"Pair style core/shell requires newton pair on");

  cut_coulsq = cut_coul * cut_coul;

// need a full neighbor list

  int irequest = neighbor->request(this);
  
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairCoreShell::write_restart(FILE *fp)
{
  write_restart_settings(fp);
  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&a[i][j],sizeof(double),1,fp);
        fwrite(&rho[i][j],sizeof(double),1,fp);
        fwrite(&c[i][j],sizeof(double),1,fp);
        fwrite(&cut_lj[i][j],sizeof(double),1,fp);
        fwrite(&kspring[i][j],sizeof(double),1,fp);
        fwrite(&kspring4[i][j],sizeof(double),1,fp);
        fwrite(&cut_spring[i][j],sizeof(double),1,fp);
        fwrite(&cut_3body[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCoreShell::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  
  allocate();
  
  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&a[i][j],sizeof(double),1,fp);
          fread(&rho[i][j],sizeof(double),1,fp);
          fread(&c[i][j],sizeof(double),1,fp);
          fread(&cut_lj[i][j],sizeof(double),1,fp);
          fread(&kspring[i][j],sizeof(double),1,fp);
          fread(&kspring4[i][j],sizeof(double),1,fp);
          fread(&cut_spring[i][j],sizeof(double),1,fp);
          fread(&cut_3body[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&a[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&rho[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&c[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_lj[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&kspring[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&kspring4[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_spring[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_3body[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairCoreShell::write_restart_settings(FILE *fp)
{
  int i;

  fwrite(&si_tag,sizeof(int),1,fp);
  fwrite(&k3body,sizeof(double),1,fp);
  fwrite(&th_theta,sizeof(double),1,fp);
  fwrite(&alf,sizeof(double),1,fp);

  for(i = 1; i <= atom->ntypes; i++){
  fwrite(&q_ion[i],sizeof(double),1,fp); 
  }

  fwrite(&cut_lj_global,sizeof(double),1,fp);
  fwrite(&cut_coul,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&tail_flag,sizeof(int),1,fp);
  
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCoreShell::read_restart_settings(FILE *fp)
{

  int i;
  
  if (comm->me == 0) {
    fread(&si_tag,sizeof(int),1,fp);
    fread(&k3body,sizeof(double),1,fp);
    fread(&th_theta,sizeof(double),1,fp);
    fread(&alf,sizeof(double),1,fp);
    memory->create(q_ion,atom->ntypes+1,"pair:q_ion");
    for(i = 1; i <= atom->ntypes; i++){
    fread(&q_ion[i],sizeof(double),1,fp);
    }	
    fread(&cut_lj_global,sizeof(double),1,fp);
    fread(&cut_coul,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
    fread(&tail_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&si_tag,1,MPI_INT,0,world);
  MPI_Bcast(&k3body,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&th_theta,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&alf,1,MPI_DOUBLE,0,world);
  for(i = 1; i <= atom->ntypes; i++){
  MPI_Bcast(&q_ion[i],1,MPI_DOUBLE,0,world);
  }
  MPI_Bcast(&cut_lj_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_coul,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&tail_flag,1,MPI_INT,0,world);
 
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairCoreShell::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g %g %g %g %g %g\n",i,a[i][i],rho[i][i],c[i][i],kspring[i][i],kspring4[i][i],cut_spring[i][i],cut_3body[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairCoreShell::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g %g %g %g %g %g\n",i,j,
              a[i][j],rho[i][j],c[i][j],kspring[i][j],kspring4[i][j],cut_spring[i][j],cut_3body[i][j],cut_lj[i][j]);
}


/* ---------------------------------------------------------------------- */
void PairCoreShell::threebody(double rsq1, double rsq2, double *delr1, double *delr2, double *fj, double *fk,int eflag, double &eng)
{
       
  double r1,rinvsq1;
  double r2,rinvsq2;
  double r3;
  double rinv12,rinv13,rinv23,cs,sn, del_theta,theta, tk;
  double fs, ft, fu;
  double minus = -1.0;
  //double th_theta = 1.910611932;
  double sign;
   double a,a11,a12,a22;
  
  
  fj[0] = 0.0;fj[1] = 0.0;fj[2] = 0.0;
  fk[0] = 0.0;fk[1] = 0.0;fk[2] = 0.0;

  r1 = sqrt(rsq1);
  rinvsq1 = 1.0/rsq1;
  
  r2 = sqrt(rsq2);
  rinvsq2 = 1.0/rsq2;
  
  rinv12 = 1.0/(r1*r2);
  
  //fprintf(screen,"\nEntry   %g \n ",th_theta);
  
   cs = (delr1[0]*delr2[0] + delr1[1]*delr2[1] + delr1[2]*delr2[2]) * rinv12;

   if (cs > 1.0) cs = 1.0;
   if (cs < -1.0) cs = -1.0;
  
   sn = sqrt(1.0 - cs*cs);
   if (sn < SMALL) sn = SMALL;
   sn = 1.0/sn;
 	
  theta = acos(cs);

	//fprintf(logfile,"i = %d;j = %d;k = %d;sn = %lf:: \t",i,j,k,sn);	
  
  del_theta = theta - th_theta; 
  
  tk = k3body*del_theta;

  if (eflag) eng = tk*del_theta;
  
  a = -2.0 * tk * sn;
  a11 = a*cs / rsq1;
  a12 = -a / (r1*r2);
  a22 = a*cs / rsq2;
  
  fj[0] = a11*delr1[0] + a12*delr2[0];
  fj[1] = a11*delr1[1] + a12*delr2[1];
  fj[2] = a11*delr1[2] + a12*delr2[2];
  fk[0] = a22*delr2[0] + a12*delr1[0];
  fk[1] = a22*delr2[1] + a12*delr1[1];
  fk[2] = a22*delr2[2] + a12*delr1[2];

}

/* ---------------------------------------------------------------------- */

void *PairCoreShell::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str,"cut_coul") == 0) return (void *) &cut_coul;
  return NULL;
}
