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

#ifdef PAIR_CLASS

PairStyle(core/shell,PairCoreShell)

#else

#ifndef LMP_PAIR_CORE_SHELL_H
#define LMP_PAIR_CORE_SHELL_H

#include "pair.h"

namespace LAMMPS_NS {

class PairCoreShell : public Pair {
 public:
  PairCoreShell(class LAMMPS *);
  virtual ~PairCoreShell();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  void write_data(FILE *);
  void write_data_all(FILE *);
  void threebody(double,double,double *,double *,double *,double *,int,double &);

  virtual void *extract(const char *, int &);

 protected:
  double cut_lj_global;
  double **cut_lj,**cut_ljsq;
  double cut_coul,cut_coulsq;
  double **a,**rho,**c;
  double **rhoinv,**buck1,**buck2,**offset;
  double *cut_respa;
  double g_ewald;
  double **kspring,**kspring4; // Spring Constant for Harmonic Spring Potential (eVA**-2)
  double **cut_spring,**cut_springsq;
  double **cut_3body,**cut_3bodysq;
  int si_tag;
  double k3body;
  double alf;	  
  double th_theta;
  
  double *q_ion;
   	
  void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

E: Pair style buck/coul/long requires atom attribute q

The atom style defined does not have these attributes.

E: Pair style requires a KSpace style

No kspace style is defined.

*/