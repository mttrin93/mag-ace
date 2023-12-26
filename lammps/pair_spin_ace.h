/* -*- c++ -*- ----------------------------------------------------------
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

PairStyle(spin/ace,PairSpinACE)

#else


#ifndef LMP_PAIR_SPIN_ACE_H
#define LMP_PAIR_SPIN_ACE_H

#include "pair_spin.h"  // IWYU pragma: export

namespace LAMMPS_NS {

class PairSpinACE : public PairSpin {
friend class FixNVESpin;
 public:
  PairSpinACE(LAMMPS *);
  virtual ~PairSpinACE();
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void *extract(const char *, int &);

  void compute(int, int);
  void compute_single_pair(int, double *);

  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);

 protected:
  // Atomic Cluster Expansion parameters
  double rcut;
  double rcut_magnetic;
  char xcfunctional[3];
  double smear;
  int nele;
  char** elelist;
  int* elenumbers;

// Param *params;                // parameter set for an I-J-K interaction
  char **elements;              // names of unique elements
  int ***elem2param;            // mapping from element triplets to parameters
  int *map;                     // mapping from atom types to elements
  double cutmax;                // max cutoff for all elements
  int nelements;                // # of unique elements
  int nparams;                  // # of stored parameter sets
  int maxparam;                 // max # of parameter sets
  int maxshort;                 // size of short neighbor list array
  int *neighshort;              // short neighbor list array

  char *potential_file_name;
  void allocate();
  
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Incorrect args in pair_spin command

Self-explanatory.

E: Spin simulations require metal unit style

Self-explanatory.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair spin requires atom attribute spin

The atom style defined does not have these attributes.

*/
