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

//#ifdef PAIR_CLASS

//PairStyle(tersoff,PairTersoff)

//#else

//#ifndef LMP_PAIR_TERSOFF_H
//#define LMP_PAIR_TERSOFF_H

//#include "pair.h"

//namespace LAMMPS_NS {

  //class PairTersoff : public Pair {
  //public:
  //PairTersoff(class LAMMPS *);
  //virtual ~PairTersoff();
  //virtual void compute(int, int);
  //void settings(int, char **);
  //void coeff(int, char **);
  //void init_style();
  //double init_one(int, int);

//protected:
/*
  struct Param {
    double lam1,lam2,lam3;
    double c,d,h;
    double gamma,powerm;
    double powern,beta;
    double biga,bigb,bigd,bigr;
    double cut,cutsq;
    double c1,c2,c3,c4;
    int ielement,jelement,kelement;
    int powermint;
    double Z_i,Z_j;              // added for TersoffZBL
    double ZBLcut,ZBLexpscale;
    double c5,ca1,ca4;           // added for TersoffMOD
    double powern_del;
  };
*/
// Param *params;                // parameter set for an I-J-K interaction
//char **elements;              // names of unique elements
//  int ***elem2param;            // mapping from element triplets to parameters
//int *map;                     // mapping from atom types to elements
//double cutmax;                // max cutoff for all elements
//int nelements;                // # of unique elements
// int nparams;                  // # of stored parameter sets
//int maxparam;                 // max # of parameter sets

// void allocate();
//virtual void read_file(char *);
//virtual void setup();
 static const double THIRD  = 1.0/3.0;
  static const double MY_PI  = 3.14159265358979323846; // pi
  static const double MY_2PI = 6.28318530717958647692; // 2pi
  static const double MY_3PI = 9.42477796076937971538; // 3pi
  static const double MY_4PI = 12.56637061435917295384; // 4pi
  static const double MY_PI2 = 1.57079632679489661923; // pi/2
  static const double MY_PI4 = 0.78539816339744830962; // pi/4
  static const double MY_PIS = 1.77245385090551602729; // sqrt(pi)
  static const double MY_SQRT2 = 1.41421356237309504880; // sqrt(2)
  static const double MY_CBRT2 = 1.25992104989487316476; // 2*(1/3)

//virtual void repulsive(Param *, double, double &, int, double &);
void repulsive(double rsq, double &fforce, double &eng,double lam1,double biga,double bigr,double bigd);
//virtual double zeta(Param *, double, double, double *, double *);
double zeta(double rsqij, double rsqik,double *delrij, double *delrik,int powermint,double lam1,double lam3,double c,double d,double h,double gamma,double bigr,double bigd);
//virtual void force_zeta(Param *, double, double, double &,double &, int, double &);
void force_zeta(double rsq, double zeta_ij,double &fforce, double &prefactor, double &eng,double c1,double c2,double c3,double c4,double powern,double beta,double bigr,double bigb,double bigd,double lam2);
void attractive(double prefactor,double rsqij, double rsqik,double *delrij, double *delrik,double *fi, double *fj, double *fk,int powermint,double lam3,double c,double d,double h,double gamma,double bigr,double bigd);
//void attractive(Param *, double, double, double, double *, double *,double *, double *, double *);
double ters_fc(double r,double bigr,double bigd);
//virtual double ters_fc(double, Param *);
double ters_fc_d(double r,double bigr,double bigd);
//virtual double ters_fc_d(double, Param *);
double ters_fa(double r,double bigr,double bigb,double bigd,double lam2);
//virtual double ters_fa(double, Param *);
double ters_fa_d(double r,double bigr,double bigd,double lam2,double bigb);
//virtual double ters_fa_d(double, Param *);
double ters_bij(double zeta, double c1,double c2,double c3,double c4,double powern,double beta);
//virtual double ters_bij(double, Param *);
double ters_bij_d(double zeta,double c1,double c2,double c3,double c4,double beta,double powern);
//virtual double ters_bij_d(double, Param *);
void ters_zetaterm_d(double prefactor,double *rij_hat, double rij,double *rik_hat, double rik,double *dri, double *drj, double *drk,int powermint,double lam3,double c,double d,double h,double gamma,double bigr,double bigd);
//virtual void ters_zetaterm_d(double, double *, double, double *, double,double *, double *, double *, Param *);
void costheta_d(double *rij_hat, double rij,double *rik_hat, double rik,double *dri, double *drj, double *drk);
void force_zeta_pot(double rsq, double zeta_ij,double &prefactor, double &eng,double c1,double c2,double c3,double c4,double powern,double beta,double bigr,double bigb,double bigd,double lam2);
void repulsive_pot(double rsq, double &eng,double lam1,double biga,double bigr,double bigd);
void calc_E_Tersoff_Pot(double *angle,int *BondNum,int Natoms,point L,BondParamsStruct *BondParams,int *Bonded,atom *particle,double *Fc,int *NBond,double &E_Tersoff_,int *NBond_once,int *Bonded_once);
//void costheta_d(double *, double, double *, double,double *, double *, double *);

  // inlined functions for efficiency
  
inline double ters_gijk(double costheta,double c,double d,double h,double gamma)
{
  double ters_c = c * c;
  double ters_d = d * d;
  double hcth   = h - costheta;
  
  return gamma*(1.0 + ters_c/ters_d - ters_c / (ters_d + hcth*hcth));
}

inline double ters_gijk_d( double costheta, double c, double d, double h, double gamma)  {
  double ters_c = c * c;
  double ters_d = d * d;
  double hcth = h - costheta;
  double numerator = -2.0 * ters_c * hcth;
  double denominator = 1.0/(ters_d + hcth*hcth);
  return gamma*numerator*denominator*denominator;
}

inline double vec3_dot( double x[3], double y[3])  {
  return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

inline void vec3_add( double x[3],  double y[3],
		      double *  z)  {
  z[0] = x[0]+y[0];  z[1] = x[1]+y[1];  z[2] = x[2]+y[2];
}

inline void vec3_scale( double k,  double x[3],double y[3])  {
  y[0] = k*x[0];  y[1] = k*x[1];  y[2] = k*x[2];
}

inline void vec3_scaleadd( double k,  double x[3],double y[3], double *z)  {
  z[0] = k*x[0]+y[0];
  z[1] = k*x[1]+y[1];
  z[2] = k*x[2]+y[2];
}
//};



//#endif
//#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair style Tersoff requires atom IDs

This is a requirement to use the Tersoff potential.

E: Pair style Tersoff requires newton pair on

See the newton command.  This is a restriction to use the Tersoff
potential.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

E: Cannot open Tersoff potential file %s

The specified potential file cannot be opened.  Check that the path
and name are correct.

E: Incorrect format in Tersoff potential file

Incorrect number of words per line in the potential file.

E: Illegal Tersoff parameter

One or more of the coefficients defined in the potential file is
invalid.

E: Potential file has duplicate entry

The potential file for a SW or Tersoff potential has more than
one entry for the same 3 ordered elements.

E: Potential file is missing an entry

The potential file for a SW or Tersoff potential does not have a
needed entry.

*/
