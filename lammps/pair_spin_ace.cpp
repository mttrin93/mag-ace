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

#include "pair_spin_ace.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "atom.h"
#include "comm.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "math_const.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "pair.h"
#include "update.h"
#include "fix_nve_spin.h"
#include "utils.h"
#include "memory.h"
#include "neigh_list.h"

using namespace LAMMPS_NS;
using namespace MathConst;

int elements_num = 104;
char const * const elements[104] = {"X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na",
                                    "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn",
                                    "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr",
                                    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb",
                                    "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
                                    "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir",
                                    "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
                                    "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"
};
int AtomicNumberByName(char* elname)
{
    int i =0;
    for(i=1; i<elements_num; i++)
        if(strcmp(elname,elements[i])==0)
            return i;
    return -1;
}

const double mub = 5.78901e-5;

// acelibini_(&nele,&rcut,&rcut_magnetic,&smear,&xcnumber,elenumbers);
extern "C" void acelibini_(int *, double *, double *, double *, int *, int *);
// acelib_(&natom,&nall,&nei,nlist,occ,occlist,nstart,nstop,eout,dlist,elist,flist,dlist_momi,elist_momi,dlist_mom,elist_mom,tlist,tlist0);
extern "C" void acelib_(int *,int *, int *,int *,int *, int *, int *, int *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);

/* ---------------------------------------------------------------------- */

// PairSpinACE::PairSpinACE(LAMMPS *lmp) : PairSpin(lmp) {
PairSpinACE::PairSpinACE(LAMMPS *lmp) : PairSpin(lmp){
  hbar = force->hplanck/MY_2PI;
  single_enable = 0;
  respa_enable = 0;
  no_virial_fdotr_compute = 1;
  lattice_flag = 0;

  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;

  nelements = 0;
  elements = NULL;
  nparams = maxparam = 0;
  //    params = NULL;
  elem2param = NULL;
  map = NULL;

  maxshort = 10;
  neighshort = NULL;
}

/* ---------------------------------------------------------------------- */

PairSpinACE::~PairSpinACE() {
    if (copymode) return;

    if (elements)
        for (int i = 0; i < nelements; i++) delete[] elements[i];
    delete[] elements;

    if (elelist)
        for (int i = 0; i < nele; i++) delete[] elelist[i];
    delete[] elelist;

    //    memory->destroy(params);
    memory->destroy(elem2param);

    if (allocated) {
        memory->destroy(setflag);
        memory->destroy(cutsq);
        memory->destroy(neighshort);
        delete[] map;
    }    
}

/* ---------------------------------------------------------------------- */

void PairSpinACE::allocate() {
  //   printf("--> PairACE::allocate\n");
   allocated = 1;
    int n = atom->ntypes;

    memory->create(setflag, n + 1, n + 1, "pair:setflag");
    memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
    memory->create(neighshort, maxshort, "pair:neighshort");
    map = new int[n + 1];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairSpinACE::settings(int narg, char **arg)
{
  //   printf("--> PairACE::settings\n");
    //pair_style 0:rcut 1:xcfunctional 2:smear 3..:elelist[1..nele]

    if (narg < 4) error->all(FLERR, "Illegal pair_style command. Correct form:\npair_style spin/ace nelements[1] rcut[4.5] rcut_magnetic[4.0] xcfunctional[pbe/lda] smear[0.1] ");

    // pair spin/* need the metal unit style
    if (strcmp(update->unit_style,"metal") != 0)
        error->all(FLERR,"Spin pair styles require metal units");
    
    nele = atof(arg[0]);
    printf("--> PairSpinACE:: nele=%d\n",nele);
    rcut = atof(arg[1]);
    printf("--> PairSpinACE:: rcut=%f\n",rcut);
    rcut_magnetic = atof(arg[2]);
    printf("--> PairSpinACE:: rcut_magnetic=%f\n",rcut_magnetic);
    strncpy(xcfunctional,arg[3],3);
    xcfunctional[3] = '\0';
    printf("--> PairSpinACE:: functional=%s\n",arg[3]);
    smear = atof(arg[4]);
    printf("--> PairSpinACE:: smear=%f\n",smear);

    cutmax = rcut;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairSpinACE::coeff(int narg, char **arg) {
     printf("--> PairSpinACE::coeff\n");
     printf("narg = %d\n",narg);
    for(int i = 0; i<narg; i++){
        printf("PairSpinACE::coeff: arg[%d] = %s\n",i,arg[i]);
    }
    //     printf("nele = %d\n",nele);
    // printf("nelements = %d\n",nelements);

    
    int i, j, n;
    int nxc;
    char buf[1024];

    if (!allocated) allocate();

    // read potential file and initialize potential parameters
//    int pot_file_name_len = (unsigned)strlen(arg[2]) + 1;
//    potential_file_name = new char[pot_file_name_len];
//    strcpy(potential_file_name, arg[2]);
//    printf("potential_file_name=%s\n",potential_file_name);

    int shift = 3;

    if (narg != shift + atom->ntypes)
        error->all(FLERR, "Incorrect args for pair coefficients");


    n =  narg - shift; //atom->ntypes; //atoi(arg[3]);

    if (n != nele)
      error->all(FLERR, "Incorrect number of atoms for coeff");
    
    elelist = new char *[nele];
    elenumbers = new int[nele];
    for (int i = 0; i<nele; i++){
      n = strlen(arg[shift+i]) + 1;
      elelist[i] = new char[n];
      strcpy(elelist[i], arg[i+shift]);
      
      int atomicnumber=AtomicNumberByName(arg[i+shift]);
      if (atomicnumber==-1){
	sprintf(buf,"%s is not a valid element name",arg[i+shift]);
	error->all(FLERR,buf);
      }
      elenumbers[i]=atomicnumber;
      printf("element: %s -> %d\n",arg[i+shift],elenumbers[i]);
    }
    

    // insure I,J args are * *
    if (strcmp(arg[0], "*") != 0 || strcmp(arg[1], "*") != 0)
        error->all(FLERR, "Incorrect args for pair coefficients");

    // read args that map atom types to elements in potential file
    // map[i] = which element the Ith atom type is, -1 if NULL
    // nelements = # of unique elements
    // elements = list of element names

    if (elelist) {
        for (i = 0; i < nelements; i++) delete[] elelist[i];
        delete[] elelist;
    }
    elelist = new char *[atom->ntypes];
    for (i = 0; i < atom->ntypes; i++) elelist[i] = NULL;

    nelements = 0;
    for (i = shift; i < narg; i++) {
        if (strcmp(arg[i], "NULL") == 0) {
            map[i - shift+1] = -1;
            continue;
        }
        for (j = 0; j < nelements; j++)
            if (strcmp(arg[i], elelist[j]) == 0) break;
        map[i - shift+1] = j;
        if (j == nelements) {
            n = strlen(arg[i]) + 1;
            elelist[j] = new char[n];
            strcpy(elelist[j], arg[i]);
            nelements++;
        }
    }

    for (int i = 0; i<nele; i++){
        printf("elelist[%d]=%s\n",i,elelist[i]);
    }

    // clear setflag since coeff() called once with I,J = * *

    n = atom->ntypes;
    for (i = 1; i <= n; i++)
        for (j = i; j <= n; j++)
            setflag[i][j] = 0;

    // set setflag i,j for type pairs where both are mapped to elements

    int count = 0;
    for (i = 1; i <= n; i++)
        for (j = i; j <= n; j++)
            if (map[i] >= 0 && map[j] >= 0) {
                setflag[i][j] = 1;
                count++;
            }

    if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");

    printf("nelements = %d\n",nelements);

    printf("xcfunctional=%s\n", xcfunctional);
    nxc = -1;
    if (strcmp(xcfunctional, "new") == 0) {
        nxc = 0;//generate reference input file only
    }    
    if (strcmp(xcfunctional, "pbe") == 0) {
        nxc = 1;//PBE
    }
    if (strcmp(xcfunctional, "lda") == 0) {
        nxc = 2; //LDA
    }
    if (nxc == -1) error->all(FLERR, "Exchange functional not recognized");
    else printf("Reference exchange functional: %s (=%d)\n", xcfunctional, nxc);

    // up 1 for index in Fortran style
    for (i = 0; i < nele; i++){
      elenumbers[i] = elenumbers[i] + 1;
    }
   
    // here everything is in place for ace call I: ini
    acelibini_(&nele, &rcut, &rcut_magnetic, &smear, &nxc, elenumbers);
    // printf("acelibini_done\n");
    
}


/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairSpinACE::init_style()
{
  PairSpin::init_style();

    if (atom->tag_enable == 0)
        error->all(FLERR, "Pair style ACE requires atom IDs");
    if (force->newton_pair == 0)
        error->all(FLERR, "Pair style ACE requires newton pair on");
    
    // need a full neighbor list

    int irequest = neighbor->request(this, instance_me);
    neighbor->requests[irequest]->half = 0;
    neighbor->requests[irequest]->full = 1;

}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairSpinACE::init_one(int i, int j) {
    if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");

    return cutmax;
}


/* ----------------------------------------------------------------------
   extract the larger cutoff
------------------------------------------------------------------------- */

void *PairSpinACE::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str,"cut") == 0) return (void *) &cutmax;
  return NULL;
}


/* ---------------------------------------------------------------------- */

void PairSpinACE::compute(int eflag, int vflag)
{
    int i,j,ii,jj,inum,jnum,itype,jtype;
    double evdwl, ecoul;
    double xi[3], eij[3], mij[3];
    double delx,dely,delz;
    double spi[3], spj[3];
    double fi[3], fmi[3];
    double local_cut2;
    double rsq, inorm;
    double msq, mnorm;
    int *ilist,*jlist,*numneigh,**firstneigh;
    int k,l;

    evdwl = 0.0;
    ev_init(eflag,vflag);

    double **x = atom->x;
    double **f = atom->f;
    double **fm = atom->fm;
    double **sp = atom->sp;
    int *type = atom->type;
        // number of atoms in cell
    int nlocal = atom->nlocal;
    int newton_pair = force->newton_pair;
    const double cutshortsq = cutmax * cutmax;

    // inum: length of the neighbrlist
    inum = list->inum;
    // ilist: list of "i" atoms for which neighbor lists exist
    ilist = list->ilist;
    //numneigh: the length of each these neigbor list
    numneigh = list->numneigh;
    // the pointer to the list of neighbors of "i"
    firstneigh = list->firstneigh;
  
    //MR new stuff
    // number of atoms including ghost atoms
    int nall = nlocal + atom->nghost;
    int n;
    int nele;
    int neicount;
    int nei;
    int natom;
    int *nlist;
    int *nstart;
    int *nstop;
    int *occlist;
    int *occ;
    double *dlist;
    double *elist;
    double *eout;
    double *flist;
    double *dlist_mom;
    double *elist_mom;
    double *dlist_momi;
    double *elist_momi;
    double *tlist;
    double *tlist0;

    //parameters read from pair_style:
    // pair_style ace nelements [1] rcut[4.5] rcut_magnetic[4.0] xcfunctional[pbe/lda] smear[0.1] nele[2] elelist[Al Co]
    // pair_style ace 4.5 4.0 pbe 0.1 1 Al

    //double rcut = 4.5;
    //double rcut_magnetic = 4.0;
    //double smear = 0.1;
    //char elelist[2];
    //char xcfunctional[]="pbe";

    double fxtmp, fytmp, fztmp;  
    double txtmp, tytmp, tztmp;
    double txtmp0, tytmp0, tztmp0;

    if (inum != nlocal) {
        printf("inum: %d nlocal: %d are different.\n", inum, nlocal);
        exit(0);
    }
    // find number of neighbors and total number of atoms, including ghost atoms
    natom = nlocal;
    nei = 0;
    for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];
        jnum = numneigh[i];
        nei = nei + jnum;
    }   
    
    // dynamical allocation for my neighborlist 
    nstart = new int [nlocal];
    nstop = new int [nlocal];
    occ = new int [nlocal];
    nlist = new int [nei];
    occlist = new int [nei];
    dlist = new double [nei];
    elist = new double [3*nei];
    flist = new double [3*nei];
    eout = new double [nlocal];
    tlist = new double [3*nei];
    tlist0 = new double [3*nlocal];
    dlist_mom = new double [nei];
    elist_mom = new double [3*nei];
    dlist_momi = new double [nlocal];
    elist_momi = new double [3*nlocal];
        
    // spin cluster expansion computation
    // loop over atoms and their neighbors
    neicount = 0;
    for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];
        k = map[type[i]];  
        occ[ii] = k;
        
        if (occ[ii] != 0) {
            printf("Error: Cannot handle multiple species.\n");
            printf("Stopping.\n");
            exit(0);
        }        
        nstart[ii] = neicount;        
//         itype = type[i];

        jlist = firstneigh[i];
        jnum = numneigh[i];
        xi[0] = x[i][0];
        xi[1] = x[i][1];
        xi[2] = x[i][2];
        spi[0] = sp[i][0];    
        spi[1] = sp[i][1];
        spi[2] = sp[i][2];
        
        dlist_momi[ii] = sp[i][3];
        elist_momi[3 * ii + 0] = spi[0];
        elist_momi[3 * ii + 1] = spi[1];
        elist_momi[3 * ii + 2] = spi[2];
                
        // loop on neighbors

        for (jj = 0; jj < jnum; jj++) {
            j = jlist[jj];
            j &= NEIGHMASK;
//             jtype = type[j];
            
            if (j >= nall) {
                printf("Error: Atoms index too large.\n");
                printf("Stopping.\n");
                exit(0);
            }

            nlist[neicount] = j;
            k = map[type[j]];
            occlist[neicount] = k;            
            
            delx = xi[0] - x[j][0];
            dely = xi[1] - x[j][1];
            delz = xi[2] - x[j][2];
            rsq = delx*delx + dely*dely + delz*delz;
            inorm = 1.0/sqrt(rsq);

            eij[0] = inorm*delx;
            eij[1] = inorm*dely;
            eij[2] = inorm*delz;
            // store distances
            dlist[neicount] = sqrt(rsq);
            if (sqrt(rsq) < 1.e-6) {
                printf("Error: Atoms with identical positions.\n");
                printf("Stopping.\n");
                exit(0);
            }
            // store direction vectors
            elist[3 * neicount + 0] = eij[0];
            elist[3 * neicount + 1] = eij[1];
            elist[3 * neicount + 2] = eij[2];
                        
            mij[0] = sp[j][0];
            mij[1] = sp[j][1];
            mij[2] = sp[j][2];
            //store magnitude moments
            dlist_mom[neicount] = sp[j][3];
            //store direction moments
            elist_mom[3 * neicount + 0] = mij[0];
            elist_mom[3 * neicount + 1] = mij[1];
            elist_mom[3 * neicount + 2] = mij[2];

            neicount++;  
        }
        nstop[ii] = neicount - 1;
    }

    if (nei != neicount) {
        printf("Something wrong with neighbor list translation: neicount\n");
        printf("Stopping.\n");
        exit(0);
    }
//    printf("neicount=%d\n",neicount);

    // up by +1 to run in Fortran style
    for (ii = 0; ii < neicount; ii++) {
        nlist[ii] = nlist[ii] + 1;
        occlist[ii] = occlist[ii] + 1;
    }
    for (ii = 0; ii < natom; ii++) {
        nstart[ii] = nstart[ii] + 1;
        nstop[ii] = nstop[ii] + 1;
        occ[ii] = occ[ii] + 1;
    }

    acelib_(&natom,&nall,&nei,nlist,occ,occlist,nstart,nstop,eout,dlist,elist,flist,dlist_momi,elist_momi,dlist_mom,elist_mom,tlist,tlist0);
    
    // check if this is necessary
    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      f[i][0] = 0.0;
      f[i][1] = 0.0;
      f[i][2] = 0.0;
      fm[i][0] = 0.0;
      fm[i][1] = 0.0;
      fm[i][2] = 0.0;
      for (jj = 0; jj < jnum; jj++) {
            jlist = firstneigh[i];
            j = jlist[jj];
            j &= NEIGHMASK;
            f[j][0] = 0.0;
            f[j][1] = 0.0;
            f[j][2] = 0.0;
            fm[j][0] = 0.0;
            fm[j][1] = 0.0;
            fm[j][2] = 0.0;
      }
    }
    
    
    // extract now the magnetic gradients that are used in the spin dynamics calculation
    k = 0;
    l = 0;
    for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];
        // xtmp = x[i][0];
        // ytmp = x[i][1];
        //ztmp = x[i][2];
        jlist = firstneigh[i];
        jnum = numneigh[i];
        
        xi[0] = x[i][0];
        xi[1] = x[i][1];
        xi[2] = x[i][2];
        
        txtmp0 = tytmp0 = tztmp0 = 0.0;
        
        txtmp0 = tlist0[l+0];
        tytmp0 = tlist0[l+1];
        tztmp0 = tlist0[l+2]; 
        
        fm[i][0] -= txtmp0;  
        fm[i][1] -= tytmp0;
        fm[i][2] -= tztmp0;
        
        l=l+3;
        
        for (jj = 0; jj < jnum; jj++) {
            j = jlist[jj];
            j &= NEIGHMASK;
                        
            delx = xi[0] - x[j][0];
            dely = xi[1] - x[j][1];
            delz = xi[2] - x[j][2];

            fxtmp = fytmp = fztmp = 0.0; 
            txtmp = tytmp = tztmp = 0.0;
            
            fxtmp = flist[k+0];
            fytmp = flist[k+1];
            fztmp = flist[k+2];
            
            txtmp=  tlist[k+0];
            tytmp = tlist[k+1];
            tztmp = tlist[k+2];
            
            //if(fxtmp!=0.0 || fytmp!=0.0 || fztmp!=0.0){
	      //                printf("Pair force (%d-%d) = (%f, %f, %f)\n",i,j,fxtmp,fytmp,fztmp);
            //}
            f[i][0] -= fxtmp;
            f[i][1] -= fytmp;
            f[i][2] -= fztmp;
            
            if (newton_pair || j < nlocal) {
                f[j][0] += fxtmp;
                f[j][1] += fytmp;
                f[j][2] += fztmp;
                
                fm[j][0] -= txtmp;  //se metti fm[i] legge tutto bene nel trajectory file
                fm[j][1] -= tytmp;
                fm[j][2] -= tztmp;
            }
                                
            k=k+3;

            if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,
                                 0.0,0.0,fxtmp, fytmp, fztmp,delx,dely,delz);
        }

        if(eflag){
            ev_tally_full(i,2.0*eout[i],0.0,0.0,0.0,0.0,0.0);
        }
//         printf("magnetic eout: %f\n",eout[i]);

	// note: incorrect force print out as still inside of loop
	//        printf("Total force II on at.%d=(%f,%f,%f)\n",i,f[i][0],f[i][1],f[i][2]);
    }

    
    delete [] nstart;
    delete [] nstop;
    delete [] occ;
    delete [] nlist;
    delete [] occlist;
    delete [] dlist;
    delete [] elist;
    delete [] flist;
    delete [] eout;
    delete [] tlist;
    delete [] tlist0;
    delete [] dlist_mom;
    delete [] elist_mom;
    delete [] dlist_momi;
    delete [] elist_momi;
    
    
//    printf("vflag_fdotr=%d\n",vflag_fdotr);
    if (vflag_fdotr) virial_fdotr_compute();
     
}



/* ----------------------------------------------------------------------
   update the pair interactions fmi acting on the spin ii
------------------------------------------------------------------------- */

void PairSpinACE::compute_single_pair(int iii, double fmi[3])
{
    int i,j,ii,jj,inum,jnum,itype,jtype,ntypes;;
    double evdwl, ecoul;
    double xi[3], eij[3], mij[3];
    double delx,dely,delz;
    double spi[3], spj[3];
    double rsq, inorm;
    double msq, mnorm;
    int *ilist,*jlist,*numneigh,**firstneigh;
    int k,locflag,l;

    double **x = atom->x;
    double **sp = atom->sp;
    int *type = atom->type;
        // number of atoms in cell
    int nlocal = atom->nlocal;
    
    // inum: length of the neighbrlist
    inum = list->inum;
    // ilist: list of "i" atoms for which neighbor lists exist
    ilist = list->ilist;
    //numneigh: the length of each these neigbor list
    numneigh = list->numneigh;
    // the pointer to the list of neighbors of "i"
    firstneigh = list->firstneigh;
  
    //MR new stuff
    // number of atoms including ghost atoms
    int nall = nlocal + atom->nghost;
    int n;
    int nele;
    int neicount;
    int nei;
    int natom;
    int *nlist;
    int *nstart;
    int *nstop;
    int *occlist;
    int *occ;
    double *dlist;
    double *elist;
    double *eout;
    double *flist;
    double *dlist_mom;
    double *elist_mom;
    double *dlist_momi;
    double *elist_momi;
    double *tlist;
    double *tlist0;

    //parameters read from pair_style:
    // pair_style ace nelements [1] rcut[8.7] xcfunctional[pbe/lda] smear[0.1] nele[2] elelist[Al Co]
    // pair_style ace 8.7 pbe 0.1 1 Al

    //double rcut = 8.7;
    //double smear = 0.1;
    //char elelist[2];
    //char xcfunctional[]="pbe";

    double txtmp, tytmp, tztmp;
    double txtmp0, tytmp0, tztmp0;

    if (inum != nlocal) {
        printf("inum: %d nlocal: %d are different.\n", inum, nlocal);
        exit(0);
    }
    // find number of neighbors and total number of atoms, including ghost atoms
    natom = nlocal;
    nei = 0;
    for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];
        jnum = numneigh[i];
        nei = nei + jnum;
    }   
    
    // dynamical allocation for my neighborlist 
    nstart = new int [nlocal];
    nstop = new int [nlocal];
    occ = new int [nlocal];
    nlist = new int [nei];
    occlist = new int [nei];
    dlist = new double [nei];
    elist = new double [3*nei];
    flist = new double [3*nei];
    eout = new double [nlocal];
    tlist = new double [3*nei];
    tlist0 = new double [3*nlocal];
    dlist_mom = new double [nei];
    elist_mom = new double [3*nei];
    dlist_momi = new double [nlocal];
    elist_momi = new double [3*nlocal];
    
    
    // check if interaction applies to type of ii

    itype = type[iii];
    ntypes = atom->ntypes;
    locflag = 0;
    k = 1;
    while (k <= ntypes) {
        if (k <= itype) {
        if (setflag[k][itype] == 1) {
            locflag =1;
            break;
        }
        k++;
        } else if (k > itype) {
        if (setflag[itype][k] == 1) {
            locflag =1;
            break;
        }
        k++;
        } else error->all(FLERR,"Wrong type number");
    }
            
    // spin cluster expansion computation
    // loop over atoms and their neighbors
    neicount = 0;
    for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];
        k = map[type[i]];  
        occ[ii] = k;
        
        if (occ[ii] != 0) {
            printf("Error: Cannot handle multiple species.\n");
            printf("Stopping.\n");
            exit(0);
        }        
        nstart[ii] = neicount;        
//         itype = type[i];

        jlist = firstneigh[i];
        jnum = numneigh[i];
        xi[0] = x[i][0];
        xi[1] = x[i][1];
        xi[2] = x[i][2];
        spi[0] = sp[i][0];    
        spi[1] = sp[i][1];
        spi[2] = sp[i][2];
        
        dlist_momi[ii] = sp[i][3];
//         mnorm = 1.0/sp[i][3];
//         elist_momi[3 * ii + 0] = mnorm*spi[0];
//         elist_momi[3 * ii + 1] = mnorm*spi[1];
//         elist_momi[3 * ii + 2] = mnorm*spi[2];
        elist_momi[3 * ii + 0] = spi[0];
        elist_momi[3 * ii + 1] = spi[1];
        elist_momi[3 * ii + 2] = spi[2];
        
        // loop on neighbors

        for (jj = 0; jj < jnum; jj++) {
            j = jlist[jj];
            j &= NEIGHMASK;
//             jtype = type[j];
            
            if (j >= nall) {
                printf("Error: Atoms index too large.\n");
                printf("Stopping.\n");
                exit(0);
            }

            nlist[neicount] = j;
            k = map[type[j]];
            occlist[neicount] = k;            
            
            delx = xi[0] - x[j][0];
            dely = xi[1] - x[j][1];
            delz = xi[2] - x[j][2];
            rsq = delx*delx + dely*dely + delz*delz;
            inorm = 1.0/sqrt(rsq);
//             eij[0] = -inorm*delx;
//             eij[1] = -inorm*dely;
//             eij[2] = -inorm*delz;
            eij[0] = inorm*delx;
            eij[1] = inorm*dely;
            eij[2] = inorm*delz;
            // store distances
            dlist[neicount] = sqrt(rsq);
            if (sqrt(rsq) < 1.e-6) {
                printf("Error: Atoms with identical positions.\n");
                printf("Stopping.\n");
                exit(0);
            }
            // store direction vectors
            elist[3 * neicount + 0] = eij[0];
            elist[3 * neicount + 1] = eij[1];
            elist[3 * neicount + 2] = eij[2];
            
//             mnorm = 1.0/sp[j][3];
//             mij[0] = mnorm*sp[j][0];
//             mij[1] = mnorm*sp[j][1];
//             mij[2] = mnorm*sp[j][2];
            mij[0] = sp[j][0];
            mij[1] = sp[j][1];
            mij[2] = sp[j][2];
            //store magnitude moments
            dlist_mom[neicount] = sp[j][3];
            //store direction moments
            elist_mom[3 * neicount + 0] = mij[0];
            elist_mom[3 * neicount + 1] = mij[1];
            elist_mom[3 * neicount + 2] = mij[2];

            neicount++;  
        }
        nstop[ii] = neicount - 1;
    }

    if (nei != neicount) {
        printf("Something wrong with neighbor list translation: neicount\n");
        printf("Stopping.\n");
        exit(0);
    }
//    printf("neicount=%d\n",neicount);

    // up by +1 to run in Fortran style
    for (ii = 0; ii < neicount; ii++) {
        nlist[ii] = nlist[ii] + 1;
        occlist[ii] = occlist[ii] + 1;
    }
    for (ii = 0; ii < natom; ii++) {
        nstart[ii] = nstart[ii] + 1;
        nstop[ii] = nstop[ii] + 1;
        occ[ii] = occ[ii] + 1;
    }

    acelib_(&natom,&nall,&nei,nlist,occ,occlist,nstart,nstop,eout,dlist,elist,flist,dlist_momi,elist_momi,dlist_mom,elist_mom,tlist,tlist0);
    
    // extract now the magnetic gradients that are used in the spin dynamics calculation
//     k = 0;
//     l = 0;
//     for (ii = 0; ii < inum; ii++) {
//         i = ilist[ii];
        
//         if (i == iii && locflag == 1) {
            // xtmp = x[i][0];
            // ytmp = x[i][1];
            //ztmp = x[i][2];
//         jlist = firstneigh[i];
//         jnum = numneigh[i];
 
    if (locflag == 1) {
      
        k = 0;
    
        xi[0] = x[iii][0];
        xi[1] = x[iii][1];
        xi[2] = x[iii][2];
        
        jlist = firstneigh[iii];
        jnum = numneigh[iii];
        
        txtmp0 = tytmp0 = tztmp0 = 0.0;
        
        l = 3*(iii-1)+1;
        
        txtmp0 = tlist0[l+0];
        tytmp0 = tlist0[l+1];
        tztmp0 = tlist0[l+2];
                        
        fmi[0] -= (mub/hbar)*txtmp0;  
        fmi[1] -= (mub/hbar)*tytmp0;
        fmi[2] -= (mub/hbar)*tztmp0;
        
        for (jj = 0; jj < jnum; jj++) {
            j = jlist[jj];
            j &= NEIGHMASK;
                        
            delx = xi[0] - x[j][0];
            dely = xi[1] - x[j][1];
            delz = xi[2] - x[j][2];

            txtmp = tytmp = tztmp = 0.0;
            
            txtmp=  tlist[k+0];
            tytmp = tlist[k+1];
            tztmp = tlist[k+2];
                            
            fmi[0] -= (mub/hbar)*txtmp;  
            fmi[1] -= (mub/hbar)*tytmp;
            fmi[2] -= (mub/hbar)*tztmp;
            
            k=k+3;

        }

//         }
//         l=l+3;
    }

    
    delete [] nstart;
    delete [] nstop;
    delete [] occ;
    delete [] nlist;
    delete [] occlist;
    delete [] dlist;
    delete [] elist;
    delete [] flist;
    delete [] eout;
    delete [] tlist;
    delete [] tlist0;
    delete [] dlist_mom;
    delete [] elist_mom;
    delete [] dlist_momi;
    delete [] elist_momi;
    
    
}



/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairSpinACE::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&cutsq[i][j],sizeof(double),1,fp);
      }
    }
  }
}

//     memory->create(setflag, n + 1, n + 1, "pair:setflag");
//     memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
//     memory->create(neighshort, maxshort, "pair:neighshort");
//     map = new int[n + 1];
    
    
/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairSpinACE::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR,&setflag[i][j],sizeof(int),1,fp,NULL,error);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
//           utils::sfread(FLERR,&J1_mag[i][j],sizeof(double),1,fp,NULL,error);
//           utils::sfread(FLERR,&J1_mech[i][j],sizeof(double),1,fp,NULL,error);
//           utils::sfread(FLERR,&J2[i][j],sizeof(double),1,fp,NULL,error);
//           utils::sfread(FLERR,&J3[i][j],sizeof(double),1,fp,NULL,error);
          utils::sfread(FLERR,&cutsq[i][j],sizeof(double),1,fp,NULL,error);
        }
//         MPI_Bcast(&J1_mag[i][j],1,MPI_DOUBLE,0,world);
//         MPI_Bcast(&J1_mech[i][j],1,MPI_DOUBLE,0,world);
//         MPI_Bcast(&J2[i][j],1,MPI_DOUBLE,0,world);
//         MPI_Bcast(&J3[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cutsq[i][j],1,MPI_DOUBLE,0,world);
      }
    }
  }
}


/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairSpinACE::write_restart_settings(FILE *fp)
{
  fwrite(&cutmax,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairSpinACE::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR,&cutmax,sizeof(double),1,fp,NULL,error);
    utils::sfread(FLERR,&offset_flag,sizeof(int),1,fp,NULL,error);
    utils::sfread(FLERR,&mix_flag,sizeof(int),1,fp,NULL,error);
  }
  MPI_Bcast(&cutmax,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

