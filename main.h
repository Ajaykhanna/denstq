#define ang2au 1.889725989;
#define au2ang 0.529177249;
#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
using namespace std;

/*****************************
 * Some variables
 * ***************************/
int natoms;
const double thresh=1.e-10;
int nelec; //total number of electrons in molecule
string statetype; //type of density to use
int nx,ny,nz;
double dx1,dx2,dx3,dy1,dy2,dy3,dz1,dz2,dz3;
double ox,oy,oz;
double ***dens, ***densrad;
double *posx,*posy,*posz;
string calctype,inputtype,donorfile,accfile;
ifstream comfile;
ofstream outfile;
const string str_Atom[] = { " X",
                            " H","He","Li","Be"," B"," C"," N"," O"," F","Ne",
                            "Na","Mg","Al","Si"," P"," S","Cl","Ar"," K","Ca",
                            "Sc","Ti"," V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
                            "Ga","Ge","As","Se","Br","Kr","Rb","Sr"," Y","Zr",
                            "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn",
                            "Sb","Te"," I","Xe","Cs","Ba","La","Ce","Pr","Nd",
														"Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb",
														"Lu","Hf","Ta"," W","Re","Os","Ir","Pt","Au","Hg",
														"Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th",
														"Pa"," U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm",
														"Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds",
														"Rg","Cn","Uut","Fl","Uup","Lv","Uus","Uuo",
                          };
class Dens {
public:
  double x,y,z;
  int atom;
  int index;
  bool xmax,ymax,zmax;
  bool xmin,ymin,zmin;
  bool surface;
  double charge;

  Dens operator= (const Dens& d);
};
/** Density copy operator **/
Dens Dens::operator=(const Dens& d) {
  this->x = d.x;
  this->y = d.y;
  this->z = d.z;
  this->atom = d.atom;
  this->index = d.index;
  this->surface = d.surface;
  this->charge = d.charge;
  return *this;
}

class Atom {
  public:
  int num;
  int atomicnum;
  string type;
  double x,y,z;
  double r;
  double dens;
  double charge;
  int denselements;
  int surfaceelements;
  double dist(int,int,int);
  double distance(Dens d);
  Dens *thisdens;
  int current;

  Atom() {dens = 0.;denselements=0;current=0;};
};

double Atom::dist(int x,int y,int z) {
  double dx = this->x-x;
  double dy = this->y-y;
  double dz = this->z-z;
  double r = sqrt(dx*dx+dy+dy+dz*dz);
  return r;
}

double Atom::distance(Dens d) {
  double dx = this->x-d.x;
  double dy = this->y-d.y;
  double dz = this->z-d.z;
  double r = sqrt(dx*dx+dy+dy+dz*dz);
  return r;
}

class Molecule {
  public:
  double *posx,*posy,*posz;
  Atom *atoms;
  int natoms;
  int nx,ny,nz;
  double dx1,dx2,dx3,dy1,dy2,dy3,dz1,dz2,dz3;
  double ox,oy,oz;
  double ***dens, ***densrad;

  void setnatoms(int n) {natoms = n;};
  void setxyz(int x, int y, int z) {nx=x;ny=y;nz=z;};
  
  Molecule();
};

Molecule::Molecule() {
 }

void projectdens(const int natoms, Atom *atoms, double ***dens, 
                double posx[], double posy[], double posz[], 
                int nx, int ny, int nz, Dens *density,
                int nelec, string type);
void calcVoronoi(int natoms,Dens *density,Atom *atoms);
Atom *collectDens(Atom *atoms,ifstream &infile);
Atom *collectDens(Molecule *mol, Atom *atoms,ifstream &infile);
double computeCoupling(Molecule *donor, Molecule *acceptor);
void calcdip(Atom *atoms); 
