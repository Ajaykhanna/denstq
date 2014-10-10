/**********************************************************************
 * Program to parse Gaussian-type 
 * cube files and project spatial densities 
 * onto atoms
 *
 * CTC 7-16-14
 *
 * Input file has the following format:
 *
 * type:   transden/coupling
 * input:  cube/dens
 * donorfile:  file.cube/.dat
 * acceptorfile: file.cube/.dat
 * DA distance:  x,dx,y,dy,z,dz (in Angstroms)
 *
 * note:  currently only cube files supported
 *        DA distance ranges from dx to x in increments of dx
 *
 * usage: dens.out inputfile
 *
 * for more info on the .cube format see
 * http://www.gaussian.com/g_tech/g_ur/u_cubegen.htm
 * *******************************************************************/
#include "main.h"

int main(int argc, char *argv[]) {
  //open input file
  comfile.open(argv[1]);

  //parse input file
  char tempc[1000];
  string temps;
  comfile.getline(tempc,1000);
  temps = strtok(tempc, ":");
  //type of calculation
  calctype = strtok(NULL, ": ");
  comfile.getline(tempc,1000);
  temps = strtok(tempc, ":");
  //input file type
  inputtype = strtok(NULL, ":");

  if (calctype == "transden") {
    cout<<calctype<<endl;
    comfile.getline(tempc,1000);
    temps = strtok(tempc, ":");
    temps = strtok(NULL, " :");
    ifstream infile;
    infile.open(temps.c_str());

    //declare atom objects and read in densities
    Atom *atoms = collectDens(atoms,infile);

    //project the density onto the nearest atom
    //writes density to "density.dat"
    outfile.open("output.dat");
    projectdens(natoms,atoms,dens,posx,posy,posz,nx,ny,nz,densrad);

    //Calculate transition dipole moment
    calcdip(atoms);
  
  } else if (calctype == "coupling") {
    cout<<calctype<<endl;

    //declare donor and acceptor objects
    Molecule donor,acceptor;

    //get donor file name
    comfile.getline(tempc,1000);
    temps = strtok(tempc, ":");
    temps = strtok(NULL, " :");
    ifstream infile;
    infile.open(temps.c_str());
    if (!infile.good()) {
      cout<<"cannot open donor file!"<<endl;
      return -1;
    }

    //declare donor atoms and read in density
    donor.atoms = collectDens(&donor,donor.atoms,infile);
    
    //project the donor density onto the nearest atom
    outfile.open("donor_dens.dat");
    projectdens(donor.natoms,donor.atoms,donor.dens,
                donor.posx,donor.posy,donor.posz,
                donor.nx,donor.ny,donor.nz,donor.densrad);
    
    //get acceptor file name
    comfile.getline(tempc,1000);
    temps = strtok(tempc, ":");
    temps = strtok(NULL, " :");
    infile.close();
    infile.open(temps.c_str());
    if(!infile.good()) {
      cout<<"cannot open acceptor file!"<<endl;
      return -1;
    }

    //declare donor atoms and read in density
    acceptor.atoms = collectDens(&acceptor,acceptor.atoms,infile);
    
    //project the acceptor density onto the nearest atom
    outfile.close();
    outfile.open("acceptor_dens.dat");
    projectdens(acceptor.natoms,acceptor.atoms,acceptor.dens,
                acceptor.posx,acceptor.posy,acceptor.posz,
                acceptor.nx,acceptor.ny,acceptor.nz,acceptor.densrad);
    outfile.close();

    //read the DA distance and increments
    comfile.getline(tempc,1000);
    temps = strtok(tempc,":");
    temps = strtok(NULL, " ,");
    double xrange = atof(temps.c_str());
    temps = strtok(NULL, " ,");
    double dx = atof(temps.c_str());
    temps = strtok(NULL, " ,");
    double yrange = atof(temps.c_str());
    temps = strtok(NULL, " ,");
    double dy = atof(temps.c_str());
    temps = strtok(NULL, " ,");
    double zrange = atof(temps.c_str());
    temps = strtok(NULL, " ,");
    double dz = atof(temps.c_str());
    int xsteps = (int)(xrange/dx);

    if (xsteps == 0) {
      cout<<"xinc=0, setting xinc=1"<<endl;
      xsteps = 1;
    }


    //return 0;

    outfile.open("coupling.dat");
    for (int j=0; j<xsteps; j++) {
      double xdisp = 0 + j*dx;
      for (int i=0; i<donor.natoms; i++) {
        donor.atoms[i].x += xdisp;
      }
      //compute coupling from two cube files
      double coupling = computeCoupling(&donor,&acceptor);
      outfile<<xdisp<<" "<<coupling<<endl;
    }
    infile.close();
 }
    
  //tidy up
  comfile.close();
  outfile.close();

  return 0;
}

/*******************************************************
 ******************************************************/

/*****************************************
 * Project the density onto nearest atom
 * ***************************************/

void projectdens(const int natoms, Atom *atoms, double ***dens, double posx[], 
  double posy[], double posz[], int nx, int ny, int nz,double ***densrad) {
  
  double vol = dx1*dy2*dz3;

  for (int i=0; i<nx; i++) {
    for (int j=0; j<ny; j++) {
      for (int k=0; k<nz; k++) {
        int closest=-1;
        double closestdist = 100.;
        for (int l=0; l<natoms; l++) {
          double dist2,closestdis;
          dist2 = ((atoms[l].x - posx[i] )*(atoms[l].x - posx[i])
                  +(atoms[l].y - posy[j] )*(atoms[l].y - posy[j])
                  +(atoms[l].z - posz[k] )*(atoms[l].z - posz[k]));
          if (dist2 < closestdist) {
            closest = l;
            closestdist = dist2;
          }
        } //end atoms
        //add density to closest atom

        atoms[closest].dens += dens[i][j][k];
      } //end z
    } //end y
  } //end x
  
  for (int i=0; i<natoms; i++)
    atoms[i].dens *= sqrt(2);
  //print atomic positions and projected densities
  for (int i=0; i<natoms; i++) {
    outfile<<atoms[i].num+1<<" "<<atoms[i].type<<" "<<atoms[i].x<<" "<<atoms[i].y<<" "<<atoms[i].z<<" "<<atoms[i].dens*vol<<endl;
  }
}

Atom *collectDens(Atom *atoms, ifstream &infile) {
//determine number of lines in file
  int nlines = -1;
  char temp[10000];
  while(infile.good()) {
    infile.getline(temp,1000);
    nlines++;
  }
  infile.clear();
  infile.seekg(0);
  infile.getline(temp,1000);
  infile.getline(temp,1000);
  infile.getline(temp,1000);

  char *tempc;
  //read first two line of data to get natoms, and origin
  natoms = atof(strtok(temp," "));
  ox = atof(strtok(NULL," "));// * au2ang;
  oy = atof(strtok(NULL," "));// * au2ang;
  oz = atof(strtok(NULL," "));// * au2ang;

  //get nx,ny,nz,dx,dy,dz
  infile.getline(temp, 1000);
  nx = atoi(strtok(temp, " "));
  dx1 = atof(strtok(NULL, " "));// * au2ang;
  dy1 = atof(strtok(NULL, " "));// / ang2au;
  dz1 = atof(strtok(NULL, " "));// / ang2au;

  infile.getline(temp, 1000);
  ny = atoi(strtok(temp, " "));
  dx2 = atof(strtok(NULL, " "));// / ang2au;
  dy2 = atof(strtok(NULL, " "));/// ang2au;
  dz2 = atof(strtok(NULL, " "));/// ang2au;
  
  infile.getline(temp, 1000);
  nz = atoi(strtok(temp, " "));
  dx3 = atof(strtok(NULL, " "));/// ang2au;
  dy3 = atof(strtok(NULL, " "));/// ang2au;
  dz3 = atof(strtok(NULL, " "));/// ang2au;
  
  atoms = new Atom[natoms];
  cout<<natoms<<" Atoms"<<endl;

//get atomic numbers, charges, and coordinates
  for (int i=0; i<natoms; i++) {
    infile.getline(temp, 1000);
    tempc = strtok(temp, " ");
    atoms[i].atomicnum = atoi(tempc);
    atoms[i].type = str_Atom[atoms[i].atomicnum];
    atoms[i].charge = atof(strtok(NULL, " "));
    atoms[i].x = atof(strtok(NULL, " "));// / ang2au;
    atoms[i].y = atof(strtok(NULL, " "));// / ang2au;
    atoms[i].z = atof(strtok(NULL, " "));// / ang2au;
    atoms[i].r = sqrt(atoms[i].x*atoms[i].x +
                  atoms[i].y*atoms[i].y + atoms[i].z*atoms[i].z);
    atoms[i].num = i;
  }

  //declare density array
  dens = new double**[nx];
  densrad = new double **[nx];
  for (int i=0; i<nx; i++) {
    dens[i] = new double*[ny];
    densrad[i] = new double*[ny];
    for (int j=0; j<ny; j++) {
      dens[i][j] = new double[nz];
      densrad[i][j] = new double[nz];
    }
  }
  posx = new double[nx];
  posy = new double[ny];
  posz = new double[nz];
  double vol = dx1*dy2*dz3;

  //import density into array and define positions
  double sum=0.;
  double dx,dy,dz;
  dx=0.; dy=0.; dz=0.;
  for (int i=0; i<nx; i++) {
    posx[i] = ox + i*dx1;
    for (int j=0; j<ny; j++) {
      posy[j] = oy + j*dy2;
      for (int k=0; k<nz; k++) {
        posz[k] = oz + k*dz3;
        infile>>dens[i][j][k];
        /* Threshold */
        //dens[i][j][k] *= vol;
        if (fabs(dens[i][j][k]) <= thresh) dens[i][j][k] = 0.;
        //dens[i][j][k] *= dx1*dy2*dz3;
        sum += dens[i][j][k];
        /*densrad[i][j][k] = sqrt(posx[i+j*nx+k*nx*ny]*posx[i+j*nx+k*nx*ny]
                               +posy[i+j*nx+k*nx*ny]*posy[i+j*nx+k*nx*ny]
                               +posz[i+j*nx+k*nx*ny]*posz[i+j*nx+k*nx*ny]);
        */
        
        dx+=dens[i][j][k]*posx[i];//*ang2au;
        dy+=dens[i][j][k]*posy[j];//*ang2au;
        dz+=dens[i][j][k]*posz[k];//*ang2au;
      }
    }
  }
  
  dx *= vol;
  dy *= vol;
  dz *= vol;

  cout<<"density dipole = "<<sqrt(dx*dx+dy*dy+dz*dz)/0.393456<<endl;
  //normalize density
  cout<<"density sum = "<<sum<<" "<<nx<<" "<<ny<<" "<<nz<<endl;
  

  return atoms;
}

Atom *collectDens(Molecule *mol, Atom *atoms, ifstream &infile) {
//determine number of lines in file
  int nlines = -1;
  char temp[10000];
  while(infile.good()) {
    infile.getline(temp,1000);
    nlines++;
  }
  infile.clear();
  infile.seekg(0);
  infile.getline(temp,1000);
  infile.getline(temp,1000);
  infile.getline(temp,1000);
  char *tempc;
  //read first two line of data to get natoms, and origin
  mol->natoms = atof(strtok(temp," "));
  mol->ox = atof(strtok(NULL," "));
  mol->oy = atof(strtok(NULL," "));
  mol->oz = atof(strtok(NULL," "));

  //get nx,ny,nz,dx,dy,dz
  infile.getline(temp, 1000);
  mol->nx = atoi(strtok(temp, " "));
  mol->dx1 = atof(strtok(NULL, " "));
  mol->dx2 = atof(strtok(NULL, " "));
  mol->dx3 = atof(strtok(NULL, " "));

  infile.getline(temp, 1000);
  mol->ny = atoi(strtok(temp, " "));
  mol->dy1 = atof(strtok(NULL, " "));
  mol->dy2 = atof(strtok(NULL, " "));
  mol->dy3 = atof(strtok(NULL, " "));
  
  infile.getline(temp, 1000);
  mol->nz = atoi(strtok(temp, " "));
  mol->dz1 = atof(strtok(NULL, " "));
  mol->dz2 = atof(strtok(NULL, " "));
  mol->dz3 = atof(strtok(NULL, " "));
  
  atoms = new Atom[mol->natoms];

//get atomic numbers, charges, and coordinates
  for (int i=0; i<mol->natoms; i++) {
    infile.getline(temp, 1000);
    tempc = strtok(temp, " ");
    atoms[i].atomicnum = atoi(tempc);
    atoms[i].type = str_Atom[atoms[i].atomicnum];
    atoms[i].charge = atof(strtok(NULL, " "));
    atoms[i].x = atof(strtok(NULL, " "));// / ang2au;
    atoms[i].y = atof(strtok(NULL, " "));// / ang2au;
    atoms[i].z = atof(strtok(NULL, " "));// / ang2au;
    atoms[i].r = sqrt(atoms[i].x*atoms[i].x +
    atoms[i].y*atoms[i].y + atoms[i].z*atoms[i].z);
    atoms[i].num = i;
  }

  //declare density array
  mol->dens = new double**[mol->nx];
  mol->densrad = new double **[mol->nx];
  for (int i=0; i<mol->nx; i++) {
    mol->dens[i] = new double*[mol->ny];
    mol->densrad[i] = new double*[mol->ny];
    for (int j=0; j<mol->ny; j++) {
      mol->dens[i][j] = new double[mol->nz];
      mol->densrad[i][j] = new double[mol->nz];
    }
  }
  mol->posx = new double[mol->nx];
  mol->posy = new double[mol->ny];
  mol->posz = new double[mol->nz];

  //import density into array and define positions
  for (int i=0; i<mol->nx; i++) {
    for (int j=0; j<mol->ny; j++) {
      for (int k=0; k<mol->nz; k++) {
        mol->posx[i] = mol->ox + i*mol->dx1 + j*mol->dx2 + k*mol->dx3;
        mol->posy[j] = mol->oy + i*mol->dy1 + j*mol->dy2 + k*mol->dy3;
        mol->posz[k] = mol->oz + i*mol->dz1 + j*mol->dz2 + k*mol->dz3;
        infile>>mol->dens[i][j][k];
        mol->densrad[i][j][k] = sqrt(mol->posx[i]*mol->posx[i]+mol->posy[j]*mol->posy[j]+mol->posz[k]*mol->posz[k]);
      }
    }
  }

  return atoms;
}

/***********************************
 * Calculate the coupling constant
 * given two molecules and their 
 * transition charges
 * *********************************/
double computeCoupling(Molecule *mold, Molecule *mola) {
  double res = 0.;
  for (int i=0; i<mold->natoms; i++) {
    for (int j=0; j<mola->natoms; j++) {
      double rda2 = (mold->atoms[i].x - mola->atoms[j].x)*(mold->atoms[i].x - mola->atoms[j].x)
                   + (mold->atoms[i].y - mola->atoms[j].y)*(mold->atoms[i].y - mola->atoms[j].y)
                   + (mold->atoms[i].z - mola->atoms[j].z)*(mold->atoms[i].z - mola->atoms[j].z);
      double rda = sqrt(rda2);
      double interaction = mold->atoms[i].dens * mola->atoms[j].dens / rda;
      res += interaction;
    }
  }
  return res;
}

/**********************************
 * Calculate transition dipole
 * moment
 * *******************************/
void calcdip(Atom *atoms) {
  double dx,dy,dz,sum,sump,summ;
  dx=0.; dy=0.; dz=0.; sum=0.; sump=0.;summ=0.;
  double vol = dx1*dy2*dz3;

  for (int i=0; i<natoms; i++) {
    dx += atoms[i].x * atoms[i].dens;
    dy += atoms[i].y * atoms[i].dens;
    dz += atoms[i].z * atoms[i].dens;
    sum += atoms[i].dens;
    if (atoms[i].dens > 0)
      sump += atoms[i].dens;
    else if (atoms[i].dens < 0)
      summ += atoms[i].dens;
  
    cout<<"atom "<<i<<", "<<atoms[i].type<<", "<<atoms[i].dens*vol<<endl;
  }
  
  dx *= vol;
  dy *= vol;
  dz *= vol;

  cout<<"net transition charge = "<<sum<<endl;
  cout<<"negative density = "<<summ<<endl;
  cout<<"positive density = "<<sump<<endl;
  cout<<"Transition dipole moment = "<<dx<<" x, "<<dy<<" y, "<<dz<<" z"<<endl;
  cout<<"Magnitude t = "<<sqrt(dx*dx+dy*dy+dz*dz)/0.393456<<" D"<<endl;
}
