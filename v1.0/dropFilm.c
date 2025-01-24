/* Title: Bouncing Droplet on thin/deep films!
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
*/

// 1 is drop

#include "axi.h"
#include "navier-stokes/centered.h"
#define FILTERED
#include "two-phaseTF.h"
#include "tension.h"
#include "adapt_wavelet_limited.h"

// Error tolerances
#define fErr (1e-3)                                 // error tolerance in VOF
#define KErr (1e-4)                                 // error tolerance in KAPPA
#define VelErr (1e-2)                            // error tolerances in velocity
#define OmegaErr (1e-2)                            // error tolerances in vorticity
#define DissErr (1e-3)

// air-water
#define Rho21 (1e-3)

// Calculations!
#define Xdist (1.040)
#define R2Drop(x,y,z) (sq(x - Xdist) + sq(y))
// domain
#define Ldomain 4                                // Dimension of the domain

// boundary conditions
u.t[left] = dirichlet(0.0);
f1[left] = 0.0;
f2[left] = 1.0;

int MAXlevel;
double tmax, We, Ohd, Bo, Ohf, hf;
#define MINlevel 3                                            // maximum level
#define tsnap (0.01)

int main(int argc, char const *argv[]) {
  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);

  if (argc < 8){
    fprintf(ferr, "Lack of command line arguments. Check! Need %d more arguments\n",8-argc);
    return 1;
  }

  MAXlevel = atoi(argv[1]);
  tmax = atof(argv[2]);
  We = atof(argv[3]); // We is 1 for 0.167 m/s <816*0.167^2*0.00075/0.017>
  Ohd = atof(argv[4]); // <0.000816/sqrt(816*0.017*0.00075) = 0.008>
  Bo = atof(argv[5]); // <816*10*0.00075^2/0.017 = 0.27>
  Ohf = atof(argv[6]);
  hf = atof(argv[7]);

  if (hf == 0){
    fprintf(ferr, "We have a problem. Wrong code. Change code or film height.\n");
    return 1;
  }
  fprintf(ferr, "Level %d tmax %g. We %g, Ohd %g, Bo %g, Ohf %3.2e, hf %3.2f\n", MAXlevel, tmax, We, Ohd, Bo, Ohf, hf);

  L0=Ldomain;
  X0=-hf; Y0=0.;
  init_grid (1 << (MINlevel));

  rho1 = 1.0; mu1 = Ohd; mu2 = Ohf;
  rho3 = Rho21; mu3 = 1e-5;

  f1.sigma = 1.0; f2.sigma = 1.0;

  run();

}

event init(t = 0){
  if(!restore (file = "dump")){
    refine((R2Drop(x,y,z) < 1.44) && (level < MAXlevel));
    fraction (f1, 1. - R2Drop(x,y,z));
    fraction (f2, -x);
    foreach () {
      u.x[] = -sqrt(We)*f1[];
      u.y[] = 0.0;
    }
    boundary((scalar *){f1, f2, u.x, u.y});
  }
}

// Gravity
event acceleration(i++) {
  face vector av = a;
  foreach_face(x){
    av.x[] -= Bo;
  }
}

int refRegion(double x, double y, double z){
  return (x > -hf && x < 0.1 && y < 2.0 ? MAXlevel+1 : MAXlevel);
}

event adapt(i++){
  scalar KAPPA1[], KAPPA2[], omega[];
  vorticity (u, omega);
  curvature(f1, KAPPA1);
  curvature(f2, KAPPA2);
  scalar D2c[];
  foreach(){
    omega[] *= (f1[]+f2[]);
    double D11 = (u.y[0,1] - u.y[0,-1])/(2*Delta);
    double D22 = (u.y[]/max(y,1e-20));
    double D33 = (u.x[1,0] - u.x[-1,0])/(2*Delta);
    double D13 = 0.5*( (u.y[1,0] - u.y[-1,0] + u.x[0,1] - u.x[0,-1])/(2*Delta) );
    double D2 = (sq(D11)+sq(D22)+sq(D33)+2.0*sq(D13));
    D2c[] = (f1[]*Ohd+f2[]*Ohf)*D2;
  }
  boundary((scalar *){D2c, omega, KAPPA1, KAPPA2});
  adapt_wavelet_limited ((scalar *){f1, f2, KAPPA1, KAPPA2, u.x, u.y, omega, D2c},
     (double[]){fErr, fErr, KErr, KErr, VelErr, VelErr, OmegaErr, DissErr},
      refRegion, MINlevel);
}
// Outputs
// static
event writingFiles (i = 0, t += tsnap; t <= tmax) {
  dump (file = "dump");
  char nameOut[80];
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump (file = nameOut);
}

event logWriting (i+=100) {
  double ke = 0., vcm = 0., wt = 0.;
  foreach (reduction(+:ke), reduction(+:vcm), reduction(+:wt)){
    ke += 2*pi*y*(0.5*rho(f1[]+f2[])*(sq(u.x[]) + sq(u.y[])))*sq(Delta);
    vcm += 2*pi*y*(f1[]*u.x[])*sq(Delta);
    wt += 2*pi*y*f1[]*sq(Delta);
  }
  vcm /= wt;
  static FILE * fp;
  if (i == 0) {
    fprintf (ferr, "i dt t ke\n");
    fp = fopen ("log", "w");
    fprintf (fp, "i dt t ke vcm\n");
    fprintf (fp, "%d %g %g %g %g\n", i, dt, t, ke, vcm);
    fclose(fp);
  } else {
    fp = fopen ("log", "a");
    fprintf (fp, "%d %g %g %g %g\n", i, dt, t, ke, vcm);
    fclose(fp);
  }
  fprintf (ferr, "%d %g %g %g %g\n", i, dt, t, ke, vcm);
}
