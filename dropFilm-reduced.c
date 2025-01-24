/* Title: Bouncing Droplet on thin/deep films!
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids

# This version is essentially same as dropFilm.c, but with reduced gravity implmentation. For the reduced gravity implementaton, please see: S. Popinet, Annu. Rev. Fluid Mech., 50:1, 49-75 (2018) (https://doi.org/10.1146/annurev-fluid-122316-045034). 

version 2.0
# Changelog
- v2.0.1 (2025-01-24)
  - Added reduced gravity implementation. To use reduced gravity, please use dropFilm-reduced.c. 
- v2.0 (2025-01-24)
  - Now, user can choose the density of the film independently of the drop.
  - Now, user can set the surface tension of the film-air interface independently. 
  - Repeating variables are: 1. Density of the drop, 2. Radius of the drop, 3. Surface tension of the drop-air interface. 
  - updated how to add default values. 
  - removed omega adaptation because this is somewhat broken in the newest basilisk version for axi cases.
  - removed adapt_wavelet_limited as it is incompatible with newest Basilisk. If you want to implement adapt_wavelet_limited, please refer to: https://github.com/comphy-lab/adapt-wavelet-limited 
  - introduced dirichlet(...) for setting boundary condition of f1 and f2. At high enough resolution, both f...[left] = ...; will work identically to f...[left] = dirichlet(...); However, when you have limited resolution, dirichlet(...) is better.
*/

// 1 is drop

#include "axi.h"
#include "navier-stokes/centered.h"
#define FILTERED
#include "two-phaseTF.h"
#include "tension.h"
#include "reduced-gravity.h"

// Error tolerances
#define fErr (1e-3)                                 // error tolerance in VOF
#define KErr (1e-4)                                 // error tolerance in KAPPA
#define VelErr (1e-2)                            // error tolerances in velocity
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
f1[left] = dirichlet(0.0);
f2[left] = dirichlet(1.0);

/*
Older version
- At high enough resolution, both f...[left] = ...; will work identically to f...[left] = dirichlet(...); However, when you have limited resolution, dirichlet(...) is better.
*/
// f1[left] = 0.0;
// f2[left] = 1.0;

int MAXlevel;
double tmax, We, Ohd, Bo, Ohf, hf, rhof, sigma23;
#define MINlevel 4                                            // maximum level
#define tsnap (0.01)

int main(int argc, char const *argv[]) {
  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);

  // Default values for all parameters
  MAXlevel = 8;  // Default maximum refinement level
  tmax = 1.0;  // Default simulation end time
  We = 4.0;   // Default Weber number
  Ohd = 1e-2; // Default Ohnesorge number (drop)
  Bo = 0.25;   // Default Bond number
  Ohf = 1e-2; // Default Ohnesorge number (film)
  hf = 0.25;    // Default film thickness
  rhof = 1.0;  // Default film density
  sigma23 = 1.0; // Default surface tension ratio

  // Print usage if help is requested
  if (argc == 2 && (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0)) {
    fprintf(stderr, "Usage: %s [MAXlevel] [tmax] [We] [Ohd] [Bo] [Ohf] [hf] [rhof] [sigma23]\n", argv[0]);
    fprintf(stderr, "Default values:\n");
    fprintf(stderr, "  MAXlevel = %d\n  tmax = %g\n  We = %g\n  Ohd = %g\n  Bo = %g\n",
            MAXlevel, tmax, We, Ohd, Bo);
    fprintf(stderr, "  Ohf = %g\n  hf = %g\n  rhof = %g\n  sigma23 = %g\n",
            Ohf, hf, rhof, sigma23);
    return 0;
  }

  // Override default values with command line arguments if provided
  if (argc > 1) MAXlevel = atoi(argv[1]);
  if (argc > 2) tmax = atof(argv[2]);
  if (argc > 3) We = atof(argv[3]);
  if (argc > 4) Ohd = atof(argv[4]);
  if (argc > 5) Bo = atof(argv[5]);
  if (argc > 6) Ohf = atof(argv[6]);
  if (argc > 7) hf = atof(argv[7]);
  if (argc > 8) rhof = atof(argv[8]);
  if (argc > 9) sigma23 = atof(argv[9]);

  // Validate critical parameters
  if (hf <= 0) {
    fprintf(stderr, "Error: Film height (hf) must be positive. Current value: %g\n", hf);
    return 1;
  }

  fprintf(stderr, "Running simulation with:\n");
  fprintf(stderr, "Level %d tmax %g\nWe %g, Ohd %g, Bo %g\nOhf %3.2e, hf %3.2f, rhof %g, sigma23 %g\n", 
          MAXlevel, tmax, We, Ohd, Bo, Ohf, hf, rhof, sigma23);

  L0=Ldomain;
  X0=-hf; Y0=0.;
  init_grid (1 << (6));

  // drop properties.
  rho1 = 1.0; mu1 = Ohd; 

  // film properties.
  rho2 = rhof; mu2 = Ohf;
  
  // air properties. Note that mu3 = 1e-5 is simply to keep the dimensionless viscosity of air negligible. You can change this value to exactly match: \eta_a/\sqrt{\rho_d\gamma R_0}. For 1 mm water drop, this will be 1.35e-4 (using silicone oil density and surface tenson in drop of 1 mm). 
  rho3 = Rho21; mu3 = 1e-5;

  // reduced gravity --> #fixME, in principle, these two will never be different. :p
  Bf1.x = -Bo;
  Bf2.x = -Bo;

  f1.sigma = 1.0; // this is always 1. 
  f2.sigma = sigma23; // ratio of surface tension coefficient of film-air to that of drop-air interface.

  /*
  **Note:** The above framework ensures non coalescence between drop and the film. As a result, the surface tesnion coefficient of film-drop interface is naturally f1.sigma+f2.sigma (so 1+sigma23). This is by design and please only change this if you know what you are doing.
  */

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
  }
}


event adapt(i++){
  scalar KAPPA1[], KAPPA2[];
  curvature(f1, KAPPA1);
  curvature(f2, KAPPA2);
  scalar D2c[];
  foreach(){
    double D11 = (u.y[0,1] - u.y[0,-1])/(2*Delta);
    double D22 = (u.y[]/max(y,1e-20));
    double D33 = (u.x[1,0] - u.x[-1,0])/(2*Delta);
    double D13 = 0.5*( (u.y[1,0] - u.y[-1,0] + u.x[0,1] - u.x[0,-1])/(2*Delta) );
    double D2 = (sq(D11)+sq(D22)+sq(D33)+2.0*sq(D13));
    D2c[] = (f1[]*Ohd+f2[]*Ohf)*D2;
  }
  adapt_wavelet ((scalar *){f1, f2, KAPPA1, KAPPA2, u.x, u.y, D2c},
     (double[]){fErr, fErr, KErr, KErr, VelErr, VelErr, DissErr},
      MAXlevel, MINlevel);
}
// Outputs
// static
event writingFiles (t = 0, t += tsnap; t <= tmax+tsnap) {
  dump (file = "dump");
  char nameOut[80];
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump (file = nameOut);
}

event logWriting (i+=100) {
  double ke = 0., vcm = 0., wt = 0.;
  foreach (reduction(+:ke), reduction(+:vcm), reduction(+:wt)){
    ke += 2*pi*y*(0.5*rho(f1[],f2[])*(sq(u.x[]) + sq(u.y[])))*sq(Delta);
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
