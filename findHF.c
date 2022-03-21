/* Title: Getting Energy
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
*/
#include "axi.h"
#include "navier-stokes/centered.h"
#include "fractions.h"

trace
double interface_energy (scalar c){
  double se = 0.;
  foreach (reduction(+:se)){
    if (c[] > 1e-6 && c[] < 1. - 1e-6) {
      coord p, n = interface_normal (point, c);
      double alpha = plane_alpha (c[], n);
      double len = line_length_center(n, alpha, &p);
      se += 2.*pi*( y + p.y*Delta )*(len*Delta); // 2*pi*\int_l (r_c)dl
    }
  }
  return se;
}

scalar f1[], f2[], *interfaces = {f1, f2};
double se1, se2;

char filename[80], nameEnergy[80];
// boundary conditions
u.t[left] = dirichlet(0.0);
f1[left] = dirichlet(0.0);
f2[left] = dirichlet(1.0);

int main(int a, char const *arguments[]) {
  sprintf (filename, "%s", arguments[1]);
  sprintf(nameEnergy, "%s", arguments[2]);

  FILE *fp;
  fp = fopen (nameEnergy, "a");
  restore (file = filename);

  f1.prolongation = fraction_refine;
  f2.prolongation = fraction_refine;

  boundary((scalar *){f1, f2, u.x, u.y});

  /*
  Do calculations start
  */
  se1 = 0., se2 = 0.;

  se1 = interface_energy (f1);
  se2 = interface_energy (f2);

  boundary((scalar *){f1, f2, u.x, u.y});

  double Zmin1 = 0.;
  foreach_boundary(bottom){
    if (f2[] > 1e-6 && f2[] < 1. - 1e-6) {
      Zmin1 = x;
      // fprintf(ferr, "%f\n", x);
    }
  }
  /*
  Do calculations end
  */

 if (t == 0){
    fprintf(ferr, "t se1 se2 Zmin\n");
    fprintf(fp, "t se1 se2 Zmin\n");
  }

  fprintf(ferr, "%f %f %f %f\n", t, se1, se2, Zmin1);
  fprintf(fp, "%f %f %f %f\n", t, se1, se2, Zmin1);
  fclose(fp);
}
