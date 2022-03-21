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
double ke, ke1, ke2, se1, se2, gpe, gpe1, gpe2, rho1, rho3, Rhor, Ohd, mu1, mu2, mu3, Ohf, eps, eps1, eps2, Bo, vcm, zcm, wt;

char filename[80], nameEnergy[80];
// boundary conditions
u.t[left] = dirichlet(0.0);
f1[left] = dirichlet(0.0);
f2[left] = dirichlet(1.0);

int main(int a, char const *arguments[]) {
  sprintf (filename, "%s", arguments[1]);
  sprintf(nameEnergy, "%s", arguments[2]);
  Rhor = atof(arguments[3]);
  Ohd = atof(arguments[4]);
  Ohf = atof(arguments[5]);
  Bo = atof(arguments[6]);

  FILE *fp;
  fp = fopen (nameEnergy, "a");
  restore (file = filename);

  rho1 = 1.0; mu1 = Ohd; mu2 = Ohf;
  rho3 = 1e-3; mu3 = 1e-5;

  f1.prolongation = fraction_refine;
  f2.prolongation = fraction_refine;

  boundary((scalar *){f1, f2, u.x, u.y});
  scalar sf1[], sf2[], *smearInterfaces = {sf1, sf2};

  foreach(){
    sf1[] = (4.*f1[] + 
    2.*(f1[0,1] + f1[0,-1] + f1[1,0] + f1[-1,0]) + 
    f1[-1,-1] + f1[1,-1] + f1[1,1] + f1[-1,1])/16.;

    sf2[] = (4.*f2[] + 
    2.*(f2[0,1] + f2[0,-1] + f2[1,0] + f2[-1,0]) + 
    f2[-1,-1] + f2[1,-1] + f2[1,1] + f2[-1,1])/16.;
  }
  for (scalar sf in smearInterfaces){
    sf.prolongation = refine_bilinear;
    boundary ({sf});
  }

  /*
  Do calculations start
  */
  ke = 0., ke1 = 0., ke2 = 0., gpe = 0., gpe1 = 0., gpe2 = 0., se1 = 0., se2 = 0., eps = 0., eps1 = 0., eps2 = 0., vcm = 0., zcm = 0., wt = 0.;

  foreach (){
    double rho = clamp(sf1[]+sf2[], 0., 1.)*(rho1 - rho3) + rho3;
    ke += 2*pi*y*(0.5*rho*(sq(u.x[]) + sq(u.y[])))*sq(Delta); // 2*pi*\int_A(0.5*rho*v^2)*r*dr*dz

    ke1 += 2*pi*y*(0.5*clamp(sf1[], 0., 1.)*rho1*(sq(u.x[]) + sq(u.y[])))*sq(Delta); // 2*pi*\int_A(0.5*rho*v^2)*r*dr*dz
    ke2 += 2*pi*y*(0.5*clamp(sf2[], 0., 1.)*rho1*(sq(u.x[]) + sq(u.y[])))*sq(Delta); // 2*pi*\int_A(0.5*rho*v^2)*r*dr*dz

    gpe += (2*pi*y)*(rho*Bo*x)*sq(Delta); // 2*pi*\int_A(rho*g*z)rdrdz
    gpe1 += (2*pi*y)*(clamp(sf1[], 0., 1.)*rho1*Bo*x)*sq(Delta); // 2*pi*\int_A(rho*g*z)rdrdz
    gpe2 += (2*pi*y)*(clamp(sf2[], 0., 1.)*rho1*Bo*x)*sq(Delta); // 2*pi*\int_A(rho*g*z)rdrdz

    zcm += (2*pi*y)*(rho1*clamp(sf1[], 0., 1.)*x)*sq(Delta);
    vcm += (2*pi*y)*(rho1*clamp(sf1[], 0., 1.)*u.x[])*sq(Delta);
    wt += (2*pi*y)*rho1*clamp(sf1[], 0., 1.)*sq(Delta);

    double D11 = (u.y[0,1] - u.y[0,-1])/(2*Delta);
    double D22 = (u.y[]/max(y,1e-20));
    double D33 = (u.x[1,0] - u.x[-1,0])/(2*Delta);
    double D13 = 0.5*( (u.y[1,0] - u.y[-1,0] + u.x[0,1] - u.x[0,-1])/(2*Delta) );
    double D2 = (sq(D11)+sq(D22)+sq(D33)+2.0*sq(D13));

    double mu = clamp(sf1[], 0., 1.)*mu1 + clamp(sf2[], 0., 1.)*mu2 + clamp(1.-sf1[]-sf2[], 0., 1.)*mu3;

    eps += (2*pi*y)*( 2*mu*D2 )*sq(Delta);
    eps1 += (2*pi*y)*( 2*clamp(sf1[], 0., 1.)*mu1*D2 )*sq(Delta);
    eps2 += (2*pi*y)*( 2*clamp(sf2[], 0., 1.)*mu2*D2 )*sq(Delta);
  }
  zcm /= wt; vcm /= wt;
  se1 = interface_energy (f1);
  se2 = interface_energy (f2);

  boundary((scalar *){f1, f2, sf1, sf2, u.x, u.y});

  double Zmin1 = 0., Zmin = 0., temp = 0.;
  foreach_boundary(bottom){
    if (f2[] > 1e-6 && f2[] < 1. - 1e-6) {
      Zmin1 = x;
      // fprintf(ferr, "%f\n", x);
    }
  }
  temp = HUGE;
  foreach_boundary(bottom){
    if (f1[] > 1e-6 && f1[] < 1. - 1e-6) {
      // fprintf(ferr, "%f\n", x);
      if (fabs(x-Zmin1) < temp){
        temp = fabs(x-Zmin1);
      }
      if (temp < 2.5*Delta){
        temp = 0.;
      }
    }
  }
  Zmin = temp;
  /*
  Do calculations end
  */

 if (t == 0){
    fprintf(ferr, "Rhor %g, Ohd %3.2e, Ohf %3.2e, Bo %g\n", Rhor, Ohd, Ohf, Bo);
    fprintf(ferr, "t ke ke1 ke2 gpe gpe1 gpe2 se1 se2 eps eps1 eps2 vcm zcm Zmin\n");
    fprintf(fp, "t ke ke1 ke2 gpe gpe1 gpe2 se1 se2 eps eps1 eps2 vcm zcm Zmin\n");
  }

  fprintf(ferr, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", t, ke, ke1, ke2, gpe, gpe1, gpe2, se1, se2, eps, eps1, eps2, vcm, zcm, Zmin);
  fprintf(fp, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", t, ke, ke1, ke2, gpe, gpe1, gpe2, se1, se2, eps, eps1, eps2, vcm, zcm, Zmin);
}
