/* Title: Interpolating data from dump files: gfs2oogl style
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
*/
#include "axi.h"
#include "navier-stokes/centered.h"

char filename[80];
int nx, ny;
double xmin, ymin, xmax, ymax;
scalar Omega[], D2[], vel[];
scalar f1[], f2[];
scalar * list = NULL;
double Ohd, Ohf, Ohs, hf, mu1, mu2, mu3;

int main(int a, char const *arguments[]){

  Ohd = atof(arguments[7]);
  Ohf = atof(arguments[8]);
  hf = atof(arguments[9]);

  L0 = 8.0;
  X0=-hf; Y0=0.;

  mu1 = Ohd; mu2 = Ohf; mu3 = 1e-5;

  // boundary conditions
  u.t[left] = dirichlet(0.0);
  f1[left] = dirichlet(0.0);
  f2[left] = dirichlet(1.0);


  sprintf (filename, "%s", arguments[1]);
  xmin = atof(arguments[2]); ymin = atof(arguments[3]);
  xmax = atof(arguments[4]); ymax = atof(arguments[5]);
  nx = atoi(arguments[6]);

  // fprintf(ferr, "Mu21 %g, mu1 %g, mu2 %g, Oh %g\n", Mu21, mu1, mu2, Oh);

  list = list_add (list, f1);
  list = list_add (list, f2);
  list = list_add (list, D2);
  list = list_add (list, Omega);
  list = list_add (list, u.x);
  list = list_add (list, u.y);
  list = list_add (list, vel);

  restore (file = filename);
  boundary((scalar *){f1, f2, u.x, u.y});

  double Uc = 0., wt = 0., Deltax, Deltay;
  vorticity(u, Omega);

  foreach(){
    double D11 = (u.y[0,1] - u.y[0,-1])/(2*Delta);
    double D22 = (u.y[]/max(y,1e-20));
    double D33 = (u.x[1,0] - u.x[-1,0])/(2*Delta);
    double D13 = 0.5*( (u.y[1,0] - u.y[-1,0] + u.x[0,1] - u.x[0,-1])/(2*Delta) );
    D2[] = (sq(D11)+sq(D22)+sq(D33)+2.0*sq(D13));

    double mu = clamp(f1[],0.,1.)*mu1 + clamp(f2[],0.,1.)*mu2 + clamp(1.-f1[]-f2[],0.,1.)*mu3;
    
    D2[] *=  2*mu;

    if (D2[] > 0.){
      D2[] = log(D2[])/log(10);
    } else {
      D2[] = -10;
    }
    Uc += (2*pi*y)*u.x[]*clamp(f1[], 0.0, 1.0)*sq(Delta);
    wt += (2*pi*y)*clamp(f1[], 0.0, 1.0)*sq(Delta);

  }

  if (wt > 0){
    Uc /= wt;
  } else {
    Uc = 0.;
  }
  foreach(){
    vel[] = u.x[] - Uc;
  }

  boundary((scalar *){vel, D2, Omega});

  FILE * fp = ferr;
  Deltax = (double)(xmax-xmin)/(nx);
  // fprintf(ferr, "%g\n", Deltax);
  ny = (int)(ymax - ymin)/Deltax;
  // fprintf(ferr, "%d\n", ny);
  Deltay = (double)(ymax-ymin)/(ny);
  // fprintf(ferr, "%g\n", Deltay);
  int len = list_len(list);
  // fprintf(ferr, "%d\n", len);
  double ** field = (double **) matrix_new (nx, ny+1, len*sizeof(double));
  for (int i = 0; i < nx; i++) {
    double x = Deltax*(i+1./2) + xmin;
    for (int j = 0; j < ny; j++) {
      double y = Deltay*(j+1./2) + ymin;
      // fprintf(ferr, "%g %g\n", x, y);
      int k = 0;
      // fprintf(ferr, "x(i) = %g(%d), y(j) = %g(%d)\n", x, i, y, j);
      for (scalar s in list){
        field[i][len*j + k++] = interpolate (s, x, y);
      }
    }
  }

  for (int i = 0; i < nx; i++) {
    double x = Deltax*i + xmin;
    for (int j = 0; j < ny; j++) {
      double y = Deltay*j + ymin;
      fprintf (fp, "%g %g", x, y);
      int k = 0;
      for (scalar s in list){
        fprintf (fp, " %g", field[i][len*j + k++]);
      }
      fputc ('\n', fp);
    }
  }
  fflush (fp);
  // fclose (fp);
  matrix_free (field);

}
