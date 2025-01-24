/* Title: Saving images with bview
# Authors: Vatsal & Youssef
# vatsalsanjay@gmail.com
# Physics of Fluids
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "view.h"

scalar f1[], f2[], omega[], vel[];
char filename[80], Imagename[80];

int main(int a, char const *arguments[]) {
  sprintf (filename, "%s", arguments[1]);
  sprintf (Imagename, "%s",arguments[2]);
  restore (file = filename);
  vorticity(u, omega);
  foreach(){
    vel[] = sqrt(sq(u.x[]) + sq(u.y[]));
  }

  // boundary conditions
  u.t[left] = dirichlet(0.0);
  f1[left] = dirichlet(0.0);
  f2[left] = dirichlet(1.0);
  f1.prolongation = fraction_refine;
  f2.prolongation = fraction_refine;
  boundary((scalar *){f1, f2, u.x, u.y, omega, vel});

  view (fov = 15.0, quat = {0,0,-0.707107,0.707107}, tx = 0.0, ty = -0.30, bg = {1,1,1}, width = 1920, height = 1016, samples = 1);
  draw_vof("f1", lw=4);
  draw_vof("f2", lw=4);
  squares ("omega", spread=9, linear=true, map = cool_warm);
  box(notics=true);
  cells(n = {0, 0, 1}, lc = {0.5,0.5,0.5}, lw=1);

  mirror ({0,1}) {
    draw_vof("f1", lw=4);
    draw_vof("f2", lw=4);
    squares ("vel", linear=true);
  }
  save (Imagename);
}
