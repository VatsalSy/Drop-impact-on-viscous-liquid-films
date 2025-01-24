/**
# Three-phase interfacial flows: two phases (f1, f2) are non-coalescing and the third (f1 = f2 = 0) is always air

# Version 2.0
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
# Last Updated: Jan 24, 2025

The interface between the fluids is tracked with a Volume-Of-Fluid
method. The volume fraction in drop is $f1=1$ and $f2=0$. In the thin film, it is $f2=1$ and $f1=0$. Air (fluid 3) is $f1 = f2 = 0$. The densities and dynamic viscosities for fluid 1 and 2 are *rho1*, *mu1*, *rho3*, *mu2*, respectively.
**Note:** The drop and the film are defined by different VoF fields, but have same properties (density and viscosity).
*/

#include "vof.h"
/**
Instead of one VoF tracer, we define two, f1 and f2.
*/
scalar f1[], f2[], *interfaces = {f1, f2};
double rho1 = 1., mu1 = 0., rho2 = 1., mu2 = 0., rho3 = 1., mu3 = 0.;
/**
Auxilliary fields are necessary to define the (variable) specific
volume $\alpha=1/\rho$ as well as the cell-centered density. */
face vector alphav[];
scalar rhov[];

event defaults (i = 0) {
  alpha = alphav;
  rho = rhov;
  /**
  If the viscosity is non-zero, we need to allocate the face-centered
  viscosity field. */
  mu = new face vector;
}

/**
The density and viscosity are defined using arithmetic averages by
default. The user can overload these definitions to use other types of
averages (i.e. harmonic). The difference comes in how we call these averages.
 * For any property A ∈ {μ,ρ}, the value is calculated as:
 * A = f₁A₁ + f₂A₂ + (1-f₁-f₂)A₃
 * where f₁, f₂ are volume fractions for phases 1 and 2,
 * A₁, A₂ are liquid properties (μ₁,μ₂ or ρ₁,ρ₂),
 * A₃ is gas property (μ₃ or ρ₃)
 */
#ifndef rho
#define rho(f1, f2)  (clamp(f1,0.,1.)*rho1 + clamp(f2,0.,1.)*rho2 + clamp(1.-f1-f2,0.,1.)*rho3)
#endif
#ifndef mu
#define mu(f1, f2)  (clamp(f1,0.,1.)*mu1 + clamp(f2,0.,1.)*mu2 + clamp(1.-f1-f2,0.,1.)*mu3)
#endif


/**
We have the option of using some "smearing" of the density/viscosity
jump. It is modified to take into account that there are two VoF tracers. */

#ifdef FILTERED
scalar sf1[], sf2[], *smearInterfaces = {sf1, sf2};
#else
#define sf1 f1
#define sf2 f2
scalar *smearInterfaces = {sf1, sf2};
#endif

event tracer_advection (i++) {

  /**
  When using smearing of the density jump, we initialise *sf* with the
  vertex-average of *f*. Introduce for loops to ensure that smearing is done properly. */
  #ifdef FILTERED
    int counter1 = 0;
    for (scalar sf in smearInterfaces){
      counter1++;
      int counter2 = 0;
      for (scalar f in interfaces){
        counter2++;
        if (counter1 == counter2){
          // fprintf(ferr, "%s %s\n", sf.name, f.name);
        #if dimension <= 2
            foreach(){
              sf[] = (4.*f[] +
          	    2.*(f[0,1] + f[0,-1] + f[1,0] + f[-1,0]) +
          	    f[-1,-1] + f[1,-1] + f[1,1] + f[-1,1])/16.;
            }
        #else // dimension == 3
            foreach(){
              sf[] = (8.*f[] +
          	    4.*(f[-1] + f[1] + f[0,1] + f[0,-1] + f[0,0,1] + f[0,0,-1]) +
          	    2.*(f[-1,1] + f[-1,0,1] + f[-1,0,-1] + f[-1,-1] +
          		f[0,1,1] + f[0,1,-1] + f[0,-1,1] + f[0,-1,-1] +
          		f[1,1] + f[1,0,1] + f[1,-1] + f[1,0,-1]) +
          	    f[1,-1,1] + f[-1,1,1] + f[-1,1,-1] + f[1,1,1] +
          	    f[1,1,-1] + f[-1,-1,-1] + f[1,-1,-1] + f[-1,-1,1])/64.;
            }
        #endif
        }
      }
    }
    #endif

  #if TREE
    for (scalar sf in smearInterfaces){
      sf.prolongation = refine_bilinear;
      sf.dirty = true; // boundary conditions need to be updated
    }
  #endif
}


event properties (i++) {

  foreach_face() {
    double ff1 = (sf1[] + sf1[-1])/2.;
    double ff2 = (sf2[] + sf2[-1])/2.;
    alphav.x[] = fm.x[]/rho(ff1, ff2);
    face vector muv = mu;
    muv.x[] = fm.x[]*mu(ff1, ff2);
  }

  foreach()
    rhov[] = cm[]*rho(sf1[], sf2[]);

#if TREE
  for (scalar sf in smearInterfaces){
    sf.prolongation = fraction_refine;
    sf.dirty = true; // boundary conditions need to be updated
  }
#endif
}
