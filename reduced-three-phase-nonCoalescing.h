/**
# Reduced gravity 

We re-express gravity in [two-phase flows](two-phase.h) as an
[interfacial force](iforce.h) as
$$
-\nabla p + \rho\mathbf{g} = 
-\nabla p' - [\rho]\mathbf{g}\cdot\mathbf{x}\mathbf{n}\delta_s
$$
with $p'= p - \rho\mathbf{g}\cdot\mathbf{x}$ the dynamic pressure and
$\rho\mathbf{g}\cdot\mathbf{x}$ the hydrostatic pressure. The corresponding 
potential is
$$
\phi = [\rho]\mathbf{G}\cdot(\mathbf{x} - \mathbf{Z})
$$
with $\mathbf{G}$ the gravity vector and $\mathbf{Z}$ an optional
reference level. */

coord Bf1 = {0.,0.,0.}, Z1 = {0.,0.,0.};
coord Bf2 = {0.,0.,0.}, Z2 = {0.,0.,0.};

/**
We need the interfacial force module as well as some
functions to compute the position of the interface. */

#include "iforce.h"
#include "curvature.h"

/**
We overload the acceleration() event to add the contribution of
gravity to the interfacial potential $\phi$.

If $\phi$ is already allocated, we add the contribution of gravity,
otherwise we allocate a new field and set it to the contribution of
gravity. */
  
event acceleration (i++)
{
  scalar phi1 = f1.phi;
  coord Gf1 = {0., 0., 0.};
  scalar phi2 = f2.phi;
  coord Gf2 = {0., 0., 0.};

  foreach_dimension(){// for non-coalescing three-phase, phase 3 always makes a precursor layer between 1 and 2
    Gf1.x = (rho3 - rho1)*Bf1.x;
    Gf2.x = (rho3 - rho2)*Bf2.x; 
  }
  
  if (phi1.i)
    position (f1, phi1, Gf1, Z1, add = true);
  else {
    phi1 = new scalar;
    position (f1, phi1, Gf1, Z1, add = false);
    f1.phi = phi1;
  }

  if (phi2.i)
    position (f2, phi2, Gf2, Z2, add = true);
  else {
    phi2 = new scalar;
    position (f2, phi2, Gf2, Z2, add = false);
    f2.phi = phi2;
  }

}