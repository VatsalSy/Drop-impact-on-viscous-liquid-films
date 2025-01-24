/**
# Diapycnal "entrainment" and "detrainment"

This implements the entrainment/detrainment between isopycnal layers
briefly described in [Hurlburt & Hogan, 2000](#hurlburt2000), section
2, equations 1), 2) and 3).

`omr` is the average entrainment velocity between layers
($\tilde{\omega}_k$ in H&H, 2000) and `Cm` is the coefficient of
additional interfacial friction associated with entrainment ($C_M$ in
H&H, 2000). They are parameters provided by the user. */

extern double omr, Cm;

/**
The maximum and minimum layer thicknesses ($h_k^+$ and $h_k^-$
respectively in H&H, 2000) are provided by the user and used to define
the weighing factors appearing in the definition of the "diapycnal
mixing velocities" $\omega_k^+$ and $\omega_k^-$ in H&H, 2000. */

extern double * hmin, * hmax;
#define wmin(h,hmin) (h >= hmin ? 0. : sq((hmin - h)/hmin))
#define wmax(h,hmax) (h <= hmax ? 0. : sq((h - hmax)/hmax))

/**
These come from the [multilayer solver](hydro.h). */

extern scalar * hl;
extern double dt, dry;

/**
This will store the initial volume of the isopycnal layers, which will
be conserved by entrainment and detrainment. */

static double * vhi = NULL;

event viscous_term (i++)
{

  /**
  If the average entrainment velocity is negative or null, we do
  nothing. */
  
  if (omr <= 0.)
    return 0;

  /**
  We compute and store the initial volume of each layer. */
  
  if (vhi == NULL) {
    vhi = malloc (sizeof(double)*nl);
    foreach_layer()
      vhi[_layer] = statsf(h).sum;
  }

  /**
  We compute the volumed-averaged entrainment/detrainment for each
  layer (bottom layer excepted). */
  
  double oma[nl], wa[nl];
  for (int l = 0; l < nl; l++)
    oma[l] = wa[l] = 0.;
  foreach(reduction(+:oma[:nl]) reduction(+:wa[:nl]))
    for (int l = nl - 1; l >= 1; l--) {
      double om = omr*(wmin((h[0,0,l]), hmin[l]) - wmax((h[0,0,l]), hmax[l]));
      if (h[0,0,l] + dt*om > dry && h[0,0,l-1] - dt*om > dry)
	oma[l] += dv()*om, wa[l] += dv();
    }

  /**
  The local entrainment/detrainment is then applied as the sum of a
  global contribution and a local contribution which depends on the
  local thickness. */
  
  foreach() {
    double hn[nl];
    coord un[nl];
    foreach_layer() {
      hn[point.l] = h[];
      foreach_dimension()
	un[point.l].x = h[]*u.x[];
    }
    for (int l = nl - 1; l >= 1; l--) {
      double om = omr*(wmin((h[0,0,l]), hmin[l]) - wmax((h[0,0,l]), hmax[l]));
      if (h[0,0,l] + dt*om > dry && h[0,0,l-1] - dt*om > dry) {
	om -= oma[l]/wa[l];
	hn[l] += dt*om;
	hn[l - 1] -= dt*om;

	/**
        These are the entrainment and detrainment terms for the
        velocity components in eqs. 1) and 2) of H&H, 2000. */
	
	int lum = om >= 0 ? l - 1 : l;
	foreach_dimension() {
	  double flux = dt*om*(u.x[0,0,lum]*(1. + Cm) - Cm*u.x[0,0,l]);
	  un[l].x     += flux;
	  un[l - 1].x -= flux;
	}
      }
    }

    /**
    We then update `h` and `u`. */
    
    foreach_layer() {
      h[] = hn[point.l];
      foreach_dimension()
	u.x[] = h[] > dry ? un[point.l].x/h[] : 0.;
    }
  }
  
  /**
  We then enforce the conservation of the volume of each layer using a
  global average flux. */

  double eh[nl];
  foreach_layer()
    eh[_layer] = vhi[_layer]/statsf(h).sum;

  foreach() {
    double hnew[nl];
    coord hu[nl];
    foreach_layer() {
      hnew[point.l] = h[]*eh[point.l];
      foreach_dimension()
	hu[point.l].x = h[]*u.x[];
    }
    for (int l = nl - 1; l >= 1; l--) {
      double dh = hnew[l] - h[0,0,l];
      if ((dh < 0. && h[0,0,l] + dh > dry && h[0,0,l-1] > dry) ||
	  (dh > 0. && h[0,0,l-1] - dh > dry)) {
	h[0,0,l] += dh;
	h[0,0,l-1] -= dh;
	int ul = dh > 0. ? l - 1 : l;
	foreach_dimension() {
	  hu[l].x += dh*u.x[0,0,ul];
	  hu[l - 1].x -= dh*u.x[0,0,ul];
	}
      }
    }
#if 0 // fixme: should be included??
    foreach_layer()
      foreach_dimension()
	u.x[] = h[] > dry ? hu[point.l].x/h[] : 0.;
#endif
  }
}

/**
We free the allocated fields to avoid memory leaks. */

event cleanup (t = end)
{
  if (vhi)
    free (vhi), vhi = NULL;
}

/**
## References

~~~bib
@article{hurlburt2000,
  title={Impact of 1/8 to 1/64 resolution on {G}ulf {S}tream model--data 
         comparisons in basin-scale subtropical {A}tlantic {O}cean models},
  author={Hurlburt, Harley E and Hogan, Patrick J},
  journal={Dynamics of Atmospheres and Oceans},
  volume={32},
  number={3-4},
  pages={283--329},
  year={2000},
  publisher={Elsevier},
  DOI={10.1016/S0377-0265(00)00050-6},
  PDF={https://apps.dtic.mil/sti/tr/pdf/ADA531039.pdf}
}
~~~
*/
