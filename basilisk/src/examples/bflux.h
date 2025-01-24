/**
# Boundary fluxes for the Gulf-Stream simulation

This implements the "ports" described in [Hurlburt & Hogan,
2000](#hulrburt2000), section 3, page 289:

*"The nonlinear simulations also include effects of the global
thermohaline circulation via ports in the northern and southern model
boundaries, including a southward DWBC and a northward upper ocean
return flow."*

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

These functions are non-zero at the locations (lon/lat) of the
southern and northern ports. These locations are not precisely
specified in H&H, 2000 (see Note b of Table 2 in H&H, 2000). */

#define northern_flux() (x > -50 && x < - 40 && y > 50.5 && val(zbs,0,0,0) < -4000)
#define southern_flux() (x > -60 && x < - 50 && y < 9.5 && val(zbs,0,0,0) < - 4000)

event viscous_term (i++)
{

  /**
  In order to distribute the fluxes within each layer, we first
  compute the volume of the northern and southern ports, for each
  layer. */
  
  double sht[nl], shb[nl];
  for (int l = 0; l < nl; l++)
    sht[l] = shb[l] = 0.;
  foreach(reduction(+:sht[:nl]) reduction(+:shb[:nl]))
    foreach_layer() {
      sht[point.l] += northern_flux()*fmax(h[] - hmin[point.l]/10., 0.)*dv();
      shb[point.l] += southern_flux()*fmax(h[] - hmin[point.l]/10., 0.)*dv();
    }

  /**
  The fluxes (in m^3^/s i.e. Sverdrups/10^6^) are given in Table 2 of
  H&H, 2000. 

  See also the 2nd paragraph page 295 which discusses in more detail
  the chosen fluxes and their control of the northward, southward
  (Deep Western Boundary Current) and upward currents. */

#if NODWBC
  static const double northern[] =
    { 0., - 0.33e6, - 2.33e6, - 4.84e6, - 6.5e6 };
  static const double southern[] =
    { 0.,   0.33e6,   2.33e6,   4.84e6,   6.5e6 };
#else
  static const double factor = 0.9;
  static const double northern[] =
    { + 14e6, - 0.33e6, - 2e6*factor - 0.33e6, - 4.5e6*factor - 0.34e6, - 6.5e6*factor };
  static const double southern[] =
    { - 13e6, 0., 2e6*factor, 4.5e6*factor, 6.5e6*factor };
#endif

  /**
  The fluxes are imposed using a thickness-weighted sum over the
  northern and southern ports. */
  
  assert (nl == 5);
  foreach()
    foreach_layer() {
      h[] += dt*northern_flux()*fmax(h[] - hmin[point.l]/10., 0.)*
	northern[point.l]/max(sht[point.l], dry);
      h[] += dt*southern_flux()*fmax(h[] - hmin[point.l]/10., 0.)*
	southern[point.l]/max(shb[point.l], dry);
      if (h[] < 0.)
	h[] = 0.;
    }
}
