/**
# Reductions on GPUs */

#include "utils.h"

int main (int argc, char * argv[])
{
  init_grid (argc > 1 ? atoi(argv[1]) : 1024);

  periodic (right);
  periodic (top);
  
  size (2.*pi);

  scalar s[];
  
  foreach (serial)
    s[] = sq (cos(2.*x)*cos(2.*y));
  
  timer t = timer_start();

  double sum;
  int iter;
  for (iter = 0; iter < 400*1024/N; iter++) {
    sum = 0.;
    foreach(reduction(max:sum))
      sum = max (sum, s[]);
  }
  
  double elapsed = timer_elapsed (t);
  printf ("N: %d elapsed: %g speed: %g\n",
	   N, elapsed, grid->tn*iter/elapsed);

  fprintf (stderr, "result: %g\n", sum);
  output_ppm (s, file = "s.png", n = 512, spread = -1);

  /**
  Check that "inout" fields work. */
  
  foreach()
    s[] += sq (0.5*sin(8.*x)*sin(8.*y));
  foreach()
    s[] = s[] + sq (0.5*sin(8.*x)*sin(8.*y));
    
  stats stat = statsf (s);
  fprintf (stderr, "min: %g max: %g\n", stat.min, stat.max);
  output_ppm (s, file = "s1.png", n = 512, spread = -1);
}
