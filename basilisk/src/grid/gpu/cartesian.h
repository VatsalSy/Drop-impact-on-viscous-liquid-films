#if DOUBLE_PRECISION
# define SINGLE_PRECISION 0
#else
# define SINGLE_PRECISION 1
#endif
#define _GPU 1
#define GRIDNAME "Cartesian (GPU)"
#define GRIDPARENT Cartesian
#define field_size() sq(N + 2)
#define grid_data() (cartesian->d)
#define field_offset(s) ((s).i*field_size())
#define depth() 0
#define GPU_CODE()							\
  "#define POINT_VARIABLES VARIABLES\n"					\
  "#define valt(s,k,l,m) _data[_index(s,m)*field_size() + (point.i + (k))*(N + 2) + point.j + (l)]\n" \
  "#define val_red_(s) _data[(s).i*field_size() + (point.i - 1)*NY + point.j - 1]\n"

#include "../cartesian.h"
@define neighborp(k,l,o) neighbor(k,l,o)
#include "gpu.h"
#include "../cartesian-common.h"

static void gpu_cartesian_methods()
{
  cartesian_methods();
  boundary_level = gpu_boundary_level;
}

#include "grid.h"
