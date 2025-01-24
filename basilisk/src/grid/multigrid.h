#if SINGLE_PRECISION
typedef float real;
#else
typedef double real;
#endif

#ifndef GRIDNAME
# define GRIDNAME "Multigrid"
#endif
#define GHOSTS 2

/* By default only one layer of ghost cells is used on the boundary to
   optimise the cost of boundary conditions. */

#ifndef BGHOSTS
@ define BGHOSTS 1
#endif

#define _I     (point.i - GHOSTS)
#define _J     (point.j - GHOSTS)
#define _K     (point.k - GHOSTS)
#define _DELTA (1./(1 << point.level))

typedef struct {
  Grid g;
  char * d;
  size_t field_size;
} Multigrid;

struct _Point {
  int i;
#if dimension > 1
  int j;
#endif
#if dimension > 2
  int k;
#endif
  int level, n;
@ifdef foreach_block
  int l;
  @define _BLOCK_INDEX , point.l
@else
  @define _BLOCK_INDEX
@endif
};
static Point last_point;

#define multigrid ((Multigrid *)grid)

// see /src/notes/index.tm
#if dimension == 1
# define _shift(l) ((1 << (l)) - 1 + 2*GHOSTS*(l))
#elif dimension == 2
# define _shift(l) (((1 << 2*(l)) - 1)/3 + 4*GHOSTS*((1 << (l)) - 1 + GHOSTS*(l)))
#elif dimension == 3
# define _shift(l) (((1 << 3*(l)) - 1)/7 + 2*GHOSTS*((1 << 2*(l)) - 1) + \
		    12*sq(GHOSTS)*((1 << (l)) - 1) + 8*cube(GHOSTS)*(l))
#endif

#define CELL(m,level,i)  (*((Cell *) &m[level][(i)*datasize]))

/***** Cartesian macros *****/
#if dimension == 1
@undef val
@def val(a,k,l,m) (((real *)multigrid->d)[point.i + (k) +
					  _shift (point.level) +
					  _index(a,m)*multigrid->field_size])
@
#elif dimension == 2
@undef val
@def val(a,k,l,m) (((real *)multigrid->d)[point.j + (l) +
					  (point.i + (k))*((1 << point.level) + 2*GHOSTS) +
					  _shift (point.level) +
					  _index(a,m)*multigrid->field_size])
@
#elif dimension == 3
@undef val
@def val(a,l,m,o) (((real *)multigrid->d)[point.k + (o) +
					  ((1 << point.level) + 2*GHOSTS)*
					  (point.j + (m) +
					   (point.i + (l))*((1 << point.level) + 2*GHOSTS)) +
					  _shift (point.level) +
					  _index(a,0)*multigrid->field_size])
@
#endif

/* low-level memory management */
#if dimension == 1
# if BGHOSTS == 1
@define allocated(...) true
# else // BGHOST != 1
@define allocated(k,l,m) (point.i+(k) >= 0 && point.i+(k) < (1 << point.level) + 2*GHOSTS)
# endif // BGHOST != 1
@def allocated_child(k,l,m) (level < depth() &&
                             point.i > 0 && point.i <= (1 << point.level) + 2)
@
#elif dimension == 2
# if BGHOSTS == 1
@define allocated(...) true
# else // BGHOST != 1
@def allocated(k,l,m) (point.i+(k) >= 0 && point.i+(k) < (1 << point.level) + 2*GHOSTS &&
		       point.j+(l) >= 0 && point.j+(l) < (1 << point.level) + 2*GHOSTS)
@
# endif // BGHOST != 1
@def allocated_child(k,l,m)  (level < depth() &&
			      point.i > 0 && point.i <= (1 << point.level) + 2 &&
			      point.j > 0 && point.j <= (1 << point.level) + 2)
@			   
#else // dimension == 3
# if BGHOSTS == 1
@define allocated(...) true
#else // BGHOST != 1
@def allocated(a,l,m) (point.i+(a) >= 0 &&
		       point.i+(a) < (1 << point.level) + 2*GHOSTS &&
		       point.j+(l) >= 0 &&
		       point.j+(l) < (1 << point.level) + 2*GHOSTS &&
		       point.k+(m) >= 0 &&
		       point.k+(m) < (1 << point.level) + 2*GHOSTS)
@
#endif // BGHOST != 1
@def allocated_child(a,l,m)  (level < depth() &&
			      point.i > 0 && point.i <= (1 << point.level) + 2 &&
			      point.j > 0 && point.j <= (1 << point.level) + 2 &&
			      point.k > 0 && point.k <= (1 << point.level) + 2)
@
#endif // dimension == 3

/***** Multigrid variables and macros *****/
@define depth()       (grid->depth)
#if dimension == 1
@def fine(a,k,l,m)
(((real *)multigrid->d)[2*point.i - GHOSTS + (k) +
			_shift (point.level + 1) +
			_index(a,m)*multigrid->field_size])
@
@def coarse(a,k,l,m)
(((real *)multigrid->d)[(point.i + GHOSTS)/2 + (k) +
			_shift (point.level - 1) +
			_index(a,m)*multigrid->field_size])
@
@def POINT_VARIABLES
  VARIABLES
  int level = point.level; NOT_UNUSED(level);
  struct { int x; } child = { 2*((point.i+GHOSTS)%2)-1 }; NOT_UNUSED(child);
  Point parent = point;	NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + GHOSTS)/2;
@
#elif dimension == 2
@def fine(a,k,l,m)
(((real *)multigrid->d)[2*point.j - GHOSTS + (l) +
			(2*point.i - GHOSTS + (k))*((1 << point.level)*2 + 2*GHOSTS) +
			_shift (point.level + 1) +
			_index(a,m)*multigrid->field_size])
@
@def coarse(a,k,l,m)
(((real *)multigrid->d)[(point.j + GHOSTS)/2 + (l) +
			((point.i + GHOSTS)/2 + (k))*((1 << point.level)/2 + 2*GHOSTS) +
			_shift (point.level - 1) +
			_index(a,m)*multigrid->field_size])
@
@def POINT_VARIABLES
  VARIABLES
  int level = point.level; NOT_UNUSED(level);
  struct { int x, y; } child = {
    2*((point.i+GHOSTS)%2)-1, 2*((point.j+GHOSTS)%2)-1
  }; NOT_UNUSED(child);
  Point parent = point;	NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + GHOSTS)/2; parent.j = (point.j + GHOSTS)/2;
@
#elif dimension == 3
@def fine(a,l,m,o)
(((real *)multigrid->d)[2*point.k - GHOSTS + (o) +
			((1 << point.level)*2 + 2*GHOSTS)*
			(2*point.j - GHOSTS + (m) +
			 (2*point.i - GHOSTS + (l))*((1 << point.level)*2 + 2*GHOSTS)) +
			_shift (point.level + 1) +
			_index(a,0)*multigrid->field_size])
@
@def coarse(a,l,m,o)
(((real *)multigrid->d)[(point.k + GHOSTS)/2 + (o) +
			((1 << point.level)/2 + 2*GHOSTS)*
			((point.j + GHOSTS)/2 + (m) +
			 ((point.i + GHOSTS)/2 + (l))*((1 << point.level)/2 + 2*GHOSTS)) +
			_shift (point.level - 1) +
			_index(a,0)*multigrid->field_size])
@
@def POINT_VARIABLES
  VARIABLES
  int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + GHOSTS)%2) - 1,
    2*((point.j + GHOSTS)%2) - 1,
    2*((point.k + GHOSTS)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point;	NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + GHOSTS)/2;
  parent.j = (point.j + GHOSTS)/2;
  parent.k = (point.k + GHOSTS)/2;
@
#endif

@def foreach_level(l)
OMP_PARALLEL() {
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
  Point point = {0};
  point.level = l; point.n = 1 << point.level;
  int _k;
  OMP(omp for schedule(static))
  for (_k = GHOSTS; _k < point.n + GHOSTS; _k++) {
    point.i = _k;
#if dimension > 1
    for (point.j = GHOSTS; point.j < point.n + GHOSTS; point.j++)
#if dimension > 2
      for (point.k = GHOSTS; point.k < point.n + GHOSTS; point.k++)
#endif
	{
#endif
          POINT_VARIABLES
@
@def end_foreach_level()
#if dimension > 1
	}
#endif
  }
}
@

@def foreach()
  OMP_PARALLEL() {
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
  Point point = {0};
  point.level = depth(); point.n = 1 << point.level;
  int _k;
  OMP(omp for schedule(static))
  for (_k = GHOSTS; _k < point.n + GHOSTS; _k++) {
    point.i = _k;
#if dimension > 1
    for (point.j = GHOSTS; point.j < point.n + GHOSTS; point.j++)
#if dimension > 2
      for (point.k = GHOSTS; point.k < point.n + GHOSTS; point.k++)
#endif
	{
#endif
          POINT_VARIABLES
@
@def end_foreach()
#if dimension > 1
	}
#endif
  }
}
@	    

@define is_active(cell) (true)
@define is_leaf(cell)   (level == depth())
@define is_local(cell)  (true)
@define leaf            2
@def refine_cell(...) do {
  fprintf (stderr, "grid depths do not match. Aborting.\n");
  assert (0);
} while (0)
@
@define tree multigrid
#include "foreach_cell.h"

@def foreach_face_generic()
  OMP_PARALLEL() {
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
  Point point = {0};
  point.level = depth(); point.n = 1 << point.level;
  int _k;
  OMP(omp for schedule(static))
  for (_k = GHOSTS; _k <= point.n + GHOSTS; _k++) {
    point.i = _k;
#if dimension > 1
    for (point.j = GHOSTS; point.j <= point.n + GHOSTS; point.j++)
#if dimension > 2
      for (point.k = GHOSTS; point.k <= point.n + GHOSTS; point.k++)
#endif
        {
#endif
	  POINT_VARIABLES
@
@def end_foreach_face_generic()
#if dimension > 1
	}
#endif
  }
}
@ 

@def foreach_vertex()
foreach_face_generic() {  
  x -= Delta/2.;
#if dimension > 1  
  y -= Delta/2.;  
#endif
#if dimension > 2
  z -= Delta/2.;  
#endif
@
@define end_foreach_vertex() } end_foreach_face_generic()

@define is_coarse() (point.level < depth())

#if dimension == 1
@define is_face_x() { int ig = -1; VARIABLES; {
@define end_is_face_x() }}

// foreach_edge?

@def foreach_child() {
  int _i = 2*point.i - GHOSTS;
  point.level++;
  point.n *= 2;
  for (int _k = 0; _k < 2; _k++) {
    point.i = _i + _k;
    POINT_VARIABLES;
@
@def end_foreach_child()
  }
  point.i = (_i + GHOSTS)/2;
  point.level--;
  point.n /= 2;
}
@
@define foreach_child_break() _k = 2

#elif dimension == 2
#define foreach_edge() foreach_face(y,x)

@define is_face_x() { int ig = -1; VARIABLES; if (point.j < point.n + GHOSTS) {
@define end_is_face_x() }}
@define is_face_y() { int jg = -1; VARIABLES; if (point.i < point.n + GHOSTS) {
@define end_is_face_y() }}
				 
@def foreach_child() {
  int _i = 2*point.i - GHOSTS, _j = 2*point.j - GHOSTS;
  point.level++;
  point.n *= 2;
  for (int _k = 0; _k < 2; _k++)
    for (int _l = 0; _l < 2; _l++) {
      point.i = _i + _k; point.j = _j + _l;
      POINT_VARIABLES;
@
@def end_foreach_child()
  }
  point.i = (_i + GHOSTS)/2; point.j = (_j + GHOSTS)/2;
  point.level--;
  point.n /= 2;
}
@
@define foreach_child_break() _k = _l = 2

#elif dimension == 3
@def foreach_vertex_aux()
foreach_vertex() {
  struct { int x, y, z; } _a = {point.i, point.j, point.k};
@
@define end_foreach_vertex_aux() } end_foreach_vertex()

#define foreach_edge()					\
  foreach_vertex_aux()					\
    foreach_dimension()					\
      if (_a.x < point.n + GHOSTS)

@define is_face_x() { int ig = -1; VARIABLES; if (point.j < point.n + GHOSTS && point.k < point.n + GHOSTS) {
@define end_is_face_x() }}
@define is_face_y() { int jg = -1; VARIABLES; if (point.i < point.n + GHOSTS && point.k < point.n + GHOSTS) {
@define end_is_face_y() }}
@define is_face_z() { int kg = -1; VARIABLES; if (point.i < point.n + GHOSTS && point.j < point.n + GHOSTS) {
@define end_is_face_z() }}
  
@def foreach_child() {
  int _i = 2*point.i - GHOSTS;
  int _j = 2*point.j - GHOSTS;
  int _k = 2*point.k - GHOSTS;
  point.level++;
  point.n *= 2;
  for (int _l = 0; _l < 2; _l++)
    for (int _m = 0; _m < 2; _m++)
      for (int _n = 0; _n < 2; _n++) {
	point.i = _i + _l; point.j = _j + _m; point.k = _k + _n;
	POINT_VARIABLES;
@
@def end_foreach_child()
  }
  point.i = (_i + GHOSTS)/2;
  point.j = (_j + GHOSTS)/2;
  point.k = (_k + GHOSTS)/2;
  point.level--;
  point.n /= 2;
}
@
@define foreach_child_break() _l = _m = _n = 2
#endif
  
@if TRASH
@ undef trash
@ define trash(list) reset(list, undefined)
@endif

#include "neighbors.h"

void reset (void * alist, double val)
{
  scalar * list = (scalar *) alist;
  for (scalar s in list)
    if (!is_constant(s))
      for (int b = 0; b < s.block; b++) {
	real * data = (real *) multigrid->d;
	data += (s.i + b)*multigrid->field_size;
	for (int i = 0; i < multigrid->field_size; i++, data++)
	  *data = val;
      }
}

// Boundaries

#if dimension == 1
@def foreach_boundary_dir(l,d)
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
  Point point = {0};
  point.level = l < 0 ? depth() : l;
  point.n = 1 << point.level;
  if (d == left) {
    point.i = GHOSTS;
    ig = -1;
  }
  else if (d == right) {
    point.i = point.n + GHOSTS - 1;
    ig = 1;
  }
  {
    POINT_VARIABLES
@
@define end_foreach_boundary_dir() }

@define neighbor(o,p,q) ((Point){point.i+o, point.level, point.n _BLOCK_INDEX})
@define is_boundary(point) (point.i < GHOSTS || point.i >= point.n + GHOSTS)

#elif dimension == 2
@def foreach_boundary_dir(l,d)
  OMP_PARALLEL() {
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
  Point point = {0};
  point.level = l < 0 ? depth() : l;
  point.n = 1 << point.level;
  int * _i = &point.j;
  if (d == left) {
    point.i = GHOSTS;
    ig = -1;
  }
  else if (d == right) {
    point.i = point.n + GHOSTS - 1;
    ig = 1;
  }
  else if (d == bottom) {
    point.j = GHOSTS;
    _i = &point.i;
    jg = -1;
  }
  else if (d == top) {
    point.j = point.n + GHOSTS - 1;
    _i = &point.i;
    jg = 1;
  }
  int _l;
  OMP(omp for schedule(static))
  for (_l = 0; _l < point.n + 2*GHOSTS; _l++) {
    *_i = _l;
    {
      POINT_VARIABLES
@
@def end_foreach_boundary_dir()
    }
  }
}
@

@def neighbor(o,p,q)
  ((Point){point.i+o, point.j+p, point.level, point.n _BLOCK_INDEX})
@
@def is_boundary(point) (point.i < GHOSTS || point.i >= point.n + GHOSTS ||
			 point.j < GHOSTS || point.j >= point.n + GHOSTS)
@

#elif dimension == 3
@def foreach_boundary_dir(l,d)
  OMP_PARALLEL() {
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
  Point point = {0};
  point.level = l < 0 ? depth() : l;
  point.n = 1 << point.level;
  int * _i = &point.j, * _j = &point.k;
  if (d == left) {
    point.i = GHOSTS;
    ig = -1;
  }
  else if (d == right) {
    point.i = point.n + GHOSTS - 1;
    ig = 1;
  }
  else if (d == bottom) {
    point.j = GHOSTS;
    _i = &point.i;
    jg = -1;
  }
  else if (d == top) {
    point.j = point.n + GHOSTS - 1;
    _i = &point.i;
    jg = 1;
  }
  else if (d == back) {
    point.k = GHOSTS;
    _i = &point.i; _j = &point.j;
    kg = -1;
  }
  else if (d == front) {
    point.k = point.n + GHOSTS - 1;
    _i = &point.i; _j = &point.j;
    kg = 1;
  }
  int _l;
  OMP(omp for schedule(static))
  for (_l = 0; _l < point.n + 2*GHOSTS; _l++) {
    *_i = _l;
    for (int _m = 0; _m < point.n + 2*GHOSTS; _m++) {
      *_j = _m;
      POINT_VARIABLES
@
@def end_foreach_boundary_dir()
    }
  }
}
@

@def neighbor(o,p,q)
  ((Point){point.i+o, point.j+p, point.k+q, point.level, point.n _BLOCK_INDEX})
@
@def is_boundary(point) (point.i < GHOSTS || point.i >= point.n + GHOSTS ||
			 point.j < GHOSTS || point.j >= point.n + GHOSTS ||
			 point.k < GHOSTS || point.k >= point.n + GHOSTS)
@

#endif // dimension == 3

@def foreach_boundary(b)
  if (default_scalar_bc[b] != periodic_bc)
    foreach_boundary_dir (depth(), b)
      if (!is_boundary(point)) {
@
@define end_foreach_boundary() } end_foreach_boundary_dir()
  
@define neighborp(k,l,o) neighbor(k,l,o)

static double periodic_bc (Point point, Point neighbor, scalar s, bool * data);

static void box_boundary_level (const Boundary * b, scalar * scalars, int l)
{
  extern double (* default_scalar_bc[]) (Point, Point, scalar, bool *);
  disable_fpe (FE_DIVBYZERO|FE_INVALID);
  for (int d = 0; d < 2*dimension; d++)
    if (default_scalar_bc[d] == periodic_bc)
      for (scalar s in scalars)
	if (!is_constant(s) && s.block > 0) {
	  if (is_vertex_scalar (s))
	    s.boundary[d] = s.boundary_homogeneous[d] = NULL;
	  else if (s.face) {
	    vector v = s.v;
	    v.x.boundary[d] = v.x.boundary_homogeneous[d] = NULL;
	  }
	}
  for (int bghost = 1; bghost <= BGHOSTS; bghost++)
    for (int d = 0; d < 2*dimension; d++) {

      scalar * list = NULL, * listb = NULL;
      for (scalar s in scalars)
	if (!is_constant(s) && s.block > 0) {
	  scalar sb = s;
#if dimension > 1
	  if (s.v.x.i >= 0) {
	    // vector component
	    int j = 0;
	    while ((&s.v.x)[j].i != s.i) j++;
	    sb = (&s.v.x)[(j - d/2 + dimension) % dimension];
	  }
#endif
	  if (sb.boundary[d] && sb.boundary[d] != periodic_bc) {
	    list = list_append (list, s);
	    listb = list_append (listb, sb);
	  }
	}
      
      if (list) {
	foreach_boundary_dir (l, d) {
	  scalar s, sb;
	  for (s,sb in list,listb) {
	    if ((s.face && sb.i == s.v.x.i) || is_vertex_scalar (s)) {
	      // normal component of face vector, or vertex scalar
	      if (bghost == 1)
		foreach_block()
		  s[(ig + 1)/2,(jg + 1)/2,(kg + 1)/2] =
		  sb.boundary[d] (point, neighborp(ig,jg,kg), s, NULL);
	    }
	    else
	      // tangential component of face vector or centered
	      foreach_block()
		s[bghost*ig,bghost*jg,bghost*kg] =
		sb.boundary[d] (neighborp((1 - bghost)*ig,
					  (1 - bghost)*jg,
					  (1 - bghost)*kg),
				neighborp(bghost*ig,bghost*jg,bghost*kg),
				s, NULL);
	  }
	}
	free (list);
	free (listb);
      }
    }
  enable_fpe (FE_DIVBYZERO|FE_INVALID);
}

/* Periodic boundaries */

#if !_MPI

#if dimension == 1

static void periodic_boundary_level_x (const Boundary * b, scalar * list, int l)
{
  scalar * list1 = NULL;
  for (scalar s in list)
    if (!is_constant(s) && s.block > 0 && s.boundary[right] == periodic_bc)
      list1 = list_add (list1, s);
  if (!list1)
    return;

  if (l == 0) {
    foreach_level(0)
      for (scalar s in list1)
	for (int b = 0; b < s.block; b++) {
	  scalar sb = {s.i + b};
	  real v = sb[];
	  foreach_neighbor()
	    sb[] = v;
	}
    free (list1);
    return;
  }

  Point point = {0};
  point.level = l < 0 ? depth() : l; point.n = 1 << point.level;
  for (int i = 0; i < GHOSTS; i++)
    for (scalar s in list1)
      for (int b = 0; b < s.block; b++) {
	scalar sb = {s.i + b};
	sb[i] = sb[i + point.n];
      }
  for (int i = point.n + GHOSTS; i < point.n + 2*GHOSTS; i++)
    for (scalar s in list1)
      for (int b = 0; b < s.block; b++) {
	scalar sb = {s.i + b};
	sb[i] = sb[i - point.n];
      }
  free (list1);
}
    
#else // dimension != 1
  
@define VT _attribute[s.i].v.y

foreach_dimension()
static void periodic_boundary_level_x (const Boundary * b, scalar * list, int l)
{
  scalar * list1 = NULL;
  for (scalar s in list)
    if (!is_constant(s) && s.block > 0) {
      if (s.face) {
	scalar vt = VT;
	if (vt.boundary[right] == periodic_bc)
	  list1 = list_add (list1, s);
      }
      else if (s.boundary[right] == periodic_bc)
	list1 = list_add (list1, s);
    }
  if (!list1)
    return;

  if (l == 0) {
    foreach_level(0)
      for (scalar s in list1)
	for (int b = 0; b < s.block; b++) {
	  scalar sb = {s.i + b};
	  real v = sb[];
	  foreach_neighbor()
	    sb[] = v;
	}
    free (list1);
    return;
  }
  
  OMP_PARALLEL() {
    Point point = {0};
    point.level = l < 0 ? depth() : l; point.n = 1 << point.level;
#if dimension == 2  
    int j;
    OMP(omp for schedule(static))
      for (j = 0; j < point.n + 2*GHOSTS; j++) {
	for (int i = 0; i < GHOSTS; i++)
	  for (scalar s in list1)
	    for (int b = 0; b < s.block; b++) {
	      scalar sb = {s.i + b};
	      sb[i,j] = sb[i + point.n,j];
	    }
	for (int i = point.n + GHOSTS; i < point.n + 2*GHOSTS; i++)
	  for (scalar s in list1)
	    for (int b = 0; b < s.block; b++) {
	      scalar sb = {s.i + b};
	      sb[i,j] = sb[i - point.n,j];
	    }
      }
#else // dimension == 3
    int j;
    OMP(omp for schedule(static))
      for (j = 0; j < point.n + 2*GHOSTS; j++)
	for (int k = 0; k < point.n + 2*GHOSTS; k++) {
	  for (int i = 0; i < GHOSTS; i++)
	    for (scalar s in list1)
	      for (int b = 0; b < s.block; b++) {
		scalar sb = {s.i + b};
		sb[i,j,k] = sb[i + point.n,j,k];
	      }
	  for (int i = point.n + GHOSTS; i < point.n + 2*GHOSTS; i++)
	    for (scalar s in list1)
	      for (int b = 0; b < s.block; b++) {
		scalar sb = {s.i + b};
		sb[i,j,k] = sb[i - point.n,j,k];
	      }
	}
#endif
  }
  free (list1);
}

@undef VT

#endif // dimension != 1  
  
#endif // !_MPI

void free_grid (void)
{
  if (!grid)
    return;
  free_boundaries();
  Multigrid * m = multigrid;
  free (m->d);
  free (m);
  grid = NULL;
}

int log_base2 (int n) {
  int m = n, r = 0;
  while (m > 1)
    m /= 2, r++;
  return (1 << r) < n ? r + 1 : r;
}
 
void init_grid (int n)
{
  free_grid();
  Multigrid * m = qmalloc (1, Multigrid);
  grid = (Grid *) m;
  grid->depth = grid->maxdepth = log_base2(n);
  N = 1 << depth();
  // mesh size
  grid->n = grid->tn = 1 << dimension*depth();
  // box boundaries
  Boundary * b = qcalloc (1, Boundary);
  b->level = box_boundary_level;
  add_boundary (b);
#if _MPI
  Boundary * mpi_boundary_new();
  mpi_boundary_new();
#else
  // periodic boundaries
  foreach_dimension() {
    Boundary * b = qcalloc (1, Boundary);
    b->level = periodic_boundary_level_x;
    add_boundary (b);
  }
#endif
  // allocate grid: this must be after mpi_boundary_new() since this modifies depth()
  m->field_size = _shift (depth() + 1);
  m->d = (char *) malloc(m->field_size*datasize);
  reset (all, 0.);
}

void realloc_scalar (int size)
{
  Multigrid * p = multigrid;
  datasize += size;
  qrealloc (p->d, p->field_size*datasize, char);
}

#if _MPI
int mpi_dims[dimension], mpi_coords[dimension];
#undef _DELTA
#undef _I
#undef _J
#undef _K
#define _DELTA (1./(1 << point.level)/mpi_dims[0])
#define _I     (point.i - GHOSTS + mpi_coords[0]*(1 << point.level))
#define _J     (point.j - GHOSTS + mpi_coords[1]*(1 << point.level))
#define _K     (point.k - GHOSTS + mpi_coords[2]*(1 << point.level))
#endif



Point locate (double xp = 0, double yp = 0, double zp = 0)
{
  Point point = {0};
  point.level = -1, point.n = 1 << depth();
#if _MPI
  point.i = (xp - X0)/L0*point.n*mpi_dims[0] + GHOSTS - mpi_coords[0]*point.n;
  if (point.i < GHOSTS || point.i >= point.n + GHOSTS)
    return point;
#if dimension >= 2
  point.j = (yp - Y0)/L0*point.n*mpi_dims[0] + GHOSTS - mpi_coords[1]*point.n;
  if (point.j < GHOSTS || point.j >= point.n + GHOSTS)
    return point;
#endif
#if dimension >= 3
  point.k = (zp - Z0)/L0*point.n*mpi_dims[0] + GHOSTS - mpi_coords[2]*point.n;
  if (point.k < GHOSTS || point.k >= point.n + GHOSTS)
    return point;
#endif  
#else // !_MPI
  point.i = (xp - X0)/L0*point.n + GHOSTS;
  if (point.i < GHOSTS || point.i >= point.n + GHOSTS)
    return point;
#if dimension >= 2
  point.j = (yp - Y0)/L0*point.n + GHOSTS;
  if (point.j < GHOSTS || point.j >= point.n + GHOSTS)
    return point;
#endif
#if dimension >= 3
  point.k = (zp - Z0)/L0*point.n + GHOSTS;
  if (point.k < GHOSTS || point.k >= point.n + GHOSTS)
    return point;
#endif  
#endif // !_MPI
  point.level = depth();
  return point;
}

#if !_GPU
# include "multigrid-common.h"
#endif

void dimensions (int nx = 0, int ny = 0, int nz = 0)
{
#if _MPI
  int p[] = {nx, ny, nz};
  for (int i = 0; i < dimension; i++)
    mpi_dims[i] = p[i];
#endif
}

#if _MPI

@if dimension == 1

@def foreach_slice_x(start, end, l) {
  Point point = {0};
  point.level = l; point.n = 1 << point.level;
  for (point.i = start; point.i < end; point.i++)
@
@define end_foreach_slice_x() }
 
@elif dimension == 2

@def foreach_slice_x(start, end, l) {
  Point point = {0};
  point.level = l; point.n = 1 << point.level;
  for (point.i = start; point.i < end; point.i++)
    for (point.j = 0; point.j < point.n + 2*GHOSTS; point.j++)
@
@define end_foreach_slice_x() }
      
@def foreach_slice_y(start, end, l) {
  Point point = {0};
  point.level = l; point.n = 1 << point.level;
  for (point.i = 0; point.i < point.n + 2*GHOSTS; point.i++)
    for (point.j = start; point.j < end; point.j++)
@
@define end_foreach_slice_y() }

@elif dimension == 3

@def foreach_slice_x(start, end, l) {
  Point point = {0};
  point.level = l; point.n = 1 << point.level;
  for (point.i = start; point.i < end; point.i++)
    for (point.j = 0; point.j < point.n + 2*GHOSTS; point.j++)
      for (point.k = 0; point.k < point.n + 2*GHOSTS; point.k++)
@
@define end_foreach_slice_x() }
      
@def foreach_slice_y(start, end, l) {
  Point point = {0};
  point.level = l; point.n = 1 << point.level;
  for (point.i = 0; point.i < point.n + 2*GHOSTS; point.i++)
    for (point.j = start; point.j < end; point.j++)
      for (point.k = 0; point.k < point.n + 2*GHOSTS; point.k++)
@
@define end_foreach_slice_y() }

@def foreach_slice_z(start, end, l) {
  Point point = {0};
  point.level = l; point.n = 1 << point.level;
  for (point.i = 0; point.i < point.n + 2*GHOSTS; point.i++)
    for (point.j = 0; point.j < point.n + 2*GHOSTS; point.j++)
      for (point.k = start; point.k < end; point.k++)
@
@define end_foreach_slice_z() }
  
@endif // dimension == 3
 
#include "multigrid-mpi.h"
#endif // _MPI
