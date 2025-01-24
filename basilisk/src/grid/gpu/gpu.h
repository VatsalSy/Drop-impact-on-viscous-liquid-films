#include <grid/gpu/glad.h>
#include <GLFW/glfw3.h>
#if DEBUG_OPENGL
#include <grid/gpu/debug.h>
#endif

#pragma autolink -L$BASILISK/grid/gpu -lglfw -lgpu -ldl

static struct {
  ///// GPU /////
  GLFWwindow * window;
  GLuint ssbo;
  bool fragment_shader;
  int current_shader;
} GPUContext = { .current_shader = -1 };

static void gpu_check_error (const char * stmt,
			     const char * fname, int line)
{
  GLenum err = glGetError();
  if (err != GL_NO_ERROR) {
    fprintf (stderr, "%s:%d: error: OpenGl %08x for '%s;'\n",
	     fname, line, err, stmt);
    abort();
  }
}

#ifdef NDEBUG
// helper macro that checks for GL errors.
#define GL_C(stmt) do {	stmt; } while (0)
#else
// helper macro that checks for GL errors.
#define GL_C(stmt) do {							\
    stmt;								\
    gpu_check_error (#stmt, __FILE__, LINENO);				\
  } while (0)
#endif

void gpu_free_solver (void)
{
  GL_C (glFinish());
  GL_C (glBindFramebuffer (GL_FRAMEBUFFER, 0));
  glDeleteBuffers (1, &GPUContext.ssbo);
  glfwTerminate();
  GPUContext.window = NULL;
}

static char * getShaderLogInfo (GLuint shader)
{
  char * infoLog = NULL;
  GLint len;
  GL_C (glGetShaderiv (shader, GL_INFO_LOG_LENGTH, &len));
  if (len > 0) {
    GLsizei actualLen;
    infoLog = malloc (len);
    GL_C (glGetShaderInfoLog (shader, len, &actualLen, infoLog));
  }
  return infoLog;
}

static char * getProgramLogInfo (GLuint program)
{
  char * infoLog = NULL;
  GLint len;
  GL_C (glGetProgramiv (program, GL_INFO_LOG_LENGTH, &len));
  if (len > 0) {
    GLsizei actualLen;
    infoLog = malloc (len);
    GL_C (glGetProgramInfoLog (program, len, &actualLen, infoLog));
  }
  return infoLog;
}

char * gpu_errors (const char * errors, const char * source, char * fout);

static GLuint createShaderFromString (const char * shaderSource,
				      const GLenum shaderType)
{
  GLuint shader;

  GL_C (shader = glCreateShader(shaderType));
  GL_C (glShaderSource (shader, 1, &shaderSource, NULL));
  GL_C (glCompileShader (shader));

  GLint compileStatus;
  GL_C (glGetShaderiv (shader, GL_COMPILE_STATUS, &compileStatus));
  if (compileStatus != GL_TRUE) {
    char * info = getShaderLogInfo (shader);
#if PRINTSHADERERROR
    fputs (shaderSource, stderr);
    fputs (info, stderr);
#endif
    char * error = gpu_errors (info, shaderSource, NULL);
    fputs (error, stderr);
    sysfree (error);
    free (info);
    glDeleteShader (shader);
    return 0;
  }

  return shader;
}

static GLuint loadNormalShader (const char * vsSource, const char * fsShader)
{
  GLuint vs = 0;
  if (vsSource) {
    vs = createShaderFromString (vsSource, GL_VERTEX_SHADER);
    if (!vs)
      return 0;
  }
  GLuint fs = createShaderFromString (fsShader, vsSource ? GL_FRAGMENT_SHADER : GL_COMPUTE_SHADER);
  if (!fs)
    return 0;
  
  GLuint shader = glCreateProgram();
  if (vs)
    glAttachShader (shader, vs);
  glAttachShader (shader, fs);
  glLinkProgram (shader);

  GLint Result;
  glGetProgramiv (shader, GL_LINK_STATUS, &Result);

  if (Result == GL_FALSE) {
    char * info = getProgramLogInfo (shader);
    fprintf (stderr, "GLSL: could not link shader \n\n%s\n%s\n%s\n",
	     info, vsSource, fsShader);
    free (info);
    glDeleteProgram (shader);
    shader = 0;
  }

  if (shader) {
    if (vs)
      glDetachShader (shader, vs);
    glDetachShader (shader, fs);
  }

  if (vs)
    glDeleteShader (vs);
  glDeleteShader (fs);
  
  return shader;
}

typedef struct {
  char * s;
  int index;
} GLString;

GLString gpu_limits_list[] = {
  {"GL_MAX_DRAW_BUFFERS", GL_MAX_DRAW_BUFFERS},
  {"GL_MAX_VERTEX_UNIFORM_COMPONENTS", GL_MAX_VERTEX_UNIFORM_COMPONENTS},
  {"GL_MAX_VERTEX_UNIFORM_BLOCKS", GL_MAX_VERTEX_UNIFORM_BLOCKS},
  {"GL_MAX_VERTEX_OUTPUT_COMPONENTS", GL_MAX_VERTEX_OUTPUT_COMPONENTS},
  {"GL_MAX_VERTEX_TEXTURE_IMAGE_UNITS", GL_MAX_VERTEX_TEXTURE_IMAGE_UNITS},
  {"GL_MAX_VERTEX_OUTPUT_COMPONENTS", GL_MAX_VERTEX_OUTPUT_COMPONENTS},
#if 1  
  {"GL_MAX_TESS_CONTROL_UNIFORM_COMPONENTS", GL_MAX_TESS_CONTROL_UNIFORM_COMPONENTS},
  {"GL_MAX_TESS_CONTROL_UNIFORM_BLOCKS", GL_MAX_TESS_CONTROL_UNIFORM_BLOCKS},
  {"GL_MAX_TESS_CONTROL_INPUT_COMPONENTS", GL_MAX_TESS_CONTROL_INPUT_COMPONENTS},
  {"GL_MAX_TESS_CONTROL_OUTPUT_COMPONENTS", GL_MAX_TESS_CONTROL_OUTPUT_COMPONENTS},
  {"GL_MAX_TESS_CONTROL_TEXTURE_IMAGE_UNITS", GL_MAX_TESS_CONTROL_TEXTURE_IMAGE_UNITS},
  {"GL_MAX_TESS_CONTROL_OUTPUT_COMPONENTS", GL_MAX_TESS_CONTROL_OUTPUT_COMPONENTS},
  {"GL_MAX_TESS_EVALUATION_UNIFORM_COMPONENTS", GL_MAX_TESS_EVALUATION_UNIFORM_COMPONENTS},
  {"GL_MAX_TESS_EVALUATION_UNIFORM_BLOCKS", GL_MAX_TESS_EVALUATION_UNIFORM_BLOCKS},
  {"GL_MAX_TESS_EVALUATION_INPUT_COMPONENTS", GL_MAX_TESS_EVALUATION_INPUT_COMPONENTS},
  {"GL_MAX_TESS_EVALUATION_OUTPUT_COMPONENTS", GL_MAX_TESS_EVALUATION_OUTPUT_COMPONENTS},
  {"GL_MAX_TESS_EVALUATION_TEXTURE_IMAGE_UNITS", GL_MAX_TESS_EVALUATION_TEXTURE_IMAGE_UNITS},
  {"GL_MAX_TESS_EVALUATION_OUTPUT_COMPONENTS", GL_MAX_TESS_EVALUATION_OUTPUT_COMPONENTS},
#endif
#if 1  
  {"GL_MAX_COMPUTE_UNIFORM_COMPONENTS", GL_MAX_COMPUTE_UNIFORM_COMPONENTS},
  {"GL_MAX_COMPUTE_UNIFORM_BLOCKS", GL_MAX_COMPUTE_UNIFORM_BLOCKS},
  {"GL_MAX_COMPUTE_TEXTURE_IMAGE_UNITS", GL_MAX_COMPUTE_TEXTURE_IMAGE_UNITS},
#endif
  {"GL_MAX_GEOMETRY_UNIFORM_COMPONENTS", GL_MAX_GEOMETRY_UNIFORM_COMPONENTS},
  {"GL_MAX_GEOMETRY_UNIFORM_BLOCKS", GL_MAX_GEOMETRY_UNIFORM_BLOCKS},
  {"GL_MAX_GEOMETRY_INPUT_COMPONENTS", GL_MAX_GEOMETRY_INPUT_COMPONENTS},
  {"GL_MAX_GEOMETRY_OUTPUT_COMPONENTS", GL_MAX_GEOMETRY_OUTPUT_COMPONENTS},
  {"GL_MAX_GEOMETRY_TEXTURE_IMAGE_UNITS", GL_MAX_GEOMETRY_TEXTURE_IMAGE_UNITS},
  {"GL_MAX_GEOMETRY_OUTPUT_COMPONENTS", GL_MAX_GEOMETRY_OUTPUT_COMPONENTS},
  {"GL_MAX_FRAGMENT_UNIFORM_COMPONENTS", GL_MAX_FRAGMENT_UNIFORM_COMPONENTS},
  {"GL_MAX_FRAGMENT_UNIFORM_BLOCKS", GL_MAX_FRAGMENT_UNIFORM_BLOCKS},
  {"GL_MAX_FRAGMENT_INPUT_COMPONENTS", GL_MAX_FRAGMENT_INPUT_COMPONENTS},
  {"GL_MAX_COMBINED_TEXTURE_IMAGE_UNITS", GL_MAX_COMBINED_TEXTURE_IMAGE_UNITS},
  {"GL_MAX_FRAGMENT_IMAGE_UNIFORMS", GL_MAX_FRAGMENT_IMAGE_UNIFORMS},
  {"GL_MAX_IMAGE_UNITS", GL_MAX_IMAGE_UNITS},
  {NULL}
};

void printWorkGroupsCapabilities()
{
  int workgroup_size[3];
  glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_SIZE, 0, &workgroup_size[0]);
  glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_SIZE, 1, &workgroup_size[1]);
  glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_SIZE, 2, &workgroup_size[2]);
  printf ("Maximal workgroup sizes:\n\tx:%u\n\ty:%u\n\tz:%u\n",
	  workgroup_size[0], workgroup_size[1], workgroup_size[2]);

  int workgroup_count[3];
  glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_COUNT, 0, &workgroup_count[0]);
  glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_COUNT, 1, &workgroup_count[1]);
  glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_COUNT, 2, &workgroup_count[2]);
  printf ("Maximum number of local invocations:\n\tx:%u\n\ty:%u\n\tz:%u\n",
	  workgroup_count[0], workgroup_count[1], workgroup_count[2]);

  int workgroup_invocations;
  glGetIntegerv (GL_MAX_COMPUTE_WORK_GROUP_INVOCATIONS, &workgroup_invocations);
  printf ("Maximum workgroup invocations:\n\t%u\n", workgroup_invocations);
}

#if TRACE == 3
# define tracing_foreach(name, file, line) tracing(name, file, line)
# define end_tracing_foreach(name, file, line) end_tracing(name, file, line)
#else
# define tracing_foreach(name, file, line)
# define end_tracing_foreach(name, file, line)
#endif

static bool _gpu_done_ = false;
@undef BEGIN_FOREACH
@def BEGIN_FOREACH if (_gpu_done_)
  _gpu_done_ = false;
 else {
   tracing_foreach ("foreach", S__FILE__, S_LINENO);
@
@undef END_FOREACH
@define END_FOREACH end_tracing_foreach ("foreach", S__FILE__, S_LINENO); }

typedef struct {
  coord p, * box, n; // region
  int level; // level
} RegionParameters;

@define DEPAREN(...) S__VA_ARGS__

@def foreach_stencil_generic(_parallel, _parameters, _externals, _kernel) {
  tracing_foreach ("foreach", S__FILE__, S_LINENO);
  static ForeachData _loop = {
    .fname = S__FILE__, .line = S_LINENO, .first = 1, .parallel = _parallel
  };
  static const char * _kernel_ = _kernel;
  External _externals_[] = _externals;
  if (baseblock) for (scalar s = baseblock[0], * i = baseblock; s.i >= 0; i++, s = *i) {
    _attribute[s.i].input = _attribute[s.i].output = false;
    _attribute[s.i].width = 0;
  }
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
  Point point = {0}; NOT_UNUSED (point);
@

@undef foreach_stencil
@def foreach_stencil(_parallel, _parameters, _externals, _kernel)
  foreach_stencil_generic(_parallel, DEPAREN(_parameters), DEPAREN(_externals), _kernel)
  RegionParameters parameters = _parameters, * _region = &parameters;
@

@undef foreach_level_stencil
@def foreach_level_stencil(_parallel, _parameters, _externals, _kernel)
  foreach_stencil_generic(_parallel, DEPAREN(_parameters), DEPAREN(_externals), _kernel)
  struct  { int level; } parameters = _parameters;
  RegionParameters _region_ = { .level = parameters.level + 1 }, * _region = &_region_;
@

@undef foreach_point_stencil
@def foreach_point_stencil(_parallel, _parameters, _externals, _kernel)
  foreach_stencil_generic(_parallel, DEPAREN(_parameters), DEPAREN(_externals), _kernel)
  RegionParameters _region_ = { .p = _parameters, .n = {1,1} }, * _region = &_region_;
@

@undef foreach_region_stencil
@def foreach_region_stencil(_parallel, _parameters, _externals, _kernel)
  foreach_stencil(_parallel, DEPAREN(_parameters), DEPAREN(_externals), _kernel)
@

@undef foreach_vertex_stencil
@def foreach_vertex_stencil(_parallel, _parameters, _externals, _kernel)
  foreach_stencil(_parallel, DEPAREN(_parameters), DEPAREN(_externals), _kernel) _loop.vertex = true;
@
 
@undef foreach_face_stencil
@def foreach_face_stencil(_parallel, _parameters, _externals, _kernel)
  foreach_stencil(_parallel, DEPAREN(_parameters), DEPAREN(_externals), _kernel)
@
 
@undef end_foreach_stencil
@def end_foreach_stencil()
#if PRINTIO
  if (baseblock) {
    fprintf (stderr, "%s:%d:", _loop.fname, _loop.line);
    for (scalar s = baseblock[0], * i = baseblock; s.i >= 0; i++, s = *i)
      if (_attribute[s.i].input || _attribute[s.i].output)
	fprintf (stderr, " %s:%d:%c:%d", _attribute[s.i].name, s.i,
		 _attribute[s.i].input && _attribute[s.i].output ? 'a' :
		 _attribute[s.i].input ? 'r' : 'w',
		 _attribute[s.i].width);
    fprintf (stderr, "\n");
  }
#endif // PRINTIO
  check_stencil (&_loop);
  _gpu_done_ = gpu_end_stencil (&_loop, _region, _externals_, _kernel_);
  _loop.first = 0;
  end_tracing_foreach ("foreach", S__FILE__, S_LINENO);
}
@

@undef end_foreach_level_stencil
@define end_foreach_level_stencil() end_foreach_stencil()

@ifndef tracing
  @ def tracing(func, file, line) do {
    if (glFinish) glFinish();
    tracing(func, file, line);
  } while(0) @
  @ def end_tracing(func, file, line) do {
    if (glFinish) glFinish();
    end_tracing(func, file, line);
  } while(0) @
@endif

bool gpu_end_stencil (ForeachData * loop, const RegionParameters * region,
		      External * externals, const char * kernel);

void realloc_ssbo()
{
  GL_C (glBindBuffer (GL_SHADER_STORAGE_BUFFER, GPUContext.ssbo));
  GL_C (glBufferData (GL_SHADER_STORAGE_BUFFER,
		      field_size()*datasize, grid_data(), GL_DYNAMIC_READ));
  GL_C (glBindBuffer (GL_SHADER_STORAGE_BUFFER, 0));
}

static void gpu_cpu_sync_scalar (scalar s, char * sep, GLenum mode);

void realloc_scalar_gpu (int size)
{
  realloc_scalar (size);
  for (scalar s in baseblock)
    if (s.gpu.stored < 0)
      gpu_cpu_sync_scalar (s, NULL, GL_MAP_READ_BIT);
  realloc_ssbo();
}

void gpu_boundary_level (scalar * list, int l)
{
  scalar * list1 = NULL;
  for (scalar s in list)
    if (s.gpu.stored > 0)
      list1 = list_prepend (list1, s);
  if (list1) {
    void cartesian_boundary_level (scalar * list, int l); 
    cartesian_boundary_level (list1, l);
    free (list1);
  }
}

#define realloc_scalar(size) realloc_scalar_gpu (size)
#define foreach_level_or_leaf(...) foreach_level(__VA_ARGS__)
#define foreach_coarse_level(...) foreach_level(__VA_ARGS__)
