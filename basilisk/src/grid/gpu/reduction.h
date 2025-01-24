real cpu_reduction (GLuint src, size_t offset, size_t nb, const char op)
{
  GL_C (glBindBuffer (GL_SHADER_STORAGE_BUFFER, src));
  GL_C (glMemoryBarrier (GL_BUFFER_UPDATE_BARRIER_BIT));
  real result = 0., * a = glMapBufferRange (GL_SHADER_STORAGE_BUFFER,
					    offset*sizeof(real), nb*sizeof(real),
					    GL_MAP_READ_BIT);
  assert (a);
  switch (op) {
  case '+':
    for (int i = 0; i < nb; i++, a++)
      result += *a;
    break;
  case 'M':
    result = *a++;
    for (int i = 1; i < nb; i++, a++)
      result = max (result, *a);
    break;
  case 'm':
    result = *a++;
    for (int i = 1; i < nb; i++, a++)
      result = min (result, *a);
    break;
  default: assert (false);
  }
  assert (glUnmapBuffer (GL_SHADER_STORAGE_BUFFER));
  GL_C (glBindBuffer (GL_SHADER_STORAGE_BUFFER, 0));
  return result;
}

real gpu_reduction (size_t offset, const char op, const RegionParameters * region, size_t nb)
{
  const int stride = 64, nwgr = 64;
  bool is_foreach_point = (region->n.x == 1 && region->n.y == 1);
  if (!is_foreach_point && nb < nwgr*stride)
    return cpu_reduction (GPUContext.ssbo, offset, nb, op);
  
  GLuint * br = gpu_grid->reduct;
  if (!br[0]) {
    GL_C (glGenBuffers (2, br));
    for (int i = 0; i < 2; i++) {
      GL_C (glBindBuffer (GL_SHADER_STORAGE_BUFFER, br[i]));
      GL_C (glBufferData (GL_SHADER_STORAGE_BUFFER,
			  (sq(N + 1)/stride + 1)*sizeof(real),
			  NULL, GL_DYNAMIC_READ));
    }
    GL_C (glBindBuffer (GL_SHADER_STORAGE_BUFFER, 0));
  }
  
  const char * start, * operation;
  static const char * startsum = "reduct = 0.;",
    * opsum = "reduct += val;";
  static const char * startmin = "reduct = val;",
    * opmin = "if (val < reduct) reduct = val;";
  static const char * startmax = "reduct = val;",
    * opmax = "if (val > reduct) reduct = val;";
  switch (op) {
  case '+': start = startsum, operation = opsum; break;
  case 'M': start = startmax, operation = opmax; break;
  case 'm': start = startmin, operation = opmin; break;
  default: // unknown reduction operation
    assert (false);
  }

  char nwgrs[20], strides[20];
  snprintf (nwgrs, 19, "%d", nwgr);
  snprintf (strides, 19, "%d", stride);
  
  char * fs =
    str_append (NULL,
		"#version 430\n", glsl_preproc,
		"layout (std430, binding = 0) readonly buffer _data_layout {"
		" real _data[]; };\n"
		"layout (std430, binding = 1) writeonly buffer _reduct_layout {"
		" real _reduct[]; };\n"
		"layout (location = 3) uniform uint offset;\n"
		"layout (location = 4) uniform uint nb;\n"
		"layout (location = 5) uniform uint nbr;\n"
		"layout (local_size_x = ", nwgrs, ") in;\n"
		"void main() {\n"
		"if (gl_GlobalInvocationID.x < nb) {\n"
		"  uint stride = ", strides, ";\n"
		"  uint index = stride*gl_GlobalInvocationID.x;\n"
		"  if (index + stride > nbr) stride = nbr - index;\n"
		"  index += offset;\n"
		"  real val = _data[index];\n"
		"  real ", start, "\n"
		"  for (uint j = 0; j < stride; j++) {\n"
		"    val = _data[index + j];\n"
		"    ", operation, "\n"
		"  }\n"
		"  _reduct[gl_GlobalInvocationID.x] = reduct;\n"
		"}}\n");
  Shader * shader;
  Adler32Hash hash;
  a32_hash_init (&hash);
  a32_hash_add (&hash, fs, strlen (fs));
  shader = load_shader (fs, a32_hash (&hash), NULL);
  assert (shader);
  if (shader->id != GPUContext.current_shader) {
    GL_C (glUseProgram (shader->id));
    GPUContext.current_shader = shader->id;
  }
  const GLint loffset = 3, lnb = 4, lnbr = 5;
  
  if (is_foreach_point) {
    real result = 0.;
    int i = (region->p.x - X0)/L0*N;
    int j = (region->p.y - Y0)/L0*N;
    if (i >= 0 && i < N && j >= 0 && j < N) {
      offset += i*N + j;
      GL_C (glUniform1ui (loffset, offset));
      GL_C (glBindBufferBase (GL_SHADER_STORAGE_BUFFER, 0, GPUContext.ssbo));
      GL_C (glBindBufferBase (GL_SHADER_STORAGE_BUFFER, 1, br[0]));
      GL_C (glUniform1ui (lnbr, 1));
      GL_C (glUniform1ui (lnb, 1));
      GL_C (glMemoryBarrier (GL_SHADER_STORAGE_BARRIER_BIT));
      GL_C (glDispatchCompute (1, 1, 1));
      GL_C (glBindBuffer (GL_SHADER_STORAGE_BUFFER, br[0]));
      GL_C (glGetBufferSubData (GL_SHADER_STORAGE_BUFFER, 0, sizeof (real), &result));
      GL_C (glBindBuffer (GL_SHADER_STORAGE_BUFFER, 0));
    }
    return result;
  }
  
  GL_C (glUniform1ui (loffset, offset));
  GL_C (glBindBufferBase (GL_SHADER_STORAGE_BUFFER, 0, GPUContext.ssbo));  
  int src = 0, dst = 1;
  while (nb >= nwgr*stride) {
    GL_C (glBindBufferBase (GL_SHADER_STORAGE_BUFFER, 1, br[dst]));      
    GL_C (glUniform1ui (lnbr, nb));
    if (nb % stride) {
      nb /= stride;
      nb++;
    }
    else
      nb /= stride;
    int ng = nb/nwgr;
    if (ng*nwgr < nb)
      ng++;
    GL_C (glUniform1ui (lnb, nb));
    GL_C (glMemoryBarrier (GL_SHADER_STORAGE_BARRIER_BIT));
    GL_C (glDispatchCompute (ng, 1, 1));
    swap (int, src, dst);
    if (offset) {
      GL_C (glUniform1ui (loffset, 0));
      offset = 0;
    }
    GL_C (glBindBufferBase (GL_SHADER_STORAGE_BUFFER, 0, br[src]));
  }

  return cpu_reduction (br[src], 0, nb, op);
}
