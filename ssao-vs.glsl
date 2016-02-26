attribute vec3 position;
varying vec3 v_position;
varying vec2 v_texcoord;

void main() {
  v_position = position;
  v_texcoord = 0.5 * (position.xy + 1.0);
  gl_Position = vec4(v_position, 1.0);
}
