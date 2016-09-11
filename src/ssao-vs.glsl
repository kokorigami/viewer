attribute vec2 position;
varying vec2 v_position;
varying vec2 v_texcoord;

void main() {
  v_position = position;
  v_texcoord = 0.5 * (position + 1.0);
  gl_Position = vec4(v_position, 0.0, 1.0);
}
