precision mediump float;

attribute vec2 position;

varying vec2 v_position;
varying vec2 v_uv;

void main() {
  gl_Position = vec4(position, 0.0, 1.0);
  v_uv = 0.5 * (position + 1.0);
  v_position = position;
}
