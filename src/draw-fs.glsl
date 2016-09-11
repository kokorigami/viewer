precision mediump float;

uniform sampler2D u_sampler;
uniform sampler2D u_origami;

varying vec2 v_uv;

void main() {
  vec4 origamiColor = texture2D(u_origami, v_uv);
  vec4 depthColor = texture2D(u_sampler, v_uv);
  gl_FragColor = origamiColor + depthColor;
}
