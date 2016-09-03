precision mediump float;

uniform sampler2D u_depthBuffer;
uniform sampler2D u_origamiBuffer;

varying vec2 uv;

void main() {
  vec4 origamiColor = texture2D(u_origamiBuffer, uv);
  vec4 depthColor = texture2D(u_depthBuffer, uv);
  gl_FragColor = origamiColor + depthColor;
}
