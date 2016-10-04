precision mediump float;

varying vec2 v_position;
varying vec2 v_uv;

uniform sampler2D u_sampler;
uniform vec2 u_screen;
uniform float u_near;
uniform float u_far;

const int numChecks = 4;

float decode_depth (float depth) {
  float adjusted = depth - 1.0;
  return adjusted * (u_far - u_near) + u_near;
}

vec3 decode_normal (vec3 encoded) {
  return encoded * 2.0 - 1.0;
}

vec4 read_depthnormal(vec2 texcoord) {
  vec4 encoded = texture2D(u_sampler, texcoord);
  vec4 decoded = vec4(decode_normal(encoded.stp), decode_depth(encoded.q));
  return decoded;
}

float read_horizon(vec2 texcoord, vec2 direction) {
  float delta = 0.0;
  vec4 depthnormal = read_depthnormal(texcoord);
  vec3 normal = depthnormal.stp;
  float depth = depthnormal.q;
  vec2 onePixel = vec2(1, 1) / u_screen;

  for (int i = 1; i < numChecks + 1; i++) {
    vec2 offset = float(i) * direction;
    vec2 offsetcoord = texcoord + onePixel * offset;
    float offsetdepth = read_depthnormal(offsetcoord).q;
    vec3 origToPoint = vec3(offsetcoord, offsetdepth) - vec3(0, 0, depth);
    delta = max(delta, dot(origToPoint, normal) * 100.0);
  }
  return delta;
}

float read_AO (vec2 texcoord) {
  float right = 0.0;
  float left = 0.0;
  float down = 0.0;
  float up = 0.0;

  right = read_horizon(texcoord, vec2(1, 0));
  left = read_horizon(texcoord, vec2(-1, 0));
  down = read_horizon(texcoord, vec2(0, 1));
  up = read_horizon(texcoord, vec2(0, -1));

  float AO = right + left + down + up;
  return AO;
}

void main() {
  float AO = read_AO(v_uv);
  gl_FragColor = vec4(0.0, 0.0, 0.0, AO);
}
