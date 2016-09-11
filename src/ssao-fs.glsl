precision mediump float;

varying vec2 v_position;
varying vec2 v_uv;

uniform sampler2D u_sampler;
uniform vec2 u_screen;
uniform float u_near;
uniform float u_far;

const int numChecks = 8;

float decode_depth (float depth) {
  return depth * (u_far - u_near) + u_near;
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
  float sinH = 0.0;
  float depth = read_depthnormal(texcoord).q;
  vec2 onePixel = vec2(1, 1) / u_screen;

  for (int i = 1; i < numChecks + 1; i++) {
    vec2 offset = float(i) * direction;
    vec2 offsetcoord = texcoord + onePixel * offset;
    float offsetdepth = read_depthnormal(offsetcoord).q;
    // float hyp = distance(vec3(0, 0, depth), vec3(offsetcoord, offsetdepth));
    // sinH = max(sinH, offsetdepth / hyp);
    sinH = max(sinH, (offsetdepth - depth));
  }
  return sinH;
}

float read_AO (vec2 texcoord) {
  vec4 decoded = read_depthnormal(texcoord);
  float sinT_X = decoded.x;
  float sinT_Y = decoded.y;
  float sinH_upX = 0.0;
  float sinH_dnX = 0.0;
  float sinH_upY = 0.0;
  float sinH_dnY = 0.0;

  // Up x
  sinH_upX = read_horizon(texcoord, vec2(1, 0));
  // Down x
  sinH_dnX = read_horizon(texcoord, vec2(-1, 0));
  // Up y
  sinH_upY = read_horizon(texcoord, vec2(0, 1));
  // Down y
  sinH_dnY = read_horizon(texcoord, vec2(0, -1));

  float AO = sinH_upX + sinH_dnX + sinH_upY + sinH_dnY;
  //float AO = (sinH_upX - sinT_X) + (sinH_dnX + sinT_X) + (sinH_upY - sinT_Y) + (sinH_dnY + sinT_Y);
  return AO;
}

void main() {
  //float AO = read_AO(v_uv);
  //gl_FragColor = vec4(AO, AO, AO, 1.0);
  gl_FragColor = vec4(1.0, 0.0, 0.0, 1.0);
}
