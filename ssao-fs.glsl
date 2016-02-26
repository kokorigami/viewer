precision mediump float;

varying vec3 v_position;
varying vec2 v_texcoord;

uniform sampler2D u_sampler;
uniform vec2 u_screen;
uniform float u_near;
uniform float u_far;

float decode_depth (float depth) {
  return depth * (u_far - u_near) + u_near;
}

vec3 decode_normal (vec3 encoded) {
  return encoded * 2.0 - 1.0;
}

vec4 read_depthnormal (vec3 position) {
  vec2 texcoord = 0.5*(position.xy + 1.0);
  vec4 encoded = texture2D(u_sampler, texcoord);
  vec4 decoded = vec4(decode_normal(encoded.stp), decode_depth(encoded.q));
  return decoded;
}

float read_horizonX (vec4 decoded0, vec3 position) {
  vec4 decoded1 = read_depthnormal(position);
  float hyp = distance(vec2(0.0, decoded0.a), vec2(position.x, decoded1.a));
  return (decoded1.a - decoded0.a) / hyp;
}

float read_horizonY (vec4 decoded0, vec3 position) {
  vec4 decoded1 = read_depthnormal(position);
  float hyp = distance(vec2(0.0, decoded0.a), vec2(position.y, decoded1.a));
  return (decoded1.a - decoded0.a) / hyp;
}

float read_AO (vec3 position) {
  vec4 decoded = read_depthnormal(position);
  float sinT_X = decoded.x;
  float sinT_Y = decoded.y;
  float sinH_upX = 0.0;
  float sinH_dnX = 0.0;
  float sinH_upY = 0.0;
  float sinH_dnY = 0.0;

  const int checkTexture = 8;
  vec3 checkPosition;

  // Up x
  for (int i = 1; i < checkTexture + 1; i++) {
    checkPosition = vec3(position.x + float(i)/u_screen.x, position.yz);
    sinH_upX = max(sinH_upX, read_horizonX(decoded, checkPosition));
  }

  // Down x
  for (int i = 1; i < checkTexture + 1; i++) {
    checkPosition = vec3(position.x - float(i)/u_screen.x, position.yz);
    sinH_dnX = max(sinH_dnX, read_horizonX(decoded, checkPosition));
  }

  // Up y
  for (int i = 1; i < checkTexture + 1; i++) {
    checkPosition = vec3(position.x, position.y + float(i)/u_screen.y, position.z);
    sinH_upY = max(sinH_upY, read_horizonX(decoded, checkPosition));
  }

  // Down y
  for (int i = 1; i < checkTexture + 1; i++) {
    checkPosition = vec3(position.x, position.y - float(i)/u_screen.y, position.z);
    sinH_dnY = max(sinH_dnY, read_horizonX(decoded, checkPosition));
  }

  float AO = (sinH_upX - sinT_X)*(sinH_dnX + sinT_X)*(sinH_upY - sinT_Y)*(sinH_dnY + sinT_Y);
  return AO;
}

void main() {
  // float AO = read_AO(v_position);
  // gl_FragColor = vec4(AO, AO, AO, 1.0);
  // gl_FragColor = vec4(1, 0, 0, 1);

  gl_FragColor = texture2D(u_sampler, v_texcoord);
}
