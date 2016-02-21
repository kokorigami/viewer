precision mediump float;
uniform sampler2D u_sampler;
uniform vec2 u_screen;
uniform float u_near;
uniform float u_far;

const int checkPixels = 16;

float decode_depth (float depth) {
  return depth * (u_far - u_near) + u_near;
}

vec3 decode_normal (vec3 encoded) {
  return encoded * 2.0 - 1.0;
}

vec2 to_position (vec3 pixel) {
  return vec2(pixel.x / u_screen.x, pixel.y / u_screen.y);
}

vec4 read_depthnormal (vec3 pixel) {
  vec2 position = to_position(pixel);
  vec4 encoded = texture2D(u_sampler, position);
  vec4 decoded = vec4(decode_normal(encoded.stp), decode_depth(encoded.q));
  return decoded;
}

float read_horizonX (vec4 decoded0, vec3 pixel) {
  vec2 position = to_position(pixel);
  vec4 decoded1 = read_depthnormal(pixel);
  float hyp = distance(vec2(0.0, decoded0.a), vec2(position.x, decoded1.a));
  return (decoded1.a - decoded0.a) / hyp;
}

float read_horizonY (vec4 decoded0, vec3 pixel) {
  vec2 position = to_position(pixel);
  vec4 decoded1 = read_depthnormal(pixel);
  float hyp = distance(vec2(0.0, decoded0.a), vec2(position.y, decoded1.a));
  return (decoded1.a - decoded0.a) / hyp;
}

float read_AO (vec3 pixel) {
  vec4 decoded = read_depthnormal(pixel);
  float sinT_X = decoded.x;
  float sinT_Y = decoded.y;
  float sinH_upX = 0.0;
  float sinH_dnX = 0.0;
  float sinH_upY = 0.0;
  float sinH_dnY = 0.0;

  // Up x
  sinH_upX = max(sinH_upX, read_horizonX(decoded, vec3(pixel.x + 1.0, pixel.yz)));
  sinH_upX = max(sinH_upX, read_horizonX(decoded, vec3(pixel.x + 2.0, pixel.yz)));
  sinH_upX = max(sinH_upX, read_horizonX(decoded, vec3(pixel.x + 3.0, pixel.yz)));
  sinH_upX = max(sinH_upX, read_horizonX(decoded, vec3(pixel.x + 4.0, pixel.yz)));
  sinH_upX = max(sinH_upX, read_horizonX(decoded, vec3(pixel.x + 5.0, pixel.yz)));
  sinH_upX = max(sinH_upX, read_horizonX(decoded, vec3(pixel.x + 6.0, pixel.yz)));
  sinH_upX = max(sinH_upX, read_horizonX(decoded, vec3(pixel.x + 7.0, pixel.yz)));
  sinH_upX = max(sinH_upX, read_horizonX(decoded, vec3(pixel.x + 8.0, pixel.yz)));

  // Down x
  sinH_dnX = max(sinH_dnX, read_horizonX(decoded, vec3(pixel.x - 1.0, pixel.yz)));
  sinH_dnX = max(sinH_dnX, read_horizonX(decoded, vec3(pixel.x - 2.0, pixel.yz)));
  sinH_dnX = max(sinH_dnX, read_horizonX(decoded, vec3(pixel.x - 3.0, pixel.yz)));
  sinH_dnX = max(sinH_dnX, read_horizonX(decoded, vec3(pixel.x - 4.0, pixel.yz)));
  sinH_dnX = max(sinH_dnX, read_horizonX(decoded, vec3(pixel.x - 5.0, pixel.yz)));
  sinH_dnX = max(sinH_dnX, read_horizonX(decoded, vec3(pixel.x - 6.0, pixel.yz)));
  sinH_dnX = max(sinH_dnX, read_horizonX(decoded, vec3(pixel.x - 7.0, pixel.yz)));
  sinH_dnX = max(sinH_dnX, read_horizonX(decoded, vec3(pixel.x - 8.0, pixel.yz)));

  // Up x
  sinH_upY = max(sinH_upY, read_horizonY(decoded, vec3(pixel.x, pixel.y + 1.0, pixel.z)));
  sinH_upY = max(sinH_upY, read_horizonY(decoded, vec3(pixel.x, pixel.y + 2.0, pixel.z)));
  sinH_upY = max(sinH_upY, read_horizonY(decoded, vec3(pixel.x, pixel.y + 3.0, pixel.z)));
  sinH_upY = max(sinH_upY, read_horizonY(decoded, vec3(pixel.x, pixel.y + 4.0, pixel.z)));
  sinH_upY = max(sinH_upY, read_horizonY(decoded, vec3(pixel.x, pixel.y + 5.0, pixel.z)));
  sinH_upY = max(sinH_upY, read_horizonY(decoded, vec3(pixel.x, pixel.y + 6.0, pixel.z)));
  sinH_upY = max(sinH_upY, read_horizonY(decoded, vec3(pixel.x, pixel.y + 7.0, pixel.z)));
  sinH_upY = max(sinH_upY, read_horizonY(decoded, vec3(pixel.x, pixel.y + 8.0, pixel.z)));

  // Down x
  sinH_dnY = max(sinH_dnY, read_horizonY(decoded, vec3(pixel.x, pixel.y - 1.0, pixel.z)));
  sinH_dnY = max(sinH_dnY, read_horizonY(decoded, vec3(pixel.x, pixel.y - 2.0, pixel.z)));
  sinH_dnY = max(sinH_dnY, read_horizonY(decoded, vec3(pixel.x, pixel.y - 3.0, pixel.z)));
  sinH_dnY = max(sinH_dnY, read_horizonY(decoded, vec3(pixel.x, pixel.y - 4.0, pixel.z)));
  sinH_dnY = max(sinH_dnY, read_horizonY(decoded, vec3(pixel.x, pixel.y - 5.0, pixel.z)));
  sinH_dnY = max(sinH_dnY, read_horizonY(decoded, vec3(pixel.x, pixel.y - 6.0, pixel.z)));
  sinH_dnY = max(sinH_dnY, read_horizonY(decoded, vec3(pixel.x, pixel.y - 7.0, pixel.z)));
  sinH_dnY = max(sinH_dnY, read_horizonY(decoded, vec3(pixel.x, pixel.y - 8.0, pixel.z)));

  float AO = (sinH_upX - sinT_X)*(sinH_dnX + sinT_X)*(sinH_upY - sinT_Y)*(sinH_dnY + sinT_Y);
  return AO;
}

void main() {
  vec3 pixel = gl_FragCoord.xyz;
  float AO = read_AO(pixel);
  gl_FragColor = vec4(AO, AO, AO, 1.0);
}
