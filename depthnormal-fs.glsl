precision mediump float;
varying vec4 v_position;
varying vec3 v_normal;

uniform mat4 u_worldViewProjection;
uniform mat4 u_worldRotation;
uniform float u_near;
uniform float u_far;

float encode_depth (vec3 position) {
  float depth = (position.z - u_near)/(u_far - u_near);
  return depth;
}

vec3 encode_normal (vec3 normal) {
  vec3 encoded = (normal + 1.0) / 2.0;
  return encoded;
}

void main() {
  vec3 position = (u_worldViewProjection * v_position).xyz; // move to vshader?
  vec3 normal = v_normal;
  gl_FragColor = vec4(encode_normal(normal), encode_depth(position));
  //float depth = encode_depth(position);
  //gl_FragColor = vec4(depth, depth, depth, 1.0);
}
