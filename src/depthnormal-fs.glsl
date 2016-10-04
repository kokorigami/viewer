precision mediump float;
varying vec4 v_position;
varying vec4 v_normal;
varying vec3 v_normalWorld;
varying vec3 v_surfaceToView;

uniform mat4 u_worldViewProjection;
uniform float u_near;
uniform float u_far;

float encode_depth (vec3 position) {
  float depth = (-position.z + u_far)/(u_far - u_near);
  return depth;
}

vec3 encode_normal (vec3 normal) {
  vec3 encoded = (normal + 1.0) / 2.0;
  return encoded;
}

void main() {
  vec3 position = (u_worldViewProjection * v_position).xyz;
  vec3 normal = normalize((u_worldViewProjection * v_normal).xyz);
  vec3 surfaceToView = normalize(v_surfaceToView);

  float faceCheck = sign(dot(normal, surfaceToView));
  vec3 reflection = vec3(faceCheck, faceCheck, faceCheck);

  gl_FragColor = vec4(encode_normal(normal * reflection), encode_depth(position));
}
