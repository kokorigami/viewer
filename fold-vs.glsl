precision mediump float;
attribute vec3 position;
attribute float foldType;
attribute float lengthSoFar;

uniform vec3 u_lightWorldPos;
uniform mat4 u_camera;
uniform mat4 u_worldViewProjection;
uniform mat4 u_world;
uniform mat4 u_worldRotation;

varying float v_foldType;
varying float v_lengthSoFar;
varying vec4 v_position;

void main() {
  v_lengthSoFar = lengthSoFar;
  v_position = vec4(position, 1);
  v_foldType = foldType;
  gl_Position = (u_worldViewProjection * v_position);
}
