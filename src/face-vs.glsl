precision mediump float;
attribute vec3 position;
attribute vec3 normal;
attribute vec3 texCoord;

// attribute vec3 color;

varying vec4 v_position;
varying vec3 v_normal;
varying vec3 v_texCoord;
varying vec3 v_surfaceToLight;
varying vec3 v_surfaceToView;

uniform vec3 u_lightWorldPos;
uniform mat4 u_camera;
uniform mat4 u_world;
uniform mat4 u_worldViewProjection;
uniform sampler2D u_texture;

void main() {
  v_normal = (u_world * vec4(normal, 0)).xyz;
  v_texCoord = texCoord;
  v_position = vec4(position, 1.0);
  v_surfaceToLight = u_lightWorldPos - (u_world * v_position).xyz;
  v_surfaceToView = (u_camera[3] - (u_world * v_position)).xyz;
  gl_Position = (u_worldViewProjection * v_position);
}
