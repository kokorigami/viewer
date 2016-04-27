precision mediump float;

uniform vec3 u_lightWorldPos;
uniform mat4 u_camera;
uniform mat4 u_worldViewProjection;
uniform mat4 u_world;
uniform mat4 u_worldRotation;

varying float v_foldType;
varying float v_lengthSoFar;
varying vec4 v_position;

const float dotSize = 0.1;
const float numDashes = 8.0;

void main() {
  float upperLimit = v_foldType + dotSize / 2.0;
  float lowerLimit = v_foldType - dotSize / 2.0;

  float section = 2.0 * fract(v_lengthSoFar * numDashes);
  float alpha = floor(section);
  if (alpha == 0.0 && section > lowerLimit && section < upperLimit)
    alpha = 1.0;

  if (alpha == 0.0)
    discard;
  else
    gl_FragColor = vec4(0, 0, 0, alpha);
}
