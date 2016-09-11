precision mediump float;
varying vec3 v_normal;
varying vec3 v_texCoord;
varying vec3 v_surfaceToLight;
varying vec3 v_surfaceToView;

uniform vec4 u_lightColor;
uniform vec4 u_ambient;
uniform float u_shininess;
uniform sampler2D u_texture0;
uniform sampler2D u_texture1;

vec4 lit(float l ,float h, float m) {
  return vec4(1.0,
              max(l, 0.0),
              (l > 0.0) ? pow(max(0.0, h), m) : 0.0,
              1.0);
}

void main() {
  vec3 normal = normalize(v_normal);
  vec3 surfaceToLight = normalize(v_surfaceToLight);
  vec3 surfaceToView = normalize(v_surfaceToView);

  float faceCheck = dot(normal, surfaceToView);
  float isFront = (sign(faceCheck) + 1.0) * 0.5;
  float isBack = (-sign(faceCheck) + 1.0) * 0.5;
  vec3 reflection = vec3(faceCheck, faceCheck, faceCheck);

  vec4 frontTexture = texture2D(u_texture0, v_texCoord.xy) * vec4(isFront, isFront, isFront, 1.0);
  vec4 backTexture = texture2D(u_texture1, v_texCoord.xy) * vec4(isBack, isBack, isBack, 1.0);
  vec4 texColor = frontTexture + backTexture;

  vec3 halfVector = normalize((surfaceToLight) + (surfaceToView * reflection));
  vec4 litR = lit(dot(normal, surfaceToLight * reflection), dot(normal, halfVector), u_shininess);

  vec4 outColor = vec4(
      (u_lightColor * (texColor * litR.y + texColor * u_ambient)).rgb,
      texColor.a);

  gl_FragColor = outColor;
}
