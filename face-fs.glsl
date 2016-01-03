precision mediump float;
varying vec4 v_color;
varying vec3 v_normal;
// varying vec2 v_texCoord;
varying vec3 v_surfaceToLight;
varying vec3 v_surfaceToView;

uniform vec4 u_lightColor;
uniform vec4 u_ambient;
uniform vec4 u_specular;
uniform float u_shininess;
uniform float u_specularFactor;
//uniform sampler2D u_texture;

vec4 lit(float l ,float h, float m) {
  return vec4(1.0,
              max(l, 0.0),
              (l > 0.0) ? pow(max(0.0, h), m) : 0.0,
              1.0);
}

void main() {
  // vec4 texColor = texture2D(u_texture, v_texCoord);
  vec4 texColor = v_color;
  vec3 normal = normalize(v_normal);
  vec3 surfaceToLight = normalize(v_surfaceToLight);
  vec3 surfaceToView = normalize(v_surfaceToView);
  vec3 halfVector = normalize(surfaceToLight + surfaceToView);

  vec4 litR = lit(dot(normal, surfaceToLight),
                dot(normal, halfVector), u_shininess);

  vec4 outColor = vec4(
      (u_lightColor * (texColor * litR.y + texColor * u_ambient)).rgb,
      texColor.a);

  gl_FragColor = outColor;
}
