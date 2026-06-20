#include "vector.h"

#include <math.h>

// Create a vector from three floats
void v_f_create_f_f_f(struct v_f &res, const double x, const double y, const double z)
{
  res.x = x;
  res.y = y;
  res.z = z;
}

/*
// Create a vector from an array of three doubles
void v_f_create_f3(struct v_f &res, const double v[3])
{
  res.x = v[0];
  res.y = v[1];
  res.z = v[2];
}

// Create a vector as a copy of another vector.
void v_f_create_v(struct v_f &res, const struct v_f src)
{
  res.x = src.x;
  res.y = src.y;
  res.z = src.z;
}
*/

// Subtract one vector from another
void v_f_sub(struct v_f &res, const struct v_f v1, const struct v_f v2)
{
  res.x = v1.x - v2.x;
  res.y = v1.y - v2.y;
  res.z = v1.z - v2.z;
}

// Add two vectors
void v_f_add(struct v_f &res, const struct v_f v1, const struct v_f v2)
{
  res.x = v1.x + v2.x;
  res.y = v1.y + v2.y;
  res.z = v1.z + v2.z;
}

// Scale a vector by a constant amount
void v_f_scale(struct v_f &res, const double c)
{
  res.x *= c;
  res.y *= c;
  res.z *= c;
}

// Compute cross-product of two vectors
void v_f_cross(struct v_f &res, const struct v_f v1, const struct v_f v2)
{
  res.x = (v1.y * v2.z) - (v2.y * v1.z);
  res.y = (v1.z * v2.x) - (v2.z * v1.x);
  res.z = (v1.x * v2.y) - (v2.x * v1.y);
}

// Compute dot-product of two vectors
double v_f_dot(const struct v_f v1, const struct v_f v2)
{
  return (v1.x * v2.x) + (v1.y * v2.y) + (v1.z * v2.z);
}

// Normalize a vector (scale to distance=1)
void v_f_norm(struct v_f &v)
{
  double dist = sqrt((v.x * v.x) + (v.y * v.y) + (v.z * v.z));
  v.x = v.x / dist;
  v.y = v.y / dist;
  v.z = v.z / dist;
}

