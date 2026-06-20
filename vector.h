// A vector structure
#ifndef VECTOR_H_
#define VECTOR_H_

struct v_f
{
  double x,y,z;
};

// Some functions for working with vectors.
void v_f_create_f_f_f(struct v_f &res, const double x, const double y, const double z);
//void v_f_create_f3(struct v_f &res, const double v[3]);
//void v_f_create_v(struct v_f &res, const struct v_f src);
void v_f_sub(struct v_f &res, const struct v_f v1, const struct v_f v2);
void v_f_add(struct v_f &res, const struct v_f v1, const struct v_f v2);
void v_f_scale(struct v_f &res, const double c);
void v_f_cross(struct v_f &res, const struct v_f v1, const struct v_f v2);
double v_f_dot(const struct v_f v1, const struct v_f v2);
void v_f_norm(struct v_f &v);

#endif


