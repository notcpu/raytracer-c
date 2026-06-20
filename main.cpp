#include <cstdio>
#include <cstdint>
#include <cmath>
#include <ctime>
#include <vector>
#include <atomic>
#include <cstddef>
#include <cstdlib>
#include <cstring>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "detail_scene.h"
#include "vector.h"

// vertex / face accessors
void get_point(struct v_f &res, unsigned int vertex) {
  const float *offset = v_mem + (vertex * 3);
  v_f_create_f_f_f(res, offset[0], offset[1], offset[2]);
}

void get_object(const unsigned int i,
                struct v_f &v1, struct v_f &v2, struct v_f &v3) {
  const uint16_t *offset = f_mem + (i * 4);
  get_point(v1, offset[0]);
  get_point(v2, offset[1]);
  get_point(v3, offset[2]);
}

// helpers

static double seconds_since(const struct timespec &start, const struct timespec &now) {
  return static_cast<double>(now.tv_sec - start.tv_sec) +
         static_cast<double>(now.tv_nsec - start.tv_nsec) / 1e9;
}

static void format_bar(char *out, size_t out_size, double progress, int bar_width) {
  if (bar_width < 4) bar_width = 4;
  if (progress < 0.0) progress = 0.0;
  if (progress > 1.0) progress = 1.0;

  int filled = static_cast<int>(progress * bar_width + 0.5);
  if (filled < 0) filled = 0;
  if (filled > bar_width) filled = bar_width;

  int pos = 0;
  if (out_size == 0) return;
  out[pos++] = '[';
  for (int i = 0; i < bar_width && pos < static_cast<int>(out_size) - 2; ++i) {
    out[pos++] = (i < filled) ? '#' : '-';
  }
  if (pos < static_cast<int>(out_size) - 1) out[pos++] = ']';
  out[pos] = '\0';
}

// moller trumbore

bool ray_tri_intersect(const struct v_f ray_direction,
                       const struct v_f ray_origin,
                       const unsigned int object,
                       double &tri_dist, double &tri_u, double &tri_v,
                       struct v_f &normal, struct v_f &point) {
  #define EPSILON 0.000001

  struct v_f v1, v2, v3;
  struct v_f e1, e2;
  struct v_f P, Q, T;
  double det, inv_det;
  double t, t_u, t_v;

  get_object(object, v1, v2, v3);

  v_f_sub(e1, v2, v1);
  v_f_sub(e2, v3, v1);

  v_f_cross(P, ray_direction, e2);
  det = v_f_dot(e1, P);

  if (fabs(det) < EPSILON) return false;
  inv_det = 1.0 / det;

  v_f_sub(T, ray_origin, v1);

  t_u = v_f_dot(T, P) * inv_det;
  if (t_u < 0.0 || t_u > 1.0) return false;

  v_f_cross(Q, T, e1);
  t_v = v_f_dot(ray_direction, Q) * inv_det;
  if (t_v < 0.0 || t_u + t_v > 1.0) return false;

  t = v_f_dot(e2, Q) * inv_det;
  if (t > tri_dist || t < EPSILON) return false;

  v_f_cross(normal, e1, e2);
  v_f_norm(normal);

  tri_u    = t_u;
  tri_v    = t_v;
  tri_dist = t;

  point.x = v1.x + tri_u * (v2.x - v1.x) + tri_v * (v3.x - v1.x);
  point.y = v1.y + tri_u * (v2.y - v1.y) + tri_v * (v3.y - v1.y);
  point.z = v1.z + tri_u * (v2.z - v1.z) + tri_v * (v3.z - v1.z);

  return true;
}

// recursive tracer

void trace(const struct v_f ray_direction,
           const struct v_f ray_origin,
           const unsigned int bounce,
           unsigned int color[3]) {

  if (bounce > max_bounce) {
    color[0] = bg[0]; color[1] = bg[1]; color[2] = bg[2];
    return;
  }

  unsigned int object = f_count;
  double dist = max_dist, t_u, t_v;
  struct v_f normal, point;

  for (unsigned int i = 0; i < f_count; i++) {
    if (ray_tri_intersect(ray_direction, ray_origin, i,
                          dist, t_u, t_v, normal, point)) {
      object = i;
    }
  }

  if (object < f_count) {
    color[0] = color[1] = color[2] = 0;

    int diffuse_value = 255;
    unsigned int idx = f_mem[object * 4 + 3];

    // Ambient
    if (mtl[idx][0]) {
      diffuse_value -= mtl[idx][0];
      double ambient_pct = mtl[idx][0] / 255.0;
      color[0] = static_cast<unsigned int>(mtl[idx][3] * ambient_pct);
      color[1] = static_cast<unsigned int>(mtl[idx][4] * ambient_pct);
      color[2] = static_cast<unsigned int>(mtl[idx][5] * ambient_pct);
    }

    // Reflective
    if (mtl[idx][1]) {
      diffuse_value -= mtl[idx][1];
      double reflect_pct = mtl[idx][1] / 255.0;
      unsigned int reflect_clr[3];

      double c = -2.0 * v_f_dot(normal, ray_direction);
      struct v_f reflect_direction;
      v_f_scale(normal, c);
      v_f_add(reflect_direction, ray_direction, normal);

      trace(reflect_direction, point, bounce + 1, reflect_clr);

      color[0] += static_cast<unsigned int>(reflect_clr[0] * reflect_pct);
      color[1] += static_cast<unsigned int>(reflect_clr[1] * reflect_pct);
      color[2] += static_cast<unsigned int>(reflect_clr[2] * reflect_pct);
    }

    // Transparency
    if (mtl[idx][2]) {
      diffuse_value -= mtl[idx][2];
      double transparent_pct = mtl[idx][2] / 255.0;
      unsigned int transparent_clr[3];

      trace(ray_direction, point, bounce + 1, transparent_clr);

      color[0] += static_cast<unsigned int>(transparent_clr[0] * transparent_pct);
      color[1] += static_cast<unsigned int>(transparent_clr[1] * transparent_pct);
      color[2] += static_cast<unsigned int>(transparent_clr[2] * transparent_pct);
    }

    // Diffuse
    if (diffuse_value) {
      double diffuse_pct = static_cast<double>(diffuse_value) / 255.0;

      struct v_f light_direction;
      get_point(light_direction, light_v);
      v_f_sub(light_direction, light_direction, point);
      v_f_norm(light_direction);

      bool lit = true;
      for (unsigned int j = 0; j < f_count; j++) {
        if (j == object) continue;
        double dist2 = max_dist;
        struct v_f junk, junk2;
        if (ray_tri_intersect(light_direction, point, j,
                               dist2, t_u, t_v, junk, junk2)) {
          lit = false;
          break;
        }
      }

      if (lit) {
        double shade = diffuse_pct * fabs(v_f_dot(light_direction, normal));
        color[0] += static_cast<unsigned int>(mtl[idx][3] * shade);
        color[1] += static_cast<unsigned int>(mtl[idx][4] * shade);
        color[2] += static_cast<unsigned int>(mtl[idx][5] * shade);
      }
    }

  } else {
    color[0] = bg[0]; color[1] = bg[1]; color[2] = bg[2];
  }

  if (color[0] > 255) color[0] = 255;
  if (color[1] > 255) color[1] = 255;
  if (color[2] > 255) color[2] = 255;
}

// main

int main(int argc, char** argv) {
  // default to scene dimensions but allow runtime overrides
  unsigned int width = defaultWidth;
  unsigned int height = defaultHeight;

  if (argc >= 2) {
    unsigned long parsed = strtoul(argv[1], nullptr, 10);
    if (parsed > 0) {
      width = static_cast<unsigned int>(parsed);
      height = width;
    }
  }

  if (argc >= 3) {
    unsigned long parsed = strtoul(argv[2], nullptr, 10);
    if (parsed > 0) {
      height = static_cast<unsigned int>(parsed);
    }
  }

  // output filename: render_NNNN.ppm where counter is stored in renders.cnt
  uint16_t renderNum = 0;
  FILE *cf = fopen("renders.cnt", "rb");
  if (cf) {
    fread(&renderNum, sizeof(renderNum), 1, cf);
    fclose(cf);
  }

  cf = fopen("renders.cnt", "wb");
  if (cf) {
    uint16_t next = renderNum + 1;
    fwrite(&next, sizeof(next), 1, cf);
    fclose(cf);
  }

  char filename[20];
  snprintf(filename, sizeof(filename), "render_%04d.ppm", renderNum);
  printf("Output file: %s\n", filename);

  FILE *imgFile = fopen(filename, "wb");
  if (!imgFile) {
    printf("Failed to open output file!\n");
    return 1;
  }

  // PPM header
  char header[64];
  int header_len = snprintf(header, sizeof(header), "P6\n%u %u\n255\n", width, height);
  if (header_len < 0) header_len = 0;
  fwrite(header, 1, static_cast<size_t>(header_len), imgFile);

  const size_t pixel_bytes = static_cast<size_t>(width) * static_cast<size_t>(height) * 3;
  const size_t total_bytes = static_cast<size_t>(header_len) + pixel_bytes;
  const double total_mib = static_cast<double>(total_bytes) / (1024.0 * 1024.0);

  // Ray setup
  double lerp_rate_x = 2.0 * tan(half_fov) / width;
  double lerp_rate_y = 2.0 * tan(half_fov) / height;
  double half_width = (width - 1) / 2.0;
  double half_height = (height - 1) / 2.0;

  struct v_f ray_origin;
  v_f_create_f_f_f(ray_origin, 0, 0, 0);

  // framebuffer for thread safe parallel render
  const size_t pixel_count = static_cast<size_t>(width) * static_cast<size_t>(height);
  std::vector<uint8_t> framebuffer(pixel_count * 3, 0);

  struct timespec t_start, t_now;
  clock_gettime(CLOCK_MONOTONIC, &t_start);
  std::atomic<unsigned int> completed_rows{0};

  printf("Rendering %ux%u...\n", width, height);
  printf("Estimated output size: %.2f MiB\n", total_mib);

#ifdef _OPENMP
  printf("OpenMP threads available: %d\n", omp_get_max_threads());
#endif

  #pragma omp parallel for schedule(dynamic)
  for (int y = 0; y < static_cast<int>(height); y++) {
    double v = half_height - y;
    for (unsigned int x = 0; x < width; x++) {
      double u = x - half_width;

      struct v_f ray;
      v_f_create_f_f_f(ray, u * lerp_rate_x, v * lerp_rate_y, -1.0);
      v_f_norm(ray);

      unsigned int color[3] = {0, 0, 0};
      trace(ray, ray_origin, 0, color);

      size_t idx = (static_cast<size_t>(y) * width + x) * 3;
      framebuffer[idx + 0] = static_cast<uint8_t>(color[0]);
      framebuffer[idx + 1] = static_cast<uint8_t>(color[1]);
      framebuffer[idx + 2] = static_cast<uint8_t>(color[2]);
    }

    unsigned int done = completed_rows.fetch_add(1, std::memory_order_relaxed) + 1;

    // single line progress display
    #pragma omp critical(progress)
    {
      clock_gettime(CLOCK_MONOTONIC, &t_now);
      double elapsed = seconds_since(t_start, t_now);
      double progress = static_cast<double>(done) / static_cast<double>(height);
      double eta = (done > 0)
        ? elapsed * (static_cast<double>(height - done) / static_cast<double>(done))
        : 0.0;

      char bar[48];
      format_bar(bar, sizeof(bar), progress, 24);

      printf("\r%6.2f%% %s %4u/%4u rows | ETA %8.2fs | elapsed %8.2fs | file %.2f MiB",
             progress * 100.0,
             bar,
             done,
             height,
             eta,
             elapsed,
             total_mib);
      fflush(stdout);
    }
  }

  printf("\n");

  // write finished image
  fwrite(framebuffer.data(), 1, framebuffer.size(), imgFile);
  fclose(imgFile);

  clock_gettime(CLOCK_MONOTONIC, &t_now);
  double total = seconds_since(t_start, t_now);
  printf("Done! render_%04d.ppm saved in %.2fs\n", renderNum, total);

  return 0;
}
