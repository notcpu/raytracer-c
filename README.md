# Ray Tracer in cpp

A small CPU ray tracer written in C++. It loads a hardcoded scene (triangles + materials), shoots rays through a virtual camera, and spits out a `.ppm` image. Supports reflection, transparency, and basic diffuse shading with shadow rays and parallelizes the render across rows with OpenMP if it's available.

## How it works

The scene geometry lives in `detail_scene.h` (vertex/face data, materials, light position, FOV, etc.) and `vector.h` handles the vector math. `main.cpp` does the actual tracing:

- **Intersection** - Möller–Trumbore ray/triangle intersection, checked against every face in the scene per ray.
- **Shading** - each material can mix ambient, reflective, transparent, and diffuse components. Reflective and transparent surfaces recurse (bounce a new ray) up to `max_bounce`. Diffuse shading casts a shadow ray toward the light and only lights the point if nothing's in the way.
- **Output** - every pixel gets traced and written into an in-memory framebuffer (thread-safe, since each row is independent), then dumped to disk as a binary PPM (`P6`) once the render finishes.
- **File naming** - output files are `render_0000.ppm`, `render_0001.ppm`, etc. The counter is tracked in `renders.cnt`, so every run gets its own number instead of overwriting the last one.

While it's running you'll get a live progress bar with ETA, elapsed time, and estimated file size, updated as rows complete.

## Usage

```bash
./raytracer [width] [height]
```

- No args → uses the scene's default resolution.
- One arg → renders a square image at that resolution (e.g. `./raytracer 800` → 800×800).
- Two args → width and height set separately (e.g. `./raytracer 1920 1080`).

Output is written to `render_NNNN.ppm` in the current directory. PPM isn't the most universally friendly format, so you'll want something like GIMP, ImageMagick, or `ffmpeg` to convert it to PNG

## Building

Two build scripts depending on which compiler you've got:

```bash
./compile.sh       # gcc
./clangcomp.sh      # clang
```

Both should produce a single binary in the project root. OpenMP support is picked up automatically at compile time if your toolchain has it. if it's missing, the tracer still works fine, just single threaded.

## Cleanup scripts

- `./remove.sh` — deletes all rendered `.ppm` files and resets `renders.cnt`, but leaves the binary alone.
- `./rmall.sh` — same as above, but also nukes the compiled binary. Use this for a full clean slate.

## Notes

- Scene data (geometry, materials, light, FOV, background color, bounce limit) is all compiled in via `detail_scene.h` — there's no scene file format or loader, so changing the scene means editing that header and recompiling.
- Render time scales with face count and resolution, since intersection testing is brute-force (no BVH or spatial acceleration structure). Its fine for simple scenes, but it will start to drag once you throw a lot of geometry at it.
- Ill try to be updating this, i dont know. It seems good enough but i might soon add other features.

## Todo

- Scene editor / friendly format for scenes for easy editing
- Realtime support (?)
- Internal PNG conversion
- More general args like exposure, angle, refraction index, etc
