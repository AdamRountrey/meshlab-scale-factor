# Validation

Validated on **6 November 2026** for the first public release, **v0.2.0**, and Windows package `meshlab-scale-factor-2025.07-win64.zip`.

## Build and runtime

- Windows x64, Release configuration, MSVC 19.44.
- MeshLab 2025.07, standard single precision; Qt 5.15.2.
- The compiled `filter_scale_factor.dll` loaded through Qt's `QPluginLoader` using an installed MeshLab 2025.07 runtime and its dependencies.

This establishes compatibility with the tested runtime. Other MeshLab versions, precision variants, operating systems, and toolchains have not been validated.

## Numerical tests

All **23 cases** in the numerical core test executable passed. Coverage includes:

- Scale, translation, reciprocal results, reordered vertices, missing floating components, and exact-position welding of STL-style duplicated vertices.
- Float-rounded coordinates, extreme scale ranges, large absolute offsets, and invalid or incompatible inputs.
- Arbitrary and 180-degree rotations; isotropic symmetric solids; planar, collinear, and thin geometry; nearly equal centroid radii.
- Rejection of asymmetric 3D reflections, changed geometry, and nonuniform scaling.
- A 20,000-point rotated and reordered cloud.

The tests use synthetic data with known transformations. They do not establish performance or reliability for every real-world mesh or every ambiguous rotation search.

## Plugin integration

The integration probe passed checks through the actual MeshLab filter interface for:

- A known factor of 2.75 with a missing floating triangle, including reciprocal and residual reporting.
- Invalid inputs, layer transforms, raw-coordinate measurement, and translation.
- Rotation disabled by default, rejection of a rotated model when disabled, and successful measurement when enabled.
- Returned proper rotation matrices and rotated layer transforms.
- Full source-to-target 4x4 transforms in all tested rotation and coordinate modes. Each matrix maps actual source vertex coordinates to corresponding reordered target coordinates, including fitted translation and scale; its last row is `[0 0 0 1]`.
- All 16 `transform_00` through `transform_33` scripting outputs and the four copyable report rows agree. The mapping checks use displayed coordinates when layer transforms are included and stored coordinates when they are excluded.
- Preservation of vertex coordinates, topology, flags, and layer transforms.

## Verification limits

The MeshLab GUI checkbox, plugin-installation flow, and report-window interaction have **not been visually verified**. The runtime integration probe exercises the plugin interface without opening those UI elements.

The algorithm's rotation search is bounded and has no global optimality guarantee. Passing these checks does not remove the geometry assumptions or numerical-precision limits described in [README.md](README.md).
