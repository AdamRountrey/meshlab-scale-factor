# MeshLab Scale Factor

Measure the uniform scale between two MeshLab layers containing the same geometry, including copies with missing small floating pieces. Optional rotation matching handles models with different orientations.

The filter reports the scale in both directions, translation, a complete 4×4 transformation matrix, an optional separate rotation matrix, and geometric residuals. It does not modify either mesh.

## Install

The release DLL supports **MeshLab 2025.07 on Windows x64, standard single precision**. MeshLab 2025.07d uses double precision and requires a separately compiled plugin.

There is **no standalone installer**. Install the DLL using MeshLab's plugin loader:

1. Open the [latest release](https://github.com/AdamRountrey/meshlab-scale-factor/releases/latest) and download the binary asset `meshlab-scale-factor-2025.07-win64.zip`.
2. Extract the ZIP to a folder.
3. In MeshLab, open **Help > Plugin Info > Load Plugins** and select the extracted `filter_scale_factor.dll`.
4. Find **Measure Scale Between Two Meshes** under **Filters > Quality Measure and Computations**. Restart MeshLab if the menu has not refreshed.

Use the plugin loader to install the DLL; you do not need to replace MeshLab's own libraries. MeshLab's official [extra-plugin installation instructions](https://github.com/cnr-isti-vclab/meshlab-extra-plugins#readme) describe this loader and the requirement to match the application's precision variant.

To update an installed version, close its report windows, remove **FilterScaleFactor** using **Uninstall** in Plugin Info, restart MeshLab, and load the new DLL. The restart releases the previous DLL before replacement.

## Measure two models

1. Import the unscaled and scaled models as separate layers.
2. Run **Filters > Quality Measure and Computations > Measure Scale Between Two Meshes**.
3. Set **Source (unscaled)** and **Target (scaled)** to the appropriate layers.
4. Leave **Compare largest connected piece only** enabled when one model may be missing floating pieces.
5. Enable **Allow rotation between models** if their orientations differ. It is off by default.
6. Click **Apply**. Read or copy the report. The result also appears in MeshLab's log.

**Multiply SOURCE dimensions by** is the source-to-target factor. A result of `2.75` means the target is 2.75 times the source's size. **Multiply TARGET dimensions by** gives the reciprocal. These are linear factors; areas and volumes change by their square and cube.

The fitted relation is:

```text
target position = scale * rotation * source position + translation
```

Translation is fitted automatically, so a different scaling pivot is supported. With rotation disabled, the rotation matrix is the identity. With it enabled, the report also includes the fitted 3×3 rotation matrix.

The copyable report and log always include the complete **source-to-target 4×4 homogeneous transform**, including when rotation is disabled. It contains scale, rotation, and translation together:

```text
T = [ s*R00  s*R01  s*R02  tx ]
    [ s*R10  s*R11  s*R12  ty ]
    [ s*R20  s*R21  s*R22  tz ]
    [   0      0      0     1 ]

[x_target, y_target, z_target, 1]^T = T * [x_source, y_source, z_source, 1]^T
```

Values are displayed row by row, using **column vectors** for multiplication. **Include layer transforms** is enabled by default: the matrix maps displayed source coordinates to displayed target coordinates, after each layer's existing transform. Disable the option to obtain a matrix mapping the stored source coordinates to the stored target coordinates. The report labels the coordinate space explicitly.

The synthetic [example models](examples/README.txt) have a known factor of **2.75**, reordered vertices, and a missing floating triangle. A second target also includes rotation and translation.

## Supported geometry

By default, the filter compares each mesh's largest connected component, measured by unique vertex count, and ignores smaller components and unused vertices. Exactly coincident vertices are welded internally, including duplicated STL triangle vertices. Vertex order may differ. Point clouds without faces use all unique points.

The retained pieces must contain the **same unique vertices**, up to positive uniform scale, translation, rounding, and optional rotation. Missing whole floating components are supported when the same main piece remains largest in both models. Equally large components are rejected as ambiguous.

Remeshing, decimation, deletion within the main piece, nonuniform scaling, and general mirror matching are outside the intended scope. Prepare matching layers if the largest-piece assumption does not apply, or disable that option when all geometry matches. Validation compares vertices, not triangle connectivity.

## Method and accuracy

With rotation disabled, sorted coordinate values along each axis determine one centered least-squares scale. With rotation enabled, scale comes from the ratio of RMS distances to the centroids. Rotation matching uses centroid-radius correspondences and a bounded search over oriented vertex pairs, with optional [Horn quaternion refinement](https://people.csail.mit.edu/bkph/papers/Absolute_Orientation.pdf).

This is a specialized implementation using classical methods, not a state-of-the-art general registration system. Horn's fit is optimal for fixed correct correspondences; the overall correspondence search has **no global optimality guarantee**. A rotation already within tolerance may be accepted without further refinement.

Both modes verify all retained points in both directions using 3D nearest-vertex distances. Unequal unique-vertex counts, excessive residuals, or an exhausted rotation search produce an error. Symmetric models may have several valid rotations, so the reported matrix need not recover the original rotation operation. Large sets with many nearly equal centroid radii may require prior alignment or a distinctive shared component.

The default **Maximum relative geometric error** is `1e-5` of the retained target's bounding-box diagonal. The report gives RMS and maximum distances in target units and relative to that diagonal. This tolerance is an acceptance threshold, not an uncertainty bound on scale.

Up to 17 significant digits are printed to preserve the computed value, **not to claim 17-digit accuracy**. MeshLab's standard build stores coordinates in single precision, and input files may already be rounded. The original factor before rounding cannot generally be recovered exactly. Small residuals support the uniform-scale model but cannot restore lost coordinate precision.

See [VALIDATION.md](VALIDATION.md) for the tested build, numerical coverage, integration checks, and GUI verification limits.

## Build from source

The repository contains the plugin, numerical core, tests, and examples. Third-party dependencies and generated binaries are not vendored.

The standalone core requires CMake 3.18 or later and a C++17 compiler. From the repository root:

```powershell
cmake -S . -B build/core -DBUILD_MESHLAB_PLUGIN=OFF
cmake --build build/core --config Release
ctest --test-dir build/core -C Release --output-on-failure
```

To build the release-compatible Windows plugin, also install Qt **5.15.2** development files and a compatible x64 MSVC toolchain, then obtain the **MeshLab-2025.07** source and its VCGLib submodule:

```powershell
git clone --branch MeshLab-2025.07 --depth 1 --recurse-submodules https://github.com/cnr-isti-vclab/meshlab.git vendor/meshlab
```

Configure with the Qt installation path for your machine, then build the plugin and checks:

```powershell
cmake -S . -B build/plugin -G "Visual Studio 17 2022" -A x64 -DBUILD_MESHLAB_PLUGIN=ON -DCMAKE_PREFIX_PATH="C:/Qt/5.15.2/msvc2019_64"
cmake --build build/plugin --config Release --target filter_scale_factor plugin_probe scale_estimator_tests
ctest --test-dir build/plugin -C Release --output-on-failure
```

`MESHLAB_SOURCE_DIR` can point to an existing matching checkout. CMake builds MeshLab's common library and this plugin without all the other filters; upstream configuration may download missing dependencies. The default build uses single precision. With the commands above, the plugin is written to `build/plugin/plugins/Release/filter_scale_factor.dll`.

`ctest` runs the numerical core suite. To check the compiled plugin interface, run `plugin_probe` with the absolute DLL path as its sole argument and matching MeshLab/Qt runtime DLLs on `PATH`. This probe does not exercise the GUI report window.

## Scripting

The filter's Python name is `get_scale_between_two_meshes`. The returned map includes `scale_source_to_target`, `scale_target_to_source`, `translation_x`, `translation_y`, `translation_z`, `rms_error`, `max_error`, `relative_rms_error`, `relative_max_error`, vertex counts, and `report`.

All sixteen scalar keys `transform_00` through `transform_33` are always present. In `transform_rc`, `r` is the zero-based row and `c` is the zero-based column of the full 4×4 source-to-target matrix above. It uses column vectors and the coordinate space selected by **Include layer transforms**.

`rotation_enabled` records the rotation option; `rotation_00` through `rotation_22` hold the separate row-major 3×3 rotation matrix, which is the identity when disabled. Saved filter scripts without `AllowRotation` retain the default of no rotation. Headless clients do not open a report window.

## License

This plugin's source code is licensed under the **GNU General Public License, version 2 or (at your option) any later version**, SPDX identifier `GPL-2.0-or-later`. It is provided without warranty. See [LICENSE](LICENSE) for the GNU GPL version 2 text. Third-party MeshLab, VCGLib, Qt, and other dependency sources retain their respective notices and terms.
