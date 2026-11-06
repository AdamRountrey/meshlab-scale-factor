// SPDX-License-Identifier: GPL-2.0-or-later
#pragma once

#include <array>
#include <cstddef>
#include <vector>

namespace meshlab_scale {

struct Vec3 {
    double x = 0.0, y = 0.0, z = 0.0;
};

struct Mesh {
    std::vector<Vec3> vertices;
    std::vector<std::array<std::size_t, 3>> faces;
};

struct Options {
    bool largestComponent = true;
    // Maximum bidirectional point distance / target bounding-box diagonal.
    double relativeTolerance = 1e-5;
    bool allowRotation = false;
};

struct Result {
    // target = scale * rotation * source + translation; rotation is row-major.
    double scale = 0.0;
    double inverseScale = 0.0;
    Vec3 translation;
    double rmsError = 0.0;
    double maxError = 0.0;
    double relativeRmsError = 0.0;
    double relativeMaxError = 0.0;
    std::size_t sourceVerticesUsed = 0, targetVerticesUsed = 0;
    // Includes exact duplicate positions and discarded components/unused vertices.
    std::size_t sourceVerticesIgnored = 0, targetVerticesIgnored = 0;
    std::array<std::array<double, 3>, 3> rotation{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
};

// Throws std::runtime_error for invalid data or a failed geometric verification.
// Both retained point sets must describe the same vertices, up to a positive
// uniform scale, translation, optional proper rotation, and coordinate rounding.
// Rotation proposals are verified against the complete retained point sets.
Result estimateScale(const Mesh& source, const Mesh& target,
                     const Options& options = {});

} // namespace meshlab_scale
