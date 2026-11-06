// SPDX-License-Identifier: GPL-2.0-or-later
#include "scale_estimator.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>

using namespace meshlab_scale;

namespace {
int passed = 0;
void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}
void close(double actual, double expected, double relative = 1e-12, double absolute = 1e-12) {
    require(std::abs(actual - expected) <= absolute + relative * std::abs(expected),
            "Unexpected numerical result: " + std::to_string(actual) + " vs " + std::to_string(expected));
}
void rejects(const std::function<void()>& call, const std::string& fragment = "") {
    try { call(); }
    catch (const std::runtime_error& e) {
        require(fragment.empty() || std::string(e.what()).find(fragment) != std::string::npos,
                std::string("Unexpected error: ") + e.what());
        return;
    }
    throw std::runtime_error("Expected rejection did not occur");
}
Mesh specimen() {
    Mesh m;
    m.vertices = {{-1.3, 2.1, 0.7}, {4.2, -0.6, 1.2}, {0.4, 0.3, -2.5},
                  {1.8, 4.3, 2.6}, {-2.7, -1.9, 3.4}, {0.01, 1.11, 2.22},
                  {2.8, 0.02, -0.19}, {-0.45, -2.8, 1.9}};
    for (std::size_t i = 1; i + 1 < m.vertices.size(); ++i) m.faces.push_back({0, i, i + 1});
    return m;
}
Mesh transform(Mesh mesh, double scale, Vec3 offset = {}) {
    for (auto& p : mesh.vertices) p = {scale * p.x + offset.x, scale * p.y + offset.y, scale * p.z + offset.z};
    return mesh;
}
using Matrix3 = std::array<std::array<double, 3>, 3>;
Matrix3 axisRotation(Vec3 axis, double angle) {
    const double length = std::sqrt(axis.x*axis.x + axis.y*axis.y + axis.z*axis.z);
    axis.x /= length; axis.y /= length; axis.z /= length;
    const double x = axis.x, y = axis.y, z = axis.z, c = std::cos(angle), s = std::sin(angle), v = 1-c;
    return {{{x*x*v+c, x*y*v-z*s, x*z*v+y*s},
             {y*x*v+z*s, y*y*v+c, y*z*v-x*s},
             {z*x*v-y*s, z*y*v+x*s, z*z*v+c}}};
}
Mesh rotateMesh(Mesh mesh, const Matrix3& r) {
    for (auto& p : mesh.vertices) p = {r[0][0]*p.x + r[0][1]*p.y + r[0][2]*p.z,
                                      r[1][0]*p.x + r[1][1]*p.y + r[1][2]*p.z,
                                      r[2][0]*p.x + r[2][1]*p.y + r[2][2]*p.z};
    return mesh;
}
void requireProper(const Matrix3& r) {
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            double value = 0;
            for (int k = 0; k < 3; ++k) value += r[i][k]*r[j][k];
            close(value, i == j ? 1.0 : 0.0);
        }
    const double determinant = r[0][0]*(r[1][1]*r[2][2]-r[1][2]*r[2][1])
        -r[0][1]*(r[1][0]*r[2][2]-r[1][2]*r[2][0])+r[0][2]*(r[1][0]*r[2][1]-r[1][1]*r[2][0]);
    close(determinant, 1.0);
}
void floatRound(Mesh& mesh) {
    for (auto& p : mesh.vertices) p = {float(p.x), float(p.y), float(p.z)};
}
void reverseOrder(Mesh& mesh) {
    std::reverse(mesh.vertices.begin(), mesh.vertices.end());
    for (auto& face : mesh.faces)
        for (auto& index : face) index = mesh.vertices.size() - 1 - index;
    std::reverse(mesh.faces.begin(), mesh.faces.end());
}
void addFloatingTriangle(Mesh& mesh) {
    const auto n = mesh.vertices.size();
    mesh.vertices.insert(mesh.vertices.end(), {{1e7, 2e7, -1e7}, {1e7 + 2, 2e7, -1e7}, {1e7, 2e7 + 3, -1e7}});
    mesh.faces.push_back({n, n + 1, n + 2});
}
void test(const char* name, const std::function<void()>& run) {
    run(); ++passed; std::cout << "PASS " << name << '\n';
}
}

int main() {
    try {
        test("scale, translation, reordered vertices and inverse", [] {
            const auto source = specimen();
            auto target = transform(source, 2.3456789012345, {12.5, -7.25, 31.125});
            reverseOrder(target);
            const auto result = estimateScale(source, target);
            close(result.scale, 2.3456789012345);
            close(result.translation.x, 12.5); close(result.translation.y, -7.25); close(result.translation.z, 31.125);
            close(result.inverseScale, 1.0 / result.scale);
            require(result.relativeMaxError < 1e-14, "Exact clone residual too high");
            const auto reverse = estimateScale(target, source);
            close(result.scale * reverse.scale, 1.0);
            close(reverse.translation.x, -12.5 / result.scale);
        });
        test("missing far floating component and unused vertex", [] {
            auto source = specimen();
            const auto target = transform(source, 0.000123456789, {5, 8, -2});
            addFloatingTriangle(source);
            source.vertices.push_back({1e25, -1e25, 1e25});
            const auto result = estimateScale(source, target);
            close(result.scale, 0.000123456789, 1e-11, 0);
            require(result.sourceVerticesUsed == 8 && result.sourceVerticesIgnored == 4, "Component accounting wrong");
            Options all; all.largestComponent = false;
            rejects([&] { estimateScale(source, target, all); }, "counts differ");
        });
        test("STL triangle soup exact-position welding", [] {
            const auto source = specimen();
            const auto scaled = transform(source, 43.123, {-12, 2, 5});
            Mesh soup;
            for (const auto& face : scaled.faces) {
                const auto n = soup.vertices.size();
                for (auto i : face) soup.vertices.push_back(scaled.vertices[i]);
                soup.faces.push_back({n, n + 1, n + 2});
            }
            const auto result = estimateScale(source, soup);
            close(result.scale, 43.123);
            require(result.targetVerticesUsed == source.vertices.size(), "Failed to weld STL vertices");
        });
        test("float32 coordinate rounding", [] {
            auto source = specimen();
            for (auto& p : source.vertices) p = {float(p.x), float(p.y), float(p.z)};
            auto target = transform(source, 123.456789, {23, -55, 150});
            for (auto& p : target.vertices) p = {float(p.x), float(p.y), float(p.z)};
            const auto result = estimateScale(source, target);
            close(result.scale, 123.456789, 1e-7, 0);
            require(result.relativeMaxError < 1e-7, "Float rounding residual too high");
            Options strict; strict.relativeTolerance = 1e-12;
            rejects([&] { estimateScale(source, target, strict); }, "verification failed");
        });
        test("extreme representable scale ranges", [] {
            const auto original = specimen();
            for (const auto& scales : {std::pair<double,double>{1e-100, 1e100}, {1e100, 1e-100}, {1e-150, 1e-120}, {1e150, 1e120}}) {
                const auto result = estimateScale(transform(original, scales.first), transform(original, scales.second));
                close(result.scale, scales.second / scales.first, 2e-14, 0);
                require(result.relativeMaxError < 1e-14, "Extreme scale residual too high");
            }
        });
        test("large absolute offsets retain local precision", [] {
            const auto source = transform(specimen(), 1, {1e10, -1e10, 1e10});
            const auto target = transform(source, 8, {-1e11, 2e11, 3e10});
            const auto result = estimateScale(source, target);
            close(result.scale, 8, 1e-6, 0);
            require(result.relativeMaxError < 1e-6, "Offset residual too high");
        });
        test("identical coordinate marginals cannot spoof 3D geometry", [] {
            Mesh a, b;
            a.vertices = {{0, 0, 0}, {1, 2, 4}, {2, 4, 1}, {4, 1, 2}};
            b.vertices = {{0, 0, 0}, {1, 4, 2}, {2, 1, 4}, {4, 2, 1}};
            rejects([&] { estimateScale(a, b); }, "verification failed");
        });
        test("nonuniform scale and changed geometry rejected", [] {
            const auto source = specimen();
            auto nonuniform = source;
            for (auto& p : nonuniform.vertices) { p.x *= 2; p.y *= 2.02; p.z *= 2; }
            rejects([&] { estimateScale(source, nonuniform); }, "verification failed");
            auto modified = transform(source, 2);
            modified.vertices[3].x += 0.1;
            rejects([&] { estimateScale(source, modified); }, "verification failed");
        });
        test("rotation rejected for nonsymmetric geometry", [] {
            const auto source = specimen();
            auto target = source;
            for (auto& p : target.vertices) { const double x = p.x; p.x = -p.y; p.y = x; }
            rejects([&] { estimateScale(source, target); }, "verification failed");
        });
        test("counts, invalid faces, nonfinite data and zero spread rejected", [] {
            const auto source = specimen();
            auto fewer = source; fewer.faces.clear(); fewer.vertices.pop_back();
            rejects([&] { estimateScale(source, fewer); }, "counts differ");
            auto invalid = source; invalid.faces[0][0] = 999;
            rejects([&] { estimateScale(source, invalid); }, "invalid face");
            invalid = source; invalid.vertices[0].x = std::numeric_limits<double>::quiet_NaN();
            rejects([&] { estimateScale(source, invalid); }, "non-finite");
            Mesh coincident; coincident.vertices = {{1, 2, 3}, {1, 2, 3}};
            rejects([&] { estimateScale(coincident, coincident); }, "distinct");
            Mesh empty; rejects([&] { estimateScale(empty, source); }, "no vertices");
            Options bad; bad.relativeTolerance = 0;
            rejects([&] { estimateScale(source, source, bad); }, "tolerance");
        });
        test("point cloud and line", [] {
            Mesh line; line.vertices = {{2, 7, -2}, {5, 7, -2}, {10, 7, -2}};
            const auto result = estimateScale(line, transform(line, 0.25, {-9, 1, 2}));
            close(result.scale, 0.25);
        });
        test("all-component option", [] {
            auto source = specimen(); addFloatingTriangle(source);
            auto target = transform(source, 0.52, {1, -2, 3}); reverseOrder(target);
            Options options; options.largestComponent = false;
            const auto result = estimateScale(source, target, options);
            close(result.scale, 0.52);
            require(result.sourceVerticesUsed == 11 && result.sourceVerticesIgnored == 0, "All-component accounting wrong");
        });
        test("ambiguous largest components rejected", [] {
            auto source = specimen();
            const auto other = transform(specimen(), 2, {30, 40, 50});
            const auto offset = source.vertices.size();
            source.vertices.insert(source.vertices.end(), other.vertices.begin(), other.vertices.end());
            for (auto face : other.faces) {
                for (auto& i : face) i += offset;
                source.faces.push_back(face);
            }
            const auto target = transform(source, 3);
            rejects([&] { estimateScale(source, target); }, "equally large");
            Options all; all.largestComponent = false;
            close(estimateScale(source, target, all).scale, 3);
        });
        test("larger shuffled point cloud", [] {
            Mesh source;
            std::mt19937 generator(42);
            std::uniform_real_distribution<double> distribution(-10, 10);
            for (int i = 0; i < 20000; ++i)
                source.vertices.push_back({distribution(generator), distribution(generator), distribution(generator)});
            auto target = transform(source, 17.123456789, {1, 2, 3});
            std::shuffle(target.vertices.begin(), target.vertices.end(), generator);
            const auto result = estimateScale(source, target);
            close(result.scale, 17.123456789);
            require(result.relativeMaxError < 1e-14, "Large point-cloud residual too high");
        });
        test("arbitrary proper rotations and full transformation output", [] {
            Options options; options.allowRotation = true;
            for (double angle : {0.0, 0.001, 0.7, 1.5707963267948966, 3.141592653589793, -2.9}) {
                const auto r = axisRotation({1, -3, 2}, angle);
                auto source = specimen();
                auto target = transform(rotateMesh(source, r), 2.3456789012345, {12, -17, 3});
                addFloatingTriangle(source); reverseOrder(target);
                const auto result = estimateScale(source, target, options);
                close(result.scale, 2.3456789012345);
                close(result.translation.x, 12); close(result.translation.y, -17); close(result.translation.z, 3);
                requireProper(result.rotation);
                for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) close(result.rotation[i][j], r[i][j]);
                const auto inverse = estimateScale(target, source, options);
                close(result.scale * inverse.scale, 1);
                require(result.relativeMaxError < 1e-13, "Rotation fit residual too high");
            }
        });
        test("rotation checkbox preserves default rejection", [] {
            const auto source = specimen();
            const auto target = rotateMesh(source, axisRotation({1,2,3}, 1.2));
            rejects([&] { estimateScale(source, target); }, "verification failed");
            Options options; options.allowRotation = true;
            close(estimateScale(source, target, options).scale, 1);
        });
        test("isotropic cube tetrahedron and octahedron", [] {
            Mesh cube, tetra, octa;
            for (double x : {-1.0, 1.0}) for (double y : {-1.0, 1.0}) for (double z : {-1.0, 1.0}) cube.vertices.push_back({x,y,z});
            tetra.vertices = {{1,1,1}, {1,-1,-1}, {-1,1,-1}, {-1,-1,1}};
            octa.vertices = {{1,0,0}, {-1,0,0}, {0,1,0}, {0,-1,0}, {0,0,1}, {0,0,-1}};
            Options options; options.allowRotation = true;
            for (auto source : {cube, tetra, octa}) {
                auto target = transform(rotateMesh(source, axisRotation({1,3,-7}, 2.7)), 0.03125, {1,-2,3});
                reverseOrder(target);
                const auto result = estimateScale(source, target, options);
                close(result.scale, 0.03125);
                requireProper(result.rotation);
                require(result.relativeMaxError < 1e-13, "Isotropic geometry fit failed");
            }
        });
        test("planar rotation out of plane and collinear geometry", [] {
            Mesh polygon, line, pair;
            polygon.vertices = {{1,0,0}, {0,1,0}, {-1,0,0}, {0,-1,0}};
            line.vertices = {{-3,0,0}, {-1,0,0}, {0,0,0}, {1,0,0}, {3,0,0}};
            pair.vertices = {{-1,0,0}, {1,0,0}};
            Options options; options.allowRotation = true;
            for (const auto& source : {polygon, line, pair}) {
                const auto target = transform(rotateMesh(source, axisRotation({3,7,1}, 2.1)), 7.25, {7,-8,9});
                const auto result = estimateScale(source, target, options);
                close(result.scale, 7.25); requireProper(result.rotation);
                require(result.relativeMaxError < 1e-13, "Planar/line residual too high");
            }
        });
        test("float-rounded rotations with repeated and near-equal radii", [] {
            Mesh shell;
            std::mt19937 generator(57);
            std::normal_distribution<double> distribution;
            for (int i = 0; i < 80; ++i) {
                Vec3 p{distribution(generator), distribution(generator), distribution(generator)};
                const double s = (1 + i*1e-8) / std::sqrt(p.x*p.x+p.y*p.y+p.z*p.z);
                p = {p.x*s, p.y*s, p.z*s};
                shell.vertices.push_back(p); shell.vertices.push_back({-p.x,-p.y,-p.z});
            }
            floatRound(shell);
            auto target = transform(rotateMesh(shell, axisRotation({2,-5,1}, 1.345)), 123.456, {23,-14,18});
            floatRound(target); reverseOrder(target);
            Options options; options.allowRotation = true;
            const auto result = estimateScale(shell, target, options);
            close(result.scale, 123.456, 2e-7, 0);
            require(result.relativeMaxError < 1e-6, "Near-tied radii residual too high");
            requireProper(result.rotation);
        });
        test("float-rounded thin and collinear source and target", [] {
            Options options; options.allowRotation = true;
            for (double thickness : {0.0, 1e-8, 1e-6}) {
                Mesh line;
                for (int i = -3; i <= 3; ++i) line.vertices.push_back({double(i), (i%2)*thickness, 0});
                auto source = transform(rotateMesh(line, axisRotation({2,1,3}, 1.1)), 1, {2,-3,4});
                floatRound(source);
                auto target = transform(rotateMesh(source, axisRotation({-3,2,1}, 2.2)), 3.125, {-1,5,8});
                floatRound(target);
                const auto result = estimateScale(source, target, options);
                close(result.scale, 3.125, 2e-7, 0);
                require(result.relativeMaxError < 1e-5, "Thin geometry residual too high");
                requireProper(result.rotation);
            }
        });
        test("chiral reflection and changed geometry rejected with rotation enabled", [] {
            const auto source = specimen();
            auto reflected = source;
            for (auto& p : reflected.vertices) p.x = -p.x;
            Options options; options.allowRotation = true;
            rejects([&] { estimateScale(source, reflected, options); }, "verification failed");
            auto changed = transform(rotateMesh(source, axisRotation({1,2,3}, 2)), 4);
            changed.vertices[2].x += 0.2;
            rejects([&] { estimateScale(source, changed, options); }, "verification failed");
            auto nonuniform = source;
            for (auto& p : nonuniform.vertices) { p.x *= 2; p.y *= 2.1; p.z *= 1.9; }
            nonuniform = rotateMesh(nonuniform, axisRotation({1,2,3}, 2));
            rejects([&] { estimateScale(source, nonuniform, options); }, "verification failed");
        });
        test("rotated extreme scales and large offsets", [] {
            Options options; options.allowRotation = true;
            const auto r = axisRotation({1,-3,2}, 0.72);
            for (const auto& scales : {std::pair<double,double>{1e-100,1e100}, {1e100,1e-100}}) {
                const auto source = transform(specimen(), scales.first);
                const auto target = transform(rotateMesh(specimen(), r), scales.second);
                close(estimateScale(source, target, options).scale, scales.second/scales.first, 2e-14, 0);
            }
            const auto source = transform(specimen(), 1, {1e8,-1e8,1e8});
            const auto target = transform(rotateMesh(source, r), 8, {-1e9,2e9,3e8});
            const auto result = estimateScale(source, target, options);
            close(result.scale, 8, 1e-7, 0);
            require(result.relativeMaxError < 1e-7, "Far-offset rotation failed");
        });
        test("large rotated reordered point cloud", [] {
            Mesh source;
            std::mt19937 generator(2026);
            std::uniform_real_distribution<double> distribution(-10,10);
            for (int i = 0; i < 20000; ++i) source.vertices.push_back({distribution(generator),distribution(generator),distribution(generator)});
            auto target = transform(rotateMesh(source, axisRotation({3,4,1}, 2.45)), 17.123456789, {-5,7,8});
            std::shuffle(target.vertices.begin(), target.vertices.end(), generator);
            Options options; options.allowRotation = true;
            const auto result = estimateScale(source, target, options);
            close(result.scale, 17.123456789);
            require(result.relativeMaxError < 1e-13, "Large rotated cloud failed");
        });
        std::cout << passed << " tests passed\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL after " << passed << " tests: " << e.what() << '\n';
        return 1;
    }
}
