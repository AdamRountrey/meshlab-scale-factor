// SPDX-License-Identifier: GPL-2.0-or-later
#include "scale_estimator.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace meshlab_scale {
namespace {

double coordinate(const Vec3& p, int axis) {
    return axis == 0 ? p.x : (axis == 1 ? p.y : p.z);
}
void setCoordinate(Vec3& p, int axis, double value) {
    if (axis == 0) p.x = value;
    else if (axis == 1) p.y = value;
    else p.z = value;
}
bool finite(const Vec3& p) {
    return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
}
Vec3 multiply(const Vec3& p, double s) {
    return {p.x * s, p.y * s, p.z * s};
}
double norm(const Vec3& p) { return std::hypot(std::hypot(p.x, p.y), p.z); }
using Matrix3 = std::array<std::array<double, 3>, 3>;
const Matrix3 identity{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
Vec3 subtract(const Vec3& a, const Vec3& b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
Vec3 rotate(const Matrix3& r, const Vec3& p) {
    return {std::fma(r[0][0], p.x, std::fma(r[0][1], p.y, r[0][2] * p.z)),
            std::fma(r[1][0], p.x, std::fma(r[1][1], p.y, r[1][2] * p.z)),
            std::fma(r[2][0], p.x, std::fma(r[2][1], p.y, r[2][2] * p.z))};
}
bool lessPoint(const Vec3& a, const Vec3& b) {
    if (a.x != b.x) return a.x < b.x;
    if (a.y != b.y) return a.y < b.y;
    return a.z < b.z;
}
bool samePoint(const Vec3& a, const Vec3& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
}

class Sum {
    double sum_ = 0.0, correction_ = 0.0;
public:
    void add(double value) {
        const double t = sum_ + value;
        if (std::abs(sum_) >= std::abs(value)) correction_ += (sum_ - t) + value;
        else correction_ += (value - t) + sum_;
        sum_ = t;
    }
    double value() const { return sum_ + correction_; }
};

class DisjointSet {
    std::vector<std::size_t> parent_, size_;
public:
    explicit DisjointSet(std::size_t count) : parent_(count), size_(count, 1) {
        std::iota(parent_.begin(), parent_.end(), 0);
    }
    std::size_t root(std::size_t i) {
        while (parent_[i] != i) {
            parent_[i] = parent_[parent_[i]];
            i = parent_[i];
        }
        return i;
    }
    void unite(std::size_t a, std::size_t b) {
        a = root(a); b = root(b);
        if (a == b) return;
        if (size_[a] < size_[b]) std::swap(a, b);
        parent_[b] = a; size_[a] += size_[b];
    }
};

std::vector<Vec3> retainedPoints(const Mesh& mesh, bool largest, const char* name) {
    if (mesh.vertices.empty()) throw std::runtime_error(std::string(name) + " mesh has no vertices.");
    for (const auto& p : mesh.vertices)
        if (!finite(p)) throw std::runtime_error(std::string(name) + " mesh contains non-finite coordinates.");
    for (const auto& face : mesh.faces)
        for (std::size_t i : face)
            if (i >= mesh.vertices.size()) throw std::runtime_error(std::string(name) + " mesh has an invalid face index.");

    std::vector<std::size_t> order(mesh.vertices.size()), remap(mesh.vertices.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        return lessPoint(mesh.vertices[a], mesh.vertices[b]);
    });
    std::vector<Vec3> unique;
    unique.reserve(mesh.vertices.size());
    for (std::size_t i : order) {
        if (unique.empty() || !samePoint(unique.back(), mesh.vertices[i])) unique.push_back(mesh.vertices[i]);
        remap[i] = unique.size() - 1;
    }
    // Exact-position welding makes component extraction work for STL triangle
    // soup and duplicated texture seams without moving any coordinates.
    if (largest && !mesh.faces.empty()) {
        DisjointSet sets(unique.size());
        std::vector<bool> used(unique.size(), false);
        for (const auto& face : mesh.faces) {
            const auto a = remap[face[0]], b = remap[face[1]], c = remap[face[2]];
            used[a] = used[b] = used[c] = true;
            sets.unite(a, b); sets.unite(a, c);
        }
        std::vector<std::size_t> counts(unique.size(), 0);
        for (std::size_t i = 0; i < unique.size(); ++i)
            if (used[i]) ++counts[sets.root(i)];
        std::size_t selected = 0;
        for (std::size_t i = 1; i < counts.size(); ++i)
            if (counts[i] > counts[selected]) selected = i;
        if (std::count(counts.begin(), counts.end(), counts[selected]) > 1)
            throw std::runtime_error(std::string(name) +
                " mesh has equally large connected components. Extract the intended component into its own mesh layer first, "
                "or disable largest-component mode when all geometry is present in both models.");
        std::vector<Vec3> component;
        component.reserve(counts[selected]);
        for (std::size_t i = 0; i < unique.size(); ++i)
            if (used[i] && sets.root(i) == selected) component.push_back(unique[i]);
        unique = std::move(component);
    }
    if (unique.size() < 2) throw std::runtime_error(std::string(name) + " mesh needs at least two distinct retained positions.");
    return unique;
}

struct Normalized {
    std::vector<Vec3> points;
    Vec3 anchor, mean;
    double extent = 0.0;
    double diagonal = 0.0;
};

Normalized normalize(std::vector<Vec3> points) {
    Normalized result;
    Vec3 low = points.front(), high = points.front();
    for (const auto& p : points) {
        for (int axis = 0; axis < 3; ++axis) {
            setCoordinate(low, axis, std::min(coordinate(low, axis), coordinate(p, axis)));
            setCoordinate(high, axis, std::max(coordinate(high, axis), coordinate(p, axis)));
        }
    }
    for (int axis = 0; axis < 3; ++axis) {
        const double lo = coordinate(low, axis), hi = coordinate(high, axis);
        const double span = hi - lo;
        // Midpoint chosen without overflow, even for opposite huge coordinates.
        setCoordinate(result.anchor, axis, std::isfinite(span) ? lo + span * 0.5 : lo * 0.5 + hi * 0.5);
    }
    for (const auto& p : points)
        for (int axis = 0; axis < 3; ++axis)
            result.extent = std::max(result.extent, std::abs(coordinate(p, axis) - coordinate(result.anchor, axis)));
    if (!(result.extent > 0.0) || !std::isfinite(result.extent))
        throw std::runtime_error("Mesh has zero or unrepresentable coordinate spread.");
    Sum means[3];
    for (auto& p : points) {
        for (int axis = 0; axis < 3; ++axis) {
            const double value = (coordinate(p, axis) - coordinate(result.anchor, axis)) / result.extent;
            setCoordinate(p, axis, value);
            means[axis].add(value);
        }
    }
    Vec3 normalizedRange;
    for (int axis = 0; axis < 3; ++axis) {
        setCoordinate(result.mean, axis, means[axis].value() / static_cast<double>(points.size()));
        setCoordinate(normalizedRange, axis,
            (coordinate(high, axis) - coordinate(result.anchor, axis)) / result.extent -
            (coordinate(low, axis) - coordinate(result.anchor, axis)) / result.extent);
    }
    result.diagonal = norm(normalizedRange);
    for (auto& p : points)
        for (int axis = 0; axis < 3; ++axis)
            setCoordinate(p, axis, coordinate(p, axis) - coordinate(result.mean, axis));
    result.points = std::move(points);
    return result;
}

// Balanced kd-tree over normalized coordinates; construction O(n log n),
// expected logarithmic queries, and linear storage. Node index is median position.
class KdTree {
    const std::vector<Vec3>& points_;
    std::vector<std::size_t> order_;
    std::vector<unsigned char> axes_;
    void build(std::size_t begin, std::size_t end) {
        if (begin == end) return;
        Vec3 lo = points_[order_[begin]], hi = lo;
        for (std::size_t i = begin + 1; i < end; ++i)
            for (int a = 0; a < 3; ++a) {
                setCoordinate(lo, a, std::min(coordinate(lo, a), coordinate(points_[order_[i]], a)));
                setCoordinate(hi, a, std::max(coordinate(hi, a), coordinate(points_[order_[i]], a)));
            }
        int axis = 0;
        for (int a = 1; a < 3; ++a)
            if (coordinate(hi, a) - coordinate(lo, a) > coordinate(hi, axis) - coordinate(lo, axis)) axis = a;
        const auto middle = begin + (end - begin) / 2;
        std::nth_element(order_.begin() + begin, order_.begin() + middle, order_.begin() + end,
            [&](std::size_t a, std::size_t b) { return coordinate(points_[a], axis) < coordinate(points_[b], axis); });
        axes_[middle] = static_cast<unsigned char>(axis);
        build(begin, middle); build(middle + 1, end);
    }
    void nearest(const Vec3& p, std::size_t begin, std::size_t end, double& best, std::size_t& index) const {
        if (begin == end) return;
        const auto middle = begin + (end - begin) / 2;
        const auto& q = points_[order_[middle]];
        const double dx = p.x - q.x, dy = p.y - q.y, dz = p.z - q.z;
        const double distance = dx * dx + dy * dy + dz * dz;
        if (distance < best) { best = distance; index = order_[middle]; }
        const double delta = coordinate(p, axes_[middle]) - coordinate(q, axes_[middle]);
        if (delta < 0.0) {
            nearest(p, begin, middle, best, index);
            if (delta * delta <= best) nearest(p, middle + 1, end, best, index);
        } else {
            nearest(p, middle + 1, end, best, index);
            if (delta * delta <= best) nearest(p, begin, middle, best, index);
        }
    }
public:
    explicit KdTree(const std::vector<Vec3>& points) : points_(points), order_(points.size()), axes_(points.size()) {
        std::iota(order_.begin(), order_.end(), 0);
        build(0, order_.size());
    }
    double squaredDistance(const Vec3& p) const {
        double best = std::numeric_limits<double>::infinity();
        std::size_t index = 0;
        nearest(p, 0, order_.size(), best, index);
        return best;
    }
    std::size_t nearestIndex(const Vec3& p) const {
        double best = std::numeric_limits<double>::infinity();
        std::size_t index = 0;
        nearest(p, 0, order_.size(), best, index);
        return index;
    }
};

// Horn's absolute orientation: the largest eigenvector of this symmetric 4x4
// matrix is a unit quaternion. Its rotation is necessarily proper (det = +1).
Matrix3 hornRotation(const std::vector<Vec3>& source, const std::vector<Vec3>& target,
                     const std::vector<std::pair<std::size_t, std::size_t>>& pairs) {
    Sum covariance[3][3];
    for (const auto& pair : pairs)
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                covariance[i][j].add(coordinate(source[pair.first], i) * coordinate(target[pair.second], j));
    double h[3][3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) h[i][j] = covariance[i][j].value();
    double n[4][4] = {
        {h[0][0] + h[1][1] + h[2][2], h[1][2] - h[2][1], h[2][0] - h[0][2], h[0][1] - h[1][0]},
        {h[1][2] - h[2][1], h[0][0] - h[1][1] - h[2][2], h[0][1] + h[1][0], h[0][2] + h[2][0]},
        {h[2][0] - h[0][2], h[0][1] + h[1][0], -h[0][0] + h[1][1] - h[2][2], h[1][2] + h[2][1]},
        {h[0][1] - h[1][0], h[0][2] + h[2][0], h[1][2] + h[2][1], -h[0][0] - h[1][1] + h[2][2]}};
    double eigenvectors[4][4] = {{1,0,0,0}, {0,1,0,0}, {0,0,1,0}, {0,0,0,1}};
    double magnitude = 0;
    for (const auto& row : n) for (double value : row) magnitude = std::max(magnitude, std::abs(value));
    for (int iteration = 0; iteration < 80; ++iteration) {
        int p = 0, q = 1;
        for (int i = 0; i < 4; ++i)
            for (int j = i + 1; j < 4; ++j)
                if (std::abs(n[i][j]) > std::abs(n[p][q])) { p = i; q = j; }
        if (std::abs(n[p][q]) <= 2e-16 * magnitude) break;
        const double angle = 0.5 * std::atan2(2 * n[p][q], n[q][q] - n[p][p]);
        const double c = std::cos(angle), s = std::sin(angle);
        const double pp = n[p][p], qq = n[q][q], pq = n[p][q];
        n[p][p] = c*c*pp - 2*s*c*pq + s*s*qq;
        n[q][q] = s*s*pp + 2*s*c*pq + c*c*qq;
        n[p][q] = n[q][p] = 0;
        for (int k = 0; k < 4; ++k) {
            if (k != p && k != q) {
                const double kp = n[k][p], kq = n[k][q];
                n[k][p] = n[p][k] = c*kp - s*kq;
                n[k][q] = n[q][k] = s*kp + c*kq;
            }
            const double vp = eigenvectors[k][p], vq = eigenvectors[k][q];
            eigenvectors[k][p] = c*vp - s*vq;
            eigenvectors[k][q] = s*vp + c*vq;
        }
    }
    int largest = 0;
    for (int i = 1; i < 4; ++i) if (n[i][i] > n[largest][largest]) largest = i;
    double w = eigenvectors[0][largest], x = eigenvectors[1][largest];
    double y = eigenvectors[2][largest], z = eigenvectors[3][largest];
    const double length = std::hypot(std::hypot(w,x), std::hypot(y,z));
    w /= length; x /= length; y /= length; z /= length;
    return {{{1 - 2*(y*y + z*z), 2*(x*y - z*w), 2*(x*z + y*w)},
             {2*(x*y + z*w), 1 - 2*(x*x + z*z), 2*(y*z - x*w)},
             {2*(x*z - y*w), 2*(y*z + x*w), 1 - 2*(x*x + y*y)}}};
}

std::array<Vec3, 3> orientedFrame(const Vec3& first, const Vec3& second) {
    const auto u = multiply(first, 1.0 / norm(first));
    const auto perpendicular = cross(u, second);
    const auto w = multiply(perpendicular, 1.0 / norm(perpendicular));
    return {u, cross(w, u), w};
}
Vec3 auxiliaryAxis(const Vec3& p) {
    int axis = 0;
    for (int i = 1; i < 3; ++i)
        if (std::abs(coordinate(p, i)) < std::abs(coordinate(p, axis))) axis = i;
    Vec3 result; setCoordinate(result, axis, 1); return result;
}
Matrix3 betweenFrames(const std::array<Vec3, 3>& source, const std::array<Vec3, 3>& target) {
    Matrix3 result{};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                result[i][j] += coordinate(target[k], i) * coordinate(source[k], j);
    return result;
}
bool properRotation(const Matrix3& r) {
    const Vec3 rows[3] = {{r[0][0], r[0][1], r[0][2]}, {r[1][0], r[1][1], r[1][2]}, {r[2][0], r[2][1], r[2][2]}};
    for (int i = 0; i < 3; ++i) {
        if (!finite(rows[i])) return false;
        for (int j = 0; j < 3; ++j)
            if (std::abs(dot(rows[i], rows[j]) - (i == j ? 1.0 : 0.0)) > 1e-10) return false;
    }
    return std::abs(dot(rows[0], cross(rows[1], rows[2])) - 1.0) < 1e-10;
}

Matrix3 findRotation(const Normalized& source, const Normalized& target, double scale, double tolerance) {
    const auto& a = source.points;
    const auto& b = target.points;
    const KdTree targetTree(b);
    const double absoluteTolerance = tolerance * target.diagonal;
    // A wider proposal threshold allows rounded anchor coordinates; acceptance
    // always uses the original user tolerance on the entire 3D geometry.
    const double proposalTolerance = std::max(4 * absoluteTolerance, 128 * std::numeric_limits<double>::epsilon());
    const double absoluteSquared = absoluteTolerance * absoluteTolerance;
    const double proposalSquared = 16 * proposalTolerance * proposalTolerance;
    std::size_t fullTrials = 0, frameTrials = 0, pairChecks = 0;
    constexpr std::size_t maxFullTrials = 16, maxFrameTrials = 4096, maxPairChecks = 2000000;
    bool searchLimited = false;
    Matrix3 accepted = identity;
    auto tryRotation = [&](Matrix3 rotation) {
        if (!properRotation(rotation)) return false;
        const auto samples = std::min<std::size_t>(64, a.size());
        for (std::size_t sample = 0; sample < samples; ++sample) {
            const auto i = sample * (a.size() - 1) / (samples - 1);
            if (targetTree.squaredDistance(multiply(rotate(rotation, a[i]), scale)) > proposalSquared) return false;
        }
        if (++fullTrials > maxFullTrials) { searchLimited = true; return false; }
        auto fullCheck = [&](const Matrix3& candidate) {
            if (!properRotation(candidate)) return false;
            std::vector<Vec3> transformed;
            transformed.reserve(a.size());
            for (const auto& p : a) {
                const auto q = multiply(rotate(candidate, p), scale);
                if (targetTree.squaredDistance(q) > absoluteSquared) return false;
                transformed.push_back(q);
            }
            const KdTree sourceTree(transformed);
            for (const auto& p : b) if (sourceTree.squaredDistance(p) > absoluteSquared) return false;
            accepted = candidate;
            return true;
        };
        // A sound frame can already satisfy the requested tolerance. Do not
        // discard it if nearest-neighbor ties would degrade a fitted proposal.
        if (fullCheck(rotation)) return true;
        std::vector<std::pair<std::size_t, std::size_t>> pairs(a.size());
        // These are local roundoff refinements of a checked rigid proposal, not
        // unconstrained registration. Scale remains the symmetric RMS ratio.
        for (int iteration = 0; iteration < 2; ++iteration) {
            for (std::size_t i = 0; i < a.size(); ++i)
                pairs[i] = {i, targetTree.nearestIndex(multiply(rotate(rotation, a[i]), scale))};
            rotation = hornRotation(a, b, pairs);
        }
        return fullCheck(rotation);
    };

    if (tryRotation(identity)) return accepted;

    using RadiusIndex = std::pair<double, std::size_t>;
    std::vector<RadiusIndex> ra, rb;
    ra.reserve(a.size()); rb.reserve(b.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
        ra.emplace_back(scale * norm(a[i]), i);
        rb.emplace_back(norm(b[i]), i);
    }
    std::sort(ra.begin(), ra.end()); std::sort(rb.begin(), rb.end());
    std::vector<std::pair<std::size_t, std::size_t>> unambiguous;
    for (std::size_t i = 0; i < ra.size(); ++i) {
        if (std::abs(ra[i].first - rb[i].first) > proposalTolerance) continue;
        if (i > 0 && (ra[i].first - ra[i-1].first <= 2*proposalTolerance || rb[i].first - rb[i-1].first <= 2*proposalTolerance)) continue;
        if (i + 1 < ra.size() && (ra[i+1].first - ra[i].first <= 2*proposalTolerance || rb[i+1].first - rb[i].first <= 2*proposalTolerance)) continue;
        unambiguous.emplace_back(ra[i].second, rb[i].second);
    }
    if (unambiguous.size() >= 3 && tryRotation(hornRotation(a, b, unambiguous))) return accepted;

    const auto first = ra.back().second;
    std::size_t second = first;
    double perpendicular = 0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        const double distance = norm(cross(a[first], a[i])) / norm(a[first]);
        if (distance > perpendicular) { perpendicular = distance; second = i; }
    }
    const bool line = perpendicular * scale <= std::max(0.1 * absoluteTolerance, 1e-13 * target.diagonal);
    const auto sourceFrame = orientedFrame(a[first], line ? auxiliaryAxis(a[first]) : a[second]);
    const double firstRadius = scale * norm(a[first]), secondRadius = scale * norm(a[second]);
    const double pairDistance = scale * norm(subtract(a[first], a[second]));
    std::vector<std::size_t> firstCandidates, secondCandidates;
    for (std::size_t i = 0; i < b.size(); ++i) {
        if (std::abs(norm(b[i]) - firstRadius) <= proposalTolerance) firstCandidates.push_back(i);
        if (std::abs(norm(b[i]) - secondRadius) <= proposalTolerance) secondCandidates.push_back(i);
    }
    for (auto i : firstCandidates) {
        if (norm(b[i]) <= std::numeric_limits<double>::epsilon()) continue;
        if (line) {
            if (++frameTrials > maxFrameTrials) { searchLimited = true; break; }
            if (tryRotation(betweenFrames(sourceFrame, orientedFrame(b[i], auxiliaryAxis(b[i]))))) return accepted;
        } else {
            for (auto j : secondCandidates) {
                if (++pairChecks > maxPairChecks) { searchLimited = true; break; }
                if (std::abs(norm(subtract(b[i], b[j])) - pairDistance) > 2*proposalTolerance) continue;
                if (norm(cross(b[i], b[j])) <= 1e-14 * norm(b[i]) * norm(b[j])) continue;
                if (++frameTrials > maxFrameTrials) { searchLimited = true; break; }
                if (tryRotation(betweenFrames(sourceFrame, orientedFrame(b[i], b[j])))) return accepted;
                if (searchLimited) break;
            }
        }
        if (searchLimited) break;
    }
    if (searchLimited)
        throw std::runtime_error("Rotation search reached its ambiguity limit (4096 frame proposals, 2000000 anchor pairs, or 16 full checks). "
            "No scale was accepted. Isolate a distinctive shared component or roughly align the models first.");
    throw std::runtime_error("Geometry verification failed: no proper 3D rotation matches the retained models at the selected tolerance. "
        "Check for changed geometry, nonuniform scale, reflection, or coordinate rounding beyond the tolerance.");
}

} // namespace

Result estimateScale(const Mesh& source, const Mesh& target, const Options& options) {
    if (!(options.relativeTolerance > 0.0) || !std::isfinite(options.relativeTolerance))
        throw std::runtime_error("Relative tolerance must be finite and positive.");
    auto sourcePoints = retainedPoints(source, options.largestComponent, "Source");
    auto targetPoints = retainedPoints(target, options.largestComponent, "Target");
    if (sourcePoints.size() != targetPoints.size()) {
        std::ostringstream message;
        message << "Retained unique vertex counts differ (" << sourcePoints.size() << " versus "
                << targetPoints.size() << "). This filter requires identical retained geometry; "
                << "resampled or partially deleted main components are not supported.";
        throw std::runtime_error(message.str());
    }
    Result result;
    result.sourceVerticesUsed = sourcePoints.size();
    result.targetVerticesUsed = targetPoints.size();
    result.sourceVerticesIgnored = source.vertices.size() - sourcePoints.size();
    result.targetVerticesIgnored = target.vertices.size() - targetPoints.size();
    const auto a = normalize(std::move(sourcePoints));
    const auto b = normalize(std::move(targetPoints));
    double normalizedScale = 0;
    if (options.allowRotation) {
        Sum sourceSquaredRadius, targetSquaredRadius;
        for (std::size_t i = 0; i < a.points.size(); ++i) {
            sourceSquaredRadius.add(dot(a.points[i], a.points[i]));
            targetSquaredRadius.add(dot(b.points[i], b.points[i]));
        }
        normalizedScale = std::sqrt(targetSquaredRadius.value() / sourceSquaredRadius.value());
        result.rotation = findRotation(a, b, normalizedScale, options.relativeTolerance);
    } else {
        Sum covariance, variance;
        std::vector<double> x(a.points.size()), y(b.points.size());
        for (int axis = 0; axis < 3; ++axis) {
            for (std::size_t i = 0; i < x.size(); ++i) {
                x[i] = coordinate(a.points[i], axis); y[i] = coordinate(b.points[i], axis);
            }
            std::sort(x.begin(), x.end()); std::sort(y.begin(), y.end());
            for (std::size_t i = 0; i < x.size(); ++i) {
                covariance.add(x[i] * y[i]); variance.add(x[i] * x[i]);
            }
        }
        normalizedScale = covariance.value() / variance.value();
    }
    result.scale = (b.extent / a.extent) * normalizedScale;
    result.inverseScale = 1.0 / result.scale;
    if (!(normalizedScale > 0.0) || !std::isfinite(normalizedScale) ||
        !(result.scale > 0.0) || !std::isfinite(result.scale) || !std::isfinite(result.inverseScale))
        throw std::runtime_error("Scale or reciprocal is not representable as a finite positive double.");

    // Keep comparisons centered and normalized. Expanding scale*p+translation
    // here would lose small details when a mesh is far from the origin.
    std::vector<Vec3> transformed;
    transformed.reserve(a.points.size());
    for (const auto& p : a.points) transformed.push_back(multiply(rotate(result.rotation, p), normalizedScale));
    const KdTree targetTree(b.points), sourceTree(transformed);
    Sum squaredErrors;
    double maxSquared = 0.0;
    for (std::size_t i = 0; i < transformed.size(); ++i) {
        const double forward = targetTree.squaredDistance(transformed[i]);
        const double backward = sourceTree.squaredDistance(b.points[i]);
        squaredErrors.add(forward); squaredErrors.add(backward);
        maxSquared = std::max(maxSquared, std::max(forward, backward));
    }
    const double normalizedRms = std::sqrt(squaredErrors.value() / (2.0 * transformed.size()));
    const double normalizedMax = std::sqrt(maxSquared);
    result.relativeRmsError = normalizedRms / b.diagonal;
    result.relativeMaxError = normalizedMax / b.diagonal;
    result.rmsError = normalizedRms * b.extent;
    result.maxError = normalizedMax * b.extent;
    if (!std::isfinite(result.relativeMaxError) || result.relativeMaxError > options.relativeTolerance) {
        std::ostringstream message;
        message.precision(8);
        message << "Geometry verification failed: maximum bidirectional vertex distance / target diagonal = "
                << result.relativeMaxError << ", tolerance = " << options.relativeTolerance
                << (options.allowRotation ? ". The retained models must have the same geometry under a proper rotation and uniform scale."
                                          : ". The retained models must have the same orientation and geometry under uniform scaling.");
        throw std::runtime_error(message.str());
    }
    const auto rotatedAnchor = rotate(result.rotation, a.anchor);
    const auto rotatedMean = rotate(result.rotation, a.mean);
    for (int axis = 0; axis < 3; ++axis) {
        // FMA reduces cancellation when source/target offsets nearly agree.
        const double anchorShift = std::fma(-result.scale, coordinate(rotatedAnchor, axis), coordinate(b.anchor, axis));
        const double meanShift = b.extent * std::fma(-normalizedScale, coordinate(rotatedMean, axis), coordinate(b.mean, axis));
        setCoordinate(result.translation, axis, anchorShift + meanShift);
    }
    if (!finite(result.translation) || !std::isfinite(result.rmsError) || !std::isfinite(result.maxError))
        throw std::runtime_error("Translation or residual is not representable as a finite double.");
    return result;
}

} // namespace meshlab_scale
