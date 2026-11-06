// SPDX-License-Identifier: GPL-2.0-or-later
// Loads the actual DLL through Qt and exercises the public MeshLab filter API.
#include <common/plugins/interfaces/filter_plugin.h>
#include <vcg/complex/allocate.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <QCoreApplication>
#include <QPluginLoader>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using Position = std::array<double, 3>;
using Rotation = std::array<std::array<double, 3>, 3>;
using Transform = std::array<std::array<double, 4>, 4>;
const std::vector<Position> points = {
    {0, 0, 0}, {2, 0, 0}, {0, 3, 0}, {0, 0, 5}, {1, 2, 1}};
const std::vector<std::array<int, 3>> faces = {
    {0, 2, 1}, {0, 1, 3}, {0, 3, 2}, {1, 2, 4}, {2, 3, 4}, {3, 1, 4}};

void require(bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(message);
}

void requireNear(double actual, double expected, const char* message)
{
    if (!std::isfinite(actual) || std::abs(actual - expected) > 1e-11 * std::max(1.0, std::abs(expected)))
        throw std::runtime_error(std::string(message) + ": got " + std::to_string(actual));
}

void requireFloatNear(double actual, double expected, const char* message)
{
    // Arbitrary rotations are rounded when stored in MeshLab's float vertices
    // and layer matrices. Keep exact-data assertions above substantially tighter.
    if (!std::isfinite(actual) || std::abs(actual - expected) > 2e-6 * std::max(1.0, std::abs(expected)))
        throw std::runtime_error(std::string(message) + ": got " + std::to_string(actual));
}

Rotation generalRotation()
{
    const double radians = std::acos(-1.0) / 180.0;
    const double cx = std::cos(17 * radians), sx = std::sin(17 * radians);
    const double cy = std::cos(-31 * radians), sy = std::sin(-31 * radians);
    const double cz = std::cos(23 * radians), sz = std::sin(23 * radians);
    // Column-vector convention: Rz(23 degrees) * Ry(-31 degrees) * Rx(17 degrees).
    return {{{cz * cy, cz * sy * sx - sz * cx, cz * sy * cx + sz * sx},
             {sz * cy, sz * sy * sx + cz * cx, sz * sy * cx - cz * sx},
             {-sy, cy * sx, cy * cx}}};
}

Position rotate(const Rotation& rotation, const Position& point)
{
    Position out{};
    for (int row = 0; row < 3; ++row)
        for (int col = 0; col < 3; ++col)
            out[row] += rotation[row][col] * point[col];
    return out;
}

void requireRotation(const std::map<std::string, QVariant>& result, const Rotation& expected)
{
    for (int row = 0; row < 3; ++row)
        for (int col = 0; col < 3; ++col) {
            const auto key = "rotation_" + std::to_string(row) + std::to_string(col);
            requireFloatNear(result.at(key).toDouble(), expected[row][col], "Incorrect fitted rotation");
        }
}

Position meshPosition(const MeshModel& model, std::size_t vertexIndex, bool useLayerTransforms)
{
    const auto& vertex = model.cm.vert.at(vertexIndex).cP();
    Position position{double(vertex[0]), double(vertex[1]), double(vertex[2])};
    if (!useLayerTransforms) return position;
    Position displayed{};
    for (int row = 0; row < 3; ++row) {
        displayed[row] = double(model.cm.Tr[row][3]);
        for (int col = 0; col < 3; ++col)
            displayed[row] += double(model.cm.Tr[row][col]) * position[col];
    }
    return displayed;
}

void requireTransformMaps(const std::map<std::string, QVariant>& result,
                          const MeshModel& source, const MeshModel& target,
                          bool useLayerTransforms, bool floatRounded = false)
{
    Transform transform{};
    for (int row = 0; row < 4; ++row)
        for (int col = 0; col < 4; ++col) {
            const auto key = "transform_" + std::to_string(row) + std::to_string(col);
            bool converted = false;
            transform[row][col] = result.at(key).toDouble(&converted);
            require(converted && std::isfinite(transform[row][col]), "Missing or invalid 4x4 transform entry");
        }
    for (int col = 0; col < 4; ++col)
        requireNear(transform[3][col], col == 3 ? 1.0 : 0.0, "Incorrect homogeneous transform last row");

    // makeMesh(..., reorder=true) places original source vertices at these
    // target indices. Use actual stored mesh positions rather than rebuilding
    // an expected matrix from the returned scale/rotation/translation values.
    // The five asymmetric points span 3D and catch transpose, missing-scale,
    // misplaced-translation and scaled-translation mistakes.
    const std::array<std::size_t, 5> targetIndices{1, 3, 4, 0, 2};
    for (std::size_t index = 0; index < targetIndices.size(); ++index) {
        const Position from = meshPosition(source, index, useLayerTransforms);
        const Position expected = meshPosition(target, targetIndices[index], useLayerTransforms);
        const std::array<double, 4> homogeneous{from[0], from[1], from[2], 1};
        std::array<double, 4> mapped{};
        for (int row = 0; row < 4; ++row)
            for (int col = 0; col < 4; ++col)
                mapped[row] += transform[row][col] * homogeneous[col];
        requireNear(mapped[3], 1, "Transform changed the homogeneous coordinate");
        for (int axis = 0; axis < 3; ++axis) {
            if (floatRounded)
                requireFloatNear(mapped[axis], expected[axis], "4x4 transform does not map source vertices to target");
            else
                requireNear(mapped[axis], expected[axis], "4x4 transform does not map source vertices to target");
        }
    }
    const QString report = result.at("report").toString();
    require(report.contains("4x4", Qt::CaseInsensitive) && report.contains("transform", Qt::CaseInsensitive),
            "Copyable report does not identify the full 4x4 transform");
    int reportRow = 0;
    for (const QString& line : report.split('\n')) {
        std::istringstream rowStream(line.toStdString());
        char opening = 0, closing = 0;
        std::array<double, 4> displayed{};
        if (rowStream >> opening >> displayed[0] >> displayed[1] >> displayed[2] >> displayed[3] >> closing) {
            if (opening != '[' || closing != ']') continue;
            require(reportRow < 4, "Copyable report has unexpected extra 4x4 matrix rows");
            for (int col = 0; col < 4; ++col)
                requireNear(displayed[col], transform[reportRow][col], "Reported 4x4 matrix differs from scripting output");
            ++reportRow;
        }
    }
    require(reportRow == 4, "Copyable report is missing rows of the full 4x4 transform");
}

MeshModel* makeMesh(MeshDocument& document, const QString& label,
                    double scale, bool floating, bool reorder)
{
    MeshModel* model = document.addNewMesh(QString(), label);
    CMeshO& mesh = model->cm;
    const std::array<int, 5> order = reorder ? std::array<int, 5>{3, 0, 4, 1, 2}
                                           : std::array<int, 5>{0, 1, 2, 3, 4};
    std::array<int, 5> inverse{};
    vcg::tri::Allocator<CMeshO>::AddVertices(mesh, floating ? 8 : 5);
    for (int i = 0; i < 5; ++i) {
        inverse[order[i]] = i;
        const Position& p = points[order[i]];
        mesh.vert[i].P() = Point3m(Scalarm(scale * p[0]), Scalarm(scale * p[1]), Scalarm(scale * p[2]));
    }
    vcg::tri::Allocator<CMeshO>::AddFaces(mesh, int(faces.size()) + (floating ? 1 : 0));
    for (std::size_t f = 0; f < faces.size(); ++f)
        for (int c = 0; c < 3; ++c)
            mesh.face[f].V(c) = &mesh.vert[inverse[faces[f][c]]];
    if (floating) {
        mesh.vert[5].P() = Point3m(Scalarm(100 * scale), Scalarm(200 * scale), Scalarm(300 * scale));
        mesh.vert[6].P() = Point3m(Scalarm(101 * scale), Scalarm(200 * scale), Scalarm(300 * scale));
        mesh.vert[7].P() = Point3m(Scalarm(100 * scale), Scalarm(201 * scale), Scalarm(300 * scale));
        for (int c = 0; c < 3; ++c) mesh.face.back().V(c) = &mesh.vert[5 + c];
    }
    mesh.Tr.SetIdentity();
    vcg::tri::UpdateBounding<CMeshO>::Box(mesh);
    return model;
}

struct Snapshot {
    std::vector<Point3m> positions;
    std::vector<int> vertexFlags;
    std::vector<int> faceFlags;
    std::vector<std::array<std::size_t, 3>> triangles;
    Matrix44m transform;
    int vertexCount;
    int faceCount;
    explicit Snapshot(const MeshModel& model) : transform(model.cm.Tr), vertexCount(model.cm.vn), faceCount(model.cm.fn)
    {
        for (const auto& v : model.cm.vert) { positions.push_back(v.cP()); vertexFlags.push_back(v.cFlags()); }
        for (const auto& f : model.cm.face) {
            faceFlags.push_back(f.cFlags());
            triangles.push_back({std::size_t(f.cV(0) - model.cm.vert.data()),
                                 std::size_t(f.cV(1) - model.cm.vert.data()),
                                 std::size_t(f.cV(2) - model.cm.vert.data())});
        }
    }
    void check(const MeshModel& model) const
    {
        require(model.cm.vn == vertexCount && model.cm.fn == faceCount, "Filter changed mesh counts");
        require(model.cm.vert.size() == positions.size() && model.cm.face.size() == faceFlags.size(), "Filter changed storage sizes");
        for (std::size_t i = 0; i < positions.size(); ++i) {
            require(model.cm.vert[i].cP() == positions[i], "Filter changed vertex coordinates");
            require(model.cm.vert[i].cFlags() == vertexFlags[i], "Filter changed vertex flags");
        }
        for (std::size_t i = 0; i < faceFlags.size(); ++i) {
            require(model.cm.face[i].cFlags() == faceFlags[i], "Filter changed face flags");
            for (int c = 0; c < 3; ++c)
                require(model.cm.face[i].cV(c) == &model.cm.vert[triangles[i][c]], "Filter changed face topology");
        }
        for (int r = 0; r < 4; ++r)
            for (int c = 0; c < 4; ++c)
                require(model.cm.Tr[r][c] == transform[r][c], "Filter changed layer transform");
    }
};

RichParameterList parameters(FilterPlugin& filter, const QAction* action,
                             const MeshDocument& document, const MeshModel& source, const MeshModel& target)
{
    RichParameterList values = filter.initParameterList(action, document);
    values.setValue("SourceMesh", IntValue(source.id()));
    values.setValue("TargetMesh", IntValue(target.id()));
    values.setValue("ShowReport", BoolValue(false));
    return values;
}

std::map<std::string, QVariant> run(FilterPlugin& filter, const QAction* action,
                                   MeshDocument& document, const RichParameterList& values)
{
    unsigned int post = MeshModel::MM_ALL;
    filter.setLog(&document.Log);
    const auto result = filter.applyFilter(action, values, document, post, nullptr);
    filter.setLog(nullptr);
    require(post == MeshModel::MM_NONE, "Read-only filter returned modifying postcondition");
    return result;
}

void mustReject(FilterPlugin& filter, const QAction* action, MeshDocument& document,
                const RichParameterList& values, const char* message)
{
    bool rejected = false;
    try { run(filter, action, document, values); }
    catch (const MLException&) { rejected = true; }
    filter.setLog(nullptr);
    require(rejected, message);
}
}

int main(int argc, char** argv)
{
    QCoreApplication application(argc, argv);
    try {
        require(argc == 2, "Usage: plugin_probe <absolute plugin DLL path>");
        QPluginLoader loader(QString::fromLocal8Bit(argv[1]));
        QObject* object = loader.instance();
        if (!object) throw std::runtime_error(loader.errorString().toStdString());
        auto* filter = qobject_cast<FilterPlugin*>(object);
        require(filter != nullptr, "DLL does not expose the MeshLab FilterPlugin interface");
        require(filter->actions().size() == 1, "Expected one scale measurement filter");
        QAction* action = filter->actions().front();
        require(filter->filterArity(action) == FilterPlugin::FIXED, "Filter must accept two explicit meshes");
        require(filter->postCondition(action) == MeshModel::MM_NONE, "Filter declares a modifying postcondition");
        std::cout << "Loaded " << filter->pluginName().toStdString() << " against MeshLab "
                  << filter->getMLVersion().first << '\n';

        {
            MeshDocument document;
            MeshModel* source = makeMesh(document, "Unscaled with floating piece", 1, true, false);
            MeshModel* target = makeMesh(document, "Scaled reordered without floating piece", 2.75, false, true);
            Snapshot beforeSource(*source), beforeTarget(*target);
            auto values = parameters(*filter, action, document, *source, *target);
            const auto result = run(*filter, action, document, values);
            requireNear(result.at("scale_source_to_target").toDouble(), 2.75, "Incorrect source-to-target scale");
            requireNear(result.at("scale_target_to_source").toDouble(), 1.0 / 2.75, "Incorrect reciprocal scale");
            requireNear(result.at("relative_max_error").toDouble(), 0, "Nonzero exact-match residual");
            requireRotation(result, Rotation{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}});
            require(result.at("source_vertices_ignored").toULongLong() == 3, "Floating piece was not ignored");
            require(result.at("source_vertices_used").toULongLong() == 5, "Wrong retained source count");
            require(!result.at("report").toString().isEmpty(), "Missing copyable result report");
            requireTransformMaps(result, *source, *target, true);
            beforeSource.check(*source);
            beforeTarget.check(*target);
            std::cout << "PASS: reordered vertices, missing floating triangle, factor 2.75, unchanged meshes\n";

            values.setValue("LargestComponent", BoolValue(false));
            mustReject(*filter, action, document, values, "All-geometry mode accepted unequal point sets");
            values.setValue("LargestComponent", BoolValue(true));
            values.setValue("TargetMesh", IntValue(source->id()));
            mustReject(*filter, action, document, values, "Accepted the same layer twice");
            values.setValue("TargetMesh", IntValue(target->id()));
            target->cm.vert[0].P()[0] += Scalarm(0.5);
            mustReject(*filter, action, document, values, "Accepted incompatible geometry");
            std::cout << "PASS: unequal components, same-layer input and incompatible geometry rejected\n";
        }
        {
            MeshDocument document;
            MeshModel* source = makeMesh(document, "Transformed source", 1, false, false);
            MeshModel* target = makeMesh(document, "Transformed target", 1, false, true);
            source->cm.Tr[0][3] = -3;
            source->cm.Tr[1][3] = 7;
            source->cm.Tr[2][3] = 8;
            for (int c = 0; c < 3; ++c) target->cm.Tr[c][c] = Scalarm(2.75);
            target->cm.Tr[0][3] = 5;
            target->cm.Tr[1][3] = -2;
            target->cm.Tr[2][3] = 7;
            Snapshot beforeSource(*source), beforeTarget(*target);
            auto values = parameters(*filter, action, document, *source, *target);
            const auto world = run(*filter, action, document, values);
            requireNear(world.at("scale_source_to_target").toDouble(), 2.75, "Layer transforms were ignored");
            requireNear(world.at("translation_x").toDouble(), 13.25, "Wrong world translation x");
            requireNear(world.at("translation_y").toDouble(), -21.25, "Wrong world translation y");
            requireNear(world.at("translation_z").toDouble(), -15, "Wrong world translation z");
            requireTransformMaps(world, *source, *target, true);
            values.setValue("UseLayerTransforms", BoolValue(false));
            const auto local = run(*filter, action, document, values);
            requireNear(local.at("scale_source_to_target").toDouble(), 1, "Raw-coordinate option ignored");
            requireTransformMaps(local, *source, *target, false);
            beforeSource.check(*source);
            beforeTarget.check(*target);
            std::cout << "PASS: layer transforms, raw-coordinate option, translations, unchanged matrices\n";
        }
        {
            MeshDocument document;
            MeshModel* source = makeMesh(document, "Unscaled with floating piece", 1, true, false);
            MeshModel* target = makeMesh(document, "Scaled, rotated and reordered", 2.75, false, true);
            const Rotation rotation = generalRotation();
            const Position offset{4, -6, 9};
            for (auto& vertex : target->cm.vert) {
                const Position p{double(vertex.P()[0]), double(vertex.P()[1]), double(vertex.P()[2])};
                const Position rotated = rotate(rotation, p);
                for (int axis = 0; axis < 3; ++axis) vertex.P()[axis] = Scalarm(rotated[axis] + offset[axis]);
            }
            vcg::tri::UpdateBounding<CMeshO>::Box(target->cm);
            Snapshot beforeSource(*source), beforeTarget(*target);
            auto values = parameters(*filter, action, document, *source, *target);
            require(values.hasParameter("AllowRotation") && !values.getBool("AllowRotation"),
                    "Rotation estimation must default to disabled");
            mustReject(*filter, action, document, values, "Rotation-disabled mode accepted rotated geometry");
            beforeSource.check(*source);
            beforeTarget.check(*target);
            values.setValue("AllowRotation", BoolValue(true));
            const auto result = run(*filter, action, document, values);
            requireFloatNear(result.at("scale_source_to_target").toDouble(), 2.75, "Incorrect rotated source-to-target scale");
            requireFloatNear(result.at("scale_target_to_source").toDouble(), 1.0 / 2.75, "Incorrect rotated reciprocal scale");
            requireFloatNear(result.at("translation_x").toDouble(), offset[0], "Incorrect rotated translation x");
            requireFloatNear(result.at("translation_y").toDouble(), offset[1], "Incorrect rotated translation y");
            requireFloatNear(result.at("translation_z").toDouble(), offset[2], "Incorrect rotated translation z");
            requireRotation(result, rotation);
            require(result.at("relative_max_error").toDouble() < 1e-6, "Rotated point-set residual too large");
            require(result.at("source_vertices_ignored").toULongLong() == 3, "Rotated comparison included floating piece");
            requireTransformMaps(result, *source, *target, true, true);
            beforeSource.check(*source);
            beforeTarget.check(*target);
            std::cout << "PASS: optional arbitrary 3D rotation, factor 2.75, translation, reordered vertices, missing piece\n";
        }
        {
            MeshDocument document;
            MeshModel* source = makeMesh(document, "Source with translated layer", 1, true, false);
            MeshModel* target = makeMesh(document, "Target with rotated layer", 1, false, true);
            const Rotation rotation = generalRotation();
            const Position sourceOffset{-3, 7, 8}, targetOffset{4, -6, 9};
            for (int row = 0; row < 3; ++row) {
                source->cm.Tr[row][3] = Scalarm(sourceOffset[row]);
                target->cm.Tr[row][3] = Scalarm(targetOffset[row]);
                for (int col = 0; col < 3; ++col)
                    target->cm.Tr[row][col] = Scalarm(2.75 * rotation[row][col]);
            }
            Snapshot beforeSource(*source), beforeTarget(*target);
            auto values = parameters(*filter, action, document, *source, *target);
            mustReject(*filter, action, document, values, "Rotation-disabled mode accepted rotated layer transforms");
            values.setValue("AllowRotation", BoolValue(true));
            const auto world = run(*filter, action, document, values);
            requireFloatNear(world.at("scale_source_to_target").toDouble(), 2.75, "Rotated layer scale was ignored");
            requireRotation(world, rotation);
            const Position rotatedSourceOffset = rotate(rotation, sourceOffset);
            requireFloatNear(world.at("translation_x").toDouble(), targetOffset[0] - 2.75 * rotatedSourceOffset[0], "Wrong rotated layer translation x");
            requireFloatNear(world.at("translation_y").toDouble(), targetOffset[1] - 2.75 * rotatedSourceOffset[1], "Wrong rotated layer translation y");
            requireFloatNear(world.at("translation_z").toDouble(), targetOffset[2] - 2.75 * rotatedSourceOffset[2], "Wrong rotated layer translation z");
            requireTransformMaps(world, *source, *target, true, true);
            values.setValue("UseLayerTransforms", BoolValue(false));
            const auto localWithRotation = run(*filter, action, document, values);
            requireFloatNear(localWithRotation.at("scale_source_to_target").toDouble(), 1, "Rotation mode ignored raw-coordinate option");
            requireTransformMaps(localWithRotation, *source, *target, false, true);
            values.setValue("AllowRotation", BoolValue(false));
            const auto localWithoutRotation = run(*filter, action, document, values);
            requireNear(localWithoutRotation.at("scale_source_to_target").toDouble(), 1, "Raw coordinates incorrectly require rotation estimation");
            requireTransformMaps(localWithoutRotation, *source, *target, false);
            beforeSource.check(*source);
            beforeTarget.check(*target);
            std::cout << "PASS: rotated live layer matrices, raw coordinates in both modes, unchanged geometry and matrices\n";
        }
        std::cout << "PASS: 4x4 transforms map source vertices to target in all coordinate and rotation modes\n";
        std::cout << "All plugin integration checks passed.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
