// SPDX-License-Identifier: GPL-2.0-or-later
#include "filter_scale_factor.h"
#include "../core/scale_estimator.h"

#include <QApplication>
#include <QClipboard>
#include <QDialog>
#include <QDialogButtonBox>
#include <QPlainTextEdit>
#include <QPushButton>
#include <QThread>
#include <QVBoxLayout>

#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

// Do not create windows in PyMeshLab, meshlabserver or other library clients.
QWidget* meshLabWindow()
{
    auto* app = qobject_cast<QApplication*>(QCoreApplication::instance());
    if (!app || QThread::currentThread() != app->thread() ||
        !QCoreApplication::applicationName().contains("MeshLab", Qt::CaseInsensitive))
        return nullptr;
    for (QWidget* widget : QApplication::topLevelWidgets()) {
        if (widget->isVisible() && widget->inherits("MainWindow"))
            return widget;
    }
    return nullptr;
}

void showReport(const QString& report, QWidget* parent)
{
    auto* dialog = new QDialog(parent);
    dialog->setAttribute(Qt::WA_DeleteOnClose);
    dialog->setWindowTitle("Scale between meshes");
    dialog->resize(760, 520);
    auto* layout = new QVBoxLayout(dialog);
    auto* text = new QPlainTextEdit(report, dialog);
    text->setReadOnly(true);
    text->setLineWrapMode(QPlainTextEdit::WidgetWidth);
    layout->addWidget(text);
    auto* buttons = new QDialogButtonBox(QDialogButtonBox::Close, dialog);
    auto* copy = buttons->addButton("Copy report", QDialogButtonBox::ActionRole);
    QObject::connect(copy, &QPushButton::clicked, dialog, [report]() {
        QApplication::clipboard()->setText(report);
    });
    QObject::connect(buttons, &QDialogButtonBox::rejected, dialog, &QDialog::close);
    layout->addWidget(buttons);
    // Modeless: no nested event loop and no blocked filter completion.
    dialog->show();
    dialog->raise();
}

meshlab_scale::Mesh copyMesh(const MeshModel& model, bool useTransform)
{
    const CMeshO& mesh = model.cm;
    const std::size_t invalid = std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> indices(mesh.vert.size(), invalid);
    meshlab_scale::Mesh out;
    out.vertices.reserve(mesh.vert.size());
    out.faces.reserve(mesh.face.size());

    if (useTransform) {
        for (int row = 0; row < 4; ++row)
            for (int col = 0; col < 4; ++col)
                if (!std::isfinite(static_cast<double>(mesh.Tr[row][col])))
                    throw std::runtime_error("Layer transform contains a nonfinite value.");
        if (mesh.Tr[3][0] != 0 || mesh.Tr[3][1] != 0 ||
            mesh.Tr[3][2] != 0 || mesh.Tr[3][3] != 1)
            throw std::runtime_error("Only affine layer transforms are supported.");
    }

    for (std::size_t i = 0; i < mesh.vert.size(); ++i) {
        const auto& vertex = mesh.vert[i];
        if (vertex.IsD())
            continue;
        const auto& p = vertex.cP();
        meshlab_scale::Vec3 position{double(p[0]), double(p[1]), double(p[2])};
        if (useTransform) {
            const double xyz[3] = {position.x, position.y, position.z};
            double transformed[3];
            for (int row = 0; row < 3; ++row) {
                transformed[row] = double(mesh.Tr[row][3]);
                for (int col = 0; col < 3; ++col)
                    transformed[row] += double(mesh.Tr[row][col]) * xyz[col];
            }
            position = {transformed[0], transformed[1], transformed[2]};
        }
        indices[i] = out.vertices.size();
        out.vertices.push_back(position);
    }

    for (const auto& face : mesh.face) {
        if (face.IsD())
            continue;
        std::array<std::size_t, 3> triangle;
        for (int corner = 0; corner < 3; ++corner) {
            const auto* vertex = face.cV(corner);
            if (!vertex)
                throw std::runtime_error("A face contains a null vertex.");
            const auto index = static_cast<std::size_t>(vertex - mesh.vert.data());
            if (index >= indices.size() || indices[index] == invalid)
                throw std::runtime_error("A face references an unavailable vertex.");
            triangle[corner] = indices[index];
        }
        out.faces.push_back(triangle);
    }
    return out;
}

QString number(double value)
{
    return QString::number(value, 'g', 17);
}

QString count(std::size_t value)
{
    return QString::number(static_cast<qulonglong>(value));
}

} // namespace

FilterScaleFactorPlugin::FilterScaleFactorPlugin()
{
    typeList = {MEASURE_SCALE_FACTOR};
    actionList.push_back(new QAction(filterName(MEASURE_SCALE_FACTOR), this));
}

QString FilterScaleFactorPlugin::pluginName() const { return "FilterScaleFactor"; }
QString FilterScaleFactorPlugin::vendor() const { return "Scale Factor Plugin"; }
QString FilterScaleFactorPlugin::filterName(ActionIDType) const
{
    return "Measure Scale Between Two Meshes";
}
QString FilterScaleFactorPlugin::pythonFilterName(ActionIDType) const
{
    return "get_scale_between_two_meshes";
}
QString FilterScaleFactorPlugin::filterInfo(ActionIDType) const
{
    return "Measure the positive uniform scale from the unscaled source to the scaled target. "
           "The meshes must have the same vertices on the retained piece; vertex order may differ. "
           "Enable Allow rotation between models if their orientations differ. "
           "By default only the largest connected piece is compared, "
           "so missing smaller floating pieces do not affect the result. Translation, including "
           "a different scaling pivot, is fitted automatically. The result is validated in 3D "
           "and includes residual error. Neither mesh is modified.";
}
FilterPlugin::FilterClass FilterScaleFactorPlugin::getClass(const QAction*) const
{
    return Measure;
}
FilterPlugin::FilterArity FilterScaleFactorPlugin::filterArity(const QAction*) const
{
    return FIXED;
}
int FilterScaleFactorPlugin::getPreConditions(const QAction*) const { return MeshModel::MM_NONE; }
int FilterScaleFactorPlugin::postCondition(const QAction*) const { return MeshModel::MM_NONE; }

RichParameterList FilterScaleFactorPlugin::initParameterList(const QAction*, const MeshDocument& md)
{
    int sourceId = md.mm() ? md.mm()->id() : -1;
    int targetId = sourceId;
    const auto begin = md.meshBegin();
    if (begin != md.meshEnd()) {
        sourceId = begin->id();
        if (targetId == sourceId) {
            auto next = begin;
            ++next;
            if (next != md.meshEnd())
                targetId = next->id();
        }
    }
    RichParameterList params;
    params.addParam(RichMesh("SourceMesh", sourceId, &md, "Source (unscaled)",
        "The factor multiplies this model's dimensions to obtain the target's dimensions."));
    params.addParam(RichMesh("TargetMesh", targetId, &md, "Target (scaled)",
        "The model whose scale is measured relative to the source."));
    params.addParam(RichBool("LargestComponent", true, "Compare largest connected piece only",
        "Ignore smaller disconnected pieces and unused vertices. Exact coincident vertices "
        "are welded internally, including STL triangle vertices. The retained pieces must "
        "contain the same unique vertices. For point clouds all unique points are compared."));
    params.addParam(RichBool("UseLayerTransforms", true, "Include layer transforms",
        "Measure the displayed coordinates, including each layer's transformation matrix. "
        "Disable to compare the stored vertex coordinates alone."));
    params.addParam(RichBool("AllowRotation", false, "Allow rotation between models",
        "Enable when one model has also been rotated. Estimate rotation, uniform scale and "
        "translation, then verify the match in 3D. The retained pieces must still contain "
        "the same unique vertices. This option can take longer; neither model is changed."));
    params.addParam(RichFloat("RelativeTolerance", 1e-5, "Maximum relative geometric error",
        "Reject the result if the maximum bidirectional nearest-vertex distance exceeds this "
        "fraction of the retained target's bounding-box diagonal. This is a geometric check, "
        "not an uncertainty bound on the estimated scale."));
    params.addParam(RichBool("ShowReport", meshLabWindow() != nullptr, "Show copyable report",
        "Open a report window in the MeshLab GUI. Results are always written to the log and "
        "returned to scripting clients. No window is opened by headless clients."));
    return params;
}

std::map<std::string, QVariant> FilterScaleFactorPlugin::applyFilter(
    const QAction* action, const RichParameterList& params, MeshDocument& md,
    unsigned int& postConditionMask, vcg::CallBackPos* cb)
{
    if (ID(action) != MEASURE_SCALE_FACTOR)
        wrongActionCalled(action);
    postConditionMask = MeshModel::MM_NONE;
    const auto* source = md.getMesh(params.getMeshId("SourceMesh"));
    const auto* target = md.getMesh(params.getMeshId("TargetMesh"));
    if (!source || !target || source == target)
        throw MLException("Load two mesh layers and choose different source and target meshes.");

    try {
        if (cb && !(*cb)(5, "Reading mesh coordinates"))
            throw std::runtime_error("Scale measurement cancelled.");
        const bool transforms = params.getBool("UseLayerTransforms");
        const auto sourceMesh = copyMesh(*source, transforms);
        const auto targetMesh = copyMesh(*target, transforms);
        meshlab_scale::Options options;
        options.largestComponent = params.getBool("LargestComponent");
        options.relativeTolerance = params.getFloat("RelativeTolerance");
        // Preserve the old default for saved filter scripts without this parameter.
        options.allowRotation = params.hasParameter("AllowRotation") && params.getBool("AllowRotation");
        if (cb && !(*cb)(25, "Measuring and validating scale"))
            throw std::runtime_error("Scale measurement cancelled.");
        const auto result = meshlab_scale::estimateScale(sourceMesh, targetMesh, options);
        // Column-vector convention: [target; 1] = transform * [source; 1].
        // This maps the coordinate space used by copyMesh, including layer
        // transforms only when the corresponding parameter is enabled.
        std::array<std::array<double, 4>, 4> transform{};
        const double translation[3] = {
            result.translation.x, result.translation.y, result.translation.z};
        for (int row = 0; row < 3; ++row) {
            for (int col = 0; col < 3; ++col)
                transform[row][col] = result.scale * result.rotation[row][col];
            transform[row][3] = translation[row];
        }
        transform[3][3] = 1.0;
        for (const auto& row : transform)
            for (double value : row)
                if (!std::isfinite(value))
                    throw std::runtime_error("The composed source-to-target transform contains an unrepresentable value.");

        QStringList lines;
        lines << "Scale between meshes"
              << "Source (unscaled): " + source->label()
              << "Target (scaled): " + target->label()
              << ""
              << "Multiply SOURCE dimensions by: " + number(result.scale)
              << "Multiply TARGET dimensions by: " + number(result.inverseScale)
              << ""
              << (options.allowRotation
                    ? "Fitted relation: target = scale * rotation * source + translation"
                    : "Fitted relation: target = scale * source + translation")
              << "Translation: (" + number(result.translation.x) + ", " +
                    number(result.translation.y) + ", " + number(result.translation.z) + ")"
              << ""
              << "Source-to-target 4x4 transform (homogeneous matrix, row-major):"
              << "Column vectors: [x_target y_target z_target 1]^T = T * [x_source y_source z_source 1]^T"
              << (transforms
                    ? "Coordinate space: displayed source to displayed target (layer transforms included)."
                    : "Coordinate space: stored source to stored target (layer transforms excluded).");
        for (const auto& row : transform)
            lines << "[ " + number(row[0]) + "  " + number(row[1]) + "  " +
                     number(row[2]) + "  " + number(row[3]) + " ]";
        lines << ""
              << "RMS geometric error (target units): " + number(result.rmsError)
              << "Maximum geometric error (target units): " + number(result.maxError)
              << "Relative RMS error: " + number(result.relativeRmsError)
              << "Relative maximum error: " + number(result.relativeMaxError)
              << "Allowed relative maximum error: " + number(options.relativeTolerance)
              << ""
              << "Source unique vertices used: " + count(result.sourceVerticesUsed) +
                    "; input vertices ignored/welded: " + count(result.sourceVerticesIgnored)
              << "Target unique vertices used: " + count(result.targetVerticesUsed) +
                    "; input vertices ignored/welded: " + count(result.targetVerticesIgnored)
              << QString("Compared: %1; layer transforms: %2")
                     .arg(options.largestComponent ? "largest connected piece (all points for point clouds)" : "all geometry")
                     .arg(transforms ? "included" : "excluded")
              << QString("Rotation fitting: %1").arg(options.allowRotation ? "enabled" : "disabled");
        if (options.allowRotation) {
            lines << "" << "Source-to-target rotation matrix (rows):";
            for (const auto& row : result.rotation)
                lines << "[ " + number(row[0]) + "  " + number(row[1]) + "  " + number(row[2]) + " ]";
            lines << "Symmetric models can have more than one equally valid rotation.";
        }
        lines << "" << "Validation passed for the compared geometry. Displayed digits preserve the computed "
                 "result; input coordinate precision limits the meaningful digits. This is not a "
                 "guarantee of the original pre-rounding scale. Neither mesh was modified.";
        const QString report = lines.join('\n');
        for (const auto& line : lines)
            log(line.toStdString());
        if (cb)
            (*cb)(100, "Scale measurement complete");
        if (params.getBool("ShowReport")) {
            if (auto* parent = meshLabWindow())
                showReport(report, parent);
        }
        std::map<std::string, QVariant> output = {
            {"scale_source_to_target", result.scale},
            {"scale_target_to_source", result.inverseScale},
            {"translation_x", result.translation.x},
            {"translation_y", result.translation.y},
            {"translation_z", result.translation.z},
            {"rms_error", result.rmsError},
            {"max_error", result.maxError},
            {"relative_rms_error", result.relativeRmsError},
            {"relative_max_error", result.relativeMaxError},
            {"source_vertices_used", QVariant::fromValue<qulonglong>(result.sourceVerticesUsed)},
            {"target_vertices_used", QVariant::fromValue<qulonglong>(result.targetVerticesUsed)},
            {"source_vertices_ignored", QVariant::fromValue<qulonglong>(result.sourceVerticesIgnored)},
            {"target_vertices_ignored", QVariant::fromValue<qulonglong>(result.targetVerticesIgnored)},
            {"report", report},
            {"rotation_enabled", options.allowRotation}
        };
        for (int row = 0; row < 3; ++row)
            for (int col = 0; col < 3; ++col)
                output["rotation_" + std::to_string(row) + std::to_string(col)] = result.rotation[row][col];
        for (int row = 0; row < 4; ++row)
            for (int col = 0; col < 4; ++col)
                output["transform_" + std::to_string(row) + std::to_string(col)] = transform[row][col];
        return output;
    } catch (const std::exception& error) {
        throw MLException(QString::fromUtf8(error.what()));
    }
}

MESHLAB_PLUGIN_NAME_EXPORTER(FilterScaleFactorPlugin)
