// SPDX-License-Identifier: GPL-2.0-or-later
#ifndef FILTER_SCALE_FACTOR_H
#define FILTER_SCALE_FACTOR_H

#include <common/plugins/interfaces/filter_plugin.h>

class FilterScaleFactorPlugin : public QObject, public FilterPlugin
{
    Q_OBJECT
    MESHLAB_PLUGIN_IID_EXPORTER(FILTER_PLUGIN_IID)
    Q_INTERFACES(FilterPlugin)

public:
    enum { MEASURE_SCALE_FACTOR };
    FilterScaleFactorPlugin();
    QString pluginName() const override;
    QString vendor() const override;
    QString filterName(ActionIDType) const override;
    QString pythonFilterName(ActionIDType) const override;
    QString filterInfo(ActionIDType) const override;
    FilterClass getClass(const QAction*) const override;
    FilterArity filterArity(const QAction*) const override;
    int getPreConditions(const QAction*) const override;
    int postCondition(const QAction*) const override;
    RichParameterList initParameterList(const QAction*, const MeshDocument&) override;
    std::map<std::string, QVariant> applyFilter(
        const QAction*, const RichParameterList&, MeshDocument&,
        unsigned int& postConditionMask, vcg::CallBackPos*) override;
};

#endif
