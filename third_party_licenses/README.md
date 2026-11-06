# License-text provenance

These files preserve the relevant license texts and selected original notices for the plugin's dependencies. Every original notice is also retained in the matching upstream source supplied with the release.

- `GPL-2.0.txt`: copied unchanged from this project's `LICENSE`, the GNU GPL version 2 text. Used with the GPL grants identified in `THIRD_PARTY_NOTICES.md`.
- `GPL-3.0.txt`: copied unchanged from `LICENSE.txt` in MeshLab commit `dc48b91ae562756a6988048c5d5c7f1d2b687256`; the pinned VCGLib source supplies the same license version.
- `MPL-2.0.txt`: the complete standard Mozilla Public License version 2.0 text, copied unchanged from MeshLab's `src/meshlabplugins/filter_cubization/CubicStylizationFiles/src/LICENSE.MPL2` at the same pinned commit. This is included for Eigen's MPL-2.0-covered headers, not as a claim that the cubization plugin is part of this binary.
- `MeshLab-NOTICES.txt`: original initial notice from the pinned MeshLab `src/common/plugins/interfaces/filter_plugin.h`.
- `VCGLib-NOTICES.txt`: original initial notice from the pinned VCGLib `vcg/space/point3.h`.
- `Eigen-NOTICES.txt`: original notices from the pinned VCGLib `eigenlib/Eigen/Core` and `eigenlib/Eigen/src/Core/util/Macros.h`.
- `Qt-NOTICES.txt`: original initial notices from Qt 5.15.2 `QtCore/qglobal.h` and `QtWidgets/qdialog.h`. The plugin distribution selects their GPL licensing option; no LGPL license-text download is needed for that choice.

The original source headers are authoritative for each file's copyright and licensing grant. These acknowledgments do not change those grants or restrict the recipients' rights under them.
