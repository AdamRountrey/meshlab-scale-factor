# Third-party notices and source availability

MeshLab Scale Factor uses the MeshLab filter interface, VCGLib and its bundled Eigen headers, and Qt. The plugin's original source is licensed under GPL-2.0-or-later. Third-party files retain their own copyright notices and licensing terms; this document does not relicense them.

The release's corresponding-source archive contains the plugin source and build scripts together with the exact MeshLab and VCGLib source revisions listed below, including Eigen and their original notices. The MeshLab and VCGLib tracked sources were unmodified when the release was built. The plugin DLL does not bundle or replace MeshLab or Qt runtime DLLs; those are supplied by the installed MeshLab application.

## MeshLab

- Version: MeshLab 2025.07, commit `dc48b91ae562756a6988048c5d5c7f1d2b687256`.
- Copyright: Visual Computing Lab, ISTI - Italian National Research Council, and the other authors identified in the upstream files.
- [Exact source](https://github.com/cnr-isti-vclab/meshlab/tree/dc48b91ae562756a6988048c5d5c7f1d2b687256).
- The filter-interface header grants GPL version 2 or any later version. The upstream root `LICENSE.txt` supplies GPL version 3. Both license texts are included: [GPL-2.0](third_party_licenses/GPL-2.0.txt), [GPL-3.0](third_party_licenses/GPL-3.0.txt).
- [Preserved interface copyright and license notice](third_party_licenses/MeshLab-NOTICES.txt).

## VCGLib

- Commit: `c94ef4e12e9ea3ae986d9af91005be8328d13719`, the submodule revision pinned by the MeshLab build.
- Copyright: Visual Computing Lab, ISTI - Italian National Research Council, and the other authors identified in the upstream files.
- [Exact source](https://github.com/cnr-isti-vclab/vcglib/tree/c94ef4e12e9ea3ae986d9af91005be8328d13719).
- The relevant template headers grant GPL version 2 or any later version; the upstream root license text is GPL version 3. The included GPL texts above apply according to each source file's grant.
- [Preserved representative header notice](third_party_licenses/VCGLib-NOTICES.txt). The complete source archive retains the notices in every upstream file.

## Eigen

- Version: 3.4.0, as identified by `Eigen/src/Core/util/Macros.h`, bundled in the pinned VCGLib source.
- Copyright: Gael Guennebaud, Benoit Jacob, and the other contributors identified in individual Eigen files.
- [Exact bundled source](https://github.com/cnr-isti-vclab/vcglib/tree/c94ef4e12e9ea3ae986d9af91005be8328d13719/eigenlib).
- The included core headers use the Mozilla Public License, version 2.0. [MPL-2.0 text](third_party_licenses/MPL-2.0.txt) and [preserved Eigen notices](third_party_licenses/Eigen-NOTICES.txt) are included. Individual files retain any additional applicable notices.

## Qt

- Version: Qt 5.15.2, official MSVC 2019 x64 development package; the plugin uses Qt Core, GUI, and Widgets interfaces.
- Copyright: The Qt Company Ltd., Intel Corporation, and the other contributors identified in the Qt sources.
- [Matching QtBase source](https://github.com/qt/qtbase/tree/v5.15.2) and [official Qt 5.15.2 source downloads](https://download.qt.io/archive/qt/5.15/5.15.2/single/).
- The relevant Qt headers offer commercial licensing, LGPL version 3, or GPL version 2 / GPL version 3 under their stated terms. This open-source plugin distribution uses the **GPL option**, with GPL version 3 available for the combined binary, consistently with the plugin source's GPL-2.0-or-later grant. It does not claim a commercial Qt license or rely on LGPL terms for distribution.
- [GPL-2.0](third_party_licenses/GPL-2.0.txt), [GPL-3.0](third_party_licenses/GPL-3.0.txt), and [preserved Qt header notices](third_party_licenses/Qt-NOTICES.txt) are included. Qt runtime libraries remain external and dynamically linked; no Qt runtime binary is included in the plugin release.

## Build and license-text provenance

The Windows binary uses Visual Studio 2022 MSVC 19.44, x64 Release, C++17, and MeshLab's standard single-precision mesh representation. The estimator computes in double precision. See [README.md](README.md) for rebuilding and [third_party_licenses/README.md](third_party_licenses/README.md) for the provenance of the copied license texts. The license copies and notices do not replace the corresponding source made available with the release.
