# SPDX-License-Identifier: GPL-2.0-or-later
# Include this directory from a configured MeshLab 2025.07 source build.
add_meshlab_plugin(filter_scale_factor
    filter_scale_factor.cpp
    filter_scale_factor.h
    ../core/scale_estimator.cpp
    ../core/scale_estimator.h)
target_compile_features(filter_scale_factor PRIVATE cxx_std_17)
target_link_libraries(filter_scale_factor PRIVATE Qt5::Widgets)
