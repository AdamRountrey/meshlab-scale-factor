Import unscaled_with_floating_piece.ply and one of the target PLY files as
separate layers in MeshLab.

Run Measure Scale Between Two Meshes. Choose:
  Source mesh: unscaled_with_floating_piece.ply
  Target mesh: scaled_2_75_without_floating_piece.ply
  Largest connected component: enabled
  Allow rotation: disabled

Expected source-to-target factor: 2.75
Expected target-to-source factor: 0.36363636363636365
Expected ignored source vertices: 3
Expected maximum geometric residual: zero, apart from floating-point roundoff

The target's vertex order differs from the source. Its floating triangle was
removed. The large distance to the floating triangle is intentional: a simple
whole-mesh bounding-box comparison would give the wrong answer.

To test rotation, instead choose target:
  scaled_2_75_rotated_without_floating_piece.ply
  Allow rotation: enabled

This target was scaled by 2.75, rotated around all three axes, and translated
by (4, -6, 9). The rotation is Rz(23 degrees) * Ry(-31 degrees) * Rx(17 degrees)
using column vectors. Its vertex order differs from the source and the floating
triangle is absent. The result should still report a factor close to 2.75;
small residuals are expected because MeshLab stores these coordinates as floats.
With Allow rotation disabled, this rotated comparison must be rejected.
