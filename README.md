# Disk Topology Stripes App

App for generating stripe patterns for knitting disc topologies. Created by Michelle Yilin Feng and Noah Barnes for CS582 using the Polyscope, Geometry Central, and Eigen libraries. Major thanks to the implementation of Knoppel et al.'s stripe pattern generator from Geometry Central.

# Compiling

Run `cmake` in the root directory of the app. Then, in the `build` directory, run `make -j10` to compile the app.

# Running

To open the app, run `./build/bin/gc_project <mesh_file_path>` from the app's root directory, where `<mesh_file_path>` is the path to the mesh you'd like to work with.

# Usage

## Selecting Boundaries

With your mesh opened in the app, click on "Pick start boundary (h=0)" to select a strip of boundary points to act as your starting boundary (same thing for "Pick end boundary (h=1)"). Once you've done this, select two boundary vertices that you would like as the endpoints of your strip of points. Selecting a vertex requires you to right-click on the vertex in the 3D display. Once you've selected your two endpoints, you must choose one last point to define the "insideness" of your strip. Without this point, it would be ambiguous what strip you want as everything is on a single boundary. For example, if you select opposite corners of a square, it would be ambiguous whether you wanted to select one pair of edges of the square or the other pair.

## Selecting Frequencies

Use the available sliders to pick a frequency you'd like for your stripes. Choose these before running "Generate Wale Curves" and "Generate Course Curves."

## Generating Stripes

Once you've picked your boundaries and frequencies, you may run "Generate H Function," "Construct gradient," "Construct Tangent," "Generate Wale Curves," and "Generate Course Curves" in succession.

If you'd like to adjust the frequency of stripes, change the sliders, then click "Generate Wale/Course Curves" again.

## Resetting

At any time, you can click "Clear boundary points" to reset your boundary selections and start over.
