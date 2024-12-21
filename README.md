# disc-topology-stripes
An algorithm to determine stripe patterns for stitching across meshes with disc topology.


 This project aims to facilitate the generation of knit structures specifically for digital meshes with disc topology. Namely, our app lets
users specify regions of their mesh as their desired start and end points during knitting. Using Geometry Central’s implementation of the approach in Kn¨oppel et al. [2015] to stripe patterns, we create and display two orthogonal stripe patterns that act as the course and wale directions to be followed by the user while knitting. The generated structures can be knitted entirely with short rows (otheriwse known as “partial rows”) with no use of increases or decreases unless necessary

To start the app, run:

gc\_project <path\_to\_mesh\_file>



in your terminal, where gc\_project is the default name of the compiled executable and <path\_to\_mesh\_file> is the file path to the mesh you'd like to work with.

With the app open, the first step is to select two strips of points for your start and end courses as mentioned previously. To do this, click on ``Pick start boundary (h=0)`` and ``Pick end boundary (h=1).`` From there, select two boundary vertices on your mesh to define the endpoints of your strip, then select one last vertex on what you want to be the ``inside'' of your strip. Without this last selection, it would be ambiguous if you want one strip of points or its complement



After this, the ``Generate H Function,`` ``Construct Gradient,`` ``Construct Tangent,`` ``Generate Wale Curves,`` and ``Generate Course Curves`` buttons can all be run in succession.

You can use the two sliders to adjust the frequency of course and wale stripes for your mesh. At any time, you can click ``Clear boundaries'' to clear your selections for the start and end courses and start over.
