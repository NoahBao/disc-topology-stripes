#include <geometrycentral/surface/manifold_surface_mesh.h>
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/direction_fields.h"

// STRIPES
#include "geometrycentral/surface/stripe_patterns.h"
#include <polyscope/polyscope.h>
#include <polyscope/point_cloud.h>
#include <queue>
#include "polyscope/curve_network.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"


// #include "include/main.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;
using T = Eigen::Triplet<double>;

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh *psMesh;

// =============FUNCTION HEADERS=============

void selectBoundary(int boundaryIndex);
void selectBoundaryPoints(int endPointIndex1, int endPointIndex2, int innerPointIndex, int boundaryId);

void constructHFunction(std::vector<int> &startVertices, std::vector<int> &endVertices);

void constructFaceGradient();
void constructVertexGradient();

void constructTangentVectors();
Vector2 projectOntoTangentPlane(Vector3 vec, Vector3 norm, Vector3 he);

// ==========================================

// =============GLOBAL VARIABLES=============

// boundaryVertices[i] = list of vertices for boundary i
// analogous for boundaryPoints
std::vector<std::vector<int>> boundaryVertices = {{}, {}};
std::vector<std::vector<Vector3>> boundaryPoints = {{}, {}};

VertexData<double> hFunction;
FaceData<Vector3> faceGradient;
VertexData<Vector3> vertexGradient;
VertexData<Vector2> courseTangentVectors;
VertexData<Vector2> waleTangentVectors;


float waleFrequency = 20.0f;
float courseFrequency = 20.0f;


bool waleStripesExists = false;
bool courseStripesExists = false;


// ===========================================================

// generates stripes using Geometry Central's stripe methods
void generateWaleStripes()
{
  // remove existing  wales curve
    if (waleStripesExists) {
        polyscope::removePointCloud("Wale Curve Points");
        polyscope::removeCurveNetwork("Wale Curve Edges");
        waleStripesExists = false;
    }

  // Generate a guiding field
  VertexData<Vector2> waleGuideField =
      geometrycentral::surface::computeSmoothestVertexDirectionField(*geometry, 2);

  // use tangent vectors instead of smoothest direction field
  waleGuideField = waleTangentVectors;

  // Compute the stripe couse pattern
  // double constantFreq = 40; // for square
  // double constantFreq = 30.; // for c shape
  // double constantFreq = 10.; // for s shape
  // double constantFreq = 3.0; // for trapozoid
  // double constantFreq = 20; // for trapozoid
  VertexData<double> frequencies(*mesh, waleFrequency);
  CornerData<double> periodicFunc;
  FaceData<int> zeroIndices;
  FaceData<int> branchIndices;
  std::tie(periodicFunc, zeroIndices, branchIndices) =
      computeStripePattern(*geometry, frequencies, waleGuideField);

  // Extract isolines
  std::vector<Vector3> isolineVerts;
  std::vector<std::array<size_t, 2>> isolineEdges;
  std::tie(isolineVerts, isolineEdges) = extractPolylinesFromStripePattern(
      *geometry, periodicFunc, zeroIndices, branchIndices, waleGuideField, false);
  polyscope::registerPointCloud("Wale Curve points", isolineVerts);
  polyscope::registerCurveNetwork("Wale Curve edges", isolineVerts, isolineEdges);
}
void generateCourseStripes()
{
    if (waleStripesExists) {
        polyscope::removePointCloud("Course Curve Points");
        polyscope::removeCurveNetwork("Course Curve Edges");
        courseStripesExists = false;
    }

    VertexData<Vector2> courseGuideField = courseTangentVectors;

    // Compute stripe pattern based on the courses guide field
    // double constantFreq = 40; // for square
    // double constantFreq = 20.;   // for c shape
    // double constantFreq = 5.; // for s shape
    // double constantFreq = 3.0; // for trapezoid shape
    // double constantFreq = 20;
    VertexData<double> frequencies(*mesh, courseFrequency);
    CornerData<double> periodicFunc;
    FaceData<int> zeroIndices;
    FaceData<int> branchIndices;
    std::tie(periodicFunc, zeroIndices, branchIndices) =
        computeStripePattern(*geometry, frequencies, courseGuideField);

    // Extract isolines
    std::vector<Vector3> isolineVerts;
    std::vector<std::array<size_t, 2>> isolineEdges;
    std::tie(isolineVerts, isolineEdges) = extractPolylinesFromStripePattern(
      *geometry, periodicFunc, zeroIndices, branchIndices, courseGuideField, false);

    polyscope::registerPointCloud("Course Curve Points", isolineVerts);
    polyscope::registerCurveNetwork("Course Curve Edges", isolineVerts, isolineEdges);

    courseStripesExists=true;
}

// select a particular boundary (start or end)
void selectBoundary(int boundaryIndex)
{
  // initialize variables
  geometry->requireVertexPositions();
  glm::vec3 pointCloudColor = boundaryIndex == 0 ? glm::vec3({1., 0.1, 1.}) : glm::vec3({1., 1., 0.1});
  std::string startOrEnd = boundaryIndex == 0 ? "Start" : "End";
  std::string pointCloudName = startOrEnd + " points (h=" + std::to_string(boundaryIndex) + ")";
  auto pc = polyscope::registerPointCloud(pointCloudName, boundaryPoints[boundaryIndex]);
  polyscope::warning("Pick two endpoints first, then choose a point to define the \"insideness\" of your selection.");

  // pick endpoint 1 and add preview for it
  int endPoint1 = psMesh->selectVertex();
  boundaryPoints[boundaryIndex].push_back(geometry->vertexPositions[endPoint1]);
  pc = polyscope::registerPointCloud(pointCloudName, boundaryPoints[boundaryIndex]);

  // same thing for endpoint 2
  int endPoint2 = psMesh->selectVertex();
  boundaryPoints[boundaryIndex].push_back(geometry->vertexPositions[endPoint2]);
  pc = polyscope::registerPointCloud(pointCloudName, boundaryPoints[boundaryIndex]);

  // pick inner point and begin search
  int innerPoint = psMesh->selectVertex();
  boundaryVertices[boundaryIndex] = {};
  boundaryPoints[boundaryIndex] = {};
  // points must be on boundary
  if (!mesh->vertex(endPoint1).isBoundary() || !mesh->vertex(endPoint2).isBoundary() || !mesh->vertex(innerPoint).isBoundary())
  {
    polyscope::warning("Please select boundary vertices only.");
  }
  else
  {
    selectBoundaryPoints(endPoint1, endPoint2, innerPoint, boundaryIndex);
    pc = polyscope::registerPointCloud(pointCloudName, boundaryPoints[boundaryIndex]);
    pc->setPointColor(pointCloudColor);
  }

  waleStripesExists = true;
}

// selects all boundary vertices that are between endpoints 1 and 2 ("in between" defined by inner point)
void selectBoundaryPoints(int endPointIndex1, int endPointIndex2, int innerPointIndex, int boundaryId)
{
  // intitialize BFS
  std::map<int, bool> isVertexVisited;
  for (Vertex vertex : mesh->vertices())
  {
    isVertexVisited[vertex.getIndex()] = false;
  }
  std::queue<int> vertices;
  vertices.push(innerPointIndex);
  isVertexVisited[innerPointIndex] = true;
  // run BFS from inner point, stop at endpoints 1 and 2
  while (!vertices.empty())
  {
    Vertex v_i = mesh->vertex(vertices.front());
    vertices.pop();

    boundaryVertices[boundaryId].push_back(v_i.getIndex());
    boundaryPoints[boundaryId].push_back(geometry->vertexPositions[v_i]);

    if (v_i.getIndex() == endPointIndex1 || v_i.getIndex() == endPointIndex2)
    {
      // stop exploring once we reach endpoints
      continue;
    }

    for (Halfedge he : v_i.incomingHalfedges())
    {
      Vertex v_j = he.tailVertex();
      if (v_j.isBoundary() && !isVertexVisited[v_j.getIndex()] && he.edge().isBoundary())
      {
        isVertexVisited[v_j.getIndex()] = true;
        vertices.push(v_j.getIndex());
      }
    }
  }
}

// uses Laplacian smoothing to construct H function based on start and stop boundaries
void constructHFunction(std::vector<int> &startVertices, std::vector<int> &endVertices)
{
  bool debug = false;
  bool useCotanWeights = true;
  bool includeDMatrix = true;
  bool doSquareLaplace = true;
  geometry->requireVertexPositions();
  geometry->requireEdgeCotanWeights();
  geometry->requireVertexDualAreas();
  geometry->requireVertexIndices();

  // =====================CONSTRUCT LAPLACE MATRIX=====================

  // construct M and D
  SparseMatrix<double> M(mesh->nVertices(), mesh->nVertices());
  SparseMatrix<double> D(mesh->nVertices(), mesh->nVertices());
  std::vector<T> m_triplets = {};
  std::vector<T> d_triplets = {};
  for (Vertex v_i : mesh->vertices())
  {
    int i = v_i.getIndex();
    double weight_sum = 0.;
    for (Halfedge edge : v_i.incomingHalfedges())
    {
      Vertex v_j = edge.tailVertex();
      int j = v_j.getIndex();
      // cotan weight : uniform
      double weight = useCotanWeights ? geometry->edgeCotanWeights[edge.edge()] : 1;
      weight_sum -= weight;
      // for all i != j
      m_triplets.push_back(T({i, j, weight}));
    }
    // for all i = j
    m_triplets.push_back(T({i, i, weight_sum}));
    double dualArea = geometry->vertexDualAreas[v_i];
    d_triplets.push_back(T({i, i, 1. / (2. * dualArea)}));
  }
  // construct L^2 = (DM)^2
  M.setFromTriplets(m_triplets.begin(), m_triplets.end());
  D.setFromTriplets(d_triplets.begin(), d_triplets.end());
  SparseMatrix<double> laplace_matrix = includeDMatrix ? D * M : M;
  SparseMatrix<double> A_mat = doSquareLaplace ? laplace_matrix * laplace_matrix : laplace_matrix;

  // ==============UPDATE MATRIX ROWS FOR FIXED VERTICES==============

  // update rows for all fixed-value vertices
  // first, transpose once to edit rows as columns
  A_mat = A_mat.transpose();
  for (int i : startVertices)
  {
    SparseMatrix<double> constrained_col = SparseMatrix<double>(mesh->nVertices(), 1);
    std::vector<T> triplets = {T({i, 0, 1})};
    constrained_col.setFromTriplets(triplets.begin(), triplets.end());
    A_mat.col(i) = constrained_col;
  }
  for (int i : endVertices)
  {
    SparseMatrix<double> constrained_col = SparseMatrix<double>(mesh->nVertices(), 1);
    std::vector<T> triplets = {T({i, 0, 1})};
    constrained_col.setFromTriplets(triplets.begin(), triplets.end());
    A_mat.col(i) = constrained_col;
  }
  // transpose again to make the edited columns the rows again
  A_mat = A_mat.transpose();

  // =========================CREATE B VECTOR=========================

  // all zero except for displaced_verts
  // those should have displacement of displacement_vector
  std::vector<T> b_triplets = {};
  for (int i : endVertices)
  {
    b_triplets.push_back(T({i, 0, 1.}));
  }
  SparseMatrix<double> b_vec = SparseMatrix<double>(mesh->nVertices(), 3);
  b_vec.setFromTriplets(b_triplets.begin(), b_triplets.end());
  if (debug)
    std::cout << "b=\n"
              << b_vec.toDense() << "\n\n";

  // ===========================SOLVE SYSTEM===========================

  Eigen::SparseLU<SparseMatrix<double>> sluSolver;
  sluSolver.compute(A_mat);
  Eigen::MatrixXd x = sluSolver.solve(b_vec);

  // ========================STORE VERTEX DATA=======================

  hFunction = VertexData<double>(*mesh);
  for (int i = 0; i < x.rows(); ++i)
  {
    hFunction[i] = x(i, 0);
  }
}

// calculates the NORMALIZED gradient of the provided H function for each face of the mesh
void constructFaceGradient()
{
  geometry->requireFaceNormals();
  geometry->requireFaceAreas();
  faceGradient = FaceData<Vector3>(*mesh);
  for (Face face : mesh->faces())
  {
    Vector3 gradientVector = {0., 0., 0.};
    Vector3 faceNormal = geometry->faceNormals[face];
    double faceArea = geometry->faceAreas[face];
    // only need to look at first 2 halfedges
    int hesProcessed = 0;
    Halfedge he = *face.adjacentHalfedges().begin();
    while (hesProcessed < 2)
    {
      Vertex tip = he.tipVertex();
      Vertex tail = he.tailVertex();
      Vertex other = he.next().tipVertex();
      Vector3 edgeVector = geometry->vertexPositions[tip] - geometry->vertexPositions[tail];
      Vector3 edgePerpVector = cross(faceNormal, edgeVector);
      double functionDiff = hesProcessed == 0 ? hFunction[other] - hFunction[tip] : hFunction[other] - hFunction[tail];
      gradientVector += functionDiff * edgePerpVector / (2. * faceArea);
      hesProcessed++;
      he = he.next();
    }
    faceGradient[face] = gradientVector.normalize();
  }
}

// calculates the gradient of the provided H function for each vertex of the mesh
// as the uniform average of the adjacent face gradients
void constructVertexGradient()
{
  vertexGradient = VertexData<Vector3>(*mesh);
  for (Vertex v : mesh->vertices())
  {
    Vector3 gradientVector = {0., 0., 0.};
    int numAdjFaces = 0;
    for (Face f : v.adjacentFaces())
    {
      gradientVector += faceGradient[f];
      numAdjFaces++;
    }
    gradientVector /= numAdjFaces;
    vertexGradient[v] = gradientVector;
  }
}

// constructs the tangent vector for each vertex based on the vertex gradient
void constructTangentVectors()
{
  geometry->requireVertexNormals();
  courseTangentVectors = VertexData<Vector2>(*mesh);
  waleTangentVectors = VertexData<Vector2>(*mesh);
  int numProcessed = 0;
  for (Vertex v : mesh->vertices())
  {
    Vector3 normalVec = geometry->vertexNormals[v];
    Halfedge he;
    for (Halfedge he_placeholder : v.outgoingHalfedges())
    {
      he = he_placeholder;
      break;
    }
    Vector3 heVec = geometry->vertexPositions[he.tipVertex()] - geometry->vertexPositions[he.tailVertex()];
    Vector2 u = projectOntoTangentPlane(vertexGradient[v], normalVec, heVec);
    double angleW = std::atan2(u.y, u.x) + M_PI / 2.0;
    std::complex<double> waleComplexDirectionField(std::complex<double>(cos(2.0 * angleW), sin(2.0 * angleW)));
    waleTangentVectors[v] = Vector2::fromComplex(waleComplexDirectionField).normalize();

    double angleC = angleW + M_PI / 2.0; 
    std::complex<double> courseComplexDirectionField(std::cos(2.0 * angleC), std::sin(2.0 * angleC));
    courseTangentVectors[v] = Vector2::fromComplex(courseComplexDirectionField).normalize();

    numProcessed++;
  }
}

Vector2 projectOntoTangentPlane(Vector3 vec, Vector3 normal, Vector3 he)
{
  Vector3 eY = cross(normal, he).normalize();
  Vector3 eX = cross(eY, normal).normalize();
  return {dot(eX, vec), dot(eY, vec)};
}

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void myCallback()
{
  if (ImGui::Button("Pick start boundary (h=0)"))
  {
    selectBoundary(0);
  }
  if (ImGui::Button("Pick end boundary (h=1)"))
  {
    selectBoundary(1);
  }
  if (ImGui::Button("Clear boundaries"))
  {
    for (int i = 0; i < boundaryVertices.size(); ++i)
    {
      boundaryVertices[i] = {};
      boundaryPoints[i] = {};
    }
    polyscope::removePointCloud("Start points (h=0)");
    polyscope::removePointCloud("End points (h=1)");
  }

  if (ImGui::Button("Generate H Function"))
  {
    if (boundaryVertices[0].size() == 0 || boundaryVertices[1].size() == 0)
    {
      polyscope::warning("Please select start and end boundaries first.");
    }
    else
    {
      constructHFunction(boundaryVertices[0], boundaryVertices[1]);
      auto hFunctionDisplay = psMesh->addVertexScalarQuantity("H Function", hFunction);
      hFunctionDisplay->setEnabled(true);
    }
  }
    ImGui::Separator();
  if (ImGui::Button("Construct Gradient"))
  {
    constructFaceGradient();
    constructVertexGradient();
    auto gradientDisplay = psMesh->addVertexVectorQuantity("Vertex Gradient", vertexGradient);
    gradientDisplay->setEnabled(true);
  }

  if (ImGui::Button("Construct Tangent Vectors"))
  {
    constructTangentVectors();
  }
      ImGui::Separator();
    bool freqChanged = false;

    freqChanged |= ImGui::SliderFloat("Wale Frequency", &waleFrequency, 0.0001f, 100.0f, "%.1f", 1.0f);
    freqChanged |= ImGui::SliderFloat("Course Frequency", &courseFrequency, 0.0001f, 100.0f, "%.1f", 1.0f);

    if (freqChanged)
    {
        generateWaleStripes();
        generateCourseStripes();
    }


  if (ImGui::Button("Generate Wale Curves"))
  {
    generateWaleStripes();
  }
  if (ImGui::Button("Generate Course Curves"))
  {
      generateCourseStripes();
  }
}

int main(int argc, char **argv)
{

  // Configure the argument parser
  args::ArgumentParser parser("geometry-central & Polyscope example project");
  args::Positional<std::string> inputFilename(parser, "mesh", "A mesh file.");

  // Parse args
  try
  {
    parser.ParseCLI(argc, argv);
  }
  catch (args::Help &h)
  {
    std::cout << parser;
    return 0;
  }
  catch (args::ParseError &e)
  {
    std::cerr << e.what() << std::endl;
    std::cerr << parser;
    return 1;
  }

  // Make sure a mesh name was given
  if (!inputFilename)
  {
    std::cerr << "Please specify a mesh file as argument" << std::endl;
    return EXIT_FAILURE;
  }

  // Initialize polyscope
  polyscope::init();

  // Set the callback function
  polyscope::state::userCallback = myCallback;

  // Load mesh
  std::tie(mesh, geometry) = readManifoldSurfaceMesh(args::get(inputFilename));

  // Register the mesh with polyscope
  psMesh = polyscope::registerSurfaceMesh(
      polyscope::guessNiceNameFromPath(args::get(inputFilename)),
      geometry->inputVertexPositions, mesh->getFaceVertexList(),
      polyscopePermutations(*mesh));

  // Set vertex tangent spaces
  geometry->requireVertexTangentBasis();
  VertexData<Vector3> vBasisX(*mesh);
  VertexData<Vector3> vBasisY(*mesh);
  for (Vertex v : mesh->vertices())
  {
    vBasisX[v] = geometry->vertexTangentBasis[v][0];
    vBasisY[v] = geometry->vertexTangentBasis[v][1];
  }

  auto vField =
      geometrycentral::surface::computeSmoothestVertexDirectionField(*geometry);
  psMesh->addVertexTangentVectorQuantity("VF", vField, vBasisX, vBasisY);

  // Give control to the polyscope gui
  polyscope::show();

  return EXIT_SUCCESS;
}