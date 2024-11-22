#include <geometrycentral/surface/manifold_surface_mesh.h>
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/direction_fields.h"

// STRIPES
#include "geometrycentral/surface/stripe_patterns.h"
#include <polyscope/polyscope.h>
#include <polyscope/point_cloud.h>
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

VertexData<double> constructHFunction(std::vector<int> &startVertices, std::vector<int> &endVertices);

void selectWholeBoundary(int rootIndex, int boundaryId);

// ==========================================

// =============GLOBAL VARIABLES=============

// boundaryVertices[i] = list of vertices for boundary i
// analogous for boundaryPoints
std::vector<std::vector<int>> boundaryVertices = {{}, {}};
std::vector<std::vector<Vector3>> boundaryPoints = {{}, {}};

// ==========================================

// =============DEPRECATED VARIABLES/FUNCTIONS=============

enum MeshShape
{
  SQUARE,
  ANNULUS_CUT
};
MeshShape inputMeshShape = SQUARE;

std::vector<std::vector<int>> determineFixedVertices();
std::vector<std::vector<int>> determineFixedVerticesSquare();
std::vector<std::vector<int>> determineFixedVerticesAnnulusCut();
int findSharedVertexId(Edge e1, Edge e2);

// ===========================================================

// generates stripes using Geometry Central's stripe methods
void generateStripes()
{
  // Generate a guiding field
  // TODO: use our own guide field based on the user's start and end boundaries
  VertexData<Vector2> guideField =
      geometrycentral::surface::computeSmoothestVertexDirectionField(*geometry, 2);

  // Compute the stripe pattern
  double constantFreq = 40.;
  VertexData<double> frequencies(*mesh, constantFreq);
  CornerData<double> periodicFunc;
  FaceData<int> zeroIndices;
  FaceData<int> branchIndices;
  std::tie(periodicFunc, zeroIndices, branchIndices) =
      computeStripePattern(*geometry, frequencies, guideField);

  // Extract isolines
  std::vector<Vector3> isolineVerts;
  std::vector<std::array<size_t, 2>> isolineEdges;
  std::tie(isolineVerts, isolineEdges) = extractPolylinesFromStripePattern(
      *geometry, periodicFunc, zeroIndices, branchIndices, guideField, false);
  polyscope::registerPointCloud("Curve points", isolineVerts);
  polyscope::registerCurveNetwork("Curve edges", isolineVerts, isolineEdges);
}

// DEPRECATED - DO NOT USE
// determines the vertices to fix based on the input mesh's shape
std::vector<std::vector<int>> determineFixedVertices()
{
  switch (inputMeshShape)
  {
  case SQUARE:
    return determineFixedVerticesSquare();
  case ANNULUS_CUT:
    return determineFixedVerticesAnnulusCut();
  }
  return {};
}

// DEPRECATED - DO NOT USE
// determines the vertices to fix for square meshes
std::vector<std::vector<int>> determineFixedVerticesSquare()
{
  geometry->requireVertexPositions();
  std::vector<int> fixedVertexIDs = {};
  double minY = 0.;
  double maxY = 0.;
  for (Vertex v_i : mesh->vertices())
  {
    // only consider boundary vertices
    if (!v_i.isBoundary())
    {
      continue;
    }
    bool isOnStartOrEnd = false;
    for (Halfedge he : v_i.incomingHalfedges())
    {
      // only consider boundary edges
      if (!he.edge().isBoundary())
      {
        continue;
      }
      Vertex v_j = he.tailVertex();
      Vector3 dispVector = geometry->vertexPositions[v_j] - geometry->vertexPositions[v_i];
      // if edge is flat
      if (dispVector.y == 0)
      {
        isOnStartOrEnd = true;
      }
    }
    if (isOnStartOrEnd)
    {
      fixedVertexIDs.push_back(v_i.getIndex());
      // std::cout << v_i.getIndex() << "\n";
      minY = std::min(minY, geometry->vertexPositions[v_i].y);
      maxY = std::max(maxY, geometry->vertexPositions[v_i].y);
    }
  }

  std::vector<int> bottomFixedVertexIds = {};
  std::vector<int> topFixedVertexIds = {};
  for (int id : fixedVertexIDs)
  {
    if (geometry->vertexPositions[id].y == minY)
    {
      bottomFixedVertexIds.push_back(id);
    }
    else
    {
      topFixedVertexIds.push_back(id);
    }
  }
  return {bottomFixedVertexIds, topFixedVertexIds};
}

// DEPRECATED - DO NOT USE
// determines the vertices to fix for annulus meshes
std::vector<std::vector<int>> determineFixedVerticesAnnulusCut()
{
  geometry->requireVertexPositions();
  std::vector<int> fixedVertexIDs = {};
  Edge prevEdge;
  bool isFirstEdge = true;
  Edge firstEdge;
  int loopNum = 0;
  for (Edge e : mesh->boundaryLoop(0).adjacentEdges())
  {
    if (isFirstEdge)
    {
      firstEdge = e;
      isFirstEdge = false;
    }
    else
    {
      Vertex v1 = mesh->vertex(findSharedVertexId(prevEdge, e));
      Vertex v0 = prevEdge.otherVertex(v1);
      Vertex v2 = e.otherVertex(v1);
      Vector3 prevEdgeVector = geometry->vertexPositions[v1] - geometry->vertexPositions[v0];
      Vector3 currEdgeVector = geometry->vertexPositions[v2] - geometry->vertexPositions[v1];
      if (std::abs(dot(prevEdgeVector.normalize(), currEdgeVector.normalize())) == 1.)
      {
        fixedVertexIDs.push_back(v1.getIndex());
        std::cout << v1.getIndex() << "\n";
      }
    }
    prevEdge = e;
  }

  std::vector<int> bottomFixedVertexIds = {};
  std::vector<int> topFixedVertexIds = {};
  for (int id : fixedVertexIDs)
  {
    bottomFixedVertexIds.push_back(id);
    topFixedVertexIds.push_back(id);
  }
  return {bottomFixedVertexIds, topFixedVertexIds};
}

// finds the ID of the vertex shared by the two edges
// returns -1 if they don't share a vertex
int findSharedVertexId(Edge e1, Edge e2)
{
  if (e1.firstVertex() == e2.firstVertex() || e1.firstVertex() == e2.secondVertex())
  {
    return e1.firstVertex().getIndex();
  }
  else if (e1.secondVertex() == e2.firstVertex() || e1.secondVertex() == e2.secondVertex())
  {
    return e1.secondVertex().getIndex();
  }
  else
  {
    return -1;
  }
}

// uses Laplacian smoothing to construct H function based on start and stop boundaries
VertexData<double> constructHFunction(std::vector<int> &startVertices, std::vector<int> &endVertices)
{
  // SparseMatrix<double> laplacianMatrix = SparseMatrix<double>(mesh->nVertices(), mesh->nVertices());
  // std::vector<Eigen::Triplet<double>> lMatTriplets = {};
  // for(Vertex v_i : mesh->vertices()) {
  //   for(Halfedge he : v_i.incomingHalfedges()) {
  //     Vertex v_j = he.tailVertex();
  //   }
  // }
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

  geometrycentral::surface::VertexData<double> hFunction(*mesh);
  for (int i = 0; i < x.rows(); ++i)
  {
    hFunction[i] = x(i, 0);
  }
  return hFunction;
}

// selects all boundary vertices that are similarly aligned to the user's selected vertex
void selectWholeBoundary(int rootIndex, int boundaryId)
{
  geometry->requireVertexPositions();
  Vector3 boundarySlope;
  for (Halfedge he : mesh->vertex(rootIndex).incomingHalfedges())
  {
    Vertex v_j = he.tailVertex();
    if (v_j.isBoundary())
    {
      boundarySlope = geometry->vertexPositions[v_j] - geometry->vertexPositions[rootIndex];
      break;
    }
  }
  boundarySlope = boundarySlope.normalize();
  // intitialize BFS
  std::map<int, bool> isVertexVisited;
  for (Vertex vertex : mesh->vertices())
  {
    isVertexVisited[vertex.getIndex()] = false;
  }
  std::queue<int> vertices;
  vertices.push(rootIndex);
  isVertexVisited[rootIndex] = true;
  // run BFS
  while (!vertices.empty())
  {
    Vertex v_i = mesh->vertex(vertices.front());
    vertices.pop();

    boundaryVertices[boundaryId].push_back(v_i.getIndex());
    boundaryPoints[boundaryId].push_back(geometry->vertexPositions[v_i]);

    for (Halfedge he : v_i.incomingHalfedges())
    {
      Vertex v_j = he.tailVertex();
      Vector3 slope = geometry->vertexPositions[v_j] - geometry->vertexPositions[v_i];
      slope = slope.normalize();
      if (v_j.isBoundary() && std::abs(dot(slope, boundarySlope)) >= 1 - 1e-6 && !isVertexVisited[v_j.getIndex()])
      {
        isVertexVisited[v_j.getIndex()] = true;
        vertices.push(v_j.getIndex());
      }
    }
  }
}

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void myCallback()
{
  if (ImGui::Button("Pick start boundary (h=0)"))
  {
    int selected = psMesh->selectVertex();
    boundaryVertices[0] = {};
    boundaryPoints[0] = {};
    selectWholeBoundary(selected, 0);
    auto pc = polyscope::registerPointCloud("Start points (h=0)", boundaryPoints[0]);
    pc->setPointColor({1., 0.1, 1.});
  }
  if (ImGui::Button("Pick end boundary (h=1)"))
  {
    int selected = psMesh->selectVertex();
    boundaryVertices[1] = {};
    boundaryPoints[1] = {};
    selectWholeBoundary(selected, 1);
    auto pc = polyscope::registerPointCloud("End points (h=1)", boundaryPoints[1]);
    pc->setPointColor({1., 1., 0.1});
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
      VertexData<double> hFunction = constructHFunction(boundaryVertices[0], boundaryVertices[1]);
      auto hFunctionDisplay = psMesh->addVertexScalarQuantity("H Function", hFunction);
      hFunctionDisplay->setEnabled(true);
    }
  }

  if (ImGui::Button("Construct Gradient"))
  {
    polyscope::warning("Unimplemented");
  }

  if (ImGui::Button("Construct Tangent Vectors"))
  {
    polyscope::warning("Unimplemented");
  }

  if (ImGui::Button("Generate Curves"))
  {
    generateStripes();
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