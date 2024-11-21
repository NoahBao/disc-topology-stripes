#include "geometrycentral/surface/manifold_surface_mesh.h"
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

using namespace geometrycentral;
using namespace geometrycentral::surface;
using T = Eigen::Triplet<double>;

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh *psMesh;

// Some algorithm parameters
float param1 = 42.0;

// Example computation function -- this one computes and registers a scalar
// quantity
void doWork()
{
  // polyscope::warning("Computing Gaussian curvature.\nalso, parameter value = " +
  //                    std::to_string(param1));

  // geometry->requireVertexGaussianCurvatures();
  // psMesh->addVertexScalarQuantity("curvature",
  //                                 geometry->vertexGaussianCurvatures,
  //                                 polyscope::DataType::SYMMETRIC);
  // Generate a guiding field
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

std::vector<int> determineFixedVertices()
{
  for (Vertex v : mesh->vertices())
  {
    if (!v.isBoundary())
    {
      continue;
    }
  }
}

VertexData<double> constructHFunction(std::vector<int> &fixedVertices)
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
  for (int i : fixedVertices)
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
  for (int i : fixedVertices)
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

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void myCallback()
{

  if (ImGui::Button("do work"))
  {
    doWork();
  }

  ImGui::SliderFloat("param", &param1, 0., 100.);
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
