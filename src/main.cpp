static char help[] = "simpleFEM multiphysics code \n\n";

#include "core/utils.hpp"
#include <eigen3/Eigen/SparseLU>
#include <petscksp.h>

int main(int argc, char **args){

  std::string meshFileName = "../rectangle_3NodeTri.vtk";
  std::string materialType = "linearElastic";
  double E=1.0, nu=0.3;
  std::string elementType = "3NodeTri";
  std::string bcType = "pressure";
  double bcValue = 1.0;
  size_t bcTag = 1;

  std::vector<size_t> dbcTag = {3};
  std::vector<int> dbcDim = {0, 1};
  std::vector<double> dbcValue = {0.0, 0.0};

  Vec         x, b;   
  Mat         A; 
  KSP         ksp;      
  PC          pc;
  PetscMPIInt    rank, size;    
  PetscFunctionBeginUser;
  const int root = 0;
  PetscInt    nlocal;
  PetscInitialize(&argc, &args, (char *)0, help);
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  PetscPrintf(PETSC_COMM_WORLD,"Number of processors = %d\n",size);
  //PetscPrintf(PETSC_COMM_SELF,"rank = %d\n",rank);

  iovtk vtkobj;
  int meshDof = 0;
  material* mat;
  element* elem = element::elementMap[elementType]();
  mesh msh(elem);
  if (rank == root){
    vtkobj = iovtk::readvtk(meshFileName);
    mat = material::materialMap[materialType]();
    mat->set_properties(E, nu);
    msh.getMeshData(vtkobj);
    Eigen::MatrixXi topo = msh.topology.data;
    int elemDof = elem->dim;
    meshDof = elemDof * msh.nodes.size;
  }
  MPI_Bcast(&meshDof , 1, MPI_INT, root, MPI_COMM_WORLD);
  VecCreate(PETSC_COMM_WORLD, &x);
  VecSetSizes(x, PETSC_DECIDE, meshDof);
  VecSetFromOptions(x);
  VecDuplicate(x, &b);
  VecGetLocalSize(x, &nlocal);

  MatCreate(PETSC_COMM_WORLD, &A);
  MatSetSizes(A, nlocal, nlocal, meshDof, meshDof);
  MatSetFromOptions(A);
  MatSetUp(A);

  std::vector<std::pair<int, double>> nodeValuePair;

  if (rank == root){
    std::cout << "Number of Nodes : " << msh.nodes.size << std::endl;
    std::cout << "Number of Elements : " << msh.topology.size << std::endl;
    Eigen::SparseMatrix<double> K(meshDof, meshDof);
    std::vector<Eigen::Triplet<double>> tripletK;
    utils::assembleMatrix(A, tripletK, msh, mat);
    domain edgeDomain(elem, 1);
    edgeDomain.getDomainData(vtkobj);
    edgeDomain.getDomainTags(vtkobj);
    bc* press = bc::bcMap[bcType](); 
    press->set_values(msh, edgeDomain, {bcTag}, bcValue);
    std::vector<bc> bcs{*press};
    Eigen::VectorXd f(meshDof); 
    f.setZero();
    utils::assembleVector(b, f, msh, bcs);

    domain dbcBoundary = edgeDomain.getSubdomain(dbcTag);
    utils::getDBCvalues(nodeValuePair, dbcBoundary, dbcDim, dbcValue, elem->dim);
    std::sort(nodeValuePair.begin(), nodeValuePair.end());

    utils::applyDirichletBC(tripletK, nodeValuePair, f);
    K.setFromTriplets(tripletK.begin(), tripletK.end());

    std::cout << "K Matrix:(" << K.rows() << "x"<< K.cols() << ") : " << std::endl << K.block(0,0,K.rows(),K.cols()) << std::endl;
    //std::cout << "f : " << f.transpose() << std::endl;

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(K);
    Eigen::VectorXd solution = solver.solve(f);
    auto topo = vtkobj.get_topology(elem->get_elemType(2));
    iovtk outvtk(vtkobj.nodes, topo, elem);
    outvtk.writeVtk("displacement.vtk", solution.data());
  }
  VecAssemblyBegin(b);
  VecAssemblyEnd(b);
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

  utils::applyDirichletBC2(A, b, nodeValuePair);

  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetOperators(ksp, A, A);
  KSPGetPC(ksp, &pc);
  PCSetType(pc, PCJACOBI);
  KSPSetTolerances(ksp, 1.e-7, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
  KSPSetFromOptions(ksp);
  KSPSolve(ksp, b, x);

  //MatZeroRowsColumns(A, numDBC, rowsDBC, 1.0, NULL, NULL);
  MatView(A, PETSC_VIEWER_STDOUT_WORLD);
  VecView(x,PETSC_VIEWER_STDOUT_WORLD);

  PetscFinalize();
  return 0;
}