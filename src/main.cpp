static char help[] = "simpleFEM multiphysics code \n\n";

#include "core/utils.hpp"
#include <eigen3/Eigen/SparseLU>
#include <petscksp.h>
#include <chrono>

int main(int argc, char **args){

  //std::string meshFileName = "../rectangle_3NodeTri_refined3.vtk";
  std::string meshFileName = args[1];
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

  iovtk vtkobj;
  int meshDof = 0;
  material* mat;
  mat = material::materialMap[materialType]();
  mat->set_properties(E, nu);

  element* elem = element::elementMap[elementType]();
  mesh msh(elem);
  domain edgeDomain(elem, 1);
  bc* press = bc::bcMap[bcType](); 

  double start_read_mesh = MPI_Wtime();
  if (rank == root){
    vtkobj = iovtk::readvtk(meshFileName);
    msh.getMeshData(vtkobj);
    int elemDof = elem->dim;
    meshDof = elemDof * msh.nodes.size;
    edgeDomain.getDomainData(vtkobj);
    edgeDomain.getDomainTags(vtkobj);
  }
  double end_read_mesh = MPI_Wtime();

  double start_comm = MPI_Wtime();
  MPI_Bcast(&meshDof , 1, MPI_INT, root, MPI_COMM_WORLD);
  MPI_Bcast(&msh.nodes.size , 1, MPI_INT, root, MPI_COMM_WORLD);
  MPI_Bcast(&msh.topology.size , 1, MPI_INT, root, MPI_COMM_WORLD);
  MPI_Bcast(&msh.topology.dim , 1, MPI_INT, root, MPI_COMM_WORLD);
  MPI_Bcast(&edgeDomain.size , 1, MPI_INT, root, MPI_COMM_WORLD);

  if(rank!=root){
    msh.topology.data.resize(msh.topology.size,elem->numNodes);
    msh.topology.data.setZero();

    msh.nodes.data.resize(msh.nodes.size,3);
    msh.nodes.data.setZero();

    edgeDomain.data.resize(edgeDomain.size, edgeDomain.elem->numNodes);
    edgeDomain.tags.resize(edgeDomain.size);
  }

  MPI_Bcast(msh.topology.data.data() , msh.topology.size *elem->numNodes, MPI_INT, root, MPI_COMM_WORLD);
  MPI_Bcast(msh.nodes.data.data() , msh.nodes.size *3, MPI_DOUBLE, root, MPI_COMM_WORLD);
  MPI_Bcast(edgeDomain.data.data() , edgeDomain.size *edgeDomain.elem->numNodes, MPI_INT, root, MPI_COMM_WORLD);
  MPI_Bcast(edgeDomain.tags.data() , edgeDomain.size, MPI_INT, root, MPI_COMM_WORLD);
  double end_comm = MPI_Wtime();

  double start_assemble = MPI_Wtime();
  VecCreate(PETSC_COMM_WORLD, &x);
  VecSetSizes(x, PETSC_DECIDE, meshDof);
  VecSetFromOptions(x);
  VecDuplicate(x, &b);
  VecGetLocalSize(x, &nlocal);

  MatCreate(PETSC_COMM_WORLD, &A);
  MatSetSizes(A, nlocal, nlocal, meshDof, meshDof);
  MatSetFromOptions(A);
  MatSetUp(A);

  utils::assembleMatrixMPI(A, msh, mat, rank, size);

  std::vector<std::pair<int, double>> nodeValuePair;
  press->set_values(msh, edgeDomain, {bcTag}, bcValue);
  std::vector<bc> bcs{*press};
  utils::assembleVectorMPI(b, msh, bcs, rank, size);

  double end_assemble = MPI_Wtime();

  double start_dbc = MPI_Wtime();
  if (rank == root){

    domain dbcBoundary = edgeDomain.getSubdomain(dbcTag);
    utils::getDBCvalues(nodeValuePair, dbcBoundary, dbcDim, dbcValue, elem->dim);
    std::sort(nodeValuePair.begin(), nodeValuePair.end());

  }
  double end_dbc = MPI_Wtime();

  double start_petsc_assemble = MPI_Wtime();
  VecAssemblyBegin(b);
  VecAssemblyEnd(b);
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  double end_petsc_assemble = MPI_Wtime();

  double start_solve = MPI_Wtime();
  utils::applyDirichletBC(A, b, nodeValuePair);
  
  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetOperators(ksp, A, A);
  KSPGetPC(ksp, &pc);
  PCSetType(pc, PCJACOBI);
  KSPSetTolerances(ksp, 1.e-7, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
  KSPSetFromOptions(ksp);
  KSPSolve(ksp, b, x);
  double end_solve = MPI_Wtime();

  //MatZeroRowsColumns(A, numDBC, rowsDBC, 1.0, NULL, NULL);
  //MatView(A, PETSC_VIEWER_STDOUT_WORLD);
  //VecView(b,PETSC_VIEWER_STDOUT_WORLD);

  double start_write_result = MPI_Wtime();
  Vec            x_seq;
  PetscScalar    *solution;
  VecScatter     ctx;
  //VecScatterCreateToAll(x,&ctx,&x_seq);
  VecScatterCreateToZero(x,&ctx,&x_seq);
  VecScatterBegin(ctx,x,x_seq,INSERT_VALUES,SCATTER_FORWARD);
  VecScatterEnd(ctx,x,x_seq,INSERT_VALUES,SCATTER_FORWARD);
  VecGetArray(x_seq,&solution);

  if(rank == root){
    auto topo = vtkobj.get_topology(elem->get_elemType(2));
    iovtk outvtk(vtkobj.nodes, topo, elem);
    outvtk.writeVtk("displacement.vtk", solution);
  }
  double end_write_result = MPI_Wtime();

  if(rank == root){
    std::cout << "Read mesh: " << end_read_mesh - start_read_mesh << " seconds" << std::endl;
    std::cout << "Communications: " << end_comm - start_comm << " seconds" << std::endl;
    std::cout << "Assemble: " << end_assemble - start_assemble << " seconds" << std::endl;
    std::cout << "dbc: " << end_dbc - start_dbc << " seconds" << std::endl;
    std::cout << "Petsc assemble: " << end_petsc_assemble - start_petsc_assemble << " seconds" << std::endl;
    std::cout << "Solve: " << end_solve - start_solve << " seconds" << std::endl;
    std::cout << "Write result: " << end_write_result - start_write_result << " seconds" << std::endl;

  }
  PetscFinalize();
  return 0;
}