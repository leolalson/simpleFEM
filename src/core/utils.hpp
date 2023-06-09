#pragma once

#include<vector>
#include<iostream>
#include <eigen3/Eigen/SparseLU>
#include <petscksp.h>

#include "bc.hpp"
#include "mesh.hpp"
#include "mat.hpp"

namespace utils{

template <typename T>
void printVector(std::vector<T> myVector, std::string name){
  std::cout << name << ":(" << myVector.size() << ") : " ;
  for(int i=0;i<myVector.size();++i){
    std::cout << myVector[i] << ", "; 
  }
  std::cout << std::endl;
}

template <typename T>
void printArray(T myVector[], int size, std::string name){
  std::cout << name << ":(" << size << ") : " ;
  for(int i=0;i<size;++i){
    std::cout << myVector[i] << ", "; 
  }
  std::cout << std::endl;
}

template <typename T>
void printVector2D(std::vector<std::vector<T>> myVector, std::string name){
  std::cout.precision(17);
  std::cout << name << ":(" << myVector.size() << ")" << std::endl;
  for(int i=0;i<myVector.size();++i){
    for(int j=0;j<myVector[i].size();++j){
      std::cout << myVector[i][j] << ", "; 
    }
    std::cout << std::endl;
  }
}

void assembleMatrix(Mat &A, std::vector<Eigen::Triplet<double>> &tripletK, mesh &msh, material* mat){
  int dim = msh.topology.dim;
  PetscInt dof = msh.elem->numNodes * dim;
  PetscInt* elemtopoDof2 = new PetscInt[dof];
  PetscScalar* value = new PetscScalar[dof * dof];
  Eigen::MatrixXd nodes2d = msh.nodes.data(Eigen::all, {0, 1});
  Eigen::MatrixXi topo = msh.topology.data;
  Eigen::MatrixXi topoDof(topo.rows(), dim * topo.cols());
  topoDof.setZero();
  for(int i=0;i<topo.cols();++i){
    for(int j=0;j<dim;++j){
      topoDof(Eigen::all, {(i*dim)+j}) = ((topo.col(i).array() * dim) + j);
    }
  }
  for(int e=0;e<msh.topology.size;++e){
    Eigen::VectorXi elemTopo = topo.row(e);
    Eigen::VectorXi elemTopoDof = topoDof.row(e);
    Eigen::MatrixXd elemNodes = nodes2d(elemTopo, Eigen::all);
    Eigen::MatrixXd K_local = msh.elem->Kmatrix(elemNodes, mat->D);

    //Eigen assembly
    for(int j=0;j<elemTopoDof.size();++j){
      for(int k=0;k<elemTopoDof.size();++k){
        tripletK.push_back(Eigen::Triplet<double> (elemTopoDof[j], elemTopoDof[k], K_local(j,k)));
      }
    }
    
    //Petsc assembly
    std::copy(elemTopoDof.data(), elemTopoDof.data() + elemTopoDof.size(), elemtopoDof2);
    value = K_local.data();
    MatSetValues(A, dof, elemtopoDof2, dof, elemtopoDof2, value, ADD_VALUES);
  }
}


void assembleVector(Vec &b, Eigen::VectorXd &f_global, mesh msh, std::vector<bc> bcs){
  bc bc0 = bcs[0];  
  int dim = msh.topology.dim;
  PetscInt dof = bc0.boundary.elem->numNodes * dim;
  PetscInt* elemtopoDof2 = new PetscInt[dof];
  PetscScalar* value = new PetscScalar[dof];
  Eigen::MatrixXd nodes2d = msh.nodes.data(Eigen::all, {0, 1});
  Eigen::MatrixXi topo = bc0.boundary.data;
  Eigen::MatrixXi topoDof(topo.rows(), dim * topo.cols());

  for(int i=0;i<topo.cols();++i){
    for(int j=0;j<dim;++j){
      topoDof(Eigen::all, {(i*dim)+j}) = ((topo.col(i).array() * dim) + j);
    }
  }

  for(int e=0;e<bc0.boundary.size;++e){
    Eigen::VectorXi elemTopo = topo.row(e);
    Eigen::VectorXi elemTopoDof = topoDof.row(e);
    Eigen::MatrixXd elemNodes = nodes2d(elemTopo, Eigen::all);
    Eigen::VectorXd f_local = bc0.boundary.elem->fvector(elemNodes, bc0.bcVector.row(e));
    for(int i=0;i<elemTopoDof.size();++i){
      f_global(elemTopoDof(i)) += f_local(i);
    }
    //Petsc assembly
    std::copy(elemTopoDof.data(), elemTopoDof.data() + elemTopoDof.size(), elemtopoDof2);
    value = f_local.data();
    VecSetValues(b, dof, elemtopoDof2, value, ADD_VALUES);
  }

}

void getDBCvalues(std::vector<std::pair<int, double>> &nodeValuePair, domain dbcBoundary, std::vector<int> dofDim, 
                    std::vector<double> dofValue, int dim){
  Eigen::MatrixXi topo = dbcBoundary.data;
  std::vector<int> topoVec(topo.data(), topo.data()+topo.cols()*topo.rows());
  std::sort( topoVec.begin(), topoVec.end());
  topoVec.erase(std::unique(topoVec.begin(), topoVec.end()), topoVec.end());

  nodeValuePair.reserve((topoVec.size() * dofDim.size()) + nodeValuePair.size());
  for(int i=0;i<dofDim.size();++i){
    for(int j=0;j<topoVec.size();++j){
      nodeValuePair.push_back(std::make_pair((topoVec[j]*dim)+dofDim[i], dofValue[i]));
    }
  }
}

void applyDirichletBC(std::vector<Eigen::Triplet<double>> &tripletK, std::vector<std::pair<int, double>> nodeValuePair,
                        Eigen::VectorXd &f_global){

  std::vector<int> nodedof;
  nodedof.reserve(nodeValuePair.size());
	std::vector<double> valuedof;
  valuedof.reserve(nodeValuePair.size());

  for(int i=0;i<nodeValuePair.size();++i){
    nodedof.push_back(nodeValuePair[i].first);	
		valuedof.push_back(nodeValuePair[i].second);
    std::cout << "nodedof: " << nodeValuePair[i].first << ", valuedof: " << nodeValuePair[i].second << std::endl;
  }

	std::vector<int>::iterator it1;
	std::vector<int>::iterator it2;
	std::vector<Eigen::Triplet<double>> tripletDBC;
	for(int index=0;index<tripletK.size();++index){
		it1 = std::find(nodedof.begin(), nodedof.end(), tripletK[index].col());
		it2 = std::find(nodedof.begin(), nodedof.end(), tripletK[index].row());
		if(it1 !=nodedof.end() || it2 !=nodedof.end()){
			tripletK.erase(tripletK.begin() + index);
			--index;
		}

	}
  for(int index=0;index<nodedof.size();++index){
    tripletDBC.push_back(Eigen::Triplet<double> (nodedof[index], nodedof[index], 1.0));
  }

	tripletK.insert(tripletK.end(), tripletDBC.begin(), tripletDBC.end());
	for(int index=0;index<nodedof.size();++index){
		f_global[nodedof[index]] = valuedof[index];
	}
	
	return;
}

void applyDirichletBC2(Mat &A, Vec &b, std::vector<std::pair<int, double>> nodeValuePair){

  PetscInt numDBC = nodeValuePair.size();
  PetscInt rowsDBC[numDBC], dof;
  Vec valuesDBC;
  PetscScalar val = 1.0;

  VecGetSize(b, &dof);
  VecCreate(PETSC_COMM_WORLD, &valuesDBC);
  VecSetSizes(valuesDBC, PETSC_DECIDE, dof);
  VecSetFromOptions(valuesDBC);

  for (PetscInt i=0;i<nodeValuePair.size();++i){
    rowsDBC[i] = nodeValuePair[i].first;
    VecSetValue(valuesDBC, nodeValuePair[i].first, nodeValuePair[i].second, INSERT_VALUES);
  }
  VecAssemblyBegin(valuesDBC);
  VecAssemblyEnd(valuesDBC);
  MatZeroRowsColumns(A, numDBC, rowsDBC, 1.0, valuesDBC, b);
	
	return;
}

} // end namespace utils