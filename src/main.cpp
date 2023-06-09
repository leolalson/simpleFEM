#include "core/utils.hpp"
#include "core/mesh.hpp"
#include "core/mat.hpp"
#include <eigen3/Eigen/SparseLU>

#include <iomanip>

Eigen::VectorXd outNormal(const Eigen::MatrixXd &Elem_Coord)
{
  Eigen::VectorXd normalVector(2);
  double length = std::sqrt(std::pow((Elem_Coord(1,0)-Elem_Coord(0,0)), 2) + std::pow((Elem_Coord(1,1)-Elem_Coord(0,1)), 2));
  normalVector(0) = (Elem_Coord(1,1)-Elem_Coord(0,1))/length;
  normalVector(1) = (Elem_Coord(0,0)-Elem_Coord(1,0))/length;
  return normalVector;
}


class bc{
  public:
  bc(){};
  ~bc(){};
  domain boundary;
  Eigen::MatrixXd bcVector;
  static std::map<std::string, bc * (*) ()> bcMap;
  virtual void set_values(mesh msh, domain bcDmn, std::vector<size_t> domainTags, std::vector<double> bcValues){};
  virtual void set_values(mesh msh, domain bcDmn, std::vector<size_t> domainTags, double bcValue){};
};

class pressure : public bc{
  public:
  pressure(){};
  ~pressure(){};
  static bc * create() { return new pressure(); }
  void set_values(mesh msh, domain bcDmn, std::vector<size_t> domainTags, std::vector<double> loadValues);
  void set_values(mesh msh, domain bcDmn, std::vector<size_t> domainTags, double bcValue);

};

std::map<std::string, bc * (*) ()> bc::bcMap = {{"pressure" , pressure::create}};

void pressure::set_values(mesh msh, domain bcDmn, std::vector<size_t> domainTags, std::vector<double> bcValues){
  this->boundary = bcDmn.getSubdomain(domainTags);
  bcVector.resize(this->boundary.size, msh.elem->dim);
  for(int i=0;i<this->boundary.size;++i){
    for(int j=0;j<msh.elem->dim;++j){
      this->bcVector(i, j) = bcValues[j];
    }
  }
}

void pressure::set_values(mesh msh, domain bcDmn, std::vector<size_t> domainTags, double bcValue){
  this->boundary = bcDmn.getSubdomain(domainTags);
  bcVector.resize(this->boundary.size, msh.elem->dim);
  Eigen::MatrixXd nodes2d = msh.nodes.data(Eigen::all, {0, 1});
  Eigen::MatrixXi topo = this->boundary.data;
  for(int i=0;i<topo.rows();++i){
    Eigen::VectorXi elemTopo = topo.row(i);
    Eigen::MatrixXd elemNodes = nodes2d(elemTopo, Eigen::all);
    Eigen::VectorXd normal = outNormal(elemNodes);
    this->bcVector.row(i) = normal * bcValue;
  }
}

void assembleMatrix(std::vector<Eigen::Triplet<double>> &tripletK, mesh &msh, material* mat){
  int dim = msh.topology.dim;
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
    for(int j=0;j<elemTopoDof.size();++j){
      for(int k=0;k<elemTopoDof.size();++k){
        tripletK.push_back(Eigen::Triplet<double> (elemTopoDof[j], elemTopoDof[k], K_local(j,k)));
      }
    }
  }
}


void assembleVector(Eigen::VectorXd &f_global, mesh msh, std::vector<bc> bcs){
  bc bc0 = bcs[0];
  int dim = msh.topology.dim;
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
    Eigen::VectorXd f_local = msh.elem->fvector(elemNodes, bc0.bcVector.row(e));
    for(int i=0;i<elemTopoDof.size();++i){
      f_global(elemTopoDof(i)) += f_local(i);
    }
  }

}

void getDBCvalues(std::vector<std::pair<int, double>> &nodeValuePair, domain dbcBoundary, std::vector<int> dofDim, 
                    std::vector<double> dofValue){
  Eigen::MatrixXi topo = dbcBoundary.data;
  std::vector<int> topoVec(topo.data(), topo.data()+topo.cols()*topo.rows());
  std::sort( topoVec.begin(), topoVec.end() );
  topoVec.erase(std::unique(topoVec.begin(), topoVec.end()), topoVec.end());

  nodeValuePair.reserve((topoVec.size() * dofDim.size()) + nodeValuePair.size());
  int dim = dbcBoundary.elem->dim;
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
  }

	std::vector<int>::iterator it1;
	std::vector<int>::iterator it2;
	std::vector<Eigen::Triplet<double>> tripletDBC;
	for(int index=0;index<tripletK.size();++index){
		it1 = std::find(nodedof.begin(), nodedof.end(), tripletK[index].col());
		it2 = std::find(nodedof.begin(), nodedof.end(), tripletK[index].row());
		if(it1 !=nodedof.end() || it2 !=nodedof.end()){
			if (tripletK[index].row() == tripletK[index].col()){
				tripletDBC.push_back(Eigen::Triplet<double> (tripletK[index].row(), tripletK[index].col(), 1.0));
			}
			tripletK.erase(tripletK.begin() + index);
			--index;
		}

	}
	tripletK.insert(tripletK.end(), tripletDBC.begin(), tripletDBC.end());
	for(int index=0;index<nodedof.size();++index){
		f_global[nodedof[index]] = valuedof[index];
	}
	
	return;
}

int main(){

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

  

  iovtk vtkobj = iovtk::readvtk(meshFileName);

  material* mat = material::materialMap[materialType]();
  mat->set_properties(E, nu);

  element* elem = element::elementMap[elementType]();
  
  mesh msh(elem);
  msh.getMeshData(vtkobj);

  int elemDof = elem->dim;
  int meshDof = elemDof * msh.nodes.size;
  std::cout << "Number of Nodes : " << msh.nodes.size << std::endl;
  std::cout << "Number of Elements : " << msh.topology.size << std::endl;

  Eigen::SparseMatrix<double> K(meshDof, meshDof);
  std::vector<Eigen::Triplet<double>> tripletK;
  assembleMatrix(tripletK, msh, mat);

  domain edgeDomain(elem, 1);
  edgeDomain.getDomainData(vtkobj);
  edgeDomain.getDomainTags(vtkobj);
  bc* press = bc::bcMap[bcType](); 
  press->set_values(msh, edgeDomain, {bcTag}, bcValue);
  std::vector<bc> bcs{*press};
  Eigen::VectorXd f(meshDof); 
  f.setZero();
  assembleVector(f, msh, bcs);

  domain dbcBoundary = edgeDomain.getSubdomain(dbcTag);

  std::vector<std::pair<int, double>> nodeValuePair;
  getDBCvalues(nodeValuePair, dbcBoundary, dbcDim, dbcValue);
  std::sort(nodeValuePair.begin(), nodeValuePair.end());

  applyDirichletBC(tripletK, nodeValuePair, f);
  K.setFromTriplets(tripletK.begin(), tripletK.end());

  //std::cout << "K Matrix:(" << K.rows() << "x"<< K.cols() << ") : " << std::endl << K.block(0,0,K.rows(),K.cols()) << std::endl;
  //std::cout << "f : " << f.transpose() << std::endl;

  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(K);
  Eigen::VectorXd solution = solver.solve(f);

  auto topo = vtkobj.get_topology(elem->get_elemType(2));
  iovtk outvtk(vtkobj.nodes, topo, elem);
  outvtk.writeVtk("displacement.vtk", solution.data());

  return 0;
}