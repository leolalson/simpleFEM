#include "bc.hpp"

Eigen::VectorXd outNormal(const Eigen::MatrixXd &Elem_Coord)
{
  Eigen::VectorXd normalVector(2);
  double length = std::sqrt(std::pow((Elem_Coord(1,0)-Elem_Coord(0,0)), 2) + std::pow((Elem_Coord(1,1)-Elem_Coord(0,1)), 2));
  normalVector(0) = (Elem_Coord(1,1)-Elem_Coord(0,1))/length;
  normalVector(1) = (Elem_Coord(0,0)-Elem_Coord(1,0))/length;
  return normalVector;
}

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