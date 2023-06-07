#include "fem.hpp"

double BMatrix(Eigen::MatrixXd &B, Eigen::MatrixXd coords, Eigen::Vector3d N, Eigen::MatrixXd dNdr){

  Eigen::MatrixXd jacobian = dNdr * coords;
  double det = jacobian.determinant();
  Eigen::MatrixXd jacobian_inv = jacobian.inverse();

  Eigen::MatrixXd dNdx = jacobian_inv * dNdr;

  B({0},{0,2,4}) = dNdx.row(0);
  B({1},{1,3,5}) = dNdx.row(1);
  B({2},{0,2,4}) = dNdx.row(1);
  B({2},{1,3,5}) = dNdx.row(0);

  return det;
}

void NMatrix(Eigen::MatrixXd &Nmat, Eigen::VectorXd N){

  Nmat({0},{0,2}) = N;
  Nmat({1},{1,3}) = N;
}

element::element(){};

element::~element(){};

std::map<std::string, element * (*) ()> element::elementMap = {{"3NodeTri" , tri3::create}};

tri3::tri3(){
  dim = 2;
  numNodes = 3;
  //elemType = 5;
  elemType = {{2, 5}, {1, 3}, {0, 1}};
  thick = 1; //thickness assumed to be 1
  dof = dim * numNodes;
  gaussPts2D = Eigen::MatrixXd{{0.1666666666667, 0.1666666666667},
                              {0.6666666666667, 0.1666666666667},
                              {0.1666666666667, 0.6666666666667}};

  gaussWts2D = Eigen::Vector3d(0.333333333333, 0.333333333333,0.333333333333);

  gaussPts1D = Eigen::MatrixXd{{-1/sqrt(3)}, {1/sqrt(3)}};
  gaussWts1D = Eigen::Vector2d(1, 1);

};

tri3::~tri3(){};




Eigen::MatrixXd tri3::Kmatrix(Eigen::MatrixXd coords, Eigen::MatrixXd D){
  Eigen::MatrixXd B(3,6);
  //Eigen::MatrixXd dNdr{{-1, 1, 0}, {-1, 0, 1}};
  Eigen::MatrixXd K_local(dof, dof);
  K_local.setZero();
  for(int i=0;i<gaussPts2D.rows();++i){
    Eigen::VectorXd rs = gaussPts2D.row(i);
    B.setZero();

    Eigen::Vector3d N(1-rs(0)-rs(1), rs(0), rs(1));
    Eigen::MatrixXd dNdr{{-1, 1, 0}, {-1, 0, 1}};

    double det = BMatrix(B, coords, N, dNdr);
    K_local += 0.5 * det * thick * gaussWts2D(i) * (B.transpose() * D * B); 
  }
  return K_local;
};

Eigen::VectorXd tri3::fvector(Eigen::MatrixXd Xe, Eigen::VectorXd F)
{
  Eigen::VectorXd f_local(4);
  f_local.setZero();
  int ngp = 2;
  double le = sqrt(pow((Xe(Xe.rows()-1, 0) - Xe(0,0)), 2) + pow((Xe(Xe.rows()-1, 1) - Xe(0,1)), 2));
  Eigen::MatrixXd Nmat(2, 4);
  // Loop over integration points of boundary element
  for (int i = 0; i < ngp; i++)
  {
    Nmat.setZero();
    Eigen::VectorXd rs = gaussPts1D.row(i);
    Eigen::Vector2d N((1-rs(0))/2, (1+rs(0))/2);
    NMatrix(Nmat, N);
    f_local += (Nmat.transpose() * F * gaussWts1D(i)) * thick * le /2;

  }
  return f_local;
};