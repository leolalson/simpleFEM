#include "mat.hpp"
#include <iostream>

material::material(){};

material::~material(){};

std::map<std::string, material * (*) ()> material::materialMap 
                              = {{"linearElastic" , linearElastic::create}}; 

//virtual void material::set_properties(double Y, double u){}


linearElastic::linearElastic(){};

linearElastic::linearElastic(double Y, double u){
  this->set_properties(Y, u);
};

linearElastic::~linearElastic(){};

void linearElastic::set_properties(double Y, double u){
  E = Y;
  nu = u;
  D = Eigen::MatrixXd::Zero(3, 3);
  double mut = E/(1 - pow(nu,2));
  D.row(0) << 1, nu, 0;
  D.row(1) << nu, 1, 0;
  D.row(2) << 0, 0, (1-nu)/2;
  D *= mut;
}