#include <map>
#include <eigen3/Eigen/Core>

class material{
  public:
  double E;
  double nu;
  material();
  ~material();
  Eigen::MatrixXd D;
  std::string type;
  static std::map<std::string, material * (*) ()> materialMap; 
  virtual void set_properties(double Y, double u){};
};

class linearElastic : public material{
  public:
  linearElastic();
  linearElastic(double E, double u);
  ~linearElastic();
  static material * create() { return new linearElastic();}
  void set_properties(double Y, double u);

};