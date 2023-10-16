#include "mesh.hpp"

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

