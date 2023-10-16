#include "mesh.hpp"

mesh::mesh(element* elem) : elem(elem){

};

void mesh::getMeshData(iovtk ioobj){
  
  std::vector<std::vector<double>> nodeData = ioobj.nodes;
  std::vector<std::vector<size_t>> topo = ioobj.get_topology(this->elem->get_elemType(2));
  Eigen::MatrixXd pointData(nodeData.size(), nodeData[0].size());
  //this->nodes.resize(nodeData.size(), nodeData[0].size());
  for(int i=0;i<nodeData.size();++i){
    for(int j=0;j<nodeData[i].size();++j){
      pointData(i,j) = nodeData[i][j];
    }
  }
  this->nodes = points(pointData);
  this->topology = domain(this->elem, elem->dim);
  this->topology.getDomainData(ioobj);
  this->topology.getDomainTags(ioobj);

}

domain::domain(element* elem, int eDim){
  this->elem = elem->getElement(elem, eDim);
  this->dim = eDim;
};

void domain::getDomainData(iovtk ioobj){
  std::vector<std::vector<size_t>> topo = ioobj.get_topology(this->elem->get_elemType(dim));
  this->size = topo.size();
  this->data.resize(topo.size(), topo[0].size());
  for(int i=0;i<topo.size();++i){
    for(int j=0;j<topo[i].size();++j){
      this->data(i,j) = topo[i][j];
    }
  }
}

std::map<size_t, size_t> domain::getDomainTags(iovtk ioobj){
  std::vector<size_t> ioTags = ioobj.get_data("CELL_DATA", this->elem->get_elemType(dim));
  size_t tagCount = 0;
  std::map<size_t, size_t> tagMap;
  tagMap[ioTags[0]] = tagCount;
  if (this->size == ioTags.size()){
    this->tags.resize(ioTags.size());
    for(int i=0;i<ioTags.size();++i){
      if(tagMap.find(ioTags[i])!=tagMap.end()){
        this->tags(i) = tagMap[ioTags[i]];
      }
      else{
        ++tagCount;
        tagMap[ioTags[i]] = tagCount;
        this->tags(i) = tagMap[ioTags[i]];
      }
    }
  }
  return tagMap;
}

domain domain::getSubdomain(std::vector<size_t> subDomainTags){
  domain* dmn = new domain();
  dmn->elem = this->elem;
  dmn->dim = this->dim;
  Eigen::MatrixXi topo = this->data;
  Eigen::VectorXi domainTags = this->tags;
  std::vector<int> pos;
  pos.reserve(domainTags.size());

  for(int i=0;i<domainTags.size();++i){
    if(std::find(subDomainTags.begin(), subDomainTags.end(), domainTags[i])!=subDomainTags.end()){
      pos.push_back(i);
    }
  }
  dmn->data = this->data(pos, Eigen::all);
  dmn->tags = this->tags(pos);
  dmn->size = pos.size();
  return *dmn;
}

