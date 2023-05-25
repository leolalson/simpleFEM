#include "io.hpp"

iovtk iovtk::readvtk(std::string fileName){
  iovtk* vtkobj = new iovtk;
  vtkobj->fileName = fileName;
  vtkobj->nodes = vtkobj->read_nodes(fileName);
  vtkobj->elements = vtkobj->read_elements(fileName);

  std::vector<std::string> scalarDataID = {"CELL_TYPES", "CELL_DATA"};
  for(int i=0;i<scalarDataID.size();++i){
    std::vector<size_t> data = vtkobj->readScalarData(fileName, scalarDataID[i]);
    if (~data.empty()){
      vtkobj->data.push_back(data);
      vtkobj->scalarDataID.push_back(scalarDataID[i]);
    }
  }
  return *vtkobj;
};

std::vector<std::vector<double>> iovtk::read_nodes(std::string fileName){

  std::ifstream fileStream(fileName);
  std::string line;
  std::vector<std::vector<double>> nodes;
  while (getline (fileStream, line)) {
    if (line.compare(0, 6, "POINTS") == 0){
      std::cout << fileName << std::endl;
      std::stringstream ss(line);
      std::string dataType;
      int dataRows;
      std::string dataFormat;
      ss >> dataType;
      ss >> dataRows;
      ss >> dataFormat;
      nodes.reserve(dataRows);
      std::stringstream nodess;
      for(int i=0;i<dataRows;++i){
        getline (fileStream, line);
        nodess.clear();
        nodess.str(line);
        std::vector<double> nodeRow(3);
        nodess >> nodeRow[0] >>  nodeRow[1] >> nodeRow[2];   
        nodes.push_back(nodeRow);   
      }
      break;
    }
  }
  fileStream.close();
  return nodes;
};

std::vector<std::vector<size_t>> iovtk::read_elements(std::string fileName){

  std::ifstream fileStream(fileName);
  std::string line;
  std::vector<std::vector<size_t>> elements;
  while (getline (fileStream, line)) {
    if(line.compare(0, 5, "CELLS") == 0){
      std::stringstream ss(line);
      std::string dataType;
      int dataRows;
      int dataCols;
      int rowCols;
      
      ss >> dataType;
      ss >> dataRows;
      ss >> dataCols;
      elements.reserve(dataRows);
      for(int i=0;i<dataRows;++i){
        getline (fileStream, line);
        std::stringstream elementss(line);
        elementss >> rowCols;
        std::vector<size_t> elementRow(rowCols);
        for (int j=0;j<rowCols;++j){
          elementss >> elementRow[j];
        }
        elements.push_back(elementRow);
      }
      break;
    }
  }
  fileStream.close();
  return elements;
};

std::vector<size_t> iovtk::readScalarData(std::string fileName, std::string dataName){
  std::ifstream fileStream(fileName);
  int skiprows = this->skiprows[dataName];
  std::string line;
  std::vector<size_t> data;
  while (getline (fileStream, line)) {
    if(line.compare(0, dataName.length(), dataName) == 0){
      std::stringstream ss(line);
      std::string dataType;
      int dataRows;
      ss >> dataType;
      ss >> dataRows;
      data.resize(dataRows);
      for(int k=0;k<skiprows;++k){
        fileStream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      }
      for(int i=0;i<dataRows;++i){
        getline(fileStream, line);
        std::stringstream edss(line);
        edss >> data[i];     
      }
    }
  }
  fileStream.close();
  return data;
};