#include "io.hpp"

iovtk iovtk::readvtk(std::string fileName){
  iovtk* vtkobj = new iovtk;
  std::ifstream fileStream(fileName);
  if (!fileStream)
	{
		std::cerr<<"Error: Invalid file name\n";
		exit(1);
	}
  fileStream.close();
  vtkobj->nodes = vtkobj->read_nodes(fileName);
  vtkobj->elements = vtkobj->read_elements(fileName);

  std::vector<std::string> dataID = {"CELL_TYPES", "CELL_DATA"};
  for(int i=0;i<dataID.size();++i){
    std::vector<size_t> data = vtkobj->readScalarData(fileName, dataID[i]);
    if (~data.empty()){
      vtkobj->dataID.push_back(dataID[i]);
      vtkobj->data[dataID[i]] = data;
      
    }
  }

  vtkobj->numNodes = vtkobj->nodes.size();
  vtkobj->numElements = vtkobj->elements.size();
  return *vtkobj;
};

std::vector<std::vector<double>> iovtk::read_nodes(std::string fileName){

  std::ifstream fileStream(fileName);
  std::string line;
  std::vector<std::vector<double>> nodes;
  while (getline (fileStream, line)) {
    if (line.compare(0, 6, "POINTS") == 0){
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

std::vector<std::vector<size_t>> iovtk::get_topology(int cellType){
  
  int type = cellTypeMap[cellType];
  std::vector<size_t> elementType = this->data["CELL_TYPES"];
  std::vector<std::vector<size_t>> elements = this->elements;
  std::vector<std::vector<size_t>> elem;
  elem.reserve(elementType.size());
  for(int i=0;i<elementType.size();++i){
    if(elementType[i] == type){
      elem.push_back(elements[i]);
    }
  }
  return elem;
}

std::vector<size_t> iovtk::get_data(std::string cellID, int cellType){
  int type = cellTypeMap[cellType];
  std::vector<size_t> typeData = this->data["CELL_TYPES"];
  std::vector<size_t> fullData = this->data[cellID];
  std::vector<size_t> dat;
  dat.reserve(fullData.size());
  for(int i=0;i<fullData.size();++i){
    if(typeData[i] == type){
      dat.push_back(fullData[i]);
    }
  }
  return dat;
}

iovtk::iovtk(std::vector<std::vector<double>> nodeData, std::vector<std::vector<size_t>> elementData, element* elem){
  this->nodes = nodeData;
  this->elements = elementData;
  this->numNodes = nodeData.size();
  this->numElements = elementData.size();

  std::vector<size_t> elemData(elementData.size(), this->cellTypeMap[elem->get_elemType(2)]);
  this->data["CELL_TYPES"] = elemData;
  this->dataID = {"CELL_TYPES"};
}

void iovtk::writeVtk(std::string fileName, const double solution[]){
    std::ofstream file_vtk;
    file_vtk.open(fileName);

    Eigen::IOFormat HeavyFmt(Eigen::FullPrecision);

    if (!file_vtk){
    	std::cerr<<"Not able to write data to vtk"<<std::endl;
    	exit (-403);
    }

    file_vtk << "# vtk DataFile Version 2.0" << "\n";
    file_vtk << "Displacement results" << "\n";
    file_vtk << "ASCII" << std::endl;
    file_vtk << "DATASET UNSTRUCTURED_GRID" << "\n";
    file_vtk << "POINTS "<< this->numNodes <<" double" << "\n";
    for (unsigned int i=0; i< this->numNodes; i++){
	    file_vtk<<std::fixed << std::setprecision(8)<< this->nodes[i][0] << " " << this->nodes[i][1] << " " << this->nodes[i][2] << "\n";
    }

    file_vtk << "CELLS "<<numElements<<" "<< numElements*4 <<"\n";

    for (unsigned int i = 0; i < this->numElements; i++){
        file_vtk << 3;
        for (int j = 0, nj = elements[i].size(); j < nj; j++){
            file_vtk << " " << elements[i][j];
        }
        file_vtk << "\n";
    }
    auto cellTypeData = data["CELL_TYPES"];
    file_vtk << "CELL_TYPES "<<this->numElements<<std::endl;
    for (unsigned int i = 0; i < this->numElements; i++){
      file_vtk << cellTypeData[i] << "\n";
    }

    file_vtk << "POINT_DATA "<< this->numNodes <<"\n";

    file_vtk << "SCALARS Displacement double 2" << "\n";
    file_vtk << "LOOKUP_TABLE default" << "\n";
    for (unsigned int i=0; i< this->numNodes; i++){
      file_vtk << solution[i*2] << " " << solution[i*2 + 1] << "\n";
    }

    file_vtk.close();
}
