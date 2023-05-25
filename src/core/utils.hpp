#include<vector>
#include<iostream>


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

} // end namespace utils