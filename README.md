# simpleFEM
To run the application, run the below commands.

```shell
mkdir build
cd build
cmake ..
make
mpirun -n 1 simpleFEM
```

Below libraries need to be installed before you build the application.

eigen3 : sudo apt install libeigen3-dev
petsc : sudo apt-get install petsc-dev