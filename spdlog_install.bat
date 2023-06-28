git clone https://github.com/gabime/spdlog.git

cd spdlog
mkdir build
cd build

cmake -G Ninja -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc ..
Ninja

cd ../../
