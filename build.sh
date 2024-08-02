rm -rf cmake_build

mkdir cmake_build
cd cmake_build

CXXFLAGS="$CXXFLAGS -Wfatal-errors -Wno-deprecated-declarations"
cmake ${CMAKE_ARGS} \
    -DBUILD_OTF_VIEWER=OFF \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX \
    ../src

make -j 4
# make install
