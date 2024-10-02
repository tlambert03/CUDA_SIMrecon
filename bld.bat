@REM rmdir /S /Q cmake_build

call "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build\vcvars64.bat"

@REM mkdir cmake_build
cd cmake_build

@REM cmake -G Ninja ^
@REM     -DBUILD_OTF_VIEWER=OFF ^
@REM     -DCMAKE_BUILD_TYPE=Release ^
@REM     -DCMAKE_INSTALL_PREFIX="%CONDA_PREFIX%/Library" ^
@REM     ../src

ninja
@REM ninja install
