@echo off

set parent_path=%~dp0
set APP=%parent_path%..\cmake_build\cudaSirecon\cudasirecon.exe

set MAKEOTF=%parent_path%..\cmake_build\otf\makeotf.exe
%MAKEOTF% %parent_path%psf.dv %parent_path%otf.otf -angle -0.804300 -ls 0.2035 -na 1.4 -nimm 1.515 -fixorigin 3 20

@REM REM test TIFF
@REM %APP% %parent_path% raw %parent_path%otf.tif -c %parent_path%config-tiff

@REM REM test MRC
%APP% %parent_path%raw.dv %parent_path%proc.dv %parent_path%otf.otf -c %parent_path%config