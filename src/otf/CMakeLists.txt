include_directories("${CMAKE_SOURCE_DIR}/cudaSirecon") # for CImg.h
add_executable(makeotf radialft.cpp)

# ######## START findFFTW: https://github.com/egpbos/findFFTW ######

configure_file(downloadFindFFTW.cmake.in findFFTW-download/CMakeLists.txt)
execute_process(
  COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-download)
if(result)
  message(FATAL_ERROR "CMake step for findFFTW failed: ${result}")
else()
  message("CMake step for findFFTW completed (${result}).")
endif()
execute_process(
  COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-download)
if(result)
  message(FATAL_ERROR "Build step for findFFTW failed: ${result}")
endif()

set(findFFTW_DIR ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-src)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${findFFTW_DIR}")

# ######## END findFFTW ##########################################

if(WIN32)
  find_package(FFTW REQUIRED COMPONENTS FLOAT_LIB)
  # include_directories( ${FFTW_INCLUDE_DIRS} )
  if(${MSVC_VERSION} GREATER 1800)
    target_link_libraries(makeotf ${LEG_STDIO})
  endif()
else()
  find_package(FFTW REQUIRED COMPONENTS FLOAT_THREADS_LIB)
endif()

message(STATUS "FFTW3 libraries: " ${FFTW_LIBRARIES})
target_link_libraries(makeotf ${FFTW_LIBRARIES} ${TIFF_LIBRARIES})

install(TARGETS makeotf RUNTIME DESTINATION bin)

if(BUILD_OTF_VIEWER)
  message(STATUS "Build OTF Viewer: TRUE")

  include_directories(${PRIISM_INCLUDE_PATH})
  link_directories(${PRIISM_LIB_PATH})
  add_executable(otfviewer otfviewer.cpp)
  target_link_libraries(otfviewer ${TIFF_LIBRARIES})

  if(WIN32)
    if(${MSVC_VERSION} GREATER 1800)
      target_link_libraries(otfviewer ${LEG_STDIO})
    endif()
  else()
    find_package(X11 REQUIRED)
    target_link_libraries(otfviewer ${X11_LIBRARIES})
  endif()

  install(TARGETS otfviewer RUNTIME DESTINATION bin)
endif()
