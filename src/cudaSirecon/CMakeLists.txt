include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  "${CMAKE_CURRENT_SOURCE_DIR}/../Buffers")

find_package(LAPACK REQUIRED)
message(STATUS "LEG_STDIO: " ${LEG_STDIO})
message(STATUS "LAPACK_FOUND: " ${LAPACK_FOUND})

# libcudasirecon shared library
add_library(cudasirecon SHARED cudaSirecon.cpp boostfs.cpp gpuFunctionsImpl.cu interface.cpp)
cuda_add_cufft_to_target(cudasirecon)
set_target_properties(cudasirecon
  PROPERTIES
  MSVC_RUNTIME_LIBRARY ""
  WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
target_link_libraries(
  cudasirecon
  Buffer
  ${TIFF_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${BOOST_LIBRARIES}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${Boost_FILESYSTEM_LIBRARY}
)
  
  
# cudasirecon executable file
add_executable(cudaSireconDriver cudaSireconDriver.cpp)
cuda_add_cufft_to_target(cudaSireconDriver)
target_link_libraries(cudaSireconDriver cudasirecon)
target_include_directories(cudaSireconDriver PUBLIC "${PROJECT_BINARY_DIR}")
set_target_properties(cudaSireconDriver
  PROPERTIES
  MSVC_RUNTIME_LIBRARY ""
  OUTPUT_NAME cudasirecon)


if(WIN32)
  if(${MSVC_VERSION} GREATER 1800)
    target_link_libraries(cudasirecon ${LEG_STDIO})
    target_link_libraries(cudaSireconDriver ${LEG_STDIO})
  endif()
else()
  target_link_libraries(cudasirecon rt)
endif(WIN32)

# added for make install to work in conda
set(HEADERS cudaSirecon.h cudaSireconImpl.h gpuFunctions.h)

install(
  TARGETS cudaSireconDriver cudasirecon
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
  
install(FILES ${HEADERS} DESTINATION include)
  