//! Requires C++11

#define cimg_use_tiff
#define cimg_use_cpp11 1
#include "CImg.h"
using namespace cimg_library;

#include <iostream>
#include <string>
#include <IMInclude.h>

int main(int argc, char *argv[])
{

  bool bRotatedFull = false;
  float gamma = 0.4;

  if (argc < 2) {
    std::cerr << "Need at least one argument" <<
      "Usage: otfviewer  [-g gamma_value] [-a] OTF_file_name\n" <<
        "         -g: gamma_value (default is 0.4)\n" <<
        "         -a: to display OTF in full F space with vertical axis being kz\n";
    return 0;
  }

  int flags, opt;
  while ((opt = getopt(argc, argv, "ag:")) != -1) {
    switch (opt) {
    case 'a':
      bRotatedFull = true;
      break;
    case 'g':
      gamma = std::stod(optarg);
      break;
    default:
      std::cerr << "\nUsage: otfviewer [-g gamma_value] [-a] OTF_file_name \n" <<
        "         -g: gamma_value (default is 0.4)\n" <<
        "         -a: to display OTF in full F space with vertical axis being kz\n";
      return 0;
    }
  }
  TIFFSetWarningHandler(NULL);
  TIFFSetErrorHandler(NULL);
  TIFF *tf = TIFFOpen(argv[optind], "r");
  CImg<> atiff;
  if (tf) {
    TIFFClose(tf);
    atiff.assign(argv[optind]);
  }
  else {
    std::cout << "Not a TIFF file; now try MRC\n";
    int istream_no = 1;
    // Suppress IM header printout; somehow only in this mode would the rest of
    // IM calls go well on Windows.
    // IMAlPrt(0);  // not needed with our dvfile implementation

    if (IMOpen(istream_no, argv[optind], "ro")) {
      std::cerr << "File " << argv[optind] << " cannot be opened.\n";
      return -1;
    }
    int ixyz[3], mxyz[3], nxyzst[3], pixeltype;
    float mmm[3];
    IW_MRC_HEADER header;
    
    IMRdHdr(istream_no, ixyz, mxyz, &pixeltype, mmm, mmm+1, mmm+2);
    IMGetHdr(istream_no, &header);

    atiff.assign(header.nx * 2, header.ny, header.nz);
    for (int z=0; z<header.nz; z++)
      IMRdSec(istream_no, atiff.data(0, 0, z));
    IMClose(istream_no);
  }

  CImg<> amplitudes(atiff.width()/2, atiff.height(), atiff.depth());
  cimg_forXYZ(amplitudes, x, y, z) {
    float r = atiff(2*x, y, z);
    float i = atiff(2*x+1, y, z);
    amplitudes(x, y, z) = pow(std::sqrt(r*r + i*i), gamma);
  }

  if (bRotatedFull) {
    CImg<> ampAlt(amplitudes.height()*2-1, amplitudes.width(), amplitudes.depth());
    cimg_forXYZ(amplitudes, x, y, z) {
      int kz = x;
      if (x > amplitudes.width()/2) kz -= amplitudes.width(); // origin at corners
      ampAlt(y+amplitudes.height(), kz+amplitudes.width()/2, z) = amplitudes(x, y, z);
      ampAlt(amplitudes.height()-y, kz+amplitudes.width()/2, z) = amplitudes(x, y, z);
    }

    amplitudes = ampAlt;
  }

  amplitudes.display(argv[optind], false);

  // Somehow doubled XY image looks the same:
  // CImg<> large = amplitudes.resize_doubleXY();
  // std::cout << amplitudes.width() << "  " << large.width() << std::endl;
  // large.display(argv[1], false);
  return 0;
}
