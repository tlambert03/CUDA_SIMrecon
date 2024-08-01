/* A new small project to radially average OTF starting from */
/* PSF dataset. Mainly for structured illumination microscopy ( */
/* 1 or 2 objectives). Also take into account the following:    */
/*  compensation for finite bead size; */
/*  the phase factor associated with the side bands; */
/*  estimate the sub-pixel position of the bead by fitting parabolas; */
/*  cleanup out-of-band noises. */

#ifdef _WIN32
#define _USE_MATH_DEFINES
#endif
#include <cmath>

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cassert>

#define cimg_use_tiff

#include <CImg.h>
using namespace cimg_library;

#include "dvfile.h"

#include <complex>
#include <fftw3.h>

#define DEFAULTNPHASES 5

 
float fitparabola(float a1, float a2, float a3);
float estimate_background(CImg<> &image, int border_size);
void determine_center(std::vector< CImg<> > &stack5phases, CImg<> &I2Mstack, float *xc, float *yc, float *zc, float *xc_i2m, float *yc_i2m, float *zc_i2m, bool two_lens, bool I2M_inc);
void apodize(int napodize,int nx, int nxExtra, int ny, float *image);
void cosapodize(int nx, int nxExtra, int ny, float *image);
void makematrix (CImg<float> &sepMatrix);
void separate(std::vector<CImg<>> &floatimage, CImg<> &sepMatrix);

void combine_reim(std::vector<std::complex<float> *> &bands, int norders, int nx, int nz, int);
void shift_center(std::complex<float> *bands, int nx, int ny, int nz, float xc, float yc, float zc);

void beadsize_compensate(std::vector<std::complex<float> *> &bands, float k0angle,
                         float linespacing, float bead_diameter, int norders, 
                         int nx, int ny, int nz, float dr, float dz);
double sphereFFT(float k, float radius);

void radialft(std::complex<float> *bands, int nx, int ny, int nz, float dr, std::complex<float> *avg, unsigned nr);
void cleanup(std::complex<float> *otfkxkz, int order, int nx, int nz, float dr, float dz,
             float linespacing, unsigned lambdanm, int icleanup, int ileave1_1, int ileave1_2,
             int ileave2, int twolens, float NA, float NIMM);
void cleanup_I2M(std::complex<float> *otfkxkz, int nx, int nz, float dkr, float dkz, int lamdanm, int icleanup, float NA, float NIMM);
void modify(std::complex<float> *otfkxkz, int nx, int nz, int *ifixz, int ifixr, int order, int twolens);
void fixorigin(std::complex<float> *otfkxkz, int nx, int nz, int kx1, int kx2);
void rescale(std::complex<float> *otfkxkz, int order, int nx, int nz, float *scalefactor, int dorescale);

void outputdata(std::string &tiffile, std::vector<std::complex<float> *> &bands, int norders,
                int nx, int ny, int nz, float dkr, float dkz, int five_bands);

void outputdata(int ostream_no, IW_MRC_HEADER *header, std::vector<std::complex<float> *> &bands,
                int norders, int nx, int ny, int nz, float dkr, float dkz, int five_bands);
void mrc_file_write(float *buffer, int nx, int ny, int nz, float rlen, float zlen, int mode, int iwave, const char *files);
int commandline(int argc, char *argv[], int * twolens, int *rescale, float *beaddiam, float *k0angle, float *linespacing, int *five_bands, int *nphases, int *interpkr, int *leavekz, int *do_compen, int *I2M_inc, std::string &I2Mfiles, float *background, int *bBgInExtHdr, int *order0gen, std::string &order0files, int *conjugate, float *na, float *nimm, int *ifixz, int *ifixr, unsigned *wavelength, float *dr, float *dz, int *bCoherentBSIM, int *bForcedPIshift, std::string &ifiles, std::string &ofiles, int *ifilein, int *ofilein);

#if defined(_MSC_VER) && (_MSC_VER >= 1900)
// for IMLIB with vs >2015
// https://stackoverflow.com/questions/30412951/unresolved-external-symbol-imp-fprintf-and-imp-iob-func-sdl2
FILE _iob[] = {*stdin, *stdout, *stderr};

extern "C" FILE *  __iob_func(void)
{
    return _iob;
}

#endif

int main(int argc, char **argv)
{
  std::string ifiles, ofiles, I2Mfiles, order0files;
  int nphases=DEFAULTNPHASES, napodize=10;
  int border_size;
  float dr=0.106f, dz=0.1f, background_dummy=-1.f, background_i2m, xcofm, ycofm, zcofm, xcofm_i2m, ycofm_i2m, zcofm_i2m;
  float bead_diameter=0.12f, scalefactor=1.0f;
  float k0angleguess=1.57f, linespacing = 0.2f ; 
  int ifixz[3], ifixr, twolens=0, dorescale=0, ifilein=0, ofilein=0, five_bands=0, do_compen=1, I2M_inc=0, Generate_band0=0, conjugate=0;
  int bBgInExtHdr = 0; /* if the background of each section is recorded in the extended header's 3rd float number (in Lin's scope) */
  int icleanup, ileavekz[3]={0,0,0}, ileave1_1, ileave1_2, ileave2, interpkr[2];
  float NA=1.4f, NIMM=1.515f;
  float user_dr=0, user_dz=0; /* will use these for dr, dz if these are greater than zero */
  int bCoherentBSIM = 0; /* In coherent Bessel-SIM, do fixorigin() differently? */
  int bForcedPIshift=0;

  IW_MRC_HEADER header, otfheader;
  IMAlPrt(0);       


  interpkr[0] = 0;
  interpkr[1] = 0;
  ifixz[0] = 1;ifixz[1]=1; ifixz[2]=1; ifixr = 0;

  unsigned wavelength = 525;

  if (!commandline(argc, argv, &twolens, &dorescale, &bead_diameter, &k0angleguess, &linespacing, &five_bands, &nphases, interpkr, ileavekz, &do_compen, &I2M_inc, I2Mfiles, &background_dummy, &bBgInExtHdr, &Generate_band0, order0files, &conjugate, &NA, &NIMM, ifixz, &ifixr, &wavelength, &user_dr, &user_dz, &bCoherentBSIM, &bForcedPIshift, ifiles, ofiles, &ifilein, &ofilein))
    return -1;

  
  printf("ileavekz[0]=%d, ifixz[0]=%d\n", ileavekz[0], ifixz[0]);

  if (!ifilein) {
    std::cout << " PSF dataset file name: ";
    getline(std::cin, ifiles);
  }

  int istream_no=1, ostream_no=2; //! For MRC file I/O
  // First try with TIFF; if unsuccessful, fall back to MRC
  bool bTIFF = false;
  TIFFSetWarningHandler(NULL);
  TIFFSetErrorHandler(NULL);
  TIFF *tf = TIFFOpen(ifiles.c_str(), "r");
  CImg <> PSFtiff;
  if (tf) {
    TIFFClose(tf);
    PSFtiff.assign(ifiles.c_str());
    bTIFF = true;
  }
  else {
    std::cout << "Not a TIFF file; now try MRC\n";
    if (IMOpen(istream_no, ifiles.c_str(), "ro")) {
      std::cerr << "File " << ifiles << " does not exist.\n";
      return -1;
    }
  }
  if (!ofilein) {
    std::cout << "Output OTF file name: ";
    getline(std::cin, ofiles);
  }
  // if (!bTIFF)
  //   if (IMOpen(ostream_no, ofiles.c_str(), "new")) {
  //     std::cerr << "File " << ofiles << " cannot be created.\n";
  //     return -1;
  //   }

  unsigned nx, nxExtra, ny, nz, nxy;
  if (bTIFF) {
    nx = PSFtiff.width();
    ny = PSFtiff.height();
    nz = PSFtiff.depth();
    if (user_dr>0)
      dr = user_dr;  // otherwise, use the default set at declaration; same for "dz"
    if (user_dz>0)
      dz = user_dz;
  }

  else {
    int ixyz[3], mxyz[3], pixeltype;
    float min, max, mean;
    IMRdHdr(istream_no, ixyz, mxyz, &pixeltype, &min, &max, &mean);
    IMGetHdr(istream_no, &header);

    nx=header.nx;
    ny=header.ny;
    nz=header.nz;
    dr = header.ylen;
    dz = header.zlen;
    wavelength = header.iwav1;
  }

  if (!I2M_inc)
    nz /= nphases;
  else 
    nz /= nphases+1;

  printf("nx=%u, ny=%u, nz=%u, dr=%f, dz=%f, wavelength=%u\n", nx, ny, nz, dr, dz, wavelength);

  float dkx = 1/(nx*dr);
  float dky = 1/(ny*dr);
  float dkz = 1/(nz*dz);

  int norders = (nphases+1)/2;

  if (nx % 2)
    nxExtra = nx+1;
  else
    nxExtra = nx+2;
  nxy = nxExtra * ny;

  std::vector<CImg<>> floatimage(nphases);
  std::vector<std::complex<float>*> bands(nphases);
  
  for(int i=0; i<nphases; i++) {
    floatimage[i].assign(nxExtra, ny, nz);
    bands[i] = (std::complex<float> *) floatimage[i].data();  /* reuse same allocation to save memory */
  }

  CImg<> I2M_image;
  if (I2M_inc)
    I2M_image.assign(nxExtra, ny, nz);

  float *background = new float[nz*nphases];
  border_size = 20;
  
  printf("Reading data...\n\n");
  int zsec = 0;  // incremented by 1 at the end of "for (phase)" loop
  CImg<> buffer;
  for(auto z=0u; z<nz; z++) {
    for(int phase=0; phase<nphases; phase++) {
      if (bTIFF)
        buffer.assign(PSFtiff.data(0, 0, zsec), nx, ny, 1, 1, true); // Reuse buffer from "PSFtiff"; check
      else {
        buffer.assign(nx, ny);
        IMRdSec(istream_no, buffer.data());
      }
      // Fix camera error at even binnings for some very old Photometric CCD
      if(buffer(nx-1, ny-1) < 0)
        buffer(nx-1, ny-1) = buffer(nx-2, ny-1);

      if (bTIFF) {
        if (background_dummy >=0)
          background[zsec] = background_dummy;
        else
          background[zsec] = estimate_background(buffer, border_size);
      }
      else if (bBgInExtHdr) {
        // Some Andor EMCCDs have hidden border pixels that can be used to report fluctuating dark
        // pixel values, which can be saved in MRC's per-exposure extra header info
        int extInts;
        float extFloats[3];
        IMRtExHdrZWT(istream_no, zsec, 0, 0, &extInts, extFloats);
        background[zsec] = extFloats[2];
      }


      for(auto i=0u; i<ny; i++) {
        // Copy from "buffer" to floatimage[phase] at section z row-by-row
        std::memcpy(floatimage[phase].data(0, i, z), buffer.data(0, i), nx*sizeof(float));
        for (auto j=0u; j<nxExtra - nx; j++)
          floatimage[phase](nx+j, i, z) = 0.;
      }

      // To-Do: why not apodize "buffer" first before copying to "floatimage[phase]"?
      float *floatsection = floatimage[phase].data(0, 0, z);
      if(napodize>=0) 
        apodize(napodize, nx, nxExtra, ny, floatsection);
      else if(napodize==-1)
        cosapodize(nx, nxExtra, ny, floatsection);
      
      zsec ++;
    } // end for(phase)

    // To-Do; but highly unlikely we'll ever need I2M OTF
    if (I2M_inc) {  // Data contains one extra section of I2M image
      if (bTIFF) {
      }
      else 
        IMRdSec(istream_no, buffer.data());
      if(buffer(nx-1, ny-1) == -1.0)  /* fix camera error at even binnings */
        buffer(nx-1, ny-1) = buffer(nx-2, ny-1);
      if (z==nz/2)
        background_i2m = estimate_background(buffer, border_size);
      for(auto i=0u; i<ny; i++) {
        // Copy from "buffer" to I2M_imag at section z row-by-row
        std::memcpy(I2M_image.data(0, i, z), buffer.data(0, i), nx*sizeof(float));
        for (auto j=0u; j<nxExtra - nx; j++)
          I2M_image(nx+j, i, z) = 0.;
      }
    }
  } // end for(z)
  
  //! Before FFT, use center band to estimate bead center position:
  determine_center(floatimage, I2M_image, &xcofm, &ycofm, &zcofm,
                   &xcofm_i2m, &ycofm_i2m, &zcofm_i2m, twolens, I2M_inc);

  printf("Center of mass is (%.3f, %.3f, %.3f)\n\n", xcofm, ycofm, zcofm);

  if (I2M_inc) {
    // printf("I2M psf's background is %.3f\n", background_i2m);
    printf("I2M psf's center of mass is (%.3f, %.3f, %.3f)\n\n", xcofm_i2m, ycofm_i2m, zcofm_i2m);
  }

#pragma omp parallel for
  for(int z=0; z<nz; z++) {
    for(int phase=0; phase<nphases; phase++) {
      CImg<> oneSection(floatimage[phase].data(0, 0, z), nxExtra, ny, 1, 1, true);
      oneSection -= background[z*nphases + phase]; // This alters floatimage because of shared buffer
    }
  }
  if (I2M_inc) 
    I2M_image -= background_i2m;

  CImg<> sepMatrix(nphases, nphases);

  if (nphases > 1) {
    makematrix(sepMatrix);
    separate(floatimage, sepMatrix);
  }

  if (Generate_band0) {
    if (bTIFF)
      floatimage[0].save_tiff(order0files.c_str());
    else
      mrc_file_write(floatimage[0], nxExtra, ny, nz, dr, dz, 0, wavelength, order0files.c_str());
  }

  printf("Before fft\n");
  unsigned nthreads = 1;
#ifndef NDEBUG
#ifdef _WIN32
  SYSTEM_INFO siSysInfo;
  GetSystemInfo(&siSysInfo);
  nthreads = siSysInfo.dwNumberOfProcessors;
#else
  nthreads = sysconf(_SC_NPROCESSORS_CONF); // Get number of processors; only works on Unix-based systems and Cygwin (but not MingW)
#endif
#endif

  printf("Number of threads used: %d\n", nthreads);

  fftwf_plan_with_nthreads(nthreads);
  fftwf_plan rfftplan3d =
    fftwf_plan_dft_r2c_3d(nz, ny, nx, floatimage[0].data(),
                          (fftwf_complex *) floatimage[0].data(), FFTW_ESTIMATE);
  for(int phase=0; phase<nphases; phase++)
    fftwf_execute_dft_r2c(rfftplan3d, floatimage[phase].data(), (fftwf_complex *) floatimage[phase].data());
  if (I2M_inc)
    fftwf_execute_dft_r2c(rfftplan3d, I2M_image.data(), (fftwf_complex *) I2M_image.data());

  fftwf_destroy_plan(rfftplan3d);
  printf("After fft\n\n");

  // Modify the phase of bands, so that it corresponds to FFT of a bead at origin
  std::cout << "Shifting center...\n";
  for (int phase=0; phase<nphases; phase++)
    shift_center(bands[phase], nx, ny, nz, xcofm, ycofm, zcofm);

  if (I2M_inc)
    shift_center((std::complex<float> *)I2M_image.data(), nx, ny, nz, xcofm_i2m, ycofm_i2m, zcofm_i2m);

  // Now compensate finite bead size for all bands
  if (do_compen)
    beadsize_compensate(bands, k0angleguess, linespacing, bead_diameter, norders, nx, ny, nz, dr, dz);

  std::vector<CImg<>> avg_output_buf(nphases);
  std::vector<std::complex<float> *> avg_output(nphases);
  //! smaller of the X-Y dimensions; to be used as the lateral dimension in avg_output
  unsigned nr = std::min(nx, ny);

  for (int phase=0; phase < nphases; phase ++) {
    avg_output_buf[phase].assign(nz, (nr/2+1)*2, 1, 1, 0.f);
    avg_output[phase] = (std::complex<float> *) avg_output_buf[phase].data();
  }

  for (auto order=0; order < norders; order ++) {
    if (order==0)
      radialft(bands[0], nx, ny, nz, dr, avg_output[0], nr);
 
    else {
      radialft(bands[2*order-1], nx, ny, nz, dr, avg_output[2*order-1], nr);
      radialft(bands[2*order], nx, ny, nz, dr, avg_output[2*order], nr);
    }
  }

  std::complex<float> * I2Mavg_output;
  if (I2M_inc) {
    I2Mavg_output = (std::complex<float> *) calloc(nz*(nr/2+1), sizeof(std::complex<float>));
    radialft((std::complex<float> *)I2M_image.data(), nx, ny, nz, dr, I2Mavg_output, nr);
  }

  icleanup = nx/2+1;
  ileave1_1 = ileavekz[0]; ileave1_2 = ileavekz[1]; ileave2 = ileavekz[2];

//  printf("ileavekz[0]=%d, ileave1_1=%d\n", ileavekz[0], ileave1_1);

  for (int phase=0; phase < nphases; phase ++) {
    modify(avg_output[phase], nx, nz, ifixz, ifixr, (phase+1)/2, twolens);
    if (bCoherentBSIM)
      {/* cleanup_CBSIM(); */}
    else if (ileave1_1 > 0) {
      cleanup(avg_output[phase], (phase+1)/2, nr, nz, dr, dz, linespacing, wavelength, icleanup, ileave1_1, ileave1_2, ileave2, twolens, NA, NIMM);
    }
    if (phase==0 && interpkr[0] > 0)
      fixorigin(avg_output[phase], nx, nz, interpkr[0], interpkr[1]);
    rescale(avg_output[phase], (phase+1)/2, nx, nz, &scalefactor, dorescale);
  }

  /* For side bands, combine bandre's and bandim's into bandplus */
  /* Shouldn't this be done later on the averaged bands? */
  if (!five_bands && nphases > 1)
    combine_reim(avg_output, norders, nx, nz, bForcedPIshift);

  if (conjugate)
    for (auto phase=0; phase<nphases; phase++)
#pragma omp parallel for
      for (auto i=0; i<(nx/2+1)*nz; i++)
        avg_output[phase][i] = std::conj(avg_output[phase][i]);

  if (bTIFF) {
    outputdata(ofiles, avg_output, norders, nx, ny, nz, dkx, dkz, five_bands);
  }
  else {
    otfheader = header;
    outputdata(ostream_no, &otfheader, avg_output, norders, nx, ny, nz, dkx, dkz, five_bands);
  }
  if (I2M_inc) {
/*     cleanup_I2M(I2Mavg_output, nx, nz, dkr, dkz, wavelength, icleanup, NA, NIMM); */
    fixorigin(I2Mavg_output, nx, nz, interpkr[0], interpkr[1]);
    rescale(I2Mavg_output, 0, nx, nz, &scalefactor, dorescale);
    mrc_file_write((float *)I2Mavg_output, nz, nx/2+1, 1, dkz, dkx, 2, wavelength, I2Mfiles.c_str());
  }

  return 0;
}

/***************************** makematrix ************************************/
/*     generates the matrix that is to be used to separate the raw indata    */
/*     into the different bands of sample information.                       */
/*****************************************************************************/
void makematrix (CImg<float> &sepMatrix)
{
  auto nphases = sepMatrix.width();
  auto norders = (nphases+1)/2;
  auto phi = 2*M_PI / nphases;
  for (auto j=0; j<nphases; j++) {
    sepMatrix(j, 0) = 1.0 / nphases;
     for (auto order=1; order<norders; order++) {
       sepMatrix(j, 2*order-1) = cos(j*order*phi)/nphases;
       sepMatrix(j, 2*order  ) = sin(j*order*phi)/nphases;
     }
  }
}

/***************************** separate **************************************/
/*     Applies image arithmetic to the image sequence floatimage[][],        */
/*     which was acquired with different phases of the illumination          */
/*     pattern, to separate the different bands of sample information.       */
/*     The coefficients are pre-stored in the matrix sepMatrix.              */
/*     The bands are returned in the same array floatimage where the         */
/*     input data was.                                                       */
/*****************************************************************************/
void separate(std::vector<CImg<>> &floatimage, CImg<> &sepMatrix)
{
  auto nphases = floatimage.size();

  //! Big temp buffer to optimize for speed:
  std::vector<CImg<>> output(nphases);
  std::vector<float *> output_ptr(nphases), floatimage_ptr(nphases);
  for (int i=0; i<nphases; i++) {
    output[i].assign(floatimage[i], "xyzc");
    output_ptr[i] = output[i].data();
    floatimage_ptr[i] = floatimage[i].data();
  }

#pragma omp parallel for
  for (int x=0; x<floatimage[0].size(); x++) {
    for (int i=0; i<nphases; i++) {
      output_ptr[i][x] = 0.0;
      for(int j=0; j<nphases; j++)
        output_ptr[i][x] += floatimage_ptr[j][x] * sepMatrix(j, i);
    }
    for(int i=0; i<nphases; i++)
      floatimage_ptr[i][x] = output_ptr[i][x];
  }
}

float estimate_background(CImg<> &image, int border_size)
{
  unsigned total=0;
  double sum = 0.;

  size_t nx = image.width();
  size_t ny = image.height();
#pragma omp parallel for shared(image) reduction(+: sum, total)
  cimg_forXY(image, x, y) {
    if (y<border_size || y>ny-border_size ||
        x<border_size || x>nx-border_size) {
      sum += image(x, y);
      total ++;
    }
  }
  return (sum/total);
}


/*  locate peak pixel to subpixel accuracy by fitting parabolas  */
void determine_center(std::vector< CImg<> > & stackNphases, CImg<> &I2M_image, float *xc, float *yc, float *zc, float *xc_i2m, float *yc_i2m, float *zc_i2m, bool bTwolens, bool bI2M)
{
  std::cout << "In determine_center()\n";
  int nx = stackNphases[0].width();
  int ny = stackNphases[0].height();
  int nz = stackNphases[0].depth();

  unsigned nphases = stackNphases.size();
  CImg<> stack3d(stackNphases[0], "xyzc", 0); // same dimension as a single-phase stack
  for (auto p=0u; p<nphases; p++)
    stack3d += stackNphases[p];
  stack3d /= nphases;

  // Search for the peak pixel
  // Be aware that stack3d is of dimension nxExtra*ny*nz
  float maxval=0.0, reval;
  int maxi=0, maxj=0, maxk=0;

  cimg_forXYZ(stack3d, x, y, z) {
    reval=stack3d(x, y, z);
    if( reval > maxval ) {
      maxval = reval;
      maxi=y; maxj=x;
      maxk=z;
    }
  }

  int iminus = maxi-1;
  int iplus = maxi+1;
  if( iminus<0 ) iminus+=ny;
  if( iplus>=ny ) iplus-=ny;
  int jminus = maxj-1;
  int jplus = maxj+1;
  if( jminus<0 ) jminus+=nx;
  if( jplus>=nx ) jplus-=nx;
  int kminus = maxk-1;
  int kplus = maxk+1;
  if( kminus<0 ) kminus+=nz;
  if( kplus>=nz ) kplus-=nz;

  float valminus,valplus; 
  if (nz > 1) {
    valminus = stack3d(maxj, maxi, kminus);
    valplus  = stack3d(maxj, maxi, kplus);
    *zc = maxk + fitparabola(valminus, maxval, valplus);
  }
  else
    *zc = maxk;
  
#ifdef ACCURATE_PEAK
  if (!twolens)
    *zc = maxk + fitparabola(valminus, maxval, valplus);
  else {
    float zcs[5], peakval[5], minval_right, minval_left;
    int n_samples = 5, kmin, kminus2, kplus2;
    float Xs[5], Ys[5];

    kminus2 = maxk-n_samples/2; kplus2 = maxk+n_samples/2;
    for (i=0; i<n_samples; i++) {
      Xs[i] = kminus2 + i;
      Ys[i] = stack3d[(kminus2+i)*nxy2 + maxi*nxExtra + maxj];
    }
    zcs[0] = fitparabola_Nsample_points(Xs, Ys, n_samples, peakval);

    minval_right = Ys[n_samples-1];
    minval_left = Ys[0];
    
    /* search forward for the next trough */
    kmin = kplus2;
    for (k=kplus2+1; k<=kplus2+5; k++)
      if (stack3d[k*nxy2+maxi*nxExtra+maxj] < minval_right) {
    minval_right = stack3d[k*nxy2+maxi*nxExtra+maxj];
    kmin = k;
      }
    for (i=0; i<n_samples; i++) {
      Xs[i] = kmin + 1 + i;
      Ys[i] = stack3d[(kmin + 1 + i)*nxy2 + maxi*nxExtra + maxj];
    }
    zcs[1] = fitparabola_Nsample_points(Xs, Ys, n_samples, peakval+1);

    /* search backward for the previous trough */
    kmin = kminus2;
    for (k=kminus2-1; k>=kminus-5; k--)
      if (stack3d[k*nxy2+maxi*nxExtra+maxj] < minval_left) {
        minval_left = stack3d[k*nxy2+maxi*nxExtra+maxj];
        kmin = k;
      }
    for (i=0; i<n_samples; i++) {
      Xs[i] = kmin - 1 - i;
      Ys[i] = stack3d[(kmin - 1 - i)*nxy2 + maxi*nxExtra + maxj];
    }
    zcs[2] = fitparabola_Nsample_points(Xs, Ys, n_samples, peakval+2);

    printf("zc0=%.3f, zc_plus=%.3f, zc_minus=%3f\n", zcs[0], zcs[1], zcs[2]);

    zcs[3] = 94; peakval[3] = 251;
    zcs[4] = 68; peakval[4] = 304;
    *zc = fitparabola_Nsample_points(zcs, peakval, 3, peakval);
 }
#endif 

  valminus = stack3d(maxj, iminus, maxk);
  valplus  = stack3d(maxj, iplus, maxk);
  *yc = maxi + fitparabola(valminus, maxval, valplus);

  valminus = stack3d(jminus, maxi, maxk);
  valplus  = stack3d(jplus, maxi, maxk);
  *xc = maxj + fitparabola(valminus, maxval, valplus);
  
  // sum = 0;
  // infocus_sec = floor(*zc);
  // for (k=0; k<nphases; k++)
  //   for (i=0; i<*yc-20; i++)
  //     for (j=0; j<nx; j++)
  //   sum += stackNphases[k][infocus_sec*nxy2 + i*nxExtra + j];
  // *background = sum / (nphases*(*yc-20)*nx);

  if (bI2M) {
    // sum = 0;
    // for (i=0; i<*yc-20; i++)
    //   for (j=0; j<nx; j++)
    // sum += I2M_image[infocus_sec*nxy2 + i*nxExtra + j];
    // *background_i2m = sum / (nx*(*yc-20));

    /* Search for the peak pixel */
    /* Be aware that I2M_image is of dimension nxExtra*ny*nz */
    maxval=0.0;
    cimg_forXYZ(I2M_image, x, y, z) {
      reval=I2M_image(x, y, z);
      if( reval > maxval ) {
        maxval = reval;
        maxi=y; maxj=x;
        maxk=z;
      }
    }

    iminus = maxi-1; iplus = maxi+1;
    if( iminus<0 ) iminus+=ny;
    if( iplus>=ny ) iplus-=ny;
    jminus = maxj-1; jplus = maxj+1;
    if( jminus<0 ) jminus+=nx;
    if( jplus>=nx ) jplus-=nx;
    kminus = maxk-1; kplus = maxk+1;
    if( kminus<0 ) kminus+=nz;
    if( kplus>=nz ) kplus-=nz;

    valminus = I2M_image(maxj, maxi, kminus);
    valplus  = I2M_image(maxj, maxi, kplus);
    *zc_i2m = maxk + fitparabola(valminus, maxval, valplus);

    valminus = I2M_image(maxj, iminus, maxk);
    valplus  = I2M_image(maxj, iplus, maxk);
    *yc_i2m = maxi + fitparabola(valminus, maxval, valplus);

    valminus = I2M_image(jminus, maxi, maxk);
    valplus  = I2M_image(jplus, maxi, maxk);
    *xc_i2m = maxj + fitparabola(valminus, maxval, valplus);
  }
}

/***************************** fitparabola **********************************/
/*     Fits a parabola to the three points (-1,a1), (0,a2), and (1,a3).     */
/*     Returns the x-value of the max (or min) of the parabola.             */
/****************************************************************************/

float fitparabola( float a1, float a2, float a3 )
{
 float slope,curve,peak;

 slope = 0.5* (a3-a1);         /* the slope at (x=0). */
 curve = (a3+a1) - 2*a2;       /* (a3-a2)-(a2-a1). The change in slope per unit of x. */
 if( curve == 0 ) 
 {
   printf("no peak: a1=%f, a2=%f, a3=%f, slope=%f, curvature=%f\n",a1,a2,a3,slope,curve);
   return( 0.0 );
 } 
 peak = -slope/curve;          /* the x value where slope = 0  */
 if( peak>1.5 || peak<-1.5 )
 {
   printf("bad peak position: a1=%f, a2=%f, a3=%f, slope=%f, curvature=%f, peak=%f\n",a1,a2,a3,slope,curve,peak);
   return( 0.0 );
 } 
 return( peak );
}

/*********************************   apodize   **********************************************/
/*   softens the edges of a singe xy section to reduce edge artifacts and improve the fits  */
/********************************************************************************************/
void apodize(int napodize,int nx, int nxExtra, int ny,float *image)
{
 float diff,fact;
 int k,l;

  for(k=0;k<nx;k++)
  {
    diff = (image[ (ny-1)*nxExtra +k ] - image[ /* 0*nxExtra+ */ k ]) / 2;
    for(l=0;l<napodize;l++)
    {
      fact = 1 - sin((((float)l+0.5)/napodize)*M_PI*0.5);
      image[l*nxExtra + k] += diff*fact;
      image[(ny-1-l)*nxExtra + k] -= diff*fact;
    }
  }
  for(l=0;l<ny;l++)
  {
    diff = (image[ l*nxExtra + nx-1 ] - image[ l*nxExtra /* +0 */ ]) / 2;
    for(k=0;k<napodize;k++)
    {
      fact = 1 - sin((((float)k+0.5)/napodize)*M_PI*0.5);
      image[l*nxExtra + k] += diff*fact;
      image[l*nxExtra + nx-1-k] -= diff*fact;
    }
  }
}

 
/*******************************   cosapodize   ********************************/
/*   softens the edges to reduce edge artifacts and improve the fits           */
/*******************************************************************************/
void cosapodize(int nx,int nxExtra, int ny,float *image)
{
 float xfact,yfact;
 int k,l;

  printf("in cosapodize\n");
  for(k=0;k<nx;k++)
  {
    xfact = sin(M_PI*((float)k+0.5)/nx);
    for(l=0;l<ny;l++)
    {
      yfact = sin(M_PI*((float)l+0.5)/ny);
      image[l*nxExtra + k] *= (xfact*yfact);
    }
  }
}

/* Find out the phase of the side bands relative to the center band    */
/* and get the real OTF of the side bands. It's based on the fact that */
/* bandre=BAND*cos(phi), and bandim=BAND*sin(-phi)                     */
void combine_reim(std::vector<std::complex<float> *> &otf, int norders, int nx, int nz, int bForcedPIshift)
{
  int order, i, nxz;
  double bandre_mag, bandim_mag, phi;

  printf("In combine_reim()\n");
  nxz = (nx/2+1)*nz;

  for (order=1; order<norders; order++) {
    bandre_mag = 0;
    bandim_mag = 0;
    for (i=0; i<nxz; i++) {
      bandre_mag += std::abs(otf[2*order-1][i]);
      bandim_mag += std::abs(otf[2*order][i]);
    }
    phi = atan(bandim_mag/bandre_mag);
    /* phi is now in the first quadrant only, since bandim_mag and bandre_mag are both positive */
    /* which quadrant phi should be in is decided by the kz=0 plane of bandre and bandim values */
    if (otf[2*order-1][nx/4*nz].real() < 0 && otf[2*order][nx/4*nz].real() >0)
      phi += M_PI;
    else if (otf[2*order-1][nx/4*nz].real() < 0 && otf[2*order][nx/4*nz].real() <0)
      phi = M_PI-phi;
    else if (otf[2*order-1][nx/4*nz].real() > 0 && otf[2*order][nx/4*nz].real() >0)
      phi = -phi;

    /* Sometimes the above logic still won't find the correct phase; user can specify additional Pi phase shift*/
    if (order==1 && bForcedPIshift)
      phi = M_PI-phi;

    printf("  phi=%f\n", phi);

    for (i=0; i<nxz; i++) {
      otf[2*order-1][i] = otf[2*order-1][i]* (float) cos(phi) + otf[2*order][i]* (float) sin(-phi);
      otf[2*order][i] = 0;  /* since bandim = band*sin(phi) and phi is essentially 0  */
    }
  }
}

void beadsize_compensate(std::vector<std::complex<float> *> &bands, float k0angle, float linespacing, float bead_diameter, int norders, int nx, int ny, int nz, float dr, float dz)
{
  printf("In beadsize_compensate()\n");
  float dkx = 1./(nx * dr);
  float dky = 1./(ny * dr);
  float dkz = 1./(nz * dz);

  int kycent = ny/2;
  int kxcent = nx/2;
  int kzcent = nz/2;
  int nxy = (nx/2+1)*ny;

  float k0mag = 1/linespacing;
  float radius = bead_diameter * 0.5;   //! the radius of the fluorescent bead

  //! The limit value of the FT of a sphere at the origin of Fourier space
  float limit_at_origin = 4*M_PI*radius*radius*radius/3;
  
  for (auto order=0; order<norders; order++) {
    std::complex<float> *bandptr=nullptr, *bandptr1=nullptr;
    if(order==0)
      bandptr = bands[order];
    else {
      bandptr  = bands[2*order-1];
      bandptr1 = bands[2*order];
    }

    float k0x = ((float)order)/(norders-1) * k0mag * cos(k0angle);
    float k0y = ((float)order)/(norders-1) * k0mag * sin(k0angle);

#pragma omp parallel for
    for (auto kin=0; kin<nz; kin++) {
      int kout = kin;
      if (kout>kzcent) kout -= nz;
      auto kz = kout * dkz;
      for (auto iin=0; iin<ny; iin++) {
        int iout = iin;
        if (iout>kycent) iout -= ny;
        auto ky = iout * dky + k0y;
        for (auto jin=0; jin<kxcent+1; jin++) {
          auto kx = jin * dkx + k0x;
          auto indin = kin*nxy+iin*(nx/2+1)+jin;
          float ratio;
          if (! (order==0 && indin==0) ) {
            //! The distance (in 1/um) in Fourier space from a point to the origin
            float rho = sqrt(kz*kz + ky*ky + kx*kx);
            ratio = sphereFFT(rho, radius) / limit_at_origin;
          }
          else
            ratio = 1;  //! the limit at the origin
          bandptr[indin] /= ratio;
          if (order > 0)
            bandptr1[indin] /= ratio;
        }
      }
    }
  }
}

/* To get rid of checkerboard effect in the OTF bands */
/* (xc, yc, zc) is the estimated center of the point source, which in most cases is the bead */
/* Converted from Fortran code. kz is treated differently than kx and ky. don't know why */
void shift_center(std::complex<float> *bands, int nx, int ny, int nz, float xc, float yc, float zc)
{
  int kycent = ny/2;
  int kxcent = nx/2;
  int kzcent = nz/2;
  int nxy = (nx/2+1)*ny;

  float dphiz = 2*M_PI*zc/nz;
  float dphiy = 2*M_PI*yc/ny;
  float dphix = 2*M_PI*xc/nx;

#pragma omp parallel for
  for (int kin=0; kin<nz; kin++) {    // the origin of Fourier space is at (0,0)
    int kz = kin;
    if (kz>kzcent) kz -= nz;
    float phi1 = dphiz*kz;      //! first part of phi
    for (int iin=0; iin<ny; iin++) {
      int ky = iin;
      if (iin>kycent) ky -= ny;
      float phi2 = dphiy*ky;   //! second part of phi */
      for (int jin=0; jin<kxcent+1; jin++) {
        int kx = jin;
        int indin = kin*nxy + iin*(nx/2+1) + jin;
        float phi = phi1 + phi2 + dphix*kx;  //! All three parts of phi
        bands[indin] *= std::complex<float>(cos(phi), sin(phi));
      }
    }
  }
}

void radialft(std::complex<float> *band, int nx, int ny, int nz, float dr,
              std::complex<float> *avg_output, unsigned nr_avg)
{
  printf("In radialft()\n");
  auto kycent = ny/2;
  auto kxcent = nx/2;
  auto kzcent = nz/2;
  auto nxy = (nx/2+1)*ny;
  auto nxz = (nx/2+1)*nz;

  int * count = (int *) calloc(nxz, sizeof(int));
  assert(count);

  float dkx = 1.f / (nx * dr);
  float dky = 1.f / (ny * dr);
  float maxK = .5f / dr; //! maximum lateral frequency afforded by the sampling rate
  float dkr = 1.f/ (nr_avg * dr);
#pragma omp parallel for
  for (auto kin=0; kin<nz; kin++) {
    int kz = kin;
    if (kin>kzcent) kz -= nz;
    for (auto iin=0; iin<ny; iin++) {
      int ky = iin;
      if (iin>kycent) ky -= ny;
      float ky_f = ky * dky;
      for (auto jin=0; jin<kxcent+1; jin++) {
        auto kx = jin;
        float kx_f = kx * dkx;
        float rdist_f = sqrt(kx_f*kx_f + ky_f*ky_f);
        if (rdist_f < maxK) {
          auto indin = kin*nxy + iin*(nx/2+1)+ jin;
          unsigned indr = rint(rdist_f / dkr);
          assert (indr < nr_avg/2 + 1); 
          unsigned indout = indr * nz + kin;
          avg_output[indout] += band[indin];
          count[indout] ++;
        }
      }
    }
  }
  
#pragma omp parallel for
  for (auto indout=0; indout<nxz; indout++) {
    if (count[indout]>0)
      avg_output[indout] /= count[indout];
  }

//! Then force real value at the origin and Hermitian symmetry everywhere else
#pragma omp parallel for
  for (auto kx=0; kx<nx/2+1; kx++) {
    auto indout0 = kx*nz+0;
    avg_output[indout0] = avg_output[indout0].real();  //! Force real value at the origin 
    for (auto kz=1; kz<=nz/2; kz++) {
      auto indout = indout0 + kz;
      auto indout_conj = indout0 + (nz-kz);
      avg_output[indout] = (avg_output[indout] + conj(avg_output[indout_conj]))/2.f;
      avg_output[indout_conj] = conj(avg_output[indout]);
    }
  }
  free(count);
}

void cleanup(std::complex<float> *otfkxkz, int order, int nr, int nz, float dr, float dz,
             float linespacing, unsigned lamdanm, int icleanup, int ileave1_1, int ileave1_2,
             int ileave2, int twolens, float NA, float NIMM)
{
  printf("In cleanup()\n");
  float lamda = lamdanm * 0.001;
  float sinalpha = NA / NIMM;
  float cosalpha = cos(asin(sinalpha));
  float krmax = 2*NA/lamda;
  float k0mag = 1.0/linespacing;

  float dkr = 1./(nr * dr);
  float dkz = 1./(nz * dz);

  if (!twolens) {
#pragma omp parallel for
    for (auto ix=0; ix<icleanup; ix++) {
      int kzstart, kzend, jotfshape;
      auto kr = ix * dkr;
      if ( kr <= krmax ) {
        auto beta = asin((NA - kr*lamda) / NIMM);
        auto kzedge = (NIMM/lamda) * ( cos(beta) - cosalpha );

        if(order==0) {
          kzstart = rint((kzedge/dkz) + 1);
          kzend = nz - kzstart;
          for (auto iz=kzstart; iz<=kzend; iz++)
            otfkxkz[ix*nz+iz] = 0;
        }
        else if (order==1) {
          jotfshape = rint((kzedge/dkz) + 0.999);
          kzend = ileave1_1 - jotfshape;
          if (kzend>=0)
            for (auto iz=0; iz<=kzend; iz++) {
              if(iz==0)
                otfkxkz[ix*nz] = 0;
              else {
                otfkxkz[ix*nz+iz] = 0;
                otfkxkz[ix*nz+nz-iz] = 0;
              }
            }
          kzstart = ileave1_2 + jotfshape;
          for (auto iz=kzstart+1; iz<nz/2+1; iz++) {
            otfkxkz[ix*nz+iz] = 0;
            otfkxkz[ix*nz+nz-iz] = 0;
          }
        }
        else { /* order == 2 */
          jotfshape = rint((kzedge/dkz) + 0.999);
          kzstart = ileave2 + jotfshape;
          for (auto iz=kzstart; iz<nz/2+1; iz++) {
            otfkxkz[ix*nz+iz] = 0;
            otfkxkz[ix*nz+nz-iz] = 0;
          }
        }
     }
      else {   /* outside of lateral resolution limit */
        for (auto iz=0; iz<nz; iz++)
          otfkxkz[ix*nz+iz] = 0;
      }
    }
  }
  else {

    float lambdaem = lamda/NIMM;
    float lambdaexc = 0.88* lambdaem;;  // 0.88 approximates a typical lambdaexc/lambdaem
    float two_over_lambdaem = 2/lambdaem;
    float two_over_lambdaexc = 2/lambdaexc;
    // alpha = asin(NA/NIMM);  //! aperture angle of objectives
    float beta = asin(k0mag/two_over_lambdaexc);  //! angle of side illumination beams

    int center_of_arc;
    if (order==0) 
      center_of_arc = ceil(two_over_lambdaexc/dkz) + 3;
    else if (order == 1) {
      center_of_arc = ceil(1. / lambdaexc * (cos(beta) + 1) / dkz) + 3;
    }
    else { // order == 2
      center_of_arc = ceil(two_over_lambdaexc * cos(beta) / dkz) + 3;
    }
    
#pragma omp parallel for
    for (auto ix=0; ix<icleanup; ix++) {
      int kzstart, kzend;
      auto kr = ix * dkr;
      if ( kr <= krmax ) {
        kzstart = ceil(center_of_arc + sqrt(two_over_lambdaem * two_over_lambdaem - kr * kr)/dkz);
        kzend = nz - kzstart;
        for (auto iz=kzstart; iz<=kzend; iz++)
          otfkxkz[ix*nz+iz] = 0;
      }
      else 
        for (auto iz=0; iz<nz; iz++)
          otfkxkz[ix*nz+iz] = 0;
    }
  }
}

void cleanup_I2M(std::complex<float> *otfkxkz, int nx, int nz, float dkr, float dkz, int lamdanm, int icleanup, float NA, float NIMM)
{
  float lamda, sinalpha, cosalpha, krmax, two_over_lambdaem, beta, lamdaem, kr;
  int ix, iz, kzstart, kzend;

  lamda = lamdanm * 0.001;
  lamdaem = lamda/NIMM;
  sinalpha = NA/NIMM;
  cosalpha = cos(asin(sinalpha));
  krmax = 2*NA/lamda;
  two_over_lambdaem = 2/lamdaem;

  for (ix=0; ix<icleanup; ix++) {
    kr = ix * dkr;
    if (kr < krmax) {
      /* clean up stuff outside the side detection bands */
      kzstart = ceil(sqrt(two_over_lambdaem * two_over_lambdaem - kr * kr)/dkz)+1;
      kzend = nz - kzstart;
      for (iz=kzstart; iz<=kzend; iz++)
        otfkxkz[ix*nz+iz] = 0;

      /* Then, clean up stuff in-between center and side band */
      beta = asin( ( NA - kr*lamdaem ) /NIMM );
      kzstart = ceil((NIMM/lamdaem) * ( cos(beta) - cosalpha )/dkz);
      kzend = floor(sqrt(two_over_lambdaem * two_over_lambdaem - krmax*krmax)/dkz)-2;
      for (iz=kzstart; iz<=kzend; iz++)
        otfkxkz[ix*nz+iz] = 0;
      for (iz=nz-kzend; iz<nz-kzstart; iz++)
        otfkxkz[ix*nz+iz] = 0;
    }
    else
      for (iz=0; iz<nz; iz++)
        otfkxkz[ix*nz+iz] = 0;
  }
}


void outputdata(std::string &tiff_filename, std::vector<std::complex<float> *> &bands,
                int norders, int nx, int ny, int nz, float dkr, float dkz, int five_bands)
{
  printf("In outputdata()\n");

  CImg<> bandtemp;
  if (five_bands)
    bandtemp.assign(nz*2, nx/2+1, 2*norders-1); /* "nz*2" because of complex type */
  else
    if (nz>1)
      bandtemp.assign(nz*2, nx/2+1, norders);
    else
      bandtemp.assign(2*(nx/2+1), 1, norders); //??

  unsigned bufsize = nz*2*(nx/2+1)*sizeof(float);
  for (auto i=0; i<norders; i++)
    if (i==0) {
      memcpy(bandtemp.data(), (float *) bands[i], bufsize);
    }
    else {
      if (five_bands) {
        memcpy(bandtemp.data(0, 0, 2*i-1), (float *) bands[2*i-1], bufsize);
        memcpy(bandtemp.data(0, 0, 2*i), (float *) bands[2*i], bufsize);
      }
      else
        memcpy(bandtemp.data(0, 0, i), bands[2*i-1], bufsize);
    }
  bandtemp.save_tiff(tiff_filename.c_str());

  // Output another TIFF file with the amplitudes square-rooted, to make it easy for checking with FIJI, say
  
  CImg<> amps(bandtemp.width()/2, bandtemp.height(), bandtemp.depth());

  cimg_forXYZ (amps, x, y, z) {
    amps(x, y, z) = sqrt(abs(std::complex<float> (bandtemp(2*x, y, z), bandtemp(2*x+1, y, z))));
  }
  std::string amptiffname(tiff_filename);
  auto pos = tiff_filename.rfind(".tif");
  if (pos == std::string::npos)
    pos = tiff_filename.rfind(".TIF");

  if (pos != std::string::npos) {
    amptiffname.insert(pos, "_OnlyForViewing");
    amps.save_tiff(amptiffname.c_str());
  }
  else
    throw -1;
}


void outputdata(int ostream_no, IW_MRC_HEADER *header,  std::vector<std::complex<float> *> &bands, int norders, int nx, int ny, int nz, float dkr, float dkz, int five_bands)
{
  header->nx = nz;
  header->ny = nx/2+1;
  if (!five_bands)
    header->nz = norders;
  else
    header->nz = norders*2-1;
  header->mode = IW_COMPLEX;
  header->xlen = dkz;
  header->ylen = dkr;
  header->zlen = 0;
  header->amin = 0;
  header->amax = 1;
  
  IMPutHdr(ostream_no, header);
  for (int i=0; i<norders; i++)
    if (i==0)
      IMWrSec(ostream_no, bands[i]);
    else {
      IMWrSec(ostream_no, bands[2*i-1]);
      if (five_bands)
        IMWrSec(ostream_no, bands[2*i]);
    }
  IMWrHdr(ostream_no, header->label, 0, header->amin, header->amax, header->amean);
  IMClose(ostream_no);
}

/* According to Mats's derivation, FT of a sphere with radius R is f(k)=(R/Pi*k^2)*g(2*Pi*k*R), where g(x)=sin(x)/x - cos(x) */
/* The limit at the origin of F space is 4*Pi*R^3/3 */
double sphereFFT(float k, float radius)
{
  double a, x;
 
  if (k>0 || k<0) {
    x = 2*M_PI*radius*k;
    a = radius / (M_PI*k*k) * (sin(x)/x - cos(x));
  }
  else
    a = 4*M_PI*pow(radius, 3)/3;

  return a;
}

void modify(std::complex<float> *otfkxkz, int nx, int nz, int *ifixz, int ifixr, int order, int twolens)
{
  int kz, kx, ind0, ind1, ind2;

  if (ifixz[0])
    for (kz=0; kz<nz; kz++)
      if (kz !=0 || order > 0)
        otfkxkz[kz] = otfkxkz[nz+kz];    /* replace the kx=0 line with kx=1 line */
  
  if (ifixr) {
    for (kx=ifixr; kx<nx/2+1; kx++) {
      ind0 = kx*nz;
      ind1 = ind0+1;
      ind2 = ind0+nz-1;
      otfkxkz[ind0] = (otfkxkz[ind1] + otfkxkz[ind2])/2.f; 
    }
  }
}

void fixorigin(std::complex<float> *otfkxkz, int nx, int nz, int kx1, int kx2)
{
  float meani, slope, avg, lineval;
  int numvals, i, j;
  double *sum, totsum=0, ysum=0, sqsum=0;

  printf("In fixorigin(), kx1=%d, kx2=%d\n", kx1, kx2);
  meani = 0.5*(kx1+kx2);
  numvals = kx2-kx1+1;
  sum = (double *) malloc((kx2+1)*sizeof(double));

  for (i=0; i<=kx2; i++) { // note "<="
    sum[i] = 0;

    if (i == 0)
      sum[i] = otfkxkz[0].real();  /* don't want to add up garbages on kz axis */
    else
      for (j=0; j<nz; j++)
        sum[i] += otfkxkz[i*nz+j].real(); /* imaginary parts will cancel each other because of symmetry */

    if (i>=kx1) {
      totsum += sum[i];
      ysum += sum[i] * (i-meani);
      sqsum += (i-meani) * (i-meani);
    }
  }
  slope = ysum / sqsum;
  avg = totsum / numvals;

  for (i=0; i<kx1; i++) {
    lineval = avg + (i-meani)*slope;
    /* otfkxkz[i*nz].re -= sum[i] - lineval; */
    otfkxkz[i*nz] = std::complex<float> (otfkxkz[i*nz].real() - (sum[i] - lineval), otfkxkz[i*nz].imag());
  }

  free(sum);
}

void rescale(std::complex<float> *otfkxkz, int order, int nx, int nz, float *scalefactor, int dorescale)
{
  int nxz, ind;

  nxz = (nx/2+1)*nz;
  if (order==0 || dorescale) { /* need to find out scalefactor */
    // for (ind=0; ind<nxz; ind++) {
    //   mag = std::abs(otfkxkz[ind]);
    //   if (mag > valmax)
    // valmax = mag;
    // }
    // *scalefactor = 1/valmax;
    *scalefactor = 1/otfkxkz[0].real();
  }
  for (ind=0; ind<nxz; ind++)
    otfkxkz[ind] *= *scalefactor;
}


std::complex<float> otfinterpolate(std::complex<float> * otf, int kr, float krscale, int kz, float kzscale, int nzotf)
/* (kr, kz) is Fourier space coords with origin at kr=kz=0 and going  betwen -nx(or ny,nz)/2 and +nx(or ny,nz)/2 */
{
  std::complex<float> otfval;

  int irindex, izindex, indices[2][2];
  float krindex, kzindex;
  float ar, az;

  krindex = kr * krscale;
  kzindex = kz * kzscale;

  irindex = floor(krindex);
  izindex = floor(kzindex);

  ar = krindex - irindex;
  az = kzindex - izindex;  // az is always 0 for 2D case, and it'll just become a 1D interp

  if (izindex == nzotf-1) {
    indices[0][0] = irindex*nzotf+izindex;
    indices[0][1] = irindex*nzotf+0;
    indices[1][0] = (irindex+1)*nzotf+izindex;
    indices[1][1] = (irindex+1)*nzotf+0;
  }
  else {
    indices[0][0] = irindex*nzotf+izindex;
    indices[0][1] = irindex*nzotf+(izindex+1);
    indices[1][0] = (irindex+1)*nzotf+izindex;
    indices[1][1] = (irindex+1)*nzotf+(izindex+1);
  }
  otfval = (1-ar)*(otf[indices[0][0]]*(1-az) + otf[indices[0][1]]*az) +
    ar*(otf[indices[1][0]]*(1-az) + otf[indices[1][1]]*az);

  return otfval;
}

void downsize_pixel(std::complex<float> *otfra[], std::complex<float> *newotf[], int norders, int nkr, int nkz, int downfactor)
/* 
   Downsize the k-space pixel size by downfactor (an integer), so that the number of pixels increases but
   the same k-space extent remains.
*/
{
  int nkr_new, nkz_new;
  int order, kr, kz;

  nkr_new = nkr*downfactor;
  nkz_new = nkz*downfactor;

  for (order=0; order <norders; order ++) {
    std::complex<float> * otforiginal = otfra[order];
    int index_new = 0;
    for (kr=0; kr<nkr_new; kr++)
      for(kz=0; kz<nkz_new; kz++) {
        newotf[order][index_new] = otfinterpolate(otforiginal, kr, 1.0/downfactor, kz, 1.0/downfactor, nkz);
        index_new ++;
      }
  }
}

void Usage()
{

   /* not complete */
  printf("Usage:\nradialft [[input file] [output file]] [Options]\n");
  printf("\nOptions:\n");
   printf("\t-nphases -- number of phases; default 5\n");
  printf("\t-beaddiam F -- the diameter of the bead in microns (default is 0.12)\n");
  printf("\t-angle F -- the k0 vector angle with which the PSF is taken (default is 0)\n");
  printf("\t-ls x -- the illumination pattern's line spacing (in microns, default is 0.2)\n");
  printf("\t-nocompen -- do not perform bead size compensation (default is to perform)\n");
  printf("\t-5bands -- to output 5 OTF bands (default is combining higher-order's real and imag bands into one output\n");
  printf("\t-fixorigin kr1 kr2 -- the starting and end pixel for interpolation along kr axis (default is 2 and 9)\n");
  printf("\t-na x -- the (effective) NA of the objective\n");
  printf("\t-nimm x -- the index of refraction of the immersion liquid\n");
  printf("\t-leavekz kz1_1 kz1_2 kz2 -- the pixels to be retained on kz axis\n");
  printf("\t-I2M I2m_otf_file_name -- data contains I2M PSF and supply the I2M otf file name\n");
  printf("\t-background bkgd_value -- use the supplied number as the background to subtract\n");
  printf("\t-help or -h -- print this message\n");
}

#include <ctype.h>  /* isdigit() */

int getinteger(char * buffer, int * ret)
{
   auto len = strlen(buffer);
   for (auto i=0; i<len; i++) 
      if (!isdigit(buffer[i])) {
     if (i==0 && buffer[i] !='+' && buffer[i] !='-' )
        return 0;
     else if (i!=0)
        return 0;
      }
   sscanf(buffer, "%d", ret);
   return 1;
}

int getintegers(char *buffer[], int *ret, int num)
{
   for (auto i=0; i<num; i++) {
      if (!getinteger(buffer[i], ret+i)) 
     return 0;
   }
   return 1;
}

int getfloat(char *buffer, float *ret)
{
   unsigned dotcount = 0;

   auto len = strlen(buffer);
   for (auto i=0; i<len; i++)
      if (!isdigit(buffer[i])) {
     if (i==0 && buffer[i] !='+' && buffer[i] !='-' && buffer[i]!='.' )
        return 0;
     else if (i!=0 && buffer[i] != '.')
        return 0;
     if (buffer[i] == '.')
        dotcount ++;
      }
   if (dotcount > 1) return 0;
   sscanf(buffer, "%f", ret);
   return 1;
   
}

int getfloats(char *buffer[], float *ret, int num)
{
   int i;

   for (i=0; i<num; i++) {
      if (!getfloat(buffer[i], ret+i)) 
     return 0;
   }
   return 1;
}

int commandline(int argc, char *argv[], int * twolens, int *rescale, float *beaddiam, float *k0angle, float *linespacing, int *five_bands, int *nphases, int *interpkr, int *leavekz, int *do_compen, int *I2M_inc,  std::string &I2Mfiles, float * background, int *bBgInExtHdr, int *order0gen, std::string &order0files, int *conjugate, float *na, float *nimm, int * ifixz, int *ifixr, unsigned *wavelength, float *dr, float *dz, int *bCoherentBSIM, int *bForcedPIshift, std::string &ifiles, std::string &ofiles, int * ifilein, int *ofilein) {
   int ncomm=1;


   while (ncomm < argc) {
     if (strcmp(argv[ncomm], "-angle") == 0) {
       if (!getfloat(argv[ncomm+1], k0angle)) {
         printf("Invalid input for switch -angle\n");
         return 0;
       }
       ncomm +=2;
     }
     else if (strcmp(argv[ncomm], "-ls") == 0) {
       if (!getfloat(argv[ncomm+1], linespacing)) {
         printf("Invalid input for switch -ls\n");
         return 0;
       }
       ncomm +=2;
     }
     else if (strcmp(argv[ncomm], "-fixorigin") == 0) {
       if (!getintegers(argv+ncomm+1, interpkr, 2)) {
         printf("Invalid input for switch -fixorigin\n");
         return 0;
       }
       ncomm +=3;
     }
     else if (strcmp(argv[ncomm], "-leavekz") == 0) {
       if (!getintegers(argv+ncomm+1, leavekz, 3)) {
         printf("Invalid input for switch -leavekz\n");
         return 0;
       }
       ncomm +=4;
     }
     else if (strcmp(argv[ncomm], "-2lenses") == 0) {
       *twolens = 1;
       ncomm ++;
     }
     else if (strcmp(argv[ncomm], "-rescale") == 0) {
       *rescale = 1;
       ncomm ++;
     }
     else if (strcmp(argv[ncomm], "-5bands") == 0) {
       *five_bands = 1;
       ncomm ++;
     }
     else if (strcmp(argv[ncomm], "-conj") == 0) {
       *conjugate = 1;
       ncomm ++;
     }
     else if (strcmp(argv[ncomm], "-nocompen") == 0) {
       *do_compen = 0;
       ncomm ++;
     }
     else if (strcmp(argv[ncomm], "-I2M") == 0) {
       *I2M_inc = 1;
       I2Mfiles = argv[ncomm+1];
       ncomm +=2;
     }
     else if (strcmp(argv[ncomm], "-gen_order0") == 0) {
       *order0gen = 1;
       order0files = argv[ncomm+1];
       ncomm +=2;
     }
     else if (strcmp(argv[ncomm], "-nphases") == 0) {
       if (!getinteger(argv[ncomm+1], nphases)) {
         printf("Invalid input for switch -nphases\n");
         return 0;
       }
       ncomm +=2;
     }
     else if (strcmp(argv[ncomm], "-beaddiam") == 0) {
       if (!getfloat(argv[ncomm+1], beaddiam)) {
         printf("Invalid input for switch -beaddiam\n");
         return 0;
       }
       ncomm += 2;
     }
     else if (strcmp(argv[ncomm], "-background") == 0) {
       if (!getfloat(argv[ncomm+1], background)) {
         printf("Invalid input for switch -background\n");
         return 0;
       }
       ncomm += 2;
     }
    else if (strcmp(argv[ncomm], "-bgInExtHdr") == 0) {
      *bBgInExtHdr = 1;
      ncomm += 1;
    }
     else if (strcmp(argv[ncomm], "-na") == 0) {
       if (!getfloat(argv[ncomm+1], na)) {
         printf("Invalid input for switch -na\n");
         return 0;
       }
       ncomm += 2;
     }
     else if (strcmp(argv[ncomm], "-nimm") == 0) {
       if (!getfloat(argv[ncomm+1], nimm)) {
         printf("Invalid input for switch -nimm\n");
         return 0;
       }
       ncomm += 2;
     }
     else if (strcmp(argv[ncomm], "-ifixkz") == 0) {
       if (!getinteger(argv[ncomm+1], ifixz)) {
         printf("Invalid input for switch -ifixkz\n");
         return 0;
       }
       ncomm += 2;
     }
     else if (strcmp(argv[ncomm], "-ifixkr") == 0) {
       if (!getinteger(argv[ncomm+1], ifixr)) {
         printf("Invalid input for switch -ifixkr\n");
         return 0;
       }
       ncomm += 2;
     }
     else if (strcmp(argv[ncomm], "-wavelength") == 0) {
       // if (!getinteger(argv[ncomm+1], wavelength)) {
       try {
         * wavelength = std::stoi(argv[ncomm+1]);
         ncomm += 2;
       }
       catch (std::exception &e) {
         std::cout << "Invalid input for switch -wavelength: " << e.what() << std::endl;
         return 0;
       }
     }
     else if (strcmp(argv[ncomm], "-xyres") == 0) {
       if (!getfloat(argv[ncomm+1], dr)) {
         printf("Invalid input for switch -xyres\n");
         return 0;
       }
       ncomm += 2;
     }
     else if (strcmp(argv[ncomm], "-zres") == 0) {
       if (!getfloat(argv[ncomm+1], dz)) {
         printf("Invalid input for switch -zres\n");
         return 0;
       }
       ncomm += 2;
     }
     else if (strcmp(argv[ncomm], "-bessel") == 0) {
       *bCoherentBSIM = 1;
       ncomm += 1;
     }
    else if (strcmp(argv[ncomm], "-PIshift") == 0) {
      *bForcedPIshift = 1;
      ncomm ++;
    }
     else if (strcmp(argv[ncomm], "-help") == 0 || strcmp(argv[ncomm], "-h") == 0) {
       Usage();
       return 0;
     }
     else if (ncomm < 3) {
       if (ncomm==1) {
         ifiles = argv[ncomm];
         *ifilein = 1;
       }
       else if (ncomm==2) {
         ofiles = argv[ncomm];
         *ofilein = 1;
       }
       ncomm ++;
     }
     else {
       printf("Invalid command line option %s\n", argv[ncomm]);
       Usage();
       return 0;
     }
   }

   return 1;
}

void mrc_file_write(float *buffer, int nx, int ny, int nz, float rlen, float zlen, int mode, int iwave, const char *files)
{
  int ostream_no=19;
  IW_MRC_HEADER header;
  int dimx, dimy, dimz, nxy, i;
  float amin=0.f, amax=1.f, amean=0.1f; // to-do

  printf("Writing output file: %s\n", files);

  // if (IMOpen(ostream_no, files, "new")) {
  //   fprintf(stderr, "File %s can not be created.\n", files);
  //   exit(-1);
  // }

  dimx = nx;
  dimy = ny;
  dimz = nz;
  nxy = nx*ny;

  switch (mode) {
  case 0: /* normal floating-point data */
    header.mode = IW_FLOAT;
    break;
  case 1: /* half Fourier space complex data */
    header.mode = IW_COMPLEX;
    dimx = nx/2 + 1;
    nxy = dimx*ny;
    break;
  case 2: /* full Fourier space complex data */
     header.mode = IW_COMPLEX;
     nxy = nx*ny*2;
     break;
  case 3: /* half Fourier floating-point data */
    header.mode = IW_FLOAT;
    dimx = nx/2+1;
    nxy = dimx * dimy;
    break;
  default:
    fprintf(stderr, "Illegal mode in mrc_file_write()\n");
    exit(-1);
  }

  header.nx = dimx;
  header.mx = dimx;
  header.ny = dimy;
  header.my = dimy;
  header.nz = dimz;
  header.mz = dimz;
  header.ylen = rlen;
  header.xlen = rlen;
  header.zlen = zlen;
  header.nxst = 0;
  header.nyst = 0;
  header.nzst = 0;
  header.num_waves = 1;
  header.num_times = 1;
  header.iwav1 = iwave;
  header.alpha = 0;
  header.beta = 0;
  header.gamma = 0;
  header.mapc = 1;
  header.mapr = 2;
  header.maps = 3;
  header.ispg = 0;
  header.nDVID = -16244;
  header.ntst = 0;
  header.inbsym = 0;
  header.nint = 0;
  header.nreal = 0;
  header.nres = 1;
  header.nzfact = 1;
  header.file_type = 0;  /* normal image type */
  header.lens = 0;
  header.interleaved = 2;
  header.nlab = 1;
  IMPutHdr(ostream_no, &header);

  for (i=0; i<nz; i++)
    IMWrSec(ostream_no, buffer + i*nxy);

  IMWrHdr(ostream_no, "Processed", 1, amin, amax, amean);
  IMClose(ostream_no);
}
