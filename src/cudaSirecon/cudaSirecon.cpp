
#include "cudaSireconImpl.h"
#include "mrc.h"
#include "SIM_reconstructor.hpp"
#include "cudasireconConfig.h"
#include <boost/filesystem.hpp>


std::string version_number = std::to_string(cudasirecon_VERSION_MAJOR) + "." +
                             std::to_string(cudasirecon_VERSION_MINOR) + "." +
                             std::to_string(cudasirecon_VERSION_PATCH);

void SetDefaultParams(ReconParams *pParams)
{
  pParams->k0startangle = 1.57193;
  pParams->linespacing = 0.177;
  pParams->na = 1.36;
  pParams->nimm = 1.515;
  pParams->ndirs = 3;
  pParams->nphases = 3;
  pParams->phaseSteps = 0;
  pParams->norders_output = 0;
  pParams->bTwolens = false;
  pParams->bFastSIM = false;
  pParams->bBessel = false;
  pParams->BesselNA = 0.45;

  pParams->otfcutoff = 0.006;
  pParams->zoomfact = 2;
  pParams->z_zoom = 1;
  pParams->nzPadTo = 0;
  pParams->explodefact= 1.0;
  pParams->bFilteroverlaps = 1;
  pParams->recalcarrays = 1; //! whether to calculate the overlaping regions between bands just once or always; used in fitk0andmodamps()
  pParams->napodize = 10;
  pParams->forceamp.assign(1, 0.f);
  // pParams->k0angles = NULL;
  pParams->bSearchforvector = 1;
  pParams->bUseTime0k0 = 1;  /* default to use time 0's k0 fit for the rest in a time series data */
  pParams->apodizeoutput = 0;  /* 0-no apodize; 1-cosine apodize; 2-triangle apodize; used in filterbands() */
  pParams->bSuppress_singularities = 1;  /* whether to dampen the OTF values near band centers; used in filterbands() */
  pParams->suppression_radius = 10;   /* if suppress_singularities is 1, the range within which suppresion is applied; used in filterbands() */
  pParams->bDampenOrder0 = false;  /* whether to dampen order 0 contribution; used in filterbands() */
  pParams->bFitallphases = 1;
  pParams->do_rescale=1; /* fading correction method: 0-no correction; 1-with correction */
  pParams->equalizez = false;
  pParams->equalizet = false;
  pParams->bNoKz0 = false;

  pParams->bRadAvgOTF = false;  /* default to use non-radially averaged OTFs */
  pParams->bOneOTFperAngle = false;  /* default to use one OTF for all SIM angles */

  pParams->bFixdrift = 0;
  pParams->drift_filter_fact = 0.0;

  pParams->constbkgd = 0.0;
  pParams->bBgInExtHdr = 0;
  pParams->bUsecorr = 0;
  pParams->electrons_per_bit = 0.6528;
  pParams->readoutNoiseVar = 32.42;  // electron^2

  pParams->bMakemodel = 0;
  pParams->bSaveAlignedRaw = 0;
  pParams->bSaveSeparated = 0;
  pParams->bSaveOverlaps = 0;
  pParams->bWriteTitle = false;

  pParams -> bTIFF = false;
}


// Functions for dealing with images

bool notGoodDimension(unsigned num)
// Good dimension is defined as one that can be fatorized into 2s, 3s, 5s, and 7s
// According to CUFFT manual, such dimension would warranty fast FFT
{
  if (num==2 || num==3 || num==5 || num==7)
    return false;
  else if (num % 2 == 0) return notGoodDimension(num / 2);
  else if (num % 3 == 0) return notGoodDimension(num / 3);
  else if (num % 5 == 0) return notGoodDimension(num / 5);
  else if (num % 7 == 0) return notGoodDimension(num / 7);
  else
    return true;
}

unsigned findOptimalDimension(unsigned inSize, int step=-1)
{
  unsigned outSize = inSize;
  while (notGoodDimension(outSize))
    outSize += step;

  return outSize;
}


void allocateOTFs(ReconParams* pParams, int sizeOTF,
                  std::vector<std::vector<GPUBuffer> >& otfs)
{
  otfs.clear();
  unsigned short nDirsOTF = 1;
  if (pParams->bOneOTFperAngle)
    nDirsOTF = pParams->ndirs;
  otfs.resize(nDirsOTF);
  for (int dir = 0; dir< nDirsOTF; dir++)
    for (int i = 0; i < pParams->norders; ++i) {
      GPUBuffer buff;
      buff.resize(sizeOTF * sizeof(cuFloatComplex));
      otfs[dir].push_back(buff);
  }
}

void allocateImageBuffers(const ReconParams& params,
    const ImageParams& imgParams, ReconData* data)
{
  data->savedBands.clear();
  data->savedBands.reserve(params.ndirs);
  for (int i = 0; i < params.ndirs; ++i) {
    data->savedBands.push_back(std::vector<GPUBuffer>());
    data->savedBands[i].reserve(params.nphases);
    for (int j = 0; j < params.nphases; ++j) {
      data->savedBands[i].push_back(GPUBuffer(
            (imgParams.nx / 2 + 1) * imgParams.ny * imgParams.nz0 *
            sizeof(cuFloatComplex), 0));
    }
  }
  for (int i = 0; i < params.ndirs; ++i) {
    for (int j = 0; j < params.nphases; ++j) {
      data->savedBands[i][j].setToZero();
    }
  }
}


void bgAndSlope(const ReconParams& myParams,
    const ImageParams& imgParams, ReconData* reconData)
{
  reconData->background.resize(sizeof(float) * imgParams.nx *
      imgParams.ny);
  reconData->slope.resize(sizeof(float) * imgParams.nx *
      imgParams.ny);

  if (myParams.bUsecorr) {
    // flatfield correction of measured data using calibration data
    printf("loading CCD calibration file\n");
    getbg_and_slope(myParams.corrfiles.c_str(),
        (float*)reconData->background.getPtr(),
        (float*)reconData->slope.getPtr(), imgParams.nx, imgParams.ny);
  } else {
    for (int i = 0; i < imgParams.nx * imgParams.ny; i++) {
      /* use the constant background value given by user */
      ((float*)(reconData->background.getPtr()))[i] = myParams.constbkgd;
      ((float*)(reconData->slope.getPtr()))[i] = 1.0;
    }
  }
  reconData->backgroundExtra = 0;
}

void getbg_and_slope(const char *corrfiles, float *background,
    float *slope, int nx, int ny)
{
  int cstream_no=10;
  int ixyz[3], mxyz[3], pixeltype;      /* variables for IMRdHdr call */
  float min, max, mean;      /* variables for IMRdHdr call */
  
  IW_MRC_HEADER header;
  if (IMOpen(cstream_no, corrfiles, "ro")) {
    fprintf(stderr, "File %s does not exist\n", corrfiles);
    exit(-1);
  }

  IMRdHdr(cstream_no, ixyz, mxyz, &pixeltype, &min, &max, &mean);
  IMGetHdr(cstream_no, &header);
  if (header.nx != nx || header.ny != ny) {
    fprintf(stderr,
        "calibration file %s has different dimension than data file",
        corrfiles);
    exit(-1);
  }
  IMAlCon(cstream_no, 0);

  IMRdSec(cstream_no, background);
  IMRdSec(cstream_no, slope);

  IMClose(cstream_no);
}

void findModulationVectorsAndPhasesForAllDirections(
    int zoffset, ReconParams* params, const ImageParams& imgParams,
    DriftParams* driftParams, ReconData* data)
{
  // Apodize (or edge softening) every 2D slice:
  apodizationDriver(zoffset, params, imgParams, driftParams, data);

  // 2D FFT every 2D slice:
  transformXYSlice(zoffset, params, imgParams, driftParams, data);

  float dkx = 1./(imgParams.dxy * imgParams.nx);
  float dky = 1./(imgParams.dxy * imgParams.ny);
  float k0magguess = 1.0 / params->linespacing;

  for (int direction = 0; direction < params->ndirs; ++direction) {

    std::vector<GPUBuffer>* rawImages = &(data->savedBands[direction]);
    std::vector<GPUBuffer>* bands = &(data->savedBands[direction]);

    std::vector<float> phaseList(params->nphases);
    if (params->phaseSteps != 0) {
      // User specified the non-ideal (i.e., not 2pi/5) phase steps for
      // each orientation. Under this circumstance, calculate the
      // non-ideal sepMatrix:
      for (int i = 0; i < params->nphases; ++i) {
        phaseList[i] = i * params->phaseSteps[direction];
      }
      makematrix(params->nphases, params->norders, direction,
          &(phaseList[0]), (&data->sepMatrix[0]),
          &(data->noiseVarFactors[0]));
    }

    // Unmixing info components in real or reciprocal space:
    separate(imgParams.nx, imgParams.ny, imgParams.nz,
        direction, params->nphases, params->norders,
        rawImages, &data->sepMatrix[0]);

#ifndef NDEBUG
    for (int phase = 0; phase < params->nphases; ++phase) {
      assert(rawImages->at(phase).hasNaNs(true) == false);
      std::cout << "Phase " << phase << "ok." << std::endl;
    }
#endif

    if (imgParams.nz > 1) {
      // 1D FFT of a stack of 2D FFTs to obtain equivalent of 3D FFT
      cufftHandle fftplan1D;
      int fftN[1] = {imgParams.nz0};
      int stride = (imgParams.nx / 2 + 1) * imgParams.ny;
      int dist = 1;
      cufftResult cuFFTErr = cufftPlanMany(&fftplan1D, 1, fftN,
                                           fftN, stride, dist,
                                           fftN, stride, dist,
                                           CUFFT_C2C,
                                           (imgParams.nx / 2 + 1) * imgParams.ny);
      if (cuFFTErr != CUFFT_SUCCESS) {
        throw std::runtime_error("CUFFT plan creation failed");
      }
      for (int i = 0; i < params->nphases; ++i) {
        cufftExecC2C(fftplan1D, (cuFloatComplex*)((*bands)[i]).getPtr(),
                     (cuFloatComplex*)((*bands)[i]).getPtr(), CUFFT_FORWARD);
      }
      cufftDestroy(fftplan1D);
    }

    if (params->bMakemodel) {
      /* Use the OTF to simulate an ideal point source; replace bands
       * with simulated data
       * DM: k0 is not initialized but has memory allocate at this point.
       * k0 initialization code is near lines 430ff in sirecon.c */
      makemodeldata(imgParams.nx, imgParams.ny, imgParams.nz0, bands,
          params->norders, data->k0[direction], imgParams.dxy, imgParams.dz,
          &data->otf[0], imgParams.wave[0], params);
    }

    /* save the separated raw if requested */
    if (params->bTIFF && params->bSaveSeparated) {
      // TO-DO
      std::cout << "No unmixed raw images were saved in TIFF mode\n";
    }
    else if (!params->bTIFF && params->bSaveSeparated) {
      CPUBuffer tmp((*rawImages)[0].getSize());
      for (int phase = 0; phase < params->nphases; ++ phase) {
        (*rawImages)[phase].set(&tmp, 0, tmp.getSize(), 0);
        for (int z = 0; z < imgParams.nz; ++z) {
          float* imgPtr = (float*)tmp.getPtr();
          IMWrSec(separated_stream_no,
              imgPtr + (z + zoffset) * (imgParams.nx/2 + 1)*2 * imgParams.ny);
        }
      }
      continue; // skip the k0 search and modamp fitting
    }


    /* After separation and FFT, the std::vector rawImages, now referred to as
     * the std::vector bands, contains the center band (bands[0])
     * , the real (bands[1], bands[3], etc.) and
     * imaginary part (bands[2], bands[4], etc.) bands
     * of the different orders */

    if (! (imgParams.ntimes > 1 && imgParams.curTimeIdx > 0
           && params->bUseTime0k0) ) {
      data->k0[direction] = data->k0guess[direction];
      printf("k0guess[direction %d] = (%f, %f) pixels\n", direction,
          data->k0guess[direction].x/dkx, data->k0guess[direction].y/dky);
    }

    // Now to fix 3D drift between dirs estimated by determinedrift_3D()
    if (direction != 0 && params->bFixdrift) {
      fixdrift_bt_dirs(bands, params->norders, driftParams->drift_bt_dirs[direction],
                       imgParams.nx, imgParams.ny, imgParams.nz0);
    }

    /* assume k0 vector not well known, so fit for it */
    cuFloatComplex amp_inv;
    cuFloatComplex amp_combo;

    // dir_ is used in upcoming calls involving data->otf, to differentiate the cases
    // of common OTF and dir-specific OTF
    int dir_=0;
    if (params->bOneOTFperAngle)
      dir_ = direction;

    if (params->bSearchforvector &&
        !(imgParams.ntimes > 1 && imgParams.curTimeIdx > 0 && params->bUseTime0k0)) {
      /* In time series, can choose to use the time 0 k0 fit for the
       * rest of the series.
       * Find initial estimate of modulation wave vector k0 by
       * cross-correlation. */
      findk0(bands, &data->overlap0, &data->overlap1, imgParams.nx,
             imgParams.ny, imgParams.nz0, params->norders,
             &(data->k0[direction]), imgParams.dxy, imgParams.dz, &(data->otf[dir_]),
             imgParams.wave[0], params);

      if (params->bSaveOverlaps) {
        // output the overlaps
        if (params->bTIFF) {
          CPUBuffer tmp0(data->overlap0.getSize()*2);
          data->overlap0.set(&tmp0, 0, data->overlap0.getSize(), 0);
          data->overlap1.set(&tmp0, 0, data->overlap1.getSize(), data->overlap0.getSize());
          CImg<> ovlp0((float* )tmp0.getPtr(), imgParams.nx*2, imgParams.ny, imgParams.nz*2, 1,
                       true);  // ovlp0 shares buffer with tmp0 (hence "true" in the final parameter)
          ovlp0.save_tiff(params->fileOverlaps.c_str());
        }
        else {
          CPUBuffer tmp0(data->overlap0.getSize());
          data->overlap0.set(&tmp0, 0, tmp0.getSize(), 0);
          CPUBuffer tmp1(data->overlap1.getSize());
          data->overlap1.set(&tmp1, 0, tmp1.getSize(), 0);
          cuFloatComplex* ol0Ptr = (cuFloatComplex*)tmp0.getPtr();
          cuFloatComplex* ol1Ptr = (cuFloatComplex*)tmp1.getPtr();
          for (int z = 0; z < imgParams.nz0; ++z) {
            IMWrSec(overlaps_stream_no, ol0Ptr + z * imgParams.nx * imgParams.ny);
            IMWrSec(overlaps_stream_no, ol1Ptr + z * imgParams.nx * imgParams.ny);
          }
        }
      }
      printf("Initial guess by findk0() of k0[direction %d] = (%f,%f) pixels\n", 
          direction, data->k0[direction].x/dkx, data->k0[direction].y/dky);

      /* refine the cross-corr estimate of k0 by a search for best k0
       * vector direction and magnitude using real space waves*/

      printf("before fitk0andmodamp\n");
      /* if (myParams.recalcarrays==0 && dist>1.0) */
      params->recalcarrays = 1;
      /* if k0 is very close to the guess, we can save time by not
       * recalculating the overlap arrays */

      fitk0andmodamps(bands, &data->overlap0, &data->overlap1, imgParams.nx,
          imgParams.ny, imgParams.nz0, params->norders, &(data->k0[direction]),
          imgParams.dxy, imgParams.dz, &(data->otf[dir_]), imgParams.wave[0],
          &data->amp[direction][0], params);

      if (imgParams.curTimeIdx == 0) {
        data->k0_time0[direction] = data->k0[direction];
      }
      /* check if the k0 vector found is reasonably close to the guess */
      vector deltak0;
      deltak0.x = data->k0[direction].x - data->k0guess[direction].x;
      deltak0.y = data->k0[direction].y - data->k0guess[direction].y;
      float dist = sqrt(deltak0.x * deltak0.x + deltak0.y * deltak0.y);
      if (dist/k0magguess > K0_WARNING_THRESH) {
        printf("WARNING: ");
      }
      printf("best fit for k0 is %.3f%% from expected value.\n", dist/k0magguess*100);

      if (imgParams.ntimes > 1 && imgParams.curTimeIdx > 0
          && dist > 2*K0_WARNING_THRESH) {
        data->k0[direction] = data->k0_time0[direction];
        printf("k0 estimate of time point 0 is used instead\n");
        for (int order = 1; order < params->norders; ++order) {
          float corr_coeff;
          if (imgParams.nz0>1)
            corr_coeff = findrealspacemodamp(bands, &data->overlap0,
              &data->overlap1, imgParams.nx, imgParams.ny, imgParams.nz0,
              0, order, data->k0[direction], imgParams.dxy, imgParams.dz,
              &(data->otf[dir_]), imgParams.wave[0], &data->amp[direction][order],
              &amp_inv, &amp_combo, 1, params);
          else
            corr_coeff = findrealspacemodamp(bands, &data->overlap0,
              &data->overlap1, imgParams.nx, imgParams.ny, imgParams.nz0,
              order-1, order, data->k0[direction], imgParams.dxy, imgParams.dz,
              &(data->otf[dir_]), imgParams.wave[0], &data->amp[direction][order],
              &amp_inv, &amp_combo, 1, params);
          printf("modamp mag=%f, phase=%f\n, correlation coeff=%f\n\n",
                 cmag(data->amp[direction][order]),
                 atan2(data->amp[direction][order].y, data->amp[direction][order].x),
                 corr_coeff);
        }
      }
    } else {
      /* assume k0 vector known, so just fit for the modulation amplitude and phase */
      printf("known k0 for direction %d = (%f, %f) pixels \n", direction, 
          data->k0[direction].x/dkx, data->k0[direction].y/dky);
      for (int order = 1; order < params->norders; ++order) {
        float corr_coeff = findrealspacemodamp(bands, &data->overlap0,
            &data->overlap1, imgParams.nx, imgParams.ny, imgParams.nz0, 
            0, order, data->k0[direction], imgParams.dxy, imgParams.dz,
            &(data->otf[dir_]), imgParams.wave[0], &data->amp[direction][order],
            &amp_inv, &amp_combo, 1, params);
        printf("modamp mag=%f, phase=%f\n",
            cmag(data->amp[direction][order]),
            atan2(data->amp[direction][order].y, data->amp[direction][order].x));
        printf("reverse modamp mag=%f, phase=%f\n", 1.0f / cmag(amp_inv),
            -atan2(amp_inv.y, amp_inv.x));
        printf("combined modamp mag=%f, phase=%f\n",
            cmag(amp_combo), atan2(amp_combo.y, amp_combo.x));
        printf("correlation coeff=%f\n\n", corr_coeff);

        if (params->bTIFF) {
          // To-DO
          std::cout << "No overlaps were saved in TIFF mode";
        }
        else if (!params->bTIFF && order == 1 && params->bSaveOverlaps) {// output the overlaps
          // output the overlaps
          CPUBuffer tmp0(data->overlap0.getSize());
          data->overlap0.set(&tmp0, 0, tmp0.getSize(), 0);
          CPUBuffer tmp1(data->overlap1.getSize());
          data->overlap1.set(&tmp1, 0, tmp1.getSize(), 0);
          float* ol0Ptr = (float*)tmp0.getPtr();
          float* ol1Ptr = (float*)tmp1.getPtr();
          for (int z = 0; z < imgParams.nz0; ++z) {
            IMWrSec(overlaps_stream_no, ol0Ptr + z * imgParams.nx * imgParams.ny);
            IMWrSec(overlaps_stream_no, ol1Ptr + z * imgParams.nx * imgParams.ny);
          }
        }
      }
    }     /* if(searchforvector) ... else ... */

    if (imgParams.nz == 1) {
      /* In 2D SIM, amp stores modamp's between each adjacent pair of
       * bands. We want to convert this to modamp w.r.t. order 0 */
      for (int order = 2; order < params->norders; ++order) {
        data->amp[direction][order] = cmul(data->amp[direction][order],
            data->amp[direction][order - 1]);
      }
    }

    if (params->forceamp[0] > 0.0) {
      /* force modamp's amplitude to be a value user provided (ideally
       * should be 1)  */
      for (int order = 1; order < params->norders; ++order) {
        float a = cmag(data->amp[direction][order]);
        if (a < params->forceamp[order-1]) {
          float ampfact = params->forceamp[order-1] / a;
          data->amp[direction][order].x *= ampfact;
          data->amp[direction][order].y *= ampfact;
          printf("modamp mag=%f, phase=%f  \n",
              cmag(data->amp[direction][order]),
              atan2(data->amp[direction][order].y, data->amp[direction][order].x));
        }
      }
    }

    /* In 2D NLSIM, we often don't trust the modamp fit between neighboring high-order components; */
    /* only if fitallphases is True do all fitted modamps get actually used; */
    /* otherwise, order 2 and above's global phase is inferred from order 1's phase.  */
    if (!params->bFitallphases) {
      float base_phase = get_phase(data->amp[direction][1]);
      cuFloatComplex expiphi;
      cuFloatComplex amplitude;
      float phi;
      for (int order = 2; order < params->norders; ++order) {
        amplitude.x = cmag(data->amp[direction][order]);
        amplitude.y = 0;
        phi = order * base_phase /*+ M_PI*(order-1)*/;
        /* sign flip for every other order happens only in saturation case? */
        expiphi.x = cos(phi);
        expiphi.y = sin(phi);
        data->amp[direction][order] = cmul(amplitude, expiphi);
      }
    }
  } /* end for (dir)  */   /* done finding the modulation vectors and phases for all directions */
}

void apodizationDriver(int zoffset, ReconParams* params,
    const ImageParams& imgParams, DriftParams* driftParams, ReconData* data)
{
  for (int direction = 0; direction < params->ndirs; ++direction) {
    /* data assumed taken with z changing fast, direction changing
     * slowly */
    std::vector<GPUBuffer>* rawImages = &(data->savedBands[direction]);
    std::vector<GPUBuffer>* bands = &(data->savedBands[direction]);

    for (int z = 0; z < imgParams.nz; ++z) {
      for (int phase = 0; phase < params->nphases; ++phase) {
        if (params->napodize > 0) {
          apodize(params->napodize, imgParams.nx, imgParams.ny,
              &(rawImages->at(phase)),
              (z + zoffset) * (imgParams.nx/2 + 1)*2 * imgParams.ny);
        } else if (params->napodize == -1) {
          cosapodize(imgParams.nx, imgParams.ny, &(*rawImages)[phase],
              (z + zoffset) * (imgParams.nx/2 +1)*2 * imgParams.ny);
        }
      } /* end for (phase), loading, flatfielding, and apodizing raw images */
    }
  } /* end for (dir)  */
}

void rescaleDriver(int it, int iw, int zoffset, ReconParams* params,
    const ImageParams& imgParams, DriftParams* driftParams, ReconData* data)
{
  for (int direction = 0; direction < params->ndirs; ++direction) {
    /* data assumed taken with z changing fast, direction changing
     * slowly */
    std::vector<GPUBuffer>* rawImages = &(data->savedBands[direction]);
    std::vector<GPUBuffer>* bands = &(data->savedBands[direction]);

    for (int z = 0; z < imgParams.nz; ++z)
      if (params->do_rescale) {
        // Goes through here
        rescale(imgParams.nx, imgParams.ny, imgParams.nz, z, zoffset,
            direction, iw, it, params->nphases, rawImages, params->equalizez,
            params->equalizet, &data->sum_dir0_phase0[0]);
      }
#ifndef NDEBUG
    std::cout<< "rescaleDriver(): " << cudaGetErrorString(cudaGetLastError()) << std::endl;
#endif
  }
}

void transformXYSlice(int zoffset, ReconParams* params,
    const ImageParams& imgParams, DriftParams* driftParams, ReconData* data)
{
  cufftHandle rfftplanGPU;
  int fftN[2] = {imgParams.ny, imgParams.nx};
  int inembed[2] = {imgParams.nx * imgParams.ny, (imgParams.nx/2 + 1)*2};
  int istride = 1;
  int idist = (imgParams.nx/2 + 1)*2 * imgParams.ny;
  int onembed[2] = {imgParams.nx * imgParams.ny, imgParams.nx /2 +1};
  int ostride = 1;
  int odist = (imgParams.nx / 2 + 1) * imgParams.ny;
  cufftResult cuFFTErr = cufftPlanMany(&rfftplanGPU, 2, &fftN[0],
      inembed, istride, idist,
      onembed, istride, odist,
      CUFFT_R2C, imgParams.nz);
  if (cuFFTErr != CUFFT_SUCCESS) {
    std::cout << "Error code: " << cuFFTErr << " at " __FILE__ << ":" << __LINE__ << std::endl;
    throw std::runtime_error("cufftPlanMany() failed.");
  }

  for (int direction = 0; direction < params->ndirs; ++direction) {
    for (int phase = 0; phase < params->nphases; ++phase) {
      cuFFTErr = cufftExecR2C(rfftplanGPU,
          (float*)data->savedBands[direction][phase].getPtr(),
          (cuFloatComplex*)data->savedBands[direction][phase].getPtr());
      if (cuFFTErr != CUFFT_SUCCESS) {
        std::cout << "Line:" << __LINE__ ;
        throw std::runtime_error("cufft failed.");
      }
    }
  }

  cufftDestroy(rfftplanGPU);
  if (cuFFTErr != CUFFT_SUCCESS) {
    std::cout << "Error code: " << cuFFTErr << " at " __FILE__ << ":" << __LINE__<< std::endl;
    throw std::runtime_error("cufftDestroy() failed.");
  }
}

/***************************** makematrix ************************************/
/*     generates the matrix that is to be used to separate the raw indata    */
/*     into the different bands of sample information.                       */
/*  Two cases:
    1. If arrPhases is NULL, we're certain that phases are equally spaced within 2Pi range;
    2. Else, user supplies a list of phases actually used (e.g, when there is drift or 
    non-ideal SLM patterns are used), and noise variance amplification factors for each 
    order and direction is calculated based on the non-ideal separation matrix
    */
/*****************************************************************************/
void makematrix(int nphases, int norders, int dir, float *arrPhases,
    float *sepMatrix, float * noiseVarFactors)
{
  std::cout << "In makematrix." << std::endl;
  /*
     norders is not necessarily nphases/2+1, in cases such as nonlinear SIM
     where the number of output orders can be smaller.
     */
  int j, order;
  float phi, sqrt_nphases, sqrt_2;

  sqrt_nphases = sqrt((float)nphases);
  sqrt_2 = sqrt(2.);

  // norders could be less than (nphases+1)/2
  if (arrPhases == 0) { /* phase values equally spaced between 0 and 2Pi */
    phi = 2.0 * M_PI / nphases;
    for (j = 0; j < nphases; ++j) {
      sepMatrix[0 * nphases + j] = 1.0;
      for (order = 1; order < norders; ++order) {
        /* With this definition, bandplus = bandre + i bandim 
         * has coefficient exp() with unit normalization */
        sepMatrix[(2 * order - 1) * nphases + j] = cos(j * order * phi);
        sepMatrix[2 * order * nphases + j] = sin(j * order * phi);
      }
    }
  }
  else {/* User supplied phase values in arrPhases */
    int *ipvt, info, i, nCols;
    float *workspace, *forwardM, *A_t_A, unit=1., zero=0.;

    /* Here we differentiate between 2 options: */
    /* 1. use direct matrix inversion;          */
    /* 2. use least-square solution.            */

    nCols = 2 * norders - 1;

    if (nCols < nphases) {  // use least-square
      ipvt = (int*)malloc(nCols * sizeof(int));
      forwardM   = (float*)malloc(nphases*nCols*sizeof(float));
      A_t_A = (float*)malloc(nCols*nCols*sizeof(float));
      /* First construct the forward matrix*/
      for (i=0; i<nphases; i++) {
        forwardM[i] = 1.0 / nphases;
        for (order=1; order<norders; order++) {
          forwardM[i+(2*order-1)*nphases] = 2 * cos(order*arrPhases[i]) / nphases;
          forwardM[i+2*order*nphases]     = 2 * sin(order*arrPhases[i]) / nphases;
        }
      }

      /*       print_matrix("forward matrix:", nphases, nCols, forwardM, nphases); */

      /* multiply transposed forward matrix with forward matrix */
      sgemm_("T", "N", &nCols, &nCols, &nphases, &unit, forwardM,  &nphases,
          forwardM, &nphases,  &zero, A_t_A,  &nCols);
      /*       print_matrix("A transpose times A:", nCols, nCols, A_t_A, nCols); */

      /* Then invert the forward matrix to form the sep matrix*/
      sgetrf_(&nCols, &nCols, A_t_A, &nCols, ipvt, &info);
      if (info==0) { // successful factorization
        workspace = (float*)malloc(nCols*sizeof(float));
        sgetri_(&nCols, A_t_A, &nCols, ipvt, workspace, &nCols, &info);
        free(workspace);
      }
      else
        printf("WARNING: sgetri() returns non 0.\n");

      if (info!=0)
        printf("WARNING: sgetrf() returns non 0.\n");

      // multiply inverse A_t_A with transposed forward matrix
      sgemm_("N", "T", &nCols, &nphases,  &nCols, &unit, A_t_A,  &nCols,
          forwardM, &nphases,  &zero, sepMatrix, &nCols);

      /*       print_matrix("Separation Matrix:", nCols, nphases, sepMatrix, nCols); */
      /* transpose back to C-style matrix */
      matrix_transpose(sepMatrix, nphases, nCols);

      free(ipvt);
      free(forwardM);
      free(A_t_A);
    }
    else {  // use direct inversion
      ipvt = (int*)malloc(nphases*sizeof(int));
      /* First construct the forward matrix (in C-style storage convention ) */
      for (i=0; i<nphases; i++) {
        sepMatrix[i*nphases] = 1.0 / nphases;
        for (order=1; order<norders; order++) {
          sepMatrix[i*nphases+(2*order-1)] = 2 * cos(order*arrPhases[i]) / nphases;
          sepMatrix[i*nphases+2*order]     = 2 * sin(order*arrPhases[i]) / nphases;
        }
      }
      /* Then invert the forward matrix to form the sep matrix*/
      sgetrf_(&nphases, &nphases, sepMatrix, &nphases, ipvt, &info);
      if (info==0) { // successful factorization
        workspace = (float*)malloc(nphases*sizeof(float));
        sgetri_(&nphases, sepMatrix, &nphases, ipvt, workspace, &nphases, &info);
        free(workspace);
      }
      else
        printf("WARNING: sgetri() returns non 0.\n");

      if (info!=0)
        printf("WARNING: sgetrf() returns non 0.\n");
      free(ipvt);
    }

    /* Report noise factors */
    for (order=0;order<norders;order++)   noiseVarFactors[dir*norders+order] = 0.0;
    for (j=0;j<nphases;j++) {
      noiseVarFactors[dir*norders+0] += pow(sepMatrix[0*nphases+j], 2);
      for (order=1;order<norders;order++)
        noiseVarFactors[dir*norders+order] += pow(sepMatrix[(2*order-1)*nphases+j], 2) + pow(sepMatrix[2*order*nphases+j], 2);
    }

    printf(" Separation noise factors: ");
    for (order=0;order<norders;order++) {
      noiseVarFactors[dir*norders+order] /= nphases;
      printf(" Order %d: %.3f,",order, sqrt(noiseVarFactors[dir*norders+order]));
    }
    printf("\n");

  }

  printf("Separation matrix:\n");
  for (j=0;j<norders*2-1;j++) {
    int i;
    for (i=0; i<nphases; i++)
      printf("%9.5f ", sepMatrix[j*nphases+i]);
    printf("\n");
  }
  printf("\n");
}

void allocSepMatrixAndNoiseVarFactors(const ReconParams& params, ReconData* reconData)
{
  reconData->sepMatrix.resize(params.nphases * (params.norders * 2 - 1));
  reconData->noiseVarFactors.resize(params.ndirs * params.norders, 1.0f);
}

void matrix_transpose(float* mat, int nRows, int nCols)
{
  int i, j;
  float* tmpmat = (float*)malloc(nRows*nCols*sizeof(float));

  for (i=0; i<nRows; i++)
    for (j=0; j<nCols; j++)
      tmpmat[j*nRows+i] = mat[i*nCols+j];
  memcpy(mat, tmpmat, nRows*nCols*sizeof(float));
  free(tmpmat);
}


float get_phase(cuFloatComplex b)
{
  return atan2(b.y, b.x);
}

float cmag(cuFloatComplex a)
{
  return (sqrt(a.x*a.x+a.y*a.y));
}

cuFloatComplex cmul(cuFloatComplex a, cuFloatComplex b)
{
  cuFloatComplex result;
  result.x = a.x * b.x - a.y * b.y;
  result.y = a.x * b.y + a.y * b.x;
  return result;
}


void dumpBands(std::vector<GPUBuffer>* bands, int nx, int ny, int nz0)
{
  int n = 0;
  for (auto i = bands->begin(); i != bands->end(); ++i) {
    std::stringstream s;

    s << "band" << n << ".tif";

    CPUBuffer buf(nx*ny*nz0*sizeof(float));
    i->set(&buf, 0, buf.getSize(), 0);
    CImg<> band_tiff((float *) buf.getPtr(), nx, ny, nz0, 1, true);
    band_tiff.save_tiff(s.str().c_str());

    std::stringstream ss;
    ss << "band" << n << ".dat";
    std::ofstream os(ss.str().c_str());
    i->dump(os, nx, 0, nx * ny * sizeof(float));
    ++n;
  }
}

template<class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
  copy(v.begin(), v.end(), std::ostream_iterator<T>(std::cout, " ")); 
  return os;
}

// IMPLEMENTATION OF CLASS SIM_RECONSTRUCTOR FOLLOWS:
SIM_Reconstructor::SIM_Reconstructor(int argc, char **argv)
{
  // Call this constructor from a command-line based program

  m_argc = argc;
  m_argv = argv;

  SetDefaultParams(&m_myParams);
  
  // define all the commandline and config file options
  setupProgramOptions();
  po::positional_options_description p;
  p.add("input-file", 1);
  p.add("output-file", 1);
  p.add("otf-file", 1);

  // parse the commandline
  store(po::command_line_parser(argc, argv).
        options(m_progopts).positional(p).run(), m_varsmap);

  if (m_varsmap.count("help") || m_argc <= 1) {
    std::cout << "cudasirecon v" + version_number + " -- Written by Lin Shao and Talley Lambert. All rights reserved.\n\n";
    std::cout << "Usage:\n cudasirecon [--options] [option-arguments] input_file(or folder) output_file(or file name pattern) otf_file\n";
    std::cout << "\n List of options:\n";
    std::cout << m_progopts << std::endl;
    exit(0);
  }
  
  if (m_varsmap.count("version")) {
    std::cout << "Version " << version_number << std::endl;
    exit(0);
  }

  notify(m_varsmap);

  if (m_config_file != "") {
    std::ifstream ifs(m_config_file.c_str());
    if (!ifs) {
      std::cout << "can not open config file: " << m_config_file << "\n";
      std::cout << "proceed without it\n";
    }
    else {
      // parse config file
      store(parse_config_file(ifs, m_progopts), m_varsmap);
      notify(m_varsmap);
    }
  }

   // fill in m_myParams fields that have not been set yet
  setParams();

  printf("nphases=%d, ndirs=%d\n", m_myParams.nphases, m_myParams.ndirs);

  // In TIFF mode, m_myParams.ifiles refers to the name of the folder raw data resides in;
  // and m_myParams.ofiles refers to a pattern in all the raw data file names.
  // To help decide if input is in TIFF or MRC format, gather all TIFF files with
  // matching names under the same folder:
  if (boost::filesystem::is_directory(m_myParams.ifiles)) {
    m_all_matching_files = gatherMatchingFiles(m_myParams.ifiles, m_myParams.ofiles);
    // If there is no name-matching TIFF files, then presumably input is MRC
    m_myParams.bTIFF = (m_all_matching_files.size() > 0);
  }

  if (m_myParams.bTIFF) {
    // Suppress "unknown field" warnings
    TIFFSetWarningHandler(NULL);
    m_imgParams.ntimes = m_all_matching_files.size();
  }

  m_OTFfile_valid = false;
  openFiles();
  // deviceMemoryUsage();

  setup();
  if (!m_myParams.bTIFF)
    ::setOutputHeader(m_myParams, m_imgParams, m_in_out_header);
  //! Load flat-field correction data
  bgAndSlope(m_myParams, m_imgParams, &m_reconData);
}


//! used by shared library
SIM_Reconstructor::SIM_Reconstructor(int nx, int ny,
                                     int nimages,      // number of raw images per volume
                                     std::string configFileName)
{
  SetDefaultParams(&m_myParams);
  // define all the commandline and config file options
  setupProgramOptions();
  std::ifstream ifs(configFileName.c_str());
  if (ifs) {
    // parse config file
    store(parse_config_file(ifs, m_progopts), m_varsmap);
    notify(m_varsmap);
  }
  // fill in m_myParams fields that have not been set yet
  setParams();

  // TODO: allow more direct access
  // will throw CImgIOException if file cannot be opened
  m_otf_tiff.assign(m_myParams.otffiles.c_str());
  m_OTFfile_valid = true;

  m_myParams.bTIFF = true;  // not really: but we don't want it to use MRC stuff

  std::cout << "zoomz: " << m_myParams.zoomfact << std::endl;

  setup(nx, ny, nimages, 1);
  bgAndSlope(m_myParams, m_imgParams, &m_reconData);
}

SIM_Reconstructor::~SIM_Reconstructor()
{
}

int SIM_Reconstructor::setupProgramOptions()
{

  // Declare a group of options that will be 
  // allowed both on command line and in a config file
  m_progopts.add_options()
    ("input-file", po::value<std::string>(), "input file name (or data folder in TIFF mode)")
    ("output-file", po::value<std::string>(), "output file name (or filename pattern in TIFF mode)")
    ("otf-file", po::value<std::string>()->required(), "OTF file name")
    ("config,c", po::value<std::string>(&m_config_file)->default_value(""),
     "name of a config file with parameters in place of command-line options")
    ("ndirs", po::value<int>(&m_myParams.ndirs)->default_value(3),
     "number of SIM  directions")
    ("nphases", po::value<int>(&m_myParams.nphases)->default_value(5),
     "number of pattern phases per SIM direction")
    ("nordersout", po::value<int>(&m_myParams.norders_output)->default_value(0),
     "number of output SIM orders (must be <= nphases//2; safe to ignore usually)")
    ("angle0", po::value<float>(&m_myParams.k0startangle)->default_value(1.648),
     "angle of the first SIM angle in radians")
    ("ls", po::value<float>(&m_myParams.linespacing)->default_value(0.172, "0.172"),
     "line spacing of SIM pattern in microns")
    ("na", po::value<float>(&m_myParams.na)->default_value(1.2, "1.2"),
     "detection objective's numerical aperture")
    ("nimm", po::value<float>(&m_myParams.nimm)->default_value(1.33, "1.33"),
     "refractive index of immersion medium")
    ("wiener", po::value<float>(&m_myParams.wiener)->default_value(0.01, "0.01"),
     "Wiener constant; lower value leads to higher resolution and noise; playing with it extensively is strongly encouraged")
    ("otfcutoff", po::value<float>(&m_myParams.otfcutoff)->default_value(0.006, "0.006"),
     "otf threshold below which it'll be considered noise and not used in makeoverlaps")
    ("zoomfact", po::value<float>(&m_myParams.zoomfact)->default_value(2.),
     "lateral zoom factor in the output over the input images; leaving it at 2 should be fine in most cases")
    ("zzoom", po::value<int>(&m_myParams.z_zoom)->default_value(1),
     "axial zoom factor; almost never needed")
    ("background", po::value<float>(&m_myParams.constbkgd)->default_value(0.),
     "camera readout background")
    ("usecorr", po::value<std::string>(), "use a flat-field correction file if provided")
    ("forcemodamp", po::value< std::vector<float> >()->multitoken(),
     "modamps to be forced to these values; useful when image quality is low and auto-estimated modamps are below, say, 0.1")
    ("k0angles", po::value< std::string >(),
     "user these pattern vector k0 angles for all directions (instead of inferring the rest agnles from angle0)")
    ("otfRA", po::bool_switch(&m_myParams.bRadAvgOTF), "using rotationally averaged OTF; otherwise using 3/2D OTF for 3/2D raw data")
    ("otfPerAngle", po::bool_switch(&m_myParams.bOneOTFperAngle), "using one OTF per SIM angle; otherwise one OTF is used for all angles, which is how it's been done traditionally")
    ("fastSI", po::bool_switch(&m_myParams.bFastSIM),
     "SIM image is organized in Z->Angle->Phase order; otherwise assuming Angle->Z->Phase image order")
    ("k0searchAll", po::value<int>(&m_myParams.bUseTime0k0)->implicit_value(false),
     "search for k0 at all time points")
    ("norescale", po::value<int>(&m_myParams.do_rescale)->implicit_value(false), "no bleach correction")
    ("equalizez", po::bool_switch(&m_myParams.equalizez), "bleach correction for z")
    ("equalizet", po::bool_switch(&m_myParams.equalizet), "bleach correction for time")
    ("dampenOrder0", po::bool_switch(&m_myParams.bDampenOrder0),
     "dampen order-0 in final assembly; do not use for 2D SIM; good choice for high-background images")
    ("nosuppress", po::value<int>(&m_myParams.bSuppress_singularities)->implicit_value(false),
     "do not suppress DC singularity in the result (good choice for 2D/TIRF data)")
    ("nokz0", po::bool_switch(&m_myParams.bNoKz0), "not using kz=0 plane of the 0th order in the final assembly (mostly for debug)")
    ("gammaApo", po::value<float>(&m_myParams.apoGamma)->default_value(1.f),
     "output apodization gamma; 1.0 means triangular apo; lower value means less dampening of high-resolution info at the tradeoff of higher noise")
    ("explodefact", po::value<float>(&m_myParams.explodefact)->default_value(1.),
     "artificially exploding the reciprocal-space distance between orders by this factor (for debug)")
    ("nofilterovlps", po::value<bool>(&m_myParams.bFilteroverlaps)->implicit_value(false),
     "not filtering the overlaping region between bands (for debug)")
    ("saveprefiltered", po::value<std::string>(),
     "save separated bands (half Fourier space) into a file and exit (for debug)")
    ("savealignedraw", po::value<std::string>(),
     "save drift-fixed raw data (half Fourier space) into a file and exit (for debug)")
    ("saveoverlaps", po::value<std::string>(),
     "save overlap0 and overlap1 (real-space complex data) into a file and exit (for debug)")
    ("2lenses", po::bool_switch(&m_myParams.bTwolens), "I5S data")
    ("bessel", po::bool_switch(&m_myParams.bBessel), "bessel-SIM data")
    ("besselExWave", po::value<float>(&m_myParams.BesselLambdaEx)->default_value(0.488, "0.488"),
     "Bessel SIM excitation wavelength in microns")
    ("besselNA", po::value<float>(&m_myParams.BesselNA)->default_value(0.144, "0.144"),
     "Bessel SIM excitation NA)")
    ("deskew", po::value<float>(&m_myParams.deskewAngle)->default_value(0.0f),
     "Deskew angle; if not 0.0 then perform deskewing before processing")
    ("deskewshift", po::value<int>(&m_myParams.extraShift)->default_value(0),
     "If deskewed, the output image's extra shift in X (positive->left)")
    ("noRecon", po::bool_switch(&m_myParams.bNoRecon)->default_value(false),
     "No reconstruction will be performed; useful when combined with --deskew")
    ("cropXY", po::value<unsigned>(&m_myParams.cropXYto)->default_value(0),
     "Crop the X-Y dimension to this number; 0 means no cropping")

    ("xyres", po::value<float>(&m_imgParams.dxy)->default_value(0.1, "0.1"),
     "x-y pixel size (only used for TIFF files)")
    ("zres", po::value<float>(&m_imgParams.dz)->default_value(0.2, "0.2"),
     "z pixel size (only used for TIFF files)")
    ("zresPSF", po::value<float>(&m_myParams.dzPSF)->default_value(0.15, "0.15"),
     "z pixel size used in PSF TIFF files)")
    ("wavelength", po::value<short>(&m_imgParams.wave[0])->default_value(530),
     "emission wavelength in nanometers (only used for TIFF files)")
    ("writeTitle", po::bool_switch(&m_myParams.bWriteTitle),
     "Write command line to MRC/DV header (may cause issues with bioformats)")
    ("help,h", "produce help message")
    ("version", "show version")
    ;

  return 0;
}

int SIM_Reconstructor::setParams()
{
  // "otf-file" is required argument

  if (m_varsmap.count("input-file")) {
    m_myParams.ifiles = m_varsmap["input-file"].as<std::string>();
  }
  
  if (m_varsmap.count("output-file")) {
    m_myParams.ofiles = m_varsmap["output-file"].as<std::string>();
  }

  if (m_varsmap.count("otf-file")) {
    m_myParams.otffiles = m_varsmap["otf-file"].as<std::string>();
  }

  if (m_varsmap.count("usecorr")) {
    m_myParams.corrfiles = m_varsmap["usecorr"].as<std::string>();
    m_myParams.bUsecorr = 1;
  }

  if (m_varsmap.count("forcemodamp")) {
    m_myParams.forceamp = m_varsmap["forcemodamp"].as< std::vector<float> >();
  }

  if (m_varsmap.count("wiener")) {
    std::cout<< "wiener=" << m_myParams.wiener << std::endl;
  }

  if (m_varsmap.count("gammaApo")) {
    // because "gammaApo" has a default value set, this block is always entered;
    // and therefore, "apodizeoutput" is always 2 (triangular or gamma apo)
    m_myParams.apodizeoutput = 2;
    std::cout << "gamma = " << m_myParams.apoGamma << std::endl;
  }

  if (m_varsmap.count("saveprefiltered")) {
    m_myParams.fileSeparated = m_varsmap["saveprefiltered"].as<std::string>();
    m_myParams.bSaveSeparated = 1;
  }

  if (m_varsmap.count("savealignedraw")) {
    m_myParams.fileRawAligned = m_varsmap["savealignedraw"].as<std::string>();
    m_myParams.bSaveAlignedRaw = 1;
  }

  if (m_varsmap.count("saveoverlaps")) {
    m_myParams.fileOverlaps = m_varsmap["saveoverlaps"].as<std::string>();
    m_myParams.bSaveOverlaps = 1;
  }

  if (m_varsmap.count("k0angles")) {
    boost::char_separator<char> sep(",");
    boost::tokenizer<boost::char_separator<char>> tokens(m_varsmap["k0angles"].as< std::string >(), sep);
    for ( auto it = tokens.begin(); it != tokens.end(); ++it)
        m_myParams.k0angles.push_back(strtod(it->c_str(), NULL));
  }

  return 0;
}

void SIM_Reconstructor::openFiles() {
  if (!m_myParams.ifiles.size() || !m_myParams.ofiles.size())
    throw std::runtime_error(
        "Calling openFiles() but no input or output files were specified");

  if (!m_myParams.bTIFF && IMOpen(istream_no, m_myParams.ifiles.c_str(), "ro"))
    // TIFF files are not opened till loadAndRescaleImage()
    throw std::runtime_error("No matching input TIFF or MRC files not found");

  // // Create output file in MRC mode;
  // // in TIFF mode, output files are not created until writeResult() is called
  // if (!m_myParams.bTIFF)
  //   if (IMOpen(ostream_no, m_myParams.ofiles.c_str(), "new")) {
  //     std::cerr << "File " << m_myParams.ofiles << " can not be created.\n";
  //     throw std::runtime_error("Error creating output file");
  //   }

  if (m_myParams.bTIFF) {
    // will throw CImgIOException if file cannot be opened
    m_otf_tiff.assign(m_myParams.otffiles.c_str());
    m_OTFfile_valid = true;
  }
  else if (IMOpen(otfstream_no, m_myParams.otffiles.c_str(), "ro"))
    throw std::runtime_error("OTF file not found");
  else
    m_OTFfile_valid = true;
}

void SIM_Reconstructor::setup(unsigned nx, unsigned ny, unsigned nImages, unsigned nChannels)
// for mem buffer inputs
{
  m_imgParams.nx = nx;
  m_imgParams.nx_raw = nx;
  m_imgParams.ny = ny;
  m_imgParams.nz = nImages;
  m_imgParams.nwaves = nChannels;
  m_imgParams.nz /= m_imgParams.nwaves * m_myParams.nphases * m_myParams.ndirs;

  setup_common();
}

//! Find out memory requirements and call setup_common allocate memory
void SIM_Reconstructor::setup()
{
  if (m_myParams.bTIFF) {
    CImg<> tiff0(m_all_matching_files[0].c_str());
    m_imgParams.nx = 0;
    m_imgParams.nx_raw = tiff0.width();
    m_imgParams.ny = tiff0.height();
    m_imgParams.nz = tiff0.depth();
    m_imgParams.nwaves = tiff0.spectrum(); // multi-color not supported yet
    m_imgParams.nz /= m_imgParams.nwaves * m_myParams.nphases * m_myParams.ndirs;
  }
  else 
    ::loadHeader(m_myParams, &m_imgParams, m_in_out_header);
  printf("nx_raw=%d, ny=%d, nz=%d\n",
         m_imgParams.nx_raw, m_imgParams.ny, m_imgParams.nz);

  setup_common();
}

void SIM_Reconstructor::setup_common()
{
  // Do we need to keep this?? -lin
  // if (m_myParams.bTIFF)
  //   m_imgParams.nz0 = ::findOptimalDimension(m_imgParams.nz);

  m_zoffset = 0;
  if (m_myParams.nzPadTo) {  // 'nzPadTo' ~never used
    m_imgParams.nz0 = m_myParams.nzPadTo;
    m_zoffset = (m_imgParams.nz0 - m_imgParams.nz) / 2;
  } else {
    m_imgParams.nz0 = m_imgParams.nz;
  }

  //! Deskewing, if requested, will be performed after this function returns;
  //  "nx" is altered to be equal to "ny" and "dz" is multiplied by sin(deskewAngle)
  if (fabs(m_myParams.deskewAngle) > 0.) {
    if (m_myParams.deskewAngle <0) m_myParams.deskewAngle += 180.;
    m_imgParams.nx = findOptimalDimension(m_imgParams.nx_raw + m_imgParams.nz0 * fabs(cos(m_myParams.deskewAngle*M_PI/180.)) * m_imgParams.dz / m_imgParams.dxy);
//    m_imgParams.nx = m_imgParams.ny;
    m_imgParams.dz_raw = m_imgParams.dz;
    m_imgParams.dz *= fabs(sin(m_myParams.deskewAngle * M_PI/180.));
  }
  else {
    m_imgParams.nx = m_imgParams.nx_raw;
    m_imgParams.dz_raw = 0;
  }

  //! If cropping is requested, change nx, ny accordingly: (does this work or make sense?)
  if (m_myParams.cropXYto > 0 && m_myParams.cropXYto < m_imgParams.nx) {
    m_imgParams.nx = m_myParams.cropXYto;
    // m_imgParams.ny = m_myParams.cropXYto;
  }

  
  printf("nx=%d, ny=%d, nz=%d, nz0 = %d, nwaves=%d\n",
         m_imgParams.nx, m_imgParams.ny, m_imgParams.nz,
         m_imgParams.nz0, m_imgParams.nwaves);
  printf("dxy=%f, dz=%f um\n", m_imgParams.dxy, m_imgParams.dz);

  // // Initialize headers for intermediate output files if requested
  // if (!m_myParams.bTIFF) {
  //   if (m_myParams.bSaveAlignedRaw) {
  //     memcpy(&aligned_header, &m_in_out_header, sizeof(m_in_out_header));
  //     IMOpen(aligned_stream_no, m_myParams.fileRawAligned.c_str(), "new");
  //     aligned_header.mode = IW_FLOAT;
  //     aligned_header.nx = m_imgParams.nx;
  //     aligned_header.inbsym = 0;
  //     IMPutHdr(aligned_stream_no, &aligned_header);
  //   }
  //   if (m_myParams.bSaveSeparated) {
  //     memcpy(&sep_header, &m_in_out_header, sizeof(m_in_out_header));
  //     IMOpen(separated_stream_no, m_myParams.fileSeparated.c_str(), "new");
  //     sep_header.nx = m_imgParams.nx/2+1;    // saved will be separated FFTs
  //     sep_header.mode = IW_COMPLEX;
  //     sep_header.inbsym = 0;
  //     IMPutHdr(separated_stream_no, &sep_header);
  //   }
  //   if (m_myParams.bSaveOverlaps) {
  //     memcpy(&overlaps_header, &m_in_out_header, sizeof(m_in_out_header));
  //     IMOpen(overlaps_stream_no, m_myParams.fileOverlaps.c_str(), "new");
  //     overlaps_header.nz = m_imgParams.nz*2*m_myParams.ndirs*m_imgParams.ntimes*m_imgParams.nwaves;
  //     overlaps_header.nx = m_imgParams.nx;
  //     overlaps_header.num_waves = 2;  // save overlap 0 and 1 as wave 0 and 1 respectively
  //     overlaps_header.interleaved = WZT_SEQUENCE;
  //     overlaps_header.mode = IW_COMPLEX;   // saved will be full-complex overlaps in real space
  //     overlaps_header.inbsym = 0;
  //     IMPutHdr(overlaps_stream_no, &overlaps_header);
  //   }
  // }

  m_myParams.norders = 0;
  if (m_myParams.norders_output != 0) {
    m_myParams.norders = m_myParams.norders_output;
  } else {
    m_myParams.norders = m_myParams.nphases / 2 + 1;
  }

  std::cout << "nphases=" << m_myParams.nphases << ", norders=" <<
    m_myParams.norders << ", ndirs=" << m_myParams.ndirs << std::endl;

  if (m_OTFfile_valid)
    setupOTFsFromFile();

  ::allocSepMatrixAndNoiseVarFactors(m_myParams, &m_reconData);
  ::makematrix(m_myParams.nphases, m_myParams.norders, 0, 0,
      &(m_reconData.sepMatrix[0]), &(m_reconData.noiseVarFactors[0]));

  ::allocateImageBuffers(m_myParams, m_imgParams, &m_reconData);

  m_imgParams.inscale = 1.0 / (m_imgParams.nx * m_imgParams.ny * m_imgParams.nz0 *
      m_myParams.zoomfact * m_myParams.zoomfact * m_myParams.z_zoom * m_myParams.ndirs);
  m_reconData.k0 = std::vector<vector>(m_myParams.ndirs);
  m_reconData.k0_time0 = std::vector<vector>(m_myParams.ndirs);
  m_reconData.k0guess = std::vector<vector>(m_myParams.ndirs);
  float delta_angle = M_PI / m_myParams.ndirs;

  float k0magguess = 1.0 / m_myParams.linespacing; // in 1/um, not assuming square images sizes
  if (m_imgParams.nz > 1) {
    int nordersIn = m_myParams.nphases / 2 + 1;
    k0magguess /= nordersIn - 1; 
  }
  for (int i = 0; i < m_myParams.ndirs; ++i) {
    float k0angleguess;
    if (m_myParams.k0angles.size() < m_myParams.ndirs) {
      k0angleguess = m_myParams.k0startangle + i * delta_angle;
    } else {
      k0angleguess = m_myParams.k0angles[i];
    }
    m_reconData.k0guess[i].x = k0magguess * cos(k0angleguess); // in 1/um
    m_reconData.k0guess[i].y = k0magguess * sin(k0angleguess); // in 1/um
  }

  m_reconData.sum_dir0_phase0 = std::vector<double>(m_imgParams.nz * m_imgParams.nwaves);
  m_reconData.amp = std::vector<std::vector<cuFloatComplex> >(
      m_myParams.ndirs, std::vector<cuFloatComplex>(m_myParams.norders));
  for (int i = 0; i < m_myParams.ndirs; ++i) {
    m_reconData.amp[i][0].x = 1.0f;
    m_reconData.amp[i][0].y = 0.0f;
  }
}


int SIM_Reconstructor::processOneVolume()
{
  // process one SIM volume (i.e., for the time point timeIdx)

  m_reconData.bigbuffer.resize(0);
  m_reconData.outbuffer.resize(0);

  m_reconData.overlap0.resize(m_imgParams.nx * m_imgParams.ny * m_imgParams.nz *
      sizeof(cuFloatComplex));
  m_reconData.overlap0.setToZero();
  m_reconData.overlap1.resize(m_imgParams.nx * m_imgParams.ny * m_imgParams.nz *
      sizeof(cuFloatComplex));
  m_reconData.overlap1.setToZero();

  findModulationVectorsAndPhasesForAllDirections(m_zoffset,
      &m_myParams, m_imgParams, &m_driftParams, &m_reconData);

  m_reconData.overlap0.resize(0);
  m_reconData.overlap1.resize(0);

  // deviceMemoryUsage();

  if (saveIntermediateDataForDebugging(m_myParams))
    return 0;

  m_reconData.bigbuffer.resize((m_myParams.zoomfact * m_imgParams.nx) *
      (m_myParams.zoomfact * m_imgParams.ny) * (m_myParams.z_zoom * m_imgParams.nz0) *
      sizeof(cuFloatComplex));
  m_reconData.bigbuffer.setToZero();
  m_reconData.outbuffer.resize((m_myParams.zoomfact * m_imgParams.nx) *
      (m_myParams.zoomfact * m_imgParams.ny) * (m_myParams.z_zoom * m_imgParams.nz0) *
      sizeof(float));
  m_reconData.outbuffer.setToZero();

  // deviceMemoryUsage();

  for (int direction = 0; direction < m_myParams.ndirs; ++direction) {
    // dir_ is used in upcoming calls involving otf, to differentiate the cases
    // of common OTF and dir-specific OTF
    int dir_ = 0;
    if (m_myParams.bOneOTFperAngle) dir_ = direction;

    filterbands(direction, &m_reconData.savedBands[direction], m_reconData.k0,
                m_myParams.ndirs, m_myParams.norders, m_reconData.otf[dir_],
                m_imgParams.dxy, m_imgParams.dz, m_reconData.amp,
                m_reconData.noiseVarFactors, m_imgParams.nx, m_imgParams.ny,
                m_imgParams.nz0, m_imgParams.wave[0], &m_myParams);
    assemblerealspacebands(
        direction, &m_reconData.outbuffer, &m_reconData.bigbuffer,
        &m_reconData.savedBands[direction], m_myParams.ndirs,
        m_myParams.norders, m_reconData.k0, m_imgParams.nx, m_imgParams.ny,
        m_imgParams.nz0, m_imgParams.dxy, m_myParams.zoomfact,
        m_myParams.z_zoom, m_myParams.explodefact);
  }
  return 1;
}

void SIM_Reconstructor::loadAndRescaleImage(int timeIdx, int waveIdx)
{
  if (m_myParams.bBessel && m_myParams.bNoRecon) return;

  ::rescaleDriver(timeIdx, waveIdx, m_zoffset, &m_myParams, m_imgParams, 
                  &m_driftParams, &m_reconData);
}

void SIM_Reconstructor::setRaw(CImg<> &input, int it, int iw)
{
  rawImage.assign(input);
  loadImageData(it, iw);
}

void SIM_Reconstructor::setFile(int it, int iw)
{
  if (m_myParams.bTIFF) {
    rawImage.assign(m_all_matching_files[it].c_str());
  }
  else {
    rawImage.assign(m_imgParams.nx_raw, m_imgParams.ny, m_imgParams.nz * m_myParams.nphases * m_myParams.ndirs);
    for (unsigned sec=0; sec<rawImage.depth(); sec++) {
      IMPosnZWT(istream_no, sec, iw, it);
      IMRdSec(istream_no, rawImage.data(0, 0, sec));
    }
  }
  loadImageData(it, iw);
}

void SIM_Reconstructor::loadImageData(int it, int iw) {
  float dskewA = m_myParams.deskewAngle;
  // if ( dskewA<0) dskewA += 180.; already done in setup_common()
  float deskewFactor = 0;
  if (fabs(dskewA) > 0.0)
    deskewFactor =
        cos(dskewA * M_PI / 180.) * m_imgParams.dz_raw / m_imgParams.dxy;

  for (int direction = 0; direction < m_myParams.ndirs; ++direction) {
    // What does "PinnedCPUBuffer" really do? Using it here messes things up.
    CPUBuffer nxExtendedBuff(sizeof(float) * (m_imgParams.nx / 2 + 1) * 2 *
                             m_imgParams.ny);

    std::vector<GPUBuffer>* rawImages = &(m_reconData.savedBands[direction]);
    int z = 0;
    for (z = 0; z < m_imgParams.nz; ++z) {
      CImg<> rawSection(m_imgParams.nx_raw, m_imgParams.ny);

      int zsec;
      // First determine which section of the raw data to load
      if (m_myParams.bFastSIM) {
        // data organized into (nz, ndirs, nphases)
        zsec = (z * m_myParams.ndirs * m_myParams.nphases +
                direction * m_myParams.nphases);
      } else {
        // data organized into (ndirs, nz, nphases)
        zsec = direction * m_imgParams.nz * m_myParams.nphases +
               z * m_myParams.nphases;
      }

      for (int phase = 0; phase < m_myParams.nphases; ++phase) {
        if (m_myParams.bBgInExtHdr) {
          /* subtract the background value of each exposure stored in MRC files'
             extended header, indexed by the section number. */
          int extInts;
          float extFloats[3];
          IMRtExHdrZWT(istream_no, zsec, iw, it, &extInts, extFloats);
          m_reconData.backgroundExtra = extFloats[2];
        }

        load_and_flatfield(rawImage, zsec, rawSection.data(),
                           (float*)m_reconData.background.getPtr(),
                           m_reconData.backgroundExtra,
                           (float*)m_reconData.slope.getPtr(),
                           m_imgParams.inscale);
        // if deskewFactor > 0, then do deskew for current z; otherwise, simply
        // transfer from "rawSection" to "nxExtendedBuff". In either case,
        // "nxExtendedBuff" holds result of current z.
        deskewOneSection(rawSection, (float*)nxExtendedBuff.getPtr(), z,
                         m_imgParams.nz, m_imgParams.nx, deskewFactor,
                         m_myParams.extraShift);
        // Transfer the data from nxExtendedBuff to device buffer previously
        // allocated (see m_reconData.savedBands) std::cout << "z=" << z;
        nxExtendedBuff.set(
            &(rawImages->at(phase)), 0,
            (m_imgParams.nx / 2 + 1) * 2 * m_imgParams.ny * sizeof(float),
            (z + m_zoffset) * (m_imgParams.nx / 2 + 1) * 2 * m_imgParams.ny *
                sizeof(float));
        // std::cout << ", phase=" << phase << std::endl;
        ++zsec;
      }  // end for (phase)
    }    // end for (z)

#ifndef NDEBUG
    for (auto i = rawImages->begin(); i != rawImages->end(); ++i)
      assert(i->hasNaNs() == false);
#endif
    //!
    if (fabs(m_myParams.deskewAngle) > 0. && m_myParams.bNoRecon) {
      // To-do: add code to off load rawImages into rawtiff_uint16
      // if (m_myParams.bNoRecon) {
      //   CImg<unsigned short> rawtiff_uint16;
      //   rawtiff_uint16.save_tiff(makeOutputFilePath(m_all_matching_files[it],
      //                                               std::string("_deskewed")).c_str());
      //   std::cout << "Deskewed images written\n";
      //   return;
      // }
    }
  }  // end for (direction)
}

void SIM_Reconstructor::setupOTFsFromFile()
{
  determine_otf_dimensions();
  ::allocateOTFs(&m_myParams, m_reconData.sizeOTF, m_reconData.otf);
  loadOTFs();
}

void SIM_Reconstructor::determine_otf_dimensions()
{
  if (m_myParams.bTIFF) {
    uint32_t nxotf, nyotf, nzotf;
    nxotf = m_otf_tiff.width()/2;  // "/2" because of real and imaginary parts
    nyotf = m_otf_tiff.height();
    nzotf = m_otf_tiff.depth();

    /* determine nzotf, nxotf, nyotf, dkrotf, dkzotf based on dataset being 2D/3D and 
       flags bRadAvgOTF and bOneOTFperAngle */

    if (m_imgParams.nz == 1) {  // 2D
      m_myParams.nxotf = nxotf;
      if (m_myParams.bRadAvgOTF)
        m_myParams.nyotf = 1;
      else
        m_myParams.nyotf = nyotf;
      m_myParams.nzotf = 1;
      // m_myParams.dkrotf = xres;  // dkrotf's unit is 1/micron
      m_myParams.dkrotf = 1/(m_imgParams.dxy * (nxotf-1)*2);  // dkrotf's unit is 1/micron
      m_myParams.dkzotf = 1;
    }
    else {   // 3D
      if (m_myParams.bRadAvgOTF) {
        m_myParams.nzotf = nxotf;
        m_myParams.nxotf = nyotf;
        m_myParams.nyotf = 1;
        m_myParams.dkzotf = 1/(m_myParams.dzPSF * m_myParams.nzotf);
        m_myParams.dkrotf = 1/(m_imgParams.dxy * (m_myParams.nxotf-1)*2);
      }
      else {
        m_myParams.nzotf = nzotf / m_myParams.norders; // each order has a 3D OTF stack (non-negative kx half of Fourier space)
        m_myParams.nxotf = nxotf;
        m_myParams.nyotf = nyotf;
        m_myParams.dkzotf = 1/(m_myParams.dzPSF * m_myParams.nzotf);
        m_myParams.dkrotf = 1/(m_imgParams.dxy * m_myParams.nxotf);
      }
    }
  }
  else {
    int ixyz[3], mxyz[3], pixeltype;
    float min, max, mean;
    IW_MRC_HEADER otfheader;
    /* Retrieve OTF file header info */
    IMRdHdr(otfstream_no, ixyz, mxyz, &pixeltype, &min, &max, &mean);
    IMGetHdr(otfstream_no, &otfheader);
    IMAlCon(otfstream_no, 0);
    /* determine nzotf, nxotf, nyotf, dkrotf, dkzotf based on dataset
     * being 2D/3D and flags bRadAvgOTF and bOneOTFperAngle */
    if (m_imgParams.nz == 1) {  // 2D, ignore bOneOTFperAngle
      m_myParams.nxotf = otfheader.nx;
      if (m_myParams.bRadAvgOTF)
        m_myParams.nyotf = 1;
      else
        m_myParams.nyotf = otfheader.ny;
      m_myParams.nzotf = 1;
      m_myParams.dkrotf = otfheader.xlen;  // dkrotf's unit is 1/micron
      m_myParams.dkzotf = 1;
    } else {   // 3D, take into account bOneOTFperAngle
      if (m_myParams.bRadAvgOTF) {
        m_myParams.nzotf = otfheader.nx;
        m_myParams.nxotf = otfheader.ny;
        m_myParams.nyotf = 1;
        m_myParams.dkzotf = otfheader.xlen;
        m_myParams.dkrotf = otfheader.ylen;
      } else {
        // each order has a 3D OTF stack (non-negative kx half of Fourier
        // space)
        m_myParams.nzotf = otfheader.nz / m_myParams.norders;
        if (m_myParams.bOneOTFperAngle)
          m_myParams.nzotf /= m_myParams.ndirs;
        m_myParams.nxotf = otfheader.nx;
        m_myParams.nyotf = otfheader.ny;
        m_myParams.dkzotf = otfheader.zlen;
        m_myParams.dkrotf = otfheader.xlen;
      }
    }
  }
  printf("nzotf=%d, dkzotf=%f, nxotf=%d, nyotf=%d, dkrotf=%f\n",
      m_myParams.nzotf, m_myParams.dkzotf, m_myParams.nxotf, m_myParams.nyotf,
      m_myParams.dkrotf);

  //! sizeOTF are determined so that correct memory can be allocated for otf
  m_reconData.sizeOTF = m_myParams.nzotf*m_myParams.nxotf*m_myParams.nyotf;
}

int SIM_Reconstructor::loadOTFs()
{
  // If one OTF per direction is used, then load all dirs of OTF
  unsigned short nDirsOTF = 1;
  if (m_myParams.bOneOTFperAngle)
    nDirsOTF = m_myParams.ndirs;

  // Load OTF data, no matter 2D, 3D, radially averaged or not.
  if (m_myParams.bTIFF) {
    size_t nzotf = m_otf_tiff.depth();
    cuFloatComplex * otfTmp;
    CPUBuffer otfTmpBuffer(m_reconData.sizeOTF * sizeof(cuFloatComplex));

    for (int dir = 0; dir < nDirsOTF; dir++) {
      for (int i=0; i<m_myParams.norders; i++) {
        // If OTF file has multiple sections, then read them into otf[i]
        if (m_imgParams.nz == 1 || m_myParams.bRadAvgOTF) {
          if (nzotf > i + dir*m_myParams.norders) {
            // each section in OTF file is OTF of one order; so load that section into otf[dir][i]
            otfTmp = (cuFloatComplex *) m_otf_tiff.data(0, 0, i);
            otfTmpBuffer.setFrom((void*) otfTmp, 0, m_reconData.sizeOTF * sizeof(cuFloatComplex), 0);
            m_reconData.otf[dir][i].setFrom(otfTmpBuffer, 0, m_reconData.sizeOTF
                                            * sizeof(cuFloatComplex), 0);
          }
          else // If there's just 1 OTF image, do not read any more and just duplicate otf[0][0] into otf[*][i]
            m_reconData.otf[0][0].set(&(m_reconData.otf[dir][i]), 0,
                                      m_reconData.otf[0][0].getSize(), 0);
        }
        else {  // non-radially averaged 3D OTF
          otfTmp = (cuFloatComplex *) m_otf_tiff.data(0, 0, i * m_myParams.nzotf);
          otfTmpBuffer.setFrom((void*) otfTmp, 0, m_reconData.sizeOTF * sizeof(cuFloatComplex), 0);
          m_reconData.otf[dir][i].setFrom(otfTmpBuffer, 0,
                                          m_reconData.sizeOTF * sizeof(cuFloatComplex),
                                          0);
        }
      }
    }
  }
  else {

    int ixyz[3], mxyz[3], pixeltype;
    float min, max, mean;
    IW_MRC_HEADER otfheader;
    /* Retrieve OTF file header info */
    IMRdHdr(otfstream_no, ixyz, mxyz, &pixeltype, &min, &max, &mean);
    IMGetHdr(otfstream_no, &otfheader);

    CPUBuffer otfTmp(m_reconData.sizeOTF * sizeof(cuFloatComplex));

    for (int dir = 0; dir< nDirsOTF; dir++) {
      for (int i = 0; i < m_myParams.norders; i++) {
        /* If OTF file has multiple sections, then read them into otf[i]; */
        if (m_imgParams.nz == 1 || m_myParams.bRadAvgOTF) {
          if (otfheader.nz > i + dir*m_myParams.norders) {
            // each section in OTF file is OTF of one order; so load that section into otf[dir][i]
            IMRdSec(otfstream_no, otfTmp.getPtr());
            otfTmp.set(&(m_reconData.otf[dir][i]), 0, otfTmp.getSize(), 0);
          } else {
            /* If it's 2D image, do not read any more and just
             * duplicate otf[0][0] into otf[*][i] */
            m_reconData.otf[0][0].set(&(m_reconData.otf[dir][i]), 0,
                                m_reconData.otf[0][0].getSize(), 0);
          }
        } else {  // non-radially averaged 3D OTF
          for (int z = 0; z < m_myParams.nzotf; ++z) {
            IMRdSec(otfstream_no, otfTmp.getPtr());
            otfTmp.set(&(m_reconData.otf[dir][i]),
                       0, m_myParams.nxotf * m_myParams.nyotf * sizeof(cuFloatComplex),
                       z * m_myParams.nxotf * m_myParams.nyotf * sizeof(cuFloatComplex));
          }
        }
      }
    }
    IMClose(otfstream_no);
  }

#ifndef NDEBUG
  for (int dir = 0; dir < nDirsOTF; dir++)
    for (auto i = m_reconData.otf[dir].begin();
         i != m_reconData.otf[dir].end(); ++i)
      assert(i->hasNaNs() == false);
#endif
  return 1;
}

// Transfer m_reconData from GPU to result buffer
void SIM_Reconstructor::getResult(float * result) {
  CPUBuffer outbufferHost((m_myParams.zoomfact * m_imgParams.nx) *
                          (m_myParams.zoomfact * m_imgParams.ny) *
                          (m_myParams.z_zoom * m_imgParams.nz0) *
                          sizeof(float));

  m_reconData.outbuffer.set(&outbufferHost, 0, outbufferHost.getSize(), 0);

  CImg<> outCimg(
      (float*)outbufferHost.getPtr(), m_myParams.zoomfact * m_imgParams.nx,
      m_myParams.zoomfact * m_imgParams.ny, m_myParams.z_zoom * m_imgParams.nz0,
      true);

  memcpy(result, outCimg.data(), outCimg.size() * sizeof(float));
}

void SIM_Reconstructor::writeResult(int it, int iw)
{
  CPUBuffer outbufferHost(
      (m_myParams.zoomfact * m_imgParams.nx) *
      (m_myParams.zoomfact * m_imgParams.ny) *
      (m_myParams.z_zoom * m_imgParams.nz0) *
      sizeof(float));
  // if (it==0)
  //   computeAminAmax(&reconData.outbuffer, m_myParams.zoomfact * m_imgParams.nx,
  //                   m_myParams.zoomfact * m_imgParams.ny, m_myParams.z_zoom * m_imgParams.nz,
  //                   &minval, &maxval);
  m_reconData.outbuffer.set(&outbufferHost, 0, outbufferHost.getSize(), 0);

#ifndef __clang__
  double t1 = omp_get_wtime();
#endif

  if (m_myParams.bTIFF) {
  CImg<> outCimg((float*) outbufferHost.getPtr(),
    m_myParams.zoomfact * m_imgParams.nx,
    m_myParams.zoomfact * m_imgParams.ny,
    m_myParams.z_zoom * m_imgParams.nz0, true);  // "true" means outCimg does not allocate host memory

  outCimg.save(makeOutputFilePath(m_all_matching_files[it], std::string("_proc")).c_str());
  }
  else {
    float maxval = -FLT_MAX;
    float minval = FLT_MAX;

    float* ptr = ((float*)outbufferHost.getPtr()) +
      (int)(m_zoffset * m_myParams.z_zoom * m_imgParams.nx * m_imgParams.ny *
            m_myParams.zoomfact * m_myParams.zoomfact);
    for (int i = 0; i < m_imgParams.nz * m_myParams.z_zoom; ++i) {
      IMWrSec(ostream_no, ptr);
      if (it == 0) {
        for (int j = 0;
             j < (int)(m_imgParams.nx * m_imgParams.ny *
                       m_myParams.zoomfact * m_myParams.zoomfact);
             ++j) {
          if (ptr[j] > maxval) {
            maxval = ptr[j];
          } else if (ptr[j] < minval) {
            minval = ptr[j];
          }
        }
      }
      ptr += (int)(m_myParams.zoomfact * m_imgParams.nx *
                   m_myParams.zoomfact * m_imgParams.ny);
    }
    if (it == 0 && iw == 0) {
      m_in_out_header.amin = minval;
      m_in_out_header.amax = maxval;
    }
  }

#ifndef __clang__
  double t2 = omp_get_wtime();
  printf("amin, amax took: %f s\n", t2 - t1);
#endif

  printf("Time point %d, wave %d done\n", it, iw);
}


void SIM_Reconstructor::closeFiles()
{
  if (m_myParams.bTIFF)
    {}
  else {
    ::IMClose(istream_no);
    ::saveCommandLineToHeader(m_argc, m_argv, m_in_out_header, m_myParams);
    ::IMClose(ostream_no);
  }
}

void load_and_flatfield(CImg<> &cimg, int section_no, float *bufDestiny,
                        float *background, float backgroundExtra,
                        float *slope, float inscale)
{
  float *buffer = cimg.data(0, 0, section_no);
  int nx = cimg.width();
  int ny = cimg.height();
  unsigned extraX = 0;
#pragma omp parallel for
  for (int l=0; l<ny; l++) {
    for (int k=0; k<nx; k++) {
      bufDestiny[l*(nx+extraX) + k] = (buffer[l*nx+k] - background[l*nx+k] - backgroundExtra) * slope[l*nx+k] * inscale;
    }
    for(int k=nx; k<nx+extraX; k++)
      bufDestiny[l*(nx+extraX) + k] = 0.0;
  }
}

void deskewOneSection(CImg<> &rawSection, float* nxp2OutBuff, int z, int nz,
                      int nx_out, float deskewFactor, int extraShift)
{
  unsigned nx_in = rawSection.width();
  unsigned ny_in = rawSection.height();
  float *in = rawSection.data();
  int nx_out_ext = (nx_out/2+1)*2;
//  std::cout << "In deskewOneSection, deskewFactor=" << deskewFactor << std::endl;
#pragma omp parallel for
  for (auto xout=0; xout<nx_out_ext; xout++) {
    float xin = xout;
    if (fabs(deskewFactor) > 0)
      xin = xout - nx_out/2. + extraShift - deskewFactor*(z-nz/2.) + nx_in/2.;

    for (int y=0; y<ny_in; y++) {
      unsigned indout = y * nx_out_ext + xout;
      if (xin >= 0 && xin < nx_in) {

        unsigned indin = y * nx_in + (unsigned int) floor(xin);

        float offset = xin - floor(xin);
        nxp2OutBuff[indout] = (1-offset) * in[indin] + offset * in[indin+1];
      }
      else
        nxp2OutBuff[indout] = 0.;
    }
  }
}
