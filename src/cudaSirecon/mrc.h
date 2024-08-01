// MRC-specific functionality
#include "cudaSirecon.h"


// MRC header parser
void loadHeader(const ReconParams& params, ImageParams* imgParams,
                IW_MRC_HEADER& header) {
  int ixyz[3];
  int mxyz[3];
  int pixeltype;
  float min;
  float max;
  float mean;
  IMRdHdr(istream_no, ixyz, mxyz, &pixeltype, &min, &max, &mean);
  IMGetHdr(istream_no, &header);
  imgParams->nx = 0;  // deferred till deskewing situation is figured out
  imgParams->nx_raw = header.nx;
  imgParams->ny = header.ny;
  imgParams->nz = header.nz / (header.num_waves * header.num_times);
  /* header.nz is defined as the total number of sections =
   * nz*nwaves*ntimes (* ndirs*nphases in our case) */
  imgParams->nz /= params.nphases * params.ndirs;
  imgParams->nwaves = header.num_waves;
  imgParams->ntimes = header.num_times;
  imgParams->wave[0] = header.iwav1;
  imgParams->wave[1] = header.iwav2;
  imgParams->wave[2] = header.iwav3;
  imgParams->wave[3] = header.iwav4;
  imgParams->wave[4] = header.iwav5;
  //! dxy: lateral pixel size; dz: axial pixel size; both in microns */
  imgParams->dxy = header.ylen;
  imgParams->dz = header.zlen;
}

void setOutputHeader(const ReconParams& myParams, const ImageParams& imgParams,
                     IW_MRC_HEADER& header) {
  header.mode = IW_FLOAT;
  header.nz =
      imgParams.nz * imgParams.nwaves * imgParams.ntimes * myParams.z_zoom;
  header.nx = imgParams.nx * myParams.zoomfact;
  header.ny *= myParams.zoomfact;
  header.xlen /= myParams.zoomfact;
  header.ylen /= myParams.zoomfact;
  header.zlen /= myParams.z_zoom;
  header.inbsym = 0;
  IMPutHdr(ostream_no, &header);
  IMAlCon(ostream_no, 0);
}

void saveCommandLineToHeader(int argc, char** argv, IW_MRC_HEADER& header,
                             const ReconParams& myParams) {
  char titles[1000];
  titles[0] = '\0';
  for (int i = 3; i < argc; ++i) {
    strcat(titles, argv[i]);
    strcat(titles, " ");
  }
  if (myParams.bWriteTitle) {
    IMAlLab(ostream_no, titles, strlen(titles) / 80 + 1);
  }
  IMWrHdr(ostream_no, header.label, 1, header.amin, header.amax, header.amean);
}

int saveIntermediateDataForDebugging(const ReconParams& params) {
  if (params.bSaveSeparated) {
    IMWrHdr(separated_stream_no, "separated bands of all directions", 1, 0, 1, 0);
    IMClose(separated_stream_no);
  }
  if (params.bSaveAlignedRaw) {
    IMWrHdr(aligned_stream_no, "drift-corrected raw images", 1,
            aligned_header.amin, aligned_header.amax, aligned_header.amean);
    IMClose(aligned_stream_no);
  }
  if (params.bSaveOverlaps) {
    IMWrHdr(overlaps_stream_no, "overlaps in real space", 1, 0, 1, 0);
    IMClose(overlaps_stream_no);
  }
  if (params.bSaveSeparated || params.bSaveAlignedRaw || params.bSaveOverlaps) {
    printf(
        "\nQuit processing because either bSaveSeparated, "
        "bSaveAlignedRaw, or bSaveOverlaps is TRUE\n");
    return 1;
  }
  return 0;
}