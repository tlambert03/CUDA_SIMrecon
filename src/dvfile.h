#pragma once

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>

/* Data Types. */
#define IW_AS_IS -1
#define IW_BYTE 0
#define IW_SHORT 1
#define IW_FLOAT 2
#define IW_COMPLEX_SHORT 3
#define IW_COMPLEX 4
#define IW_EMTOM 5
#define IW_USHORT 6
#define IW_LONG 7

/* image sequence definitions */
#define ZTW_SEQUENCE 0 /* "non-interleaved", by defintion       */
#define WZT_SEQUENCE 1 /* "interleaved", from R3D and others    */
#define ZWT_SEQUENCE 2 /* new sequence. Unsupported as of 11/97 */

static const std::unordered_map<int, size_t> pixelTypeSizes = {
    {IW_BYTE, sizeof(uint8_t)},      {IW_SHORT, sizeof(int16_t)},
    {IW_FLOAT, sizeof(float)},       {IW_COMPLEX_SHORT, 2 * sizeof(int16_t)},
    {IW_COMPLEX, 2 * sizeof(float)}, {IW_EMTOM, sizeof(int16_t)},
    {IW_USHORT, sizeof(uint16_t)},   {IW_LONG, sizeof(int32_t)}};

typedef struct IW_MRC_Header {
  int32_t nx, ny, nz;        // nz : nplanes*nwave*ntime
  int32_t mode;              // data type
  int32_t nxst, nyst, nzst;  // index of the first col/row/section
  int32_t mx, my, mz;        // number of intervals in x/y/z
  float xlen, ylen, zlen;    // pixel spacing for x/y/z
  float alpha, beta, gamma;  // cell angles
  int32_t mapc, mapr, maps;  // column/row/section axis
  float amin, amax, amean;   // min/max/mean intensity
  int32_t ispg, inbsym;      // space group number, number of bytes in extended header
  int16_t nDVID, nblank;     // ID value, unused
  int32_t ntst;              // starting time index (used for time series data)
  char ibyte[24];            // 24 bytes of blank space
  int16_t nint, nreal;       // number of integers/floats in extended header per section
  int16_t nres, nzfact;      // number of sub-resolution data sets, reduction quotient for z axis
  float min2, max2, min3, max3, min4, max4;  // min/max intensity for 2nd, 3rd, 4th wavelengths
  int16_t file_type, lens, n1, n2, v1, v2;   // file type, lens ID, n1, n2, v1, v2
  float min5, max5;                          // min/max intensity for 5th wavelength
  int16_t num_times;                         // number of time points
  int16_t interleaved;                       // (0 = ZTW, 1 = WZT, 2 = ZWT)
  float tilt_x, tilt_y, tilt_z;              // x/y/z axis tilt angles
  int16_t num_waves, iwav1, iwav2, iwav3, iwav4, iwav5;  // number & values of wavelengths
  float zorig, xorig, yorig;                             // z/x/y origin
  int32_t nlab;                                          // number of titles
  char label[800];

  std::string sequence_order() const {
    switch (interleaved) {
      case 0: return "CTZ";
      case 1: return "TZC";
      case 2: return "TCZ";
      default: return "CTZ";
    }
  }

  int num_planes() const { return nz / (num_waves ? num_waves : 1) / (num_times ? num_times : 1); }

  std::string image_type() const {
    switch (file_type) {
      case 0:
      case 100: return "NORMAL";
      case 1: return "TILT_SERIES";
      case 2: return "STEREO_TILT_SERIES";
      case 3: return "AVERAGED_IMAGES";
      case 4: return "AVERAGED_STEREO_PAIRS";
      case 5: return "EM_TILT_SERIES";
      case 20: return "MULTIPOSITION";
      case 8000: return "PUPIL_FUNCTION";
      default: return "UNKNOWN";
    }
  }

  void print() {
    std::cout << "Header:" << std::endl;
    std::cout << "  Dimensions: " << ny << "x" << nx << "x" << num_planes() << std::endl;
    std::cout << "  Number of wavelengths: " << num_waves << std::endl;
    std::cout << "  Number of time points: " << num_times << std::endl;
    std::cout << "  Pixel size: " << mode << std::endl;
    std::cout << "  Pixel spacing: " << xlen << "x" << ylen << "x" << zlen << std::endl;
    std::cout << "  mxyz: " << mx << "x" << my << "x" << mz << std::endl;
    std::cout << "  Cell angles: " << alpha << "x" << beta << "x" << gamma << std::endl;
    std::cout << "  Min/Max/Mean: " << amin << "/" << amax << "/" << amean << std::endl;
    std::cout << "  Image type: " << image_type() << std::endl;
    std::cout << "  Sequence order: " << sequence_order() << std::endl;
  }

} IW_MRC_HEADER, *IW_MRC_HEADER_PTR;

class DVFile {
 private:
  std::unique_ptr<std::fstream> _file;
  std::string _path;
  bool _big_endian;
  IW_MRC_Header hdr;
  bool closed = true;

  // Private default constructor
  DVFile() = default;

 public:
  // Constructor for opening an existing file
  DVFile(const std::string& path) {
    _path = path;
    _file = std::make_unique<std::fstream>(path, std::ios::in | std::ios::out | std::ios::binary);
    if (!_file->is_open()) {
      throw std::runtime_error("Failed to open file");
    }

    // Determine byte order
    _file->seekg(24 * 4);
    char dvid[2];
    _file->read(dvid, 2);
    if (dvid[0] == (char)0xA0 && dvid[1] == (char)0xC0) {
      _big_endian = false;
    } else if (dvid[0] == (char)0xC0 && dvid[1] == (char)0xA0) {
      _big_endian = true;
    } else {
      throw std::runtime_error(path + " is not a recognized DV file.");
    }

    // Read header
    _file->seekg(0);
    _file->read(reinterpret_cast<char*>(&hdr), sizeof(IW_MRC_Header));
    closed = false;
  }

  static std::unique_ptr<DVFile> createNew(const std::string& path) {
    std::unique_ptr<DVFile> dvfile(new DVFile());  // Use the private default constructor
    dvfile->_path = path;
    dvfile->_file = std::make_unique<std::fstream>(
        path, std::ios::in | std::ios::out | std::ios::binary | std::ios::trunc);

    if (!dvfile->_file->is_open()) {
      throw std::runtime_error("Failed to create file");
    }

    dvfile->closed = false;
    return dvfile;
  }

  // this is only here for the IVE API
  void setCurrentZWT(int z, int w, int t) {
    _validateZWT(z, w, t);

    size_t frame_size = hdr.ny * hdr.nx * getPixelSize();
    int header_size = 1024 + hdr.inbsym;
    int section_offset = (t * hdr.num_waves * hdr.num_planes() + w * hdr.num_planes() + z);
    _file->seekg(header_size + section_offset * frame_size);
  }

  void readSec(void* array) {
    if (closed) {
      throw std::runtime_error("Cannot read from closed file. Please reopen with .open()");
    }
    size_t frame_size = hdr.ny * hdr.nx * getPixelSize();
    _file->read(reinterpret_cast<char*>(array), frame_size);
  }

  void readSec(void* array, int t, int w, int z) {
    setCurrentZWT(z, w, t);
    readSec(array);
  }

  void writeSection(const void* array) {
    if (closed) {
      throw std::runtime_error("Cannot write to closed file. Please reopen with .open()");
    }
    size_t frame_size = hdr.ny * hdr.nx * getPixelSize();
    _file->write(reinterpret_cast<const char*>(array), frame_size);
  }

  size_t getPixelSize() { return pixelTypeSizes.at(hdr.mode); }

  void open() {
    if (closed) {
      _file->open(_path, std::ios::binary);
      if (!_file->is_open()) {
        throw std::runtime_error("Failed to open file");
      }
      closed = false;
    }
  }

  void close() {
    if (!closed) {
      _file->close();
      closed = true;
    }
  }

  ~DVFile() { close(); }

  std::string getPath() const { return _path; }

  IW_MRC_Header getHeader() const { return hdr; }

  void putHeader(const IW_MRC_Header& header) {
    if (closed) {
      throw std::runtime_error("Cannot write to closed file. Please reopen with .open()");
    }
    _file->seekp(0);
    _file->write(reinterpret_cast<const char*>(&header), sizeof(IW_MRC_Header));
    hdr = header;
  }

  bool isClosed() const { return closed; }

  std::map<std::string, int> sizes() {
    int num_real_z =
        hdr.nz / (hdr.num_waves ? hdr.num_waves : 1) / (hdr.num_times ? hdr.num_times : 1);
    std::map<std::string, int> d = {{"T", hdr.num_times},
                                    {"C", hdr.num_waves},
                                    {"Z", num_real_z},
                                    {"Y", hdr.ny},
                                    {"X", hdr.nx}};
    std::string axes = hdr.sequence_order() + "YX";
    std::map<std::string, int> sizes;
    for (char c : axes) {
      std::string key(1, c);  // Convert char to std::string
      sizes[key] = d[key];
    }
    return sizes;
  }

 private:
  void _validateZWT(int z, int w, int t) {
    if (t >= hdr.num_times) {
      throw std::runtime_error("Time index out of range");
    }
    if (w >= hdr.num_waves) {
      throw std::runtime_error("Wavelength index out of range");
    }
    if (z >= hdr.num_planes()) {
      throw std::runtime_error("Section index out of range");
    }
  }
};

//////////////////////////////////////////////////////////////////////////////
// IVE API
//////////////////////////////////////////////////////////////////////////////

// stream -> DVFile
static std::map<int, std::unique_ptr<DVFile>> dvfile_map;

inline DVFile& getDVFile(int istream) {
  auto it = dvfile_map.find(istream);
  if (it == dvfile_map.end()) {
    throw std::runtime_error("Stream not found: " + std::to_string(istream));
  }
  return *(it->second);
}

/**
 * @brief Open an image file and attach it to a stream.
 *
 * @param istream The input stream to be used for the operation.
 * @param name The name of the file to be opened.
 * @param attrib The file mode.
 *  - "ro" Opens an existing file or image window for reading only.
 *  - "new" Creates a file or image window and opens it for reading and writing.
 * @return int 0 if successful and 1 if not.
 */
inline int IMOpen(int istream, const char* name, const char* attrib) {
  // Check if the stream identifier is already in use and close it if necessary
  if (dvfile_map.find(istream) != dvfile_map.end()) {
    dvfile_map[istream]->close();
    dvfile_map.erase(istream);
    std::cerr << "Warning: Reusing stream identifier " << istream << ". Previous stream closed."
              << std::endl;
  }

  if (std::string(attrib) == "ro") {
    try {
      dvfile_map[istream] = std::make_unique<DVFile>(name);
    } catch (const std::exception& e) {
      std::cerr << "Error: " << e.what() << std::endl;
      return 1;  // Return non-zero to indicate failure
    }
  } else if (std::string(attrib) == "new") {
    try {
      dvfile_map[istream] = DVFile::createNew(name);
    } catch (const std::exception& e) {
      std::cerr << "Error: " << e.what() << std::endl;
      return 1;  // Return non-zero to indicate failure
    }
  } else {
    std::cerr << "Unknown file mode: " << attrib << std::endl;
    return 1;  // Return non-zero to indicate failure
  }
  return 0;  // Return 0 to indicate success
}

inline void IMClose(int istream) {
  // call destructor of DVFile
  dvfile_map.erase(istream);
}

inline void IMGetHdr(int istream, IW_MRC_HEADER* header) {
  *header = getDVFile(istream).getHeader();
}

inline void IMRdHdr(int istream, int ixyz[3], int mxyz[3], int* imode, float* min, float* max,
                    float* mean) {
  IW_MRC_HEADER header;
  IMGetHdr(istream, &header);
  ixyz[0] = header.nx;
  ixyz[1] = header.ny;
  ixyz[2] = header.nz;
  mxyz[0] = header.mx;
  mxyz[1] = header.my;
  mxyz[2] = header.mz;
  *imode = header.mode;
  *min = header.amin;
  *max = header.amax;
  *mean = header.amean;
}

/**
 * @brief Set the image conversion mode during read/write operations from image storage.
 *
 * By default in IVE, images that are read from image storage are converted to
 * 4-byte floating-point data. Similarly, when images are written to image
 * storage they are converted to the data type indicated by the image data type
 * associated with the corresponding stream (see IMAlMode). The default in IVE
 * is ConversionFlag=TRUE.
 * We, however, don't ever convert the data type of the image data. So for now,
 * this is a no-op.
 *
 * @param istream The input stream to be used for the operation.
 * @param flag The flag indicating the type of operation to be performed.
 */
inline void IMAlCon(int istream, int flag) {
  // if flag is 1, warn:
  if (flag == 1) {
    std::cerr << "Warning: IMAlCon is not implemented. ConversionFlag=TRUE is not supported."
              << std::endl;
  }
}

/**
 * @brief Change the image titles.
 *
 * @param istream The input stream to be used for the operation.
 * @param titles The titles to be changed.  Contains at least NumTitles title strings, each of
 * which must contain exactly 80 characters.
 * @param num_titles The number of titles to be changed.
 */
inline void IMAlLab(int istream, const char* labels, int nl) {
  std::cerr << "Warning: IMAlLab is not implemented." << std::endl;
}

/**
 * @brief  Enable or disable printing to standard output ("stdout").
 *
 * Certain IM functions will print information to stdout if Format=TRUE, which is the default. To
 * disable printing, set flag=FALSE.
 *
 * @param flag The flag indicating the type of operation to be performed.
 */
inline void IMAlPrt(int flag) {
  if (flag == 1) {
    std::cerr << "Warning: IMAlPrt is not implemented." << std::endl;
  }
}

/**
 * @brief Position the read/write point at a particular Z, W, T section.
 *
 * If the stream points to a scratch window, IMPosnZWT can only change the destination wavelength.
 *
 * @param istream The input stream to be used for the operation.
 * @param iz The Z section number.
 * @param iw The wavelength number.
 * @param it The time-point number.
 * @return int 0 if successful and 1 if not.
 */
inline int IMPosnZWT(int istream, int iz, int iw, int it) {
  DVFile& dvfile = getDVFile(istream);

  try {
    dvfile.setCurrentZWT(iz, iw, it);
    return 0;
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}

/**
 * @brief Read the next section.
 *
 * Reads the next section into ImgBuffer and advances the file pointer to the
 * section after that. The results are undefined if ImgBuffer does not have at
 * least nx * ny elements or the file pointer does not point to the beginning of
 * a section.
 * In most cases, ImgBuffer will contain floating-point data. When
 * image conversion is off, however, the data type of ImgBuffer should
 * correspond to whatever data type is actually stored. For example, if the
 * image data are stored as 16-bit integers, then ImgBuffer should point to a
 * 16-bit buffer. See IMAlCon.
 *
 * @param istream The input stream to be used for the operation.
 * @param ImgBuffer The image buffer to store the data.
 */
inline void IMRdSec(int istream, void* ImgBuffer) {
  try {
    getDVFile(istream).readSec(ImgBuffer);
  } catch (const std::runtime_error& e) {
    std::cerr << "Error reading section: " << e.what() << std::endl;
    // Handle the error appropriately, e.g., by rethrowing or returning an error code
    throw;  // Rethrow the exception to propagate it further
  }
}

inline void IMWrSec(int istream, const void* array) {
  try {
    getDVFile(istream).writeSection(array);
  } catch (const std::runtime_error& e) {
    std::cerr << "Error writing section: " << e.what() << std::endl;
    // Handle the error appropriately, e.g., by rethrowing or returning an error code
    throw;  // Rethrow the exception to propagate it further
  }
}

/**
 * @brief Put an entire header into a stream.
 *
 * Header should point to a memory location that contains a complete image header.
 *
 * @param istream The input stream to be used for the operation.
 * @param header The header to be saved in the stream.
 */
inline void IMPutHdr(int istream, const IW_MRC_HEADER* header) {
  getDVFile(istream).putHeader(*header);
}

/**
 * @brief Write the image header to the storage device.
 *
 * Write image header associated with StreamNum to the storage device. Use this
 * function to save the results of all IMAl functions. Header modifications are
 * not saved until IMWrHdr is used! The contents of Title will be saved in the
 * header according to the method indicated by ntflag. The minimum, maximum, and
 * mean intensity - Min, Max, and Mean, respectively - of wavelength 0 are also
 * saved in the header every time IMWrHdr is used.
 *
 * @param istream The input stream to be used for the operation.
 * @param title The title to be saved in the header.
 * @param ntflag The method to be used to save the title.
 *  - 0: use Title as the only title
 *  - 1: add Title to the end of the list
 */
inline void IMWrHdr(int istream, const char title[80], int ntflag, float dmin, float dmax,
                    float dmean) {
  DVFile& dvfile = getDVFile(istream);
  IW_MRC_HEADER header = dvfile.getHeader();
  header.amin = dmin;
  header.amax = dmax;
  header.amean = dmean;
  if (ntflag == 0) {
    // use Title as the only title
    strncpy(header.label, title, 80);
  } else if (ntflag == 1) {
    // FIXME  this is wrong.
    // Append title to the end of the list
    std::string new_title = title;
    new_title += " ";
    new_title += header.label;
    strncpy(header.label, new_title.c_str(), 80);
  } else {
    throw std::runtime_error("Invalid ntflag: " + std::to_string(ntflag));
  }

  dvfile.putHeader(header);
}

/**
 * @brief Return extended header values for a particular Z section, wavelength,
 * and time-point.
 *
 * The integer and floating-point values for the requested Z section (ZSecNum),
 * wavelength (WaveNum), and time-point (TimeNum) are returned in IntValues and
 * FloatValues, respectively.
 *
 */
inline void IMRtExHdrZWT(int istream, int iz, int iw, int it, int ival[], float rval[]) {
  std::cerr << "Warning: IMRtExHdrZWT is not implemented." << std::endl;
}
