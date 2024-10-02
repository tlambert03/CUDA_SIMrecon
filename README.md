# cudasirecon

Mats Gustafsson & Lin Shao's 3-beam SIM reconstruction software, with CUDA acceleration.

Algorithm as described in [Gustafsson et al (2008) *Biophys. J.* **94(12)**: 4957–4970. doi: 10.1529/biophysj.107.120345](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2397368/)


## Installation

Packages for Linux and Windows on conda-forge:

```bash
$ conda install -c conda-forge cudasirecon

# you may wish to put it in a new environment ("sim" here):
$ conda create -n sim -y -c conda-forge cudasirecon
$ conda activate sim
```

In the examples shown above and below, *Shell* command line is used. The same works mostly in a Anaconda Prompt window on Windows, except for using`\`in place of`/`as the folder name dividers

## Usage
### For`cudasirecon`

```bash
# the binary will be available as
$ cudasirecon

# for command line help
$ cudasirecon --help
# or:
$ cudasirecon -h

# a typical call for tiff files
$ cudasirecon some/folder file_pattern otf.tif -c config_file

# a typical call for dv/mrc files
$ cudasirecon data.dv data-proc.dv otf.dv -c config_file

```

#### Notes
- `cudasirecon`should be supplied (unless to just display help message) with at least three arguments, which are interpreted differently depending on whether the input files are in TIFF or MRC/DV format:
  - For MRC/DV format, the three arguments represent`input_file_name`,`output_file_name`, and `otf_file_name`in that order
  - For TIFF format,  the three arguments represent `input_files_folder`,`pattern_shared_among_the_input_files`, and`otf_file_name`in that order; the output files are saved in a subfolder named `GPUsirecon` (subject to change) under `input_files_folder`
  - The order of these three arguments, to be called *positional arguments* hereafter, is essential because their position determines how the program interprets them (unless a specific option names preceeds them; more on that below).
  - The reason for such difference is that MRC/DV format usually stores a full time-series of 3D or 2D images in a single file, while in TIFF format typically each time point is saved in a separate TIFF file with a slightly different names from other files in the same series.
- Besides these three arguments, `cudasirecon` is typically supplied with options or flags (a detailed list is below). For example `--ndirs 3` informs the program the number of SIM pattern angles; or `--bessel` to inform that the data was acquired in lattice-light-sheet mode. All options that require an argument have a default setting if that option is not specified (for example one can skip the `--ndirs 3` option because that's the default). There is no rules on ordering of the options/flags and they can appear either before or after or even in-between the three positional arguments.
- The `-c` or `--config` option allows using a config file (as demo'ed in the above Shell scipt) to specify multiple options, with one`option_name=option_value`pair per each line and the choices for`option_name`are the same as the command-line option names (the long form preceded with`--`). An example 3D SIM config file may look like this:

```bash
nimm=1.515
background=90
wiener=0.001
# angles of illumination in radians
k0angles=-0.804300,-1.8555,0.238800
ls=0.2035
ndirs=3
nphases=5
na=1.42
otfRA=1
dampenOrder0=1
```

- OTF files are converted from PSF files using`makeotf`program also included in this package; more on its usage below.


### Full list of options/flags for`cudasirecon`

```txt
$ cudasirecon --help
Usage:
 cudasirecon [--options] [option-arguments] input_file(or folder) output_file(or file name pattern) otf_file

 List of options:
  --input-file arg            input file name (or data folder in TIFF mode)
  --output-file arg           output file name (or filename pattern in TIFF mode)
  --otf-file arg              OTF file name
  -c [ --config ] arg         name of a config file with parameters in place of command-line options
  --ndirs arg (=3)            number of SIM  directions
  --nphases arg (=5)          number of pattern phases per SIM direction
  --nordersout arg (=0)       number of output SIM orders (must be <= nphases//2; safe to ignore usually)
  --angle0 arg (=1.648)       angle of the first SIM angle in radians
  --ls arg (=0.172)           line spacing of SIM pattern in microns
  --na arg (=1.2)             detection objective's numerical aperture
  --nimm arg (=1.33)          refractive index of immersion medium
  --wiener arg (=0.01)        Wiener constant; lower value leads to higher resolution and noise;
                              playing with it extensively is strongly encouraged
  --zoomfact arg (=2)         lateral zoom factor in the output over the input images;
                              leaving it at 2 should be fine in most cases
  --zzoom arg (=1)            axial zoom factor; almost never needed
  --background arg (=0)       camera readout background
  --usecorr arg               use a flat-field correction file if provided
  --forcemodamp arg           modamps to be forced to these values; useful when
                              image quality is low and auto-estimated modamps 
                              are below, say, 0.1
  --k0angles arg              user these pattern vector k0 angles for all 
                              directions (instead of inferring the rest agnles 
                              from angle0)
  --otfRA                     using rotationally averaged OTF; otherwise using 
                              3/2D OTF for 3/2D raw data
  --otfPerAngle               using one OTF per SIM angle; otherwise one OTF is
                              used for all angles, which is how it's been done 
                              traditionally
  --fastSI                    SIM image is organized in Z->Angle->Phase order; 
                              otherwise assuming Angle->Z->Phase image order
  --k0searchAll [=arg(=0)]    search for k0 at all time points
  --norescale [=arg(=0)]      no bleach correction
  --equalizez                 bleach correction for z
  --equalizet                 bleach correction for time
  --dampenOrder0              dampen order-0 in final assembly; do not use for 2D SIM; good choice for high-background images
  --nosuppress [=arg(=0)]     do not suppress DC singularity in the result (good choice for 2D/TIRF data)
  --nokz0                     not using kz=0 plane of the 0th order in the final assembly (mostly for debug)
  --gammaApo arg (=1)         output apodization gamma; 1.0 means triangular 
                              apo; lower value means less dampening of 
                              high-resolution info at the tradeoff of higher 
                              noise
  --explodefact arg (=1)      artificially exploding the reciprocal-space 
                              distance between orders by this factor (for debug)
  --nofilterovlps [=arg(=0)]  not filtering the overlaping region between bands (for debug)
  --saveprefiltered arg       save separated bands (half Fourier space) into a file and exit (for debug)
  --savealignedraw arg        save drift-fixed raw data (half Fourier space) into a file and exit (for debug)
  --saveoverlaps arg          save overlap0 and overlap1 (real-space complex data) into a file and exit (for debug)
  --2lenses                   I5S data
  --bessel                    bessel-SIM data
  --besselExWave arg (=0.488) Bessel SIM excitation wavelength in microns
  --besselNA arg (=0.144)     Bessel SIM excitation NA)
  --deskew arg (=0)           Deskew angle; if not 0.0 then perform deskewing before processing
  --deskewshift arg (=0)      If deskewed, the output image's extra shift in X (positive->left)
  --noRecon                   No reconstruction will be performed; useful when combined with --deskew
  --cropXY arg (=0)           Crop the X-Y dimension to this number; 0 means no cropping
  --xyres arg (=0.1)          x-y pixel size (only used for TIFF files)
  --zres arg (=0.2)           z pixel size (only used for TIFF files)
  --zresPSF arg (=0.15)       z pixel size used in PSF TIFF files)
  --wavelength arg (=530)     emission wavelength in nanometers (only used for TIFF files)
  --writeTitle                Write command line to MRC/DV header (may cause issues with bioformats)
  -h [ --help ]               produce help message
  --version                   show version

```

#### Notes on the essential flags
- `--ndirs`and `--nphases` informs the program of the number of SIM angles and phases within each angle, respectively.
- `--ls`, `--angle0` (or `--k0angles`) inform the program of a good estimation on the SIM pattern's line spacing (in &#x00B5;m) and angles (in radians). If using `--angle0`, then the rest of the angles are successively incremented by &#x03C0;/`ndirs`; if using `--k0angles` then list all angles in order and separated by commas without any space.
- `--wiener`specifies the Wiener constant that balances between resolution and amplified noise artifact. The lower it is, the higher the resolution but also more amplified-noise-related artifacts and vice versa. 
- `--background`informs the estimated background level of the input images. The camera's dark image's mean intensity would be a good guess to start with. Too high of a background value supplied can be fatal to the reconstruction result.
- `--otfRA` is needed if rotationally averaged OTF is used
- `--fastSI`is needed if Z->Angle->Phase acquisition order is used
- For lattice-light-sheet SIM, one would need the `--bessel` flag and might want to specify the following flags: `--besselNA`, `--besselNA`, `--deskew`, `--deskewshift`

### For`makeotf`
The utility of `makeotf` is to generate rotationally averaged OTFs for all SIM orders. The input to it is a single-direction SIM data set acquire with one point-like object as the sample. In general this OTF can be used for reconstructing multi-angle SIM data.

```
# And example of how to call makeotf:

$ makeotf /path/to/psf.dv otf.dv -angle -1.855500 -ls 0.2075 -na 1.4 -nimm 1.515 -fixorigin 3 20 -leavekz 7 11 3

# The example above shows practically all flags that are needed for generating a 3D SIM OTF; if 2D SIM or LLSM-SIM you should skip the last two flags.
```
#### Essential flags

- `-fixorigin`takes two integers (>0) as arguments. Purpose of this flag is to suppress the singularity of order-0 (i.e., wide-field) 3D OTF at the origin by extrapolating the OTF's amplitudes between those pixels along the k<sub>r</sub> axis toward the origin and use that value in the end result. 3 and 20 are in general good for PSFs obtained with at least 256x256 pixels
 - **Do not** supply this option to process LLSM-SIM PSFs, as the issue with origin singularity is much less severe.
- `-leavekz`takes three integers as arguments. Purpose of this flag is to zero out the region outside of the OTF support. The first two numbers correspond to the two pixels on the positive k<sub>z</sub> axis of the order-1 OTF, between which the order-1 OTF support intersects the k<sub>z</sub> axis. The third number corresponds to a pixel on the positive k<sub>z</sub> axis, between which and the origin the order-2 OTF support intersects the k<sub>z</sub> axis. With these other related numbers such as NA and wavelengths, the program could calculate what's inside and what's outside the OTF supports and then zero out the outside parts.
 - To get these numbers right, it usually requires trial and error approach and examining the OTF output to see if the carving-out makes sense.
 - In TIFF mode, along with the OTF file the program also outputs a companion TIFF file, whose name contains `OnlyForViewing`, that can be easily opened to examine the OTF.
 - **Do not** supply this option to process LLSM-SIM PSFs, as the OTF supports are not well defined geometrically.


  
## GPU requirements

This software requires a CUDA-compatible NVIDIA GPU.

The libraries available on conda-forge have been compiled against different
versions of the CUDA toolkit.  The required CUDA libraries are bundled in the
conda distributions so you don't need to install the CUDA toolkit separately.
If desired, you can pick which version of CUDA you'd like based on your needs,
but please note that different versions of the CUDA toolkit have different GPU
driver requirements:

To specify a specific cudatoolkit version, install as follows (for instance, to use
`cudatoolkit=10.2`):

```sh
# NOTE: conda-forge coming soon... not available yet.
conda install -c conda-forge cudasirecon cudatoolkit=10.2
```

| CUDA  | Linux driver | Win driver |
| ----- | ------------ | ---------- |
| 10.2  | ≥ 440.33     | ≥ 441.22   |
| 11.0  | ≥ 450.36.06  | ≥ 451.22   |
| 11.1  | ≥ 455.23     | ≥ 456.38   |
| 11.2  | ≥ 460.27.03  | ≥ 460.82   |


If your CUDA Driver version is too low for the version of cudasirecon that you have installed, you may get an error that looks like: `!!Error occurred: cudaSetDevice failed`

If you run into trouble, feel free to [open an issue](https://github.com/scopetools/cudasirecon/issues) and describe your setup.


## Multichannel reconstruction

`cudasirecon` does not currently accept multi-channel files.  So it is necessary
to temporarily pull out each channel into a new file prior to reconstruction.
The provided `recon.py` script is an example of how to use the
[`mrc`](https://github.com/tlambert03/mrc) package to extract individual
channels from a .dv file, reconstruct them, and merge them back (and clean up
the intermediate file).  It is used as follows (note, `mrc`, `numpy`, and
`cudasirecon` must be in your path):

```
python recon.py /path/to/raw_data.dv
```

`recon.py` will also accept any key value pairs that `cudasirecon` also accepts
(to see that full list, type `cudasirecon -h` at the command prompt).  For
instance, to override just a couple of the reconstruction parameters, you could
do something like this:

```
python recon.py /path/to/raw_data.dv wiener 0.001 background 150
```

There are a couple of hard-coded filepaths in `recon.py`.  Specifically, it
currently expects to find the OTFs and config files in the same directory as the
recon.py script.  You can change that by putting in an absolute directory to
some other folder for the variables at the top of the file:

```python
# path to your otf directory.  Defaults to the same as the recon.py file
OTF_DIR = os.path.abspath(os.path.dirname(__file__))
# path to your config directory.  Defaults to the same as the recon.py file
CONFIG_DIR = os.path.abspath(os.path.dirname(__file__))
```

Note also, that the config and OTF files must be named with the emission
wavelength in the filenames. For example: `config528.txt` and `otf528.otf`


## Compiling locally from source

#### on linux:
```sh
conda env create -f environment-linux.yml
conda activate simbuild
./build.sh  # see build.sh for cmake details

# test it:
./test_data/run.sh
```

#### on windows:
```sh
conda env create -f environment-windows.yml
conda activate simbuild
bld.bat  # see bld.bat for cmake details
```
