This is a program for setting up exponential disk galaxies in Modified Newtonian Dynamics using the algebraic MOND approximation (implemented in the force functions at the end of the dice_vel.c file in the src folder). The mass of each particle is varied to maintain a similar mass in the central region and a similar radial resolution in the outer region. This maximises the chance that the disk is stable. The outer radial limit can be altered but should exceed about 10 scale lengths. It is also easy to set the code up to use equal mass particles, in this case search dice_structure.c in the src folder for the two comments with "equal mass" and follow the instructions. Within each component, the relative mass of each particle is specified by its gal->w property, but particles in different components with the same gal->w could have different masses. The normalisation is done using the fix_masses routine developed by Roy Truelove and Indranil Banik. The cylindrical radius of each particle is calculated from a random number. The Markov Chain Monte Carlo procedure in the original algorithm has been disabled. Instead, you must be able to tell the computer how to go from some random numbers to a position such that the particles follow the correct density distribution, when weighted according to the weights you specify. Example: cylindrical radius is uniform out to 20 scale lengths, the weight is u*exp(-u) where u = r_cyl/scale_length, thus yielding an exponential disk. This is just an example and would lead to an unstable disk, but other decompositions are of course possible between the particle distribution and relative weights. Their product determines the actual density distribution.

To change file name:
Search for Galaxy_name in dice.h, need to provide names if galaxy has only one component or if there are multiple components.
If longer file names are needed: go to dice_io.c and dice_init.c and allocate more space to outfilename.

To change parameters:
Change examples:test_m33.config to point to appropriate params file
Change examples:params_files:testM33.params
No new variables have been added for each component of each galaxy
The most important parameters are m200 (total mass) for the whole galaxy and, for each component, mass_frac, scale_length and flatz. The last is the aspect ratio = scale_height/scale_length.
Note the algorithm only supports exponential radial profiles with sech_sq vertical profiles.
The outer radial limit of each component is set in the headers file dice.h, search for r_outer_limit_comp.
To have equal mass particles, go to dice_structure.c, search for equal mass particles and follow the instructions. Two changes are required. In this case, r_outer_limit_comp is not used and the particle positions are just allocated according to the exponential surface density law. Please ensure the size of the grid is large enough to include all the particles (the parameters called boxsize).
By changing these blocks of code, it is possible to set up other density profiles.
The code also has sections for spherical components which can be adapted to your chosen spherical density profile. For these cases, please set flatz to 1.0.

It is recommended to only change the component scale lengths, outer limits (except for equal mass particles), mass fractions, aspect ratios and particle numbers as well as m200 for the total galaxy mass.


Running instructions:
The results are stored in build:bin.
Navigate to build:bin directory to start.
Repeatable mode (can keep repeating this, including changes to the code at the end of each run):
cd ..
make
make install
cd bin
./dice ../../example/test_m33.config




Installation instructions:
Disclaimer:
This Wiki page intends to guide you in the process of compiling and installing the DICE software.
Before trying to install DICE, please check that the following software/libraries are installed on your system: CMake, GSL, FFTW3

Installing dependencies:
It is preferable to install your local version of GSL and FFTW, especially if your are running on a cluster.
First, make sure you have a local folder in your $HOME:
cd ~
mkdir local

Add to your .bashrc(linux systems) or .profile (darwin systems):
export PATH=$HOME/local/bin:$PATH

To perform a local installation of the libraries, execute the following command lines:
wget http://mirror.switch.ch/ftp/mirror/gnu/gsl/gsl-2.4.tar.gz
tar -zxvf gsl-2.4.tar.gz
cd gsl-2.4
./configure --prefix=$HOME/local --enable-shared
make
make install
cd ..
wget http://www.fftw.org/fftw-3.3.5.tar.gz
tar -zxvf fftw-3.3.5.tar.gz
cd fftw-3.3.5
./configure --prefix=$HOME/local --enable-threads --enable-shared
make
make install
cd ..

Compile & Install
You can checkout the latest version in the Git repository ny typing:

git clone https://github.com/GFThomas/MOND/init_conditions/disc

The DICE package comes with the CMake cross-platform build system. So technically, you don’t have to worry so much about the compilation. Make sure you have cmake installed by typing:

cmake --version
If your version is older than cmake 2.6, you will have to update your system to a more recent version of cmake. To generate the makefile, type:

cd disc
mkdir build
cd build
cmake ..
If no errors are detected through this step, compile the code and install the code like this:
make
make install

By default, DICE is installed in $HOME/local/disc, but you can specify a different installation directory using the flag -DCMAKE_INSTALL_PREFIX=/install/path.
A cmake macro is implemented to locate standard installations of the required libraries. Nevertheless, if you installed them in a different way, use the following keywords to help cmake locating these libraries:
cmake .. -DGSL_PATH=/path/to/gsl -DFFTW3_PATH=/path/to/fftw3 -DFFTW3_THREADS_PATH=/path/to/fftw3_threads

Furthermore, it is possible to specify the compiler flags to cmake. For example, debug options for the icc compiler:
cmake .. -DCMAKE_C_FLAGS="-g -O2 -traceback -check=uninit -debug all -ftz -ftrapuv"

Once installed, think about adding to the PATH environment variable the location of the DICE binary. If you use bash and under standard installation process, add to your .bashrc (or .profile for OSX systems) :
export PATH=$HOME/local/bin:$PATH

Optionally, if you want the code to run faster, you can use the threading abilities of the FFTW library for the gravitational potential computation and some openMP optimisations. In this case, and assuming that you have installed the treaded version of FFTW3, then replace the previous commands by:

cmake .. -DENABLE_THREADS=ON
make
make install
CMake will automatically look for the fftw3_threads library, and link it. Compilation flags can be specified with the keyword -DCMAKE_C_FLAGS

If you have some errors referring to unknown reference _intel_fast_memset_, it probably means that FFTW was compiled with intel icc compiler with specific optimisation flags. Thus, you will need icc to compile DICE.

If you want to get the latest version in the Git repository, just type:
cd /path/to/dice (which should be home/disc)
git pull




Reference:
http://adsabs.harvard.edu/abs/2016ascl.soft07002P

The QUMOND analogue of the Toomre disk stability condition is obtained from https://arxiv.org/abs/1808.10545.
The algorithm is based on applying QUMOND via the algebraic MOND approximation. The critical MOND adjustments are in the dice_vel.c file.

For queries, please contact Indranil Banik at indranilbanik1992@gmail.com.
