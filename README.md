# gromacs-tutorials

## Introduction

GROMACS is a versatile molecular dynamics (MD) package primarily used to perform MD simulations on biochemical molecules such as proteins, lipids and nucleic acids that have a lot of complicated bonded interactions. However, because it is extremely fast at calculating the non-bonded interactions that usually dominate simulations, it is increasingly used in research on non-biological systems such as polymers and fluid dynamics.

### GROMACS Installation
1. All commands are for Debian-based distros
2. Check C and C++ compilers versions (`gcc --version` and `g++ --version`)
3. If these are not up to date with the latest version, install latest versions via following commands
* `sudo apt update`
* `sudo apt install gcc g++`
3. Check CMake version via `cmake --version`. If it is not version 3.18.4 or later, update via:
* `sudo apt update`
* `sudo apt install cmake`
4. Download [latest GROMACS tarball](https://manual.gromacs.org/current/download.html)
5. Move GROMACS tarball to directory you want to install in
* `mv gromacs-2024.2.tar.gz <directory>`
6. Extract and decompress GROMACS tarball via:
* `tar xfz gromacs-2024.2.tar.gz`
7. Install GROMACS via following commands:
* `cd gromacs-2024.2 && mkdir build && cd build && cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON && make && make check && sudo make install && source /usr/local/gromacs/bin/GMXRC`
8. Install GROMACS tools command on your system via command:
* `sudo apt install gromacs`
