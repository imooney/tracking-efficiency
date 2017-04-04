# PYTHIA configuration file.
# Generated on Wed Sep 7 14:15:32 EDT 2016 with the user supplied options:
# --enable-shared

# Install directory prefixes.
PREFIX_BIN=/Users/imoo/pythia8219/bin
PREFIX_INCLUDE=/Users/imoo/pythia8219/include
PREFIX_LIB=/Users/imoo/pythia8219/lib
PREFIX_SHARE=/Users/imoo/pythia8219/share/Pythia8

# Compilation flags (see ./configure --help for further documentation).
ENABLE_SHARED=true
CXX=g++
CXX_COMMON=-O2 -ansi -pedantic -W -Wall -Wshadow -fPIC
CXX_SHARED=-dynamiclib -flat_namespace -undefined suppress
CXX_SONAME=-Wl,-dylib_install_name
LIB_SUFFIX=.dylib

# EVTGEN configuration.
EVTGEN_USE=false
EVTGEN_BIN=
EVTGEN_INCLUDE=./
EVTGEN_LIB=./

# FASTJET3 configuration.
FASTJET3_USE= true
FASTJET3_BIN= /Users/imoo/fastjet-install/bin/
FASTJET3_INCLUDE= /Users/imoo/fastjet-install/include/
FASTJET3_LIB= /Users/imoo/fastjet-install/lib/

# HEPMC2 configuration.
HEPMC2_USE=false
HEPMC2_BIN=
HEPMC2_INCLUDE=./
HEPMC2_LIB=./

# HEPMC3 configuration.
HEPMC3_USE=false
HEPMC3_BIN=
HEPMC3_INCLUDE=./
HEPMC3_LIB=./

# LHAPDF5 configuration.
LHAPDF5_USE=false
LHAPDF5_BIN=
LHAPDF5_INCLUDE=./
LHAPDF5_LIB=./
LHAPDF5_PLUGIN=LHAPDF5.h

# LHAPDF6 configuration.
LHAPDF6_USE=false
LHAPDF6_BIN=
LHAPDF6_INCLUDE=./
LHAPDF6_LIB=./
LHAPDF6_PLUGIN=LHAPDF5.h

# POWHEG configuration.
POWHEG_USE=false
POWHEG_BIN=
POWHEG_INCLUDE=./
POWHEG_LIB=./

# PROMC configuration.
PROMC_USE=false
PROMC_BIN=
PROMC_INCLUDE=./
PROMC_LIB=./

# ROOT configuration.
ROOT_USE= true
ROOT_BIN= /Users/imoo/root/bin/
ROOT_INCLUDE= /Users/imoo/root/include/
ROOT_LIB= /Users/imoo/root/lib/

# GZIP configuration.
GZIP_USE=false
GZIP_BIN=
GZIP_INCLUDE=./
GZIP_LIB=./

# BOOST configuration.
BOOST_USE=false
BOOST_BIN=
BOOST_INCLUDE=./
BOOST_LIB=./

# PYTHON configuration.
PYTHON_USE=false
PYTHON_BIN=
PYTHON_INCLUDE=./
PYTHON_LIB=./
