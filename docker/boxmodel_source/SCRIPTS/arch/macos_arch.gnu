# Architecture file for macOS (Auto-generated)
export my_netcdf_fortran_lib='unknown option: --libdir
Usage: nf-config [OPTION]

Available values for OPTION include:

  --help        display this help message and exit
  --all         display all options
  --cc          C compiler
  --fc          Fortran compiler
  --cflags      pre-processor and compiler flags
  --fflags      flags needed to compile a Fortran program
  --has-dap     whether OPeNDAP is enabled in this build
  --has-nc2     whether NetCDF-2 API is enabled
  --has-nc4     whether NetCDF-4/HDF-5 is enabled in this build
  --has-f90     whether Fortran 90 API is enabled in this build
  --has-f03     whether Fortran 2003 API is enabled in this build
  --flibs       libraries needed to link a Fortran program
  --prefix      Install prefix
  --includedir  Include directory
  --version     Library version'
export my_netcdf_fortran_inc='/opt/homebrew/Cellar/netcdf-fortran/4.6.2/include'
export my_netcdf_fortran_bin='/opt/homebrew/Cellar/netcdf-fortran/4.6.2/bin'

export my_netcdf_c_lib='/opt/homebrew/Cellar/netcdf/4.9.3_1/lib'
export my_netcdf_c_inc='/opt/homebrew/Cellar/netcdf/4.9.3_1/include'
export my_netcdf_c_bin='/opt/homebrew/Cellar/netcdf/4.9.3_1/bin'

export LD_LIBRARY_PATH=${my_netcdf_c_lib}:${my_netcdf_fortran_lib}:$LD_LIBRARY_PATH
export my_hdr=Makefile.hdr.gnu
export FC=gfortran
