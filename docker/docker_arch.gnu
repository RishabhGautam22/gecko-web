#!/bin/bash
# Architecture file for Docker (Debian/Ubuntu)

# NetCDF paths for Debian/Ubuntu (installed via apt-get)
export my_netcdf_fortran_lib='/usr/lib/x86_64-linux-gnu'
export my_netcdf_fortran_inc='/usr/include'
export my_netcdf_fortran_bin='/usr/bin'

export my_netcdf_c_lib='/usr/lib/x86_64-linux-gnu'
export my_netcdf_c_inc='/usr/include'
export my_netcdf_c_bin='/usr/bin'

# HDF5 (usually not explicitly needed if linked via netcdf, but good to have)
export my_hdf5_lib='/usr/lib/x86_64-linux-gnu'

export LD_LIBRARY_PATH=${my_netcdf_c_lib}:${my_netcdf_fortran_lib}:$LD_LIBRARY_PATH

# Use the GNU makefile header
export my_hdr=Makefile.hdr.gnu
export FC=gfortran
