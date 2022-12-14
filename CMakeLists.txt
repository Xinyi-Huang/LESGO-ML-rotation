cmake_minimum_required(VERSION 2.8.9)

################################################################################
#                               User Inputs
################################################################################

# Specify Compiler
# Options: Intel and GNU
set(CMAKE_Fortran_COMPILER_ID "Intel")

# Set the name of the computer you are compiling on.
# Supported computers: panoramix, njord, marcc, darwin-brew
# Otherwise, default settings are used
set(hostname "aci-psu")

#--Set global VERBOSE flag;
set(VERBOSE false)

# Set extra safety measures 
set(USE_SAFETYMODE true)

# Output additional quantities to file. May be useful
# for diagnosing issues
set(OUTPUT_EXTRA false)

# Enable parallelization support
set(USE_MPI true)

# Read and write binary files as big_endian or little_endian
# Options: DEFAULT, LITTLE, or BIG
set(WRITE_ENDIAN "DEFAULT")
set(READ_ENDIAN "DEFAULT")

# Select a turbine model
# USE_TURBINE: Actuator Disk Model without Rotation
# USE_ATM: Actuator Line/Sector Model 
set(USE_TURBINES false)
set(USE_ATM false)

# Select the data output
set(USE_CGNS false)

# Select level set
set(USE_LVLSET false)

# Concurrent precursor simulation 
# Requires configuration of two builds: upstream and downstream
set(USE_CPS false)

# Homogeneous isotropic turbulence input
set(USE_HIT false)

# Shift the domain in the y (spanwise) direction to eliminate
# streaks effects on averaging
set(USE_STREAKS false)

################################################################################
#                         Set up build (do not change)
################################################################################
# Set compiler versions for Cmake
#  Intel
if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    add_definitions(-DPPIFORT)
    set(compiler_name "ifort")
#  GNU
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    add_definitions(-DPPGFORTRAN)
    set(compiler_name "gfortran")
endif ()
find_program(WHICH NAMES which)
exec_program(${WHICH} ARGS ${compiler_name} OUTPUT_VARIABLE compiler_path)
set(CMAKE_Fortran_COMPILER ${compiler_path})

# Project details
enable_language (Fortran)
project (lesgo)
set(exec_name "lesgo")

# Library links
link_libraries(fftw3)    

# Common include and link directories
include_directories(/usr/local/include)
link_directories(/usr/local/lib)
if ( NOT(${hostname} STREQUAL "darwin-brew") )
    include_directories($ENV{HOME}/include /usr/include)
    link_directories($ENV{HOME}/lib /usr/lib64)
endif ()

# Include host-specific paths
#  panoramix
if (${hostname} STREQUAL "panoramix")
    include_directories(/opt/cgns-3.2.1-hdf5-mpi/include 
    /opt/hdf5-1.8.12-mpi/include /usr/local/hdf5/include /usr/local/cgns/include 
    /act/fftw3/mpich/gcc-4.7.2/include /act/hdf5-1.8.15/mpich/gcc-4.7.2/include 
    /act/cgns-3.3.0/mpich/gcc-4.7.2/include)
    link_directories(/opt/hdf5-1.8.12-mpi/lib  /opt/cgns-3.2.1-hdf5-mpi/lib 
    /usr/local/hdf5/lib  /usr/local/cgns/lib /act/fftw3/mpich/gcc-4.7.2/lib 
    /act/hdf5-1.8.15/mpich/gcc-4.7.2/lib /act/cgns-3.3.0/mpich/gcc-4.7.2/lib)
#  marcc
elseif (${hostname} STREQUAL "marcc")
    # Intel
    if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
        include_directories(/cm/shared/apps/Intel/fftw/3.3.4/include)
        link_directories(/cm/shared/apps/Intel/fftw/3.3.4/lib)
    # GNU -- These have to be installed locally for now.
    elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
        if (USE_CGNS)
            include_directories($ENV{HOME}/cgns/include/)
            link_directories($ENV{HOME}/hdf5/lib/ $ENV{HOME}/cgns/lib/)
        endif ()
    endif ()
#  njord
elseif (${hostname} STREQUAL "njord")
    # nothing yet, but can add if needed
#  all other systems
elseif (${hostname} STREQUAL "certainty")
        include_directories(/share/apps/lib/fftw/3.3.4/mvapich2-2.0rc1-intel-14/include)
        link_directories(/share/apps/lib/fftw/3.3.4/mvapich2-2.0rc1-intel-14/lib)
elseif (${hostname} STREQUAL "generic")
    # nothing yet, but can add if needed
    if (USE_CGNS)
        include_directories($ENV{HOME}/cgns/include/)
        link_directories($ENV{HOME}/hdf5/lib/ $ENV{HOME}/cgns/lib/)
    endif ()
endif ()

# List common sources
set (Sources cfl_util.f90 clocks.f90 convec.f90 derivatives.f90 divstress_uv.f90
     divstress_w.f90 dns_stress.f90 emul_complex.f90 fft.f90 finalize.f90 
     forcing.f90 fringe_util.f90 functions.f90 grid.f90 ic.f90 ic_dns.f90 
     initial.f90 initialize.f90 input_util.f90 interpolag_Sdep.f90 
     interpolag_Ssim.f90 io.f90 lagrange_Sdep.f90 lagrange_Ssim.f90 main.f90 
     messages.f90 open_file.f90 padd.f90 param.f90 param_output.f90 
     press_stag_array.f90 rmsdiv.f90 scaledep_dynamic.f90 
     sgs_param.f90 sgs_stag_util.f90 sim_param.f90 stat_defs.f90 
     std_dynamic.f90 string_util.f90 test_filtermodule.f90 tridag_array.f90 
     types.f90 unpadd.f90 wallstress.f90 wallstress_dns.f90 init_random_seed.f90)

if (USE_MPI)
    add_definitions(-DPPMPI)
    set (CMAKE_Fortran_COMPILER mpiifort)
    set(Sources ${Sources} mpi_transpose_mod.f90 mpi_defs.f90)
    set(exec_name "${exec_name}-mpi")
endif (USE_MPI)

if (USE_CPS)
    add_definitions(-DPPCPS)
    set(Sources ${Sources} concurrent_precursor.f90)
    set(exec_name "${exec_name}-cps")
endif (USE_CPS)

if (USE_HIT)
    add_definitions(-DPPHIT)
    set(Sources ${Sources} hit_inflow.f90)
    set(exec_name "${exec_name}-HIT")
endif (USE_HIT)

if (USE_STREAKS)
    add_definitions(-DPPSTREAKS)
    set(exec_name "${exec_name}-streaks")
endif (USE_STREAKS)

if (USE_LVLSET)
    add_definitions(-DPPLVLSET)
    set(Sources ${Sources} level_set_base.f90 level_set.f90 linear_simple.f90 
        trees_pre_ls.f90 trees_base_ls.f90 trees_setup_ls.f90 trees_io_ls.f90 
        trees_global_fmask_ls.f90)
    set(exec_name "${exec_name}-ls")
endif (USE_LVLSET)

if (USE_TURBINES)
    add_definitions(-DPPTURBINES)
    set(Sources ${Sources} turbines.f90 turbines_base.f90 util.f90)
    set(exec_name "${exec_name}-turbines")
endif (USE_TURBINES)

if (USE_ATM)
    add_definitions(-DPPATM)
    set(Sources ${Sources} atm_base.f90 atm_input_util.f90 
        actuator_turbine_model.f90 atm_lesgo_interface.f90)
    set(exec_name "${exec_name}-ATM")
endif (USE_ATM)

if (OUTPUT_EXTRA)
    add_definitions(-DPPOUTPUT_EXTRA)
    set(exec_name "${exec_name}-exout")
endif (OUTPUT_EXTRA)

if (USE_DYN_TN)
    add_definitions(-DPPDYN_TN)
    set(exec_name "${exec_name}-dyntn")
endif (USE_DYN_TN)

if (USE_SAFETYMODE)
    add_definitions(-DPPSAFETYMODE)
else (USE_SAFETYMODE)
    set(exec_name "${exec_name}-safety_off")
endif (USE_SAFETYMODE)

if (USE_CGNS)
    add_definitions(-DPPCGNS)
    link_libraries(cgns hdf5)   
    set(exec_name "${exec_name}-cgns") 
endif (USE_CGNS)

if (VERBOSE)
    add_definitions(-DPPVERBOSE)
endif (VERBOSE)

if (WRITE_ENDIAN STREQUAL "LITTLE")
    add_definitions(-DPPWRITE_LITTLE_ENDIAN)
endif (WRITE_ENDIAN STREQUAL "LITTLE")

if (WRITE_ENDIAN STREQUAL "BIG")
    add_definitions(-DPPWRITE_BIG_ENDIAN)
endif (WRITE_ENDIAN STREQUAL "BIG")

if (READ_ENDIAN STREQUAL "LITTLE")
    add_definitions(-DPPREAD_LITTLE_ENDIAN)
endif (READ_ENDIAN STREQUAL "LITTLE")

if (READ_ENDIAN STREQUAL "BIG")
    add_definitions(-DPPREAD_BIG_ENDIAN)
endif (READ_ENDIAN STREQUAL "BIG")

# Intel
if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    set (CMAKE_Fortran_FLAGS "-cpp -O2 -funroll-loops -no-wrap-margin -assume byterecl")
# GNU
else (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    set (CMAKE_Fortran_FLAGS "-cpp -O2 -ffast-math -funroll-loops -fall-intrinsics")
endif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")

# main lesgo executable
add_executable(${exec_name} ${Sources})

