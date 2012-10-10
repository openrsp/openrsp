if(NOT DEFINED DEFAULT_Fortran_FLAGS_SET)

if(CMAKE_Fortran_COMPILER_ID MATCHES GNU) # this is gfortran
    set(CMAKE_Fortran_FLAGS         "-DVAR_GFORTRAN -DGFORTRAN=445 -g -fbacktrace")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0")
    set(CMAKE_Fortran_FLAGS_RELEASE "-Wuninitialized -O3 -ffast-math -funroll-loops -ftree-vectorize")
    if(ENABLE_STATIC_LINKING)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -static"
            )
    endif()
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -fdefault-integer-8"
            )
    endif()
    if(ENABLE_BOUNDS_CHECK)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -fbounds-check"
            )
    endif()
    if(ENABLE_CODE_COVERAGE)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -ftest-coverage"
            )
    endif()

#   radovan: for magnus, this is how you can set compiler flags for one source file only
#   set_source_files_properties(openrsp/cgto_diff_eri.f90
#       PROPERTIES COMPILE_FLAGS
#       "-raboof"
#       )
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES G95)
    message(FATAL_ERROR "g95 is not supported")
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    add_definitions(-DVAR_IFORT)
    set(CMAKE_Fortran_FLAGS         "-g -w -fpp -assume byterecl -traceback")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -xW -ip")
    if(ENABLE_STATIC_LINKING)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -static-libgcc -static-intel"
            )
    endif()
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -i8"
            )
    endif()
    if(ENABLE_BOUNDS_CHECK)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -check bounds -fpstkchk -check pointers -check uninit -check output_conversion -ftrapuv"
            )
    endif()

    if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        message("--Switch off warnings due to incompatibility XCode 4 and Intel 11 on OsX 10.6")
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -Qoption,ld,-w"
            )
    endif()
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
    set(CMAKE_Fortran_FLAGS         "-DVAR_PGF90")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-g -O0 -Mframe")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -mcmodel=medium -fast -Munroll")
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -m64 -i8"
            )
    else()
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -m32"
            )
    endif()
    if(ENABLE_BOUNDS_CHECK)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} "
            )
    endif()
    if(ENABLE_CODE_COVERAGE)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} "
            )
    endif()
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set(CMAKE_Fortran_FLAGS         "-qzerosize -qextname")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-g")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -q64"
            )
    endif()

    set_source_files_properties(${FREE_FORTRAN_SOURCES}
        PROPERTIES COMPILE_FLAGS
        "-qfree"
        )
    set_source_files_properties(${FIXED_FORTRAN_SOURCES}
        PROPERTIES COMPILE_FLAGS
        "-qfixed"
        )
endif()

save_compiler_flags(Fortran)
endif()
