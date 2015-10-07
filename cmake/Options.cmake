# Options for making the library
option(OPENRSP_USER_CONTEXT "Enable user context in callback functions." OFF)
if(OPENRSP_USER_CONTEXT)
    add_definitions(-DOPENRSP_C_USER_CONTEXT)
endif()
set(OPENRSP_PERT_LABEL_BIT "10"
    CACHE STRING "Number of bits for a perturbation label.")
add_definitions(-DOPENRSP_PERT_LABEL_BIT=${OPENRSP_PERT_LABEL_BIT})
option(OPENRSP_BUILD_WEB "Build OpenRSP from WEB files." OFF)
option(OPENRSP_TEST_EXECUTABLE "Build test suite as excutables." ON)
option(OPENRSP_FORTRAN_API "Build Fortran 2003 APIs." OFF)
# OpenRSP built on top of the QcMatrix library
set(QCMATRIX_HEADER_DIR None
    CACHE STRING "Directory of header files of QcMatrix library.")
include_directories(${QCMATRIX_HEADER_DIR})
if(OPENRSP_FORTRAN_API)
    set(QCMATRIX_MODULE_DIR None
        CACHE STRING "Directory of Fortran modules of QcMatrix library.")
    include_directories(${QCMATRIX_MODULE_DIR})
    # OpenRSP Fortran APIs need to enable user context of OpenRSP C parts
    add_definitions(-DOPENRSP_C_USER_CONTEXT)
    if(OPENRSP_USER_CONTEXT)
        add_definitions(-DOPENRSP_F_USER_CONTEXT)
    endif()
endif()
set(QCMATRIX_LIB None
    CACHE STRING "Name of QcMatrix library with absolute path.")
