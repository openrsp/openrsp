# Options for making the library
set(OPENRSP_PERT_LABEL_BIT "10"
    CACHE STRING "Number of bits for a perturbation label.")
add_definitions(-DOPENRSP_PERT_LABEL_BIT=${OPENRSP_PERT_LABEL_BIT})
option(OPENRSP_ZERO_BASED "Zero-based numbering." ON)
if(OPENRSP_ZERO_BASED)
    add_definitions(-DOPENRSP_ZERO_BASED)
endif()
option(OPENRSP_BUILD_WEB "Build OpenRSP from WEB files (only useful for developers)." OFF)
# OpenRSP built on top of the QcMatrix library
set(QCMATRIX_HEADER_DIR None
    CACHE STRING "Directory of header files of QcMatrix library.")
include_directories(${QCMATRIX_HEADER_DIR})
set(QCMATRIX_LIB None
    CACHE STRING "Name of QcMatrix library with absolute path.")
