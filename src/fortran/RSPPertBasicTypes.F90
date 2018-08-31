! basic data types
#include "api/qcmatrix_c_type.h"

module RSPPertBasicTypes_f
    use, intrinsic :: iso_c_binding
    use qcmatrix_f, only: QINT

    implicit none

    ! <datatype name='QcPertInt'>
    !   Data type of integers to represent perturbation labels
    ! </datatype>
    ! <datatype name='C_QCPERTINT'>
    !   Integers of perturbation labels to interoperates with C code
    ! </datatype>
    ! <constant name='QCPERTINT_MAX'>
    !   Maximal value of an object of the <QcPertInt> type
    ! </constant>
    !integer(kind=4), parameter, public :: QcPertInt = 8
    !integer, parameter, public :: C_QCPERTINT = C_LONG
    !integer(kind=4), parameter, public :: QcPertInt = 4
    integer, parameter, public :: QcPertInt = QINT
    !integer, parameter, public :: C_QCPERTINT = C_INT
    integer, parameter, public :: C_QCPERTINT = C_QINT
end module RSPPertBasicTypes_f

