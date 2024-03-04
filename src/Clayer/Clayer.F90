!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module provides the interface to the rest of the C code in the 
! C layer. This layer is used for easier interface between fortran and 
! C (e.g. avoiding function pointers and other non-trivial things to 
! link to in Fortran). Interoperability is enabled by using the 
! iso_c_binding intrinsic module. 

! Notes
!======
! Note 1: this module serves as a wrapper module to include all 
! possible different C modules immediately. Currently, this is limited 
! to only CSparseF.F90. 

module Clayer 

    use CSparseF 

end module 
