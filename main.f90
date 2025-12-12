program main

! use, intrinsic :: omp_lib

 use common_arrays
! use utility
 use initialization
! use molecular_sites
 use vdw_ccpol23plus 
! use statistics 
 

  implicit none

  integer :: i,j,is,js,k
  real(4) iniTime, deltaTime, tStart, tFinish
  real(4) iniTimeTotal, deltaTimeTotal, tStartTotal, tFinishTotal
  real(8) :: sumvir_test
  real(8) :: imx,imy,imz,iqx,iqy,iqz


  open(100, file='OUTPUT',STATUS='UNKNOWN', POSITION='APPEND')

  call start

  write(100,'(A250)')'------------------------------------------------------------------------------------------&
         --------------------------------------------------------------------------------------------------------------------&                   
         ---------------------------------------------------------------------------------------------------------------------'
  write(100,'(A10, A13, A17, A16, A17, A16, A15, A15,A19, A15, A15, A15, A15,A18,A18,A15)') &                                                       
                     "STEP","VDW(kcal/mol)","IND(kcal/mol)","3B(kcal/mol)","TOTPE(kcal/mol)"
  write(100,'(A250)')'------------------------------------------------------------------------------------------&
         --------------------------------------------------------------------------------------------------------------------&
         ---------------------------------------------------------------------------------------------------------------------'

  stpvdwpe = vdwpe
  stpindpe = indpe
  stp3bpe = fs2pe+fs3pe
  stppe = stpvdwpe+stpcpe+stpindpe+stp3bpe

  write(100,'(1x,i8,1p,15e16.8)') step,stpvdwpe,stpindpe,stp3bpe,stppe


close(100)

end program main
