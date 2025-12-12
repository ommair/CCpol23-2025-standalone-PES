module vdw_ccpol23plus

  use, intrinsic :: iso_fortran_env
  use common_arrays

  
  implicit none
  private

  public :: ccpol23plus

contains

  subroutine ccpol23plus(xx,yy,zz,xxs,yys,zzs,qvdwe)
    implicit none
    real(8) :: xx(nom),yy(nom),zz(nom)
    real(8) :: xxs(nsite,nom),yys(nsite,nom),zzs(nsite,nom)
    real(8) :: boxi,qqvdwe

    integer :: i,j,is,js,k
    real(8) :: wij,qis,qjs
    real(8) :: xcc,ycc,zcc
    real(8) :: dx,dy,dz
    real(8) :: rccsq
    real(8) :: rssq,rsx,rsy,rsz
    real(8) :: vfij
    real(8) :: sqr
    real(8) :: eks,u_exp,du_exp
    real(8) :: u_ele, f1,df1, du_ele
    real(8) :: u_ind,f6,f8,f10,df6,df8,df10
    real(8) :: r6,r8,r10,du_ind
    real(8) :: fdamp
    real(8) :: fxsij,fysij,fzsij
    ! Shared
    real(8), allocatable :: lfxs(:, :), lfys(:, :), lfzs(:, :)

    real(8) :: qvdwe, qvirial, qvir_vdw
    real(8) :: lqvdwe, lqvirial,lvdwlrc,lvirlrc
    real(8), parameter :: gamma = 20.d0
    real(8), parameter :: sqr0 = 1.4d0 

    qvdwe=0.d0
    lqvdwe=0.d0
   

    do i=1,nm-1

       do j=i+1, nm

          do is=1,ns

             do js=1,ns
                 
                rsx = xxs(is,i)-xxs(js,j)
                rsy = yys(is,i)-yys(js,j)
                rsz = zzs(is,i)-zzs(js,j) 
                  
                qis = charge(is)
                qjs = charge(js)

                  rssq = rsx**2 + rsy**2 + rsz**2
            
                  sqr = sqrt(rssq) ! square root of r^2

                  do k=1,npi

                     if( (an(is) .eq. an1(k) .and. an(js) .eq. an2(k)) .or. &
                         (an(js) .eq. an1(k) .and. an(is) .eq. an2(k)) )then         
 
                      ! ccpol8s Energy, Virial and Forces                    

                      ! exponential (excahnge - repulsion) 
                       eks = exp(-beta(k)*sqr)  !  exp(-beta r)
                       u_exp = eks*(c0(k)+c1(k)*sqr+c2(k)*sqr**2+c3(k)*sqr**3)

                       fdamp = 1.d0/( 1.d0+exp(-gamma*(sqr-sqr0))) ! short range damping function
                       u_exp = u_exp * fdamp   ! u_exp in atomic units       

                       ! short range coulumb contribution

                       f1 = tt(1, d1(k), sqr)
                       u_ele = f1*qis*qjs/sqr

    !                   write(*,*) f1,qis,qjs,sqr

                       ! induction-dispersion

                       f6 = tt(6, d6(k), sqr) 
                       f8 = tt(8, d8(k), sqr)
                       f10 = tt(10, d10(k), sqr) 
                       r6 = rssq*rssq*rssq
                       r8 = rssq*r6
                       r10 = rssq*r8     
 
                       u_ind = - f6*c6(k)/r6 - f8*c8(k)/r8 - f10*c10(k)/r10 ! u_ind in atomic units
                           
                       lqvdwe = lqvdwe + u_exp + u_ele + u_ind 

                       !write(*,*) u_exp, u_ele, u_ind

                     end if

                  end do

             end do
           
          end do

       end do

    enddo

    qvdwe = qvdwe + lqvdwe

    qvdwe = qvdwe*h2kcal

  end subroutine ccpol23plus

  function tt(n, delta, r)
    implicit none
    integer :: i,ncn
    integer, intent(in) :: n
    real(8), intent(in) :: delta, r
    real(8) :: tt
    real(8) :: term,ssum,deltar

    deltar = delta * r
    term = 1.0d0
    ssum = 1.0d0
    ncn = n
    do i=1, ncn
       term = term*deltar/i
       ssum = ssum + term
    end do
    tt = 1.0 - exp(-deltar)*ssum

  end function tt

end module vdw_ccpol23plus 
