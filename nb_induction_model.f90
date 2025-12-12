module  nb_induction_model

  use, intrinsic :: iso_fortran_env
  use common_arrays

  
  implicit none
  private

  public :: nb_induction

  public :: rspace_efield_stat

  public :: rspace_efield_idm

  public :: iteridm

contains

  subroutine nb_induction(xx,yy,zz,xxs,yys,zzs,qindpe)

    implicit none

    integer :: is,js,i,j,k,isteps
    real(8) :: xx(nom),yy(nom),zz(nom)
    real(8) :: xxs(nsite,nom),yys(nsite,nom),zzs(nsite,nom)
    real(8) :: boxi,idmlength
    real(8) :: qindpe    

    alphad = alpha

    call rspace_efield_stat(xx,yy,zz,xxs,yys,zzs)

    call iteridm(xx,yy,zz,xxs,yys,zzs,qindpe)  ! it includes all contribution from real, recip, self 

  end subroutine nb_induction

  subroutine rspace_efield_stat(xx,yy,zz,xxs,yys,zzs)
    implicit none
    real(8) :: xx(nom),yy(nom),zz(nom)
    real(8) :: xxs(nsite,nom),yys(nsite,nom),zzs(nsite,nom)
    real(8) :: alp,boxi

    integer :: is,js,i,j,k,isteps
    real(8) :: dist2,rccsq,dist,doti,dotj
    real(8) :: qis,qjs  ! charges on sites
    real(8) :: xcc,ycc,zcc,dx,dy,dz,ddx,ddy,ddz,d1i,d2i,d3i,d5i,d7i
    real(8) :: expar2,f2,f3

    do i=1,nm
       do is=1,ns
          rE0x(is,i) = 0.d0 ! initiating permanent charges fields
          rE0y(is,i) = 0.d0
          rE0z(is,i) = 0.d0
       end do
    end do


      do i=1,nm 

         do j=1,nm    ! loop over molecule j   starts

          if (j .ne. i) then

            do is=1,ns

               do js=1,ns

                  qis = charge(is)
                  qjs = charge(js)

                  ddx =  (xxs(is,i) - xxs(js,j))  ! separation between polarization
                  ddy =  (yys(is,i) - yys(js,j))  ! center of A and charged sites of B
                  ddz =  (zzs(is,i) - zzs(js,j))

                  dist2 = ddx**2 + ddy**2 + ddz**2

                  dist = sqrt(dist2)
                  d1i = 1.d0/dist
                  d2i = d1i**2
                  d3i = d2i*d1i
    
                  f2 = tt(2,delta2,dist)
                 
                  rE0x(is,i) = rE0x(is,i) + f2*qjs*ddx*d3i
                  rE0y(is,i) = rE0y(is,i) + f2*qjs*ddy*d3i
                  rE0z(is,i) = rE0z(is,i) + f2*qjs*ddz*d3i
 

               end do

            end do

          endif 
 
        end do        ! loop over molecule j ends

    end do           ! loop over molecule i ends

  end subroutine rspace_efield_stat


  subroutine iteridm(xx,yy,zz,xxs,yys,zzs,qindpe)
    implicit none
    real(8) :: xx(nom),yy(nom),zz(nom)
    real(8) :: xxs(nsite,nom),yys(nsite,nom),zzs(nsite,nom)
    real(8) :: boxi,qindpe

    integer :: is,js,i,j,k,isteps
!    real(8) :: Etx(nsite,nom),Ety(nsite,nom),Etz(nsite,nom)  ! total efield 
!    real(8) :: shtEtx(nsite,nom),shtEty(nsite,nom),shtEtz(nsite,nom)  ! short range total efield
    real(8), allocatable :: Etx(:, :), Ety(:, :), Etz(:, :)  ! total efield 
    real(8), allocatable :: shtEtx(:, :), shtEty(:, :), shtEtz(:, :)  ! short range total efield
    real(8) :: efldxi,efldyi,efldzi,efldxj,efldyj,efldzj
    real(8) :: xpolc,ypolc,zpolc !separation between dipole centers
    real(8) :: dist2,rccsq,dist,doti,dotj
    real(8) :: qis,qjs  ! charges on sites
    real(8) :: xcc,ycc,zcc,dx,dy,dz,ddx,ddy,ddz,d1i,d2i,d3i,d5i,d7i
    real(8) :: polr,d3ii,d5ii,expar2,selffac
    real(8) :: thr_iter,change
    real(8) :: energy,engtst
    real(8), parameter :: maxit = 500
    real(8) :: polE1x,polE1y,polE1z
    real(8) lchange,lenergy

    allocate(Etx(ns,nm))
    allocate(Ety(ns,nm))
    allocate(Etz(ns,nm))
    allocate(shtEtx(ns,nm))
    allocate(shtEty(ns,nm))
    allocate(shtEtz(ns,nm))

    do i=1,nm
       do is=1,ns

       idmx(is,i) = 0.d0 ! initiating induced dipoles
       idmy(is,i) = 0.d0
       idmz(is,i) = 0.d0

       E0x(is,i) = 0.d0 ! initiating permanent charges fields
       E0y(is,i) = 0.d0
       E0z(is,i) = 0.d0

       Eidmx(is,i) = 0.d0 ! initiating permanent charges fields
       Eidmy(is,i) = 0.d0
       Eidmz(is,i) = 0.d0  

       Etx(is,i) = 0.d0 ! initiating induced dipoles
       Ety(is,i) = 0.d0
       Etz(is,i) = 0.d0

       end do
    end do

!**************************************************************************
! Now iterate to calculate converged polarization energy        starts here

    thr_iter = 1.d-20
    change = 10.d0
    isteps = 0

    do while (change.gt.thr_iter.and.isteps.lt.maxit)  ! while loop starts

    energy = 0.0d0
    change = 0.d0 


    call rspace_efield_idm(xx,yy,zz,xxs,yys,zzs)              

    lchange = 0.0D0
    lenergy = 0.0D0

      do i=1,nm    ! loop over ith molecule starts

         do is=1,ns 

            E0x(is,i) = rE0x(is,i)
            E0y(is,i) = rE0y(is,i)
            E0z(is,i) = rE0z(is,i)  

            Eidmx(is,i) = rEidmx(is,i)
            Eidmy(is,i) = rEidmy(is,i)
            Eidmz(is,i) = rEidmz(is,i)

            Etx(is,i) =  Eidmx(is,i) + E0x(is,i) ! total efield from Ewald sum 
            Ety(is,i) =  Eidmy(is,i) + E0y(is,i) 
            Etz(is,i) =  Eidmz(is,i) + E0z(is,i) 

            polr = polarizability(is) !apol(is,i)

            polE1x = polr* Etx(is,i) 
            polE1y = polr* Ety(is,i) 
            polE1z = polr* Etz(is,i)   
 
            lchange = (idmx(is,i)-polE1x)**2 + &
                      (idmy(is,i)-polE1y)**2 + &
                      (idmz(is,i)-polE1z)**2 + lchange

            idmx(is,i) = polE1x
            idmy(is,i) = polE1y
            idmz(is,i) = polE1z

            lenergy = -0.5d0*polr*(Etx(is,i)*E0x(is,i) + &
                                  Ety(is,i)*E0y(is,i) + &
                                  Etz(is,i)*E0z(is,i)) +  lenergy

!            energy = -0.5d0*(idmx(is,i)*E0x(is,i) + &
!                             idmy(is,i)*E0y(is,i) + &
!                             idmz(is,i)*E0z(is,i)) +  energy

         end do   ! loop over is 

      end do  ! loop over ith molecule ends  

      change = change + lchange
      energy = energy + lenergy

      isteps = isteps + 1


    end do       ! while loop ends

    if (isteps.ge.maxit) then
        write (*,*) 'No convergence in indN_iter'
        write (*,'(a,g12.3)') 'energy change=',change
        write (*,'(a,g12.3)') 'thr_iter=',thr_iter
        stop
    end if

    qindpe = energy*h2kcal

!**************************************************************************
! calculation of converged polarization energy        ends here

    deallocate(Etx)
    deallocate(Ety)
    deallocate(Etz)
    deallocate(shtEtx)
    deallocate(shtEty)
    deallocate(shtEtz)

  end subroutine iteridm

  subroutine rspace_efield_idm(xx,yy,zz,xxs,yys,zzs)
    implicit none
    real(8) :: xx(nom),yy(nom),zz(nom)
    real(8) :: xxs(nsite,nom),yys(nsite,nom),zzs(nsite,nom)
    real(8) :: alp,boxi,f2,f3

    integer :: is,js,i,j,k,isteps
    real(8) :: efldxi,efldyi,efldzi,efldxj,efldyj,efldzj
    real(8) :: dist2,rccsq,dist,doti,dotj,d3ii,d5ii,drr,expar2
    real(8) :: qis,qjs  ! charges on sites
    real(8) :: xcc,ycc,zcc,dx,dy,dz,ddx,ddy,ddz,d1i,d2i,d3i,d5i,d7i

    do i=1,nm
       do is=1,ns

       rEidmx(is,i) = 0.d0 ! initiating permanent charges fields
       rEidmy(is,i) = 0.d0
       rEidmz(is,i) = 0.d0

       end do
    end do

    do i=1,nm   ! loop over ith molecule starts

       do j=1,nm   ! loop over jth molecule starts

          if (j .NE. i) then

             do is=1,ns

                do js=1,ns    ! loop over sites of molecule j  Starts

                   ddx =  xxs(is,i) - xxs(js,j)
                   ddy =  yys(is,i) - yys(js,j)
                   ddz =  zzs(is,i) - zzs(js,j)

                   dist2 = ddx**2 + ddy**2 + ddz**2

                   dist = sqrt(dist2)
                   d1i = 1.d0/dist
                   d2i = d1i**2
                   d3i = d2i*d1i
                   d5i = d3i*d2i

                   f3 =tt(3,delta3,dist)

                      ! dot product of induced dipole and separation vector
                   doti = idmx(is,i)*ddx+idmy(is,i)*ddy+idmz(is,i)*ddz
                   dotj = idmx(js,j)*ddx+idmy(js,j)*ddy+idmz(js,j)*ddz 

                   ! calculate the  efield of inducd dipole of molecule j at
                   ! induced dipole of molecule i

                   efldxi = f3*(3.d0*d5i*dotj*ddx - idmx(js,j)*d3i)
                   efldyi = f3*(3.d0*d5i*dotj*ddy - idmy(js,j)*d3i)
                   efldzi = f3*(3.d0*d5i*dotj*ddz - idmz(js,j)*d3i)

                      ! calculate total efield at induced dipole at molecule i
                      ! because of parmenant charges
                      !  and induced dipoles at molecule j

                      rEidmx(is,i) = rEidmx(is,i) + efldxi
                      rEidmy(is,i) = rEidmy(is,i) + efldyi
                      rEidmz(is,i) = rEidmz(is,i) + efldzi

!                       end if

                end do     ! loop over sites of molecule j  ends

             end do

          endif

       end do      ! loop over jth molecule ends

    end do  ! loop over ith molecule ends

  end subroutine rspace_efield_idm


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

end module nb_induction_model
