module initialization

  use, intrinsic :: iso_fortran_env
  use common_arrays
  use vdw_ccpol23plus
  use nb_induction_model
  use threebody_potential 

  implicit none  

  private
  public :: start, control_parameters, molecule

contains

  subroutine start
       
    implicit none
  
    integer :: i,j,is,js,k,m,ics,is3,nlin0,ii3,ii,iis

    call molecule               ! reads molecular site's masses,charges,polarizabilites

    call read_geo              !  reads only atomic corrds 
    call generate_offatomic_sites   ! generates off-atomics sites using ratios

    do i=1,nm
       do is=1,ns

          x(i) = x(i) + massites(is)*xs(is,i)/totm
          y(i) = y(i) + massites(is)*ys(is,i)/totm
          z(i) = z(i) + massites(is)*zs(is,i)/totm

       end do
    end do


    call ccpol23plus(x,y,z,xs,ys,zs,vdwpe)

    call nb_induction(x,y,z,xs,ys,zs,indpe) 

    call threebody_term(x,y,z,xs,ys,zs,fs2pe,fs3pe) 
 

  write(*,'(2x,a,2g20.12)')'vdw         (kcal/mol)' , vdwpe                                                                                                          
  write(*,'(2x,a,2g20.12)')'NB(ind)     (kcal/mol)' , indpe                                                                                                          
  write(*,'(2x,a,2g20.12)')'3B(Fs2+Fs3) (kcal/mol)' , fs2pe+fs3pe                                                                                                    
  write(*,'(2x,a,2g20.12)')'CCpol23+    (kcal/mol)' , vdwpe+indpe+fs2pe+fs3pe   

  end subroutine start

  subroutine control_parameters

    implicit none
  
    open(unit=10, file='CONTROL', form='formatted')

 
    read(10,*) pot_name   ! name of potential required 
 

    close(10)

  end subroutine control_parameters

  subroutine molecule

    implicit none

    integer :: i,j,is,js,k,ics 
    integer :: is3,nlin0,ii3

    open(unit=20, file='FIELD',form='formatted')
    read(20,*) system_name
    read(20,*) ns 
    read(20,*)

    do is=1,ns
       read(20,*)an(is),massites(is),charge(is),polarizability(is)
    end do

    read(20,*)
    read(20,*)
    do is=1,ns
       read(20,*)an(is),coff(is,1),coff(is,2),coff(is,3)
    end do

    read(20,*)
    read(20,*) npi
    read(20,*)

       do k=1,npi
          !read(20,*) an1(k),an2(k),expalp(k),beta(k),a1(k),a2(k),a3(k)
          read(20,*) an1(k),an2(k),beta(k),c0(k),c1(k),c2(k),c3(k), &
                c6(k),c8(k),c10(k),d1(k),d6(k),d8(k),d10(k),qa(k),qb(k)
       end do  

    close(20)

! reading 3body parameters starts 

!*************  short range ************
    open(unit=40, file='data3b_full_SR_5p5', form='formatted')
!    open(unit=40, file='data3b_full', form='formatted')
    read(40,*)

    do is=1,17
    do js=1,17
       read(40,*) bS3_s(is,js),gS3_s(is,js),r0S3_s(is,js),buf1(is),buf2(js)
       !write(*,*) bS3(is,js),gS3(is,js),r0S3(is,js)
    end do
    end do


!   reading non-linear parameters fs3
    read(40,*)
    read(40,*) nlin0
    do is3=1,nlin0
       read(40,*) ii3, cs3_s(is3)
    end do

    close(40)

!*************  Long range ************
    open(unit=50, file='data3b_full_LR_4p5', form='formatted')
!!    open(unit=50, file='data3b_full_LR_4p5_filter0', form='formatted')
    read(50,*)

    do is=1,17
    do js=1,17
       read(50,*) bS3_l(is,js),gS3_l(is,js),r0S3_l(is,js),buf1(is),buf2(js)
       !write(*,*) bS3(is,js),gS3(is,js),r0S3(is,js)
    end do
    end do


!   reading non-linear parameters fs3
    read(50,*)
    read(50,*) nlin0
    do is3=1,nlin0
       read(50,*) ii3, cs3_l(is3)
    end do

    close(50)


! reading 3body parameters ends  

    totm = 0.d0     ! initaiting total mass of the molecule
    natom = 0      ! initiating number of massive centers

    do i=1,ns
       if (massites(i).ne.0.d0) then
         natom = natom + 1            ! total number of massive sites
         totm  = totm + massites(i)    ! total mass of a molecule
       end if
    end do

!    write(*,*) totm

  end subroutine molecule


 subroutine read_geo

  implicit none

  integer :: i,is,ii,iis,ii3,is3,js,j,k,nt,ks,ls
  real(8) :: xt(256),yt(256),zt(256),e1(256),e2(256),e3(256)

!  reading configuration file (geometry)
  open (10,file='input_geo')
   read (10,*) nm
!   read (10,*) ns

    do i=1,nm
       do is=1,3 !ns
          read (10,*) atom(is),xstmp(is,i),ystmp(is,i),zstmp(is,i)  ! read geo of atomic sites in Ang

          xst(is,i) = xstmp(is,i)/bohr2a
          yst(is,i) = ystmp(is,i)/bohr2a
          zst(is,i) = zstmp(is,i)/bohr2a
       end do
    end do

  close(10)

  end subroutine read_geo


  subroutine generate_offatomic_sites

  implicit none

  integer :: i,is,ii,iis,ii3,is3,js,j,k,nt,ks,ls
  real(8) :: xt(nsite,nom),yt(nsite,nom),zt(nsite,nom)
  real(8) :: x12,y12,z12,x13,y13,z13,cross(3)

!  rs = r1 + c1*r21 + c2*r31 + c3*(r21 crossproduct r31)

  do i=1,nm

        x12 = xst(1,i)-xst(2,i)
        y12 = yst(1,i)-yst(2,i)
        z12 = zst(1,i)-zst(2,i)

        x13 = xst(1,i)-xst(3,i)
        y13 = yst(1,i)-yst(3,i)
        z13 = zst(1,i)-zst(3,i)

        cross(1) = y12 * z13 - z12 * y13
        cross(2) = z12 * x13 - x12 * z13
        cross(3) = x12 * y13 - y12 * x13

        do is=1,ns

        if (is .lt. 4) then

        xt(is,i) = xst(is,i) !+coff(is,1)*x12+coff(is,2)*x13+coff(is,3)*cross(1)
        yt(is,i) = yst(is,i) !+coff(is,1)*y12+coff(is,2)*y13+coff(is,3)*cross(2)
        zt(is,i) = zst(is,i) !+coff(is,1)*z12+coff(is,2)*z13+coff(is,3)*cross(3)

        xs(is,i)= xt(is,i)
        ys(is,i)= yt(is,i)
        zs(is,i)= zt(is,i)

        else

        xt(is,i) = xst(1,i)+coff(is,1)*x12+coff(is,2)*x13+coff(is,3)*cross(1)
        yt(is,i) = yst(1,i)+coff(is,1)*y12+coff(is,2)*y13+coff(is,3)*cross(2)
        zt(is,i) = zst(1,i)+coff(is,1)*z12+coff(is,2)*z13+coff(is,3)*cross(3)

        xs(is,i)= xt(is,i)
        ys(is,i)= yt(is,i)
        zs(is,i)= zt(is,i)

        end if

        end do

  end do


  end subroutine generate_offatomic_sites

end module initialization
