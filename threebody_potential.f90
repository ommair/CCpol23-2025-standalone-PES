module  threebody_potential

  use, intrinsic :: iso_fortran_env
  use common_arrays

  implicit none
  private

  public :: threebody_term
!  public :: FS2_term, FS3_term

contains  

  subroutine threebody_term(xxx,yyy,zzz,xxxs,yyys,zzzs,qfs2,qfs3)

    implicit none
    real(8) :: xxx(nom),yyy(nom),zzz(nom)
    real(8) :: xxxs(nsite,nom),yyys(nsite,nom),zzzs(nsite,nom)
    real(8) :: boxii
    real(8) :: qfs2,qfs3

    integer :: is,js,i,j,k,isteps
!    real(8) :: xx(nom),yy(nom),zz(nom)
!    real(8) :: xxs(nsite,nom),yys(nsite,nom),zzs(nsite,nom)
    real(8), allocatable :: xx(:), yy(:), zz(:)
    real(8), allocatable :: xxs(:, :), yys(:, :), zzs(:, :)
    real(8) :: boxi

    real(8) :: tol3,eprec3

   allocate(xx(nom));
   allocate(yy(nom));
   allocate(zz(nom));
   allocate(xxs(nsite,nom));
   allocate(yys(nsite,nom));
   allocate(zzs(nsite,nom));

    do i=1,nm

       xx(i) = xxx(i)
       yy(i) = yyy(i)
       zz(i) = zzz(i)

       do is=1,ns
          xxs(is,i) = (xxxs(is,i))
          yys(is,i) = (yyys(is,i))
          zzs(is,i) = (zzzs(is,i))
       end do


    end do

    boxi = boxii/bohr2a
   
!    call FS2_term(xx,yy,zz,xxs,yys,zzs,boxi,qfs2)

    call FS3_term(xx,yy,zz,xxs,yys,zzs,qfs3)


    deallocate(xx);
    deallocate(yy);
    deallocate(zz);
    deallocate(xxs);
    deallocate(yys);
    deallocate(zzs);


  end subroutine threebody_term

  subroutine FS2_term(xx,yy,zz,xxs,yys,zzs,qfs2)
    implicit none
    real(8) :: xx(nom),yy(nom),zz(nom)
    real(8) :: xxs(nsite,nom),yys(nsite,nom),zzs(nsite,nom)
    real(8) :: boxi
    real(8) :: alp3,qfs2,lqvirial,qvirial,fs2_vir

!    real(8) :: xxx(nom),yyy(nom),zzz(nom)
!    real(8) :: xxxs(nsite,nom),yyys(nsite,nom),zzzs(nsite,nom)
!    real(8) :: charge3b(nsite)
    real(8), allocatable :: xxx(:), yyy(:), zzz(:);
    real(8), allocatable :: xxxs(:, :), yyys(:, :), zzzs(:, :);
    real(8), allocatable :: charge3b(:);
    real(8) :: boxii  ! charges on sites
   
    integer :: k
    real(8) qqfs2,qqfs20,efs2
    integer :: is,js,ks,i,j,isteps,l
    real(8) :: dsij2,dski2,dskj2,dsij,dski,dskj
    real(8) :: dsij1i,rc2_3b

    real(8) :: dij2,dki2,dkj2
    real(8) :: dij,dij1i,dij3i
    real(8) :: dkj,dkj1i,dkj3i
    real(8) :: dki,dki1i,dki3i

!    real(8) :: qis,qjs,qks,boxi  ! charges on sites
    real(8) :: qis,qjs,qks  ! charges on sites
    real(8) :: eqij,eqijsum 

    real(8) :: dxij,dyij,dzij,dxki,dyki,dzki,dxkj,dykj,dzkj
    real(8) :: xij,yij,zij,xki,yki,zki,xkj,ykj,zkj

    real(8) :: xsij,ysij,zsij
    real(8) :: xqi,yqi,zqi,xqj,yqj,zqj          ! position of additional charges w.r.t com
    real(8) :: dxkqi,dykqi,dzkqi,dxkqj,dykqj,dzkqj     ! sep of site k and xqi ...
    real(8) :: dxkci,dykci,dzkci,dxkcj,dykcj,dzkcj     ! sep of site k and com i ...

    real(8) :: dkqi2,dkqj2
    real(8) :: dkqi,dkqj,dkqi1i,dkqj1i,dkqi2i,dkqi3i,dkqj2i,dkqj3i

    real(8) :: dkci2,dkcj2
    real(8) :: dkci,dkcj,dkci1i,dkcj1i,dkci2i,dkci3i,dkcj2i,dkcj3i
    real(8) :: dotkqi_ij,dotkqj_ij             
    real(8) :: term0,termxy1,termxy2    
    real(8) :: gradxqij,gradyqij,gradzqij
    real(8) :: torqxija,torqyija,torqzija
    real(8) :: torqxijb,torqyijb,torqzijb  

    real(8) :: fackci,fackcj,fackqi,fackqj

    real(8) :: fxi1,fyi1,fzi1,fxi2,fyi2,fzi2,fxi3,fyi3,fzi3,fxi4,fyi4,fzi4
    real(8) :: fxj1,fyj1,fzj1,fxj2,fyj2,fzj2,fxj3,fyj3,fzj3,fxj4,fyj4,fzj4
    real(8) :: fxk,fyk,fzk

    real(8) :: vpxi,vpyi,vpzi 
    real(8) :: vpxj,vpyj,vpzj

    real(8) :: qfxt, qfyt, qfzt

    real(8) :: erf_ki_kj, derf_ki,derf_kj

    !   long range damping as used by Robert starts
    real(8), parameter :: epsprim = 1.d-10
    real(8), parameter :: delta33 = 10.d0 !log(1.d0/epsprim)/log(2.d0)
    real(8) :: r0dmp,flrki,flrkj,dflrki,dflrkj,dflrkiv,dflrkjv
    real(8) :: pom,pomkj,pomki
    real(8) :: termx1,termy1,termz1,termx2,termy2,termz2
    !   long range damping as used by Robert ends

    real(8), parameter :: eta3 =0.48056977605478400d0  !0.693231401521011**2.d0  ! eta OO  HH OH
    real(8), parameter :: beta3 = 0.85836606579445000d0 !/bohr2a !2.d0*0.429183032897225  ! beta OO HH OH
    real(8), parameter :: g3 = 0.400000000000000 !*bohr2a 


    allocate(charge3b(nsite));

    alp3= alpha3
    rc2_3b = rcutsd3
    r0dmp = r0_tilda

    do i=1,ns
       charge3b(i) = 0.d0
    end do

    charge3b(1) = 0.258504607653328389d0
    charge3b(2) = 0.564050362569283648d0
    charge3b(3) = 0.564050362569283648d0
    charge3b(4) = -0.693302666395947953d0
    charge3b(5) = -0.693302666395947953d0
 
    qqfs2 = 0.d0
    qqfs20 = 0.d0

    allocate(xxx(nom));
    allocate(yyy(nom));
    allocate(zzz(nom));
    allocate(xxxs(nsite, nom));
    allocate(yyys(nsite, nom));
    allocate(zzzs(nsite, nom));

    do k=1,nm

       do i=1,nm-1

          if (k .NE. i) then

             do j=i+1,nm

                if (j .NE. k) then

                   dxij = xx(i)-xx(j)   ! com sep of  i - j
                   dyij = yy(i)-yy(j)
                   dzij = zz(i)-zz(j)

                   dxki = xx(k)-xx(i)   ! com sep of  k - i
                   dyki = yy(k)-yy(i)
                   dzki = zz(k)-zz(i)

                   dxkj = xx(k)-xx(j)   ! com sep of  k - j
                   dykj = yy(k)-yy(j)
                   dzkj = zz(k)-zz(j)
 
                   xij = dxij 
                   yij = dyij 
                   zij = dzij 

                   xki = dxki 
                   yki = dyki 
                   zki = dzki 

                   xkj = dxkj 
                   ykj = dykj 
                   zkj = dzkj 

                   dij2 = xij**2 + yij**2 + zij**2   ! rij2
                   dki2 = xki**2 + yki**2 + zki**2   ! rki2
                   dkj2 = xkj**2 + ykj**2 + zkj**2   ! rkj2

                    ! long range damping ****  dkj2 dki2
                    pom = exp(alp3*(sqrt(dki2)-r0dmp))  ! ks and com i+
                    pomki = 1.d0/(1.d0+pom)
                    flrki = pomki**delta33
                    dflrki = -alp3*delta33*pom*pomki*flrki
                    dflrkiv = -alp3*delta33*pom*pomki*sqrt(dki2)

                    pom = exp(alp3*(sqrt(dkj2)-r0dmp))  ! ks and com j+
                    pomkj = 1.d0/(1.d0+pom)
                    flrkj = pomkj**delta33
                    dflrkj = -alp3*delta33*pom*pomkj*flrkj 
                    dflrkjv = -alp3*delta33*pom*pomkj*sqrt(dkj2)

                    !write(*,*) flrki,flrkj,dflrki,dflrkj,alp3,r0dmp,pom,delta33,pomki
 
                    do ks=1, 5  ! loop ks

                       qks= charge3b(ks) !chgs(ks,k)

!                       if (qks .ne. 0.d0 ) then

                      ! Calculate the "exchange charge", i.e., eta*sum(exp(-beta*r_ab))
                    
                          do is=1,3    ! loop is         only for OO, OH, HH

                             do js=1,3 ! loop js
                           
                                xsij =  (xxs(is,i) - xxs(js,j))
                                ysij =  (yys(is,i) - yys(js,j))
                                zsij =  (zzs(is,i) - zzs(js,j))     

                                dsij2 = xsij**2 + ysij**2 + zsij**2 

                                dsij = sqrt(dsij2)
                                dsij1i = 1.d0/dsij 
                               
                                eqij = eta3*exp(-beta3*dsij)

                                ! Calculate the positions of exchange charges 

                                dij  = sqrt(dij2) ! separation between i-j
                                dij1i = 1.d0/dij  ! inverse separation between i-j
                                dij3i = dij1i**3

                                dki  = sqrt(dki2) ! separation between k-i
                                dki1i = 1.d0/dki  ! inverse separation between k-i
                                dki3i = dki1i**3 

                                dkj  = sqrt(dkj2) ! separation between k-j
                                dkj1i = 1.d0/dkj  ! inverse separation between k-j
                                dkj3i = dkj1i**3

                                xqi = xx(i) + g3*xij*dij1i  ! comp of pos of i- charges
                                yqi = yy(i) + g3*yij*dij1i
                                zqi = zz(i) + g3*zij*dij1i

                                xqj = xx(j) - g3*xij*dij1i  ! comp of pos of i- charges
                                yqj = yy(j) - g3*yij*dij1i
                                zqj = zz(j) - g3*zij*dij1i
                                
                               ! separation of site ks and additional charges i- and j-
  
                                dxkqi = (xx(k)+xxs(ks,k)) - xqi       ! sep b/w charge at ks and i-
                                dykqi = (yy(k)+yys(ks,k)) - yqi
                                dzkqi = (zz(k)+zzs(ks,k)) - zqi

                                dkqi2 = dxkqi**2 + dykqi**2 + dzkqi**2
                                dkqi = sqrt(dkqi2)
                                dkqi1i = 1.d0/dkqi
                                dkqi2i = dkqi1i**2
                                dkqi3i = dkqi1i**3

                                dxkqj = (xx(k)+xxs(ks,k)) - xqj       ! sep b/w charge at ks and j-
                                dykqj = (yy(k)+yys(ks,k)) - yqj
                                dzkqj = (zz(k)+zzs(ks,k)) - zqj

                                dkqj2 = dxkqj**2 + dykqj**2 + dzkqj**2
                                dkqj = sqrt(dkqj2)
                                dkqj1i = 1.d0/dkqj
                                dkqj2i = dkqj1i**2
                                dkqj3i = dkqj1i**3

                                ! separation of site ks and additional charges i+ and j+

                                dxkci = (xx(k)+xxs(ks,k)) - xx(i)     ! sep ks and com i+
                                dykci = (yy(k)+yys(ks,k)) - yy(i)
                                dzkci = (zz(k)+zzs(ks,k)) - zz(i)

                                dkci2 = dxkci**2 + dykci**2 + dzkci**2
                                dkci = sqrt(dkci2)
                                dkci1i = 1.d0/dkci
                                dkci2i = dkci1i**2
                                dkci3i = dkci1i**3

                                dxkcj = (xx(k)+xxs(ks,k)) - xx(j)     ! sep ks and com j+
                                dykcj = (yy(k)+yys(ks,k)) - yy(j)
                                dzkcj = (zz(k)+zzs(ks,k)) - zz(j)

                                dkcj2 = dxkcj**2 + dykcj**2 + dzkcj**2
                                dkcj = sqrt(dkcj2)
                                dkcj1i = 1.d0/dkcj
                                dkcj2i = dkcj1i**2
                                dkcj3i = dkcj1i**3

                               ! *** FS2 energy contribution *** 
                         
                               term0 = (dkci1i + dkcj1i - dkqi1i - dkqj1i) 
                               qqfs20 = qqfs20 + qks*eqij*term0
                               qqfs2 = qqfs2 + qks*eqij*term0*flrki*flrkj
 
                             end do ! loop js
                          end do    ! loop is

                    end do  ! loop ks
 
                endif 

             end do  

          endif 

       end do 

    end do    

    deallocate(xxx)
    deallocate(yyy)
    deallocate(zzz)
    deallocate(xxxs)
    deallocate(yyys)
    deallocate(zzzs)
    deallocate(charge3b)

    qfs2 = qqfs2*h2kcal

end subroutine FS2_term

subroutine FS3_term(xx,yy,zz,xxs,yys,zzs,qfs3)
    implicit none
    real(8) :: xx(nom),yy(nom),zz(nom),RAB,RBC,RAC,Ravg
    real(8) :: xxs(nsite,nom),yys(nsite,nom),zzs(nsite,nom)
    real(8) :: boxi,rmin,rmax,smfun,cfact,cf,smfun2
    real(8) :: qfs3,lqvirial,qvirial

!    real(8) :: xxx(nom),yyy(nom),zzz(nom)
!    real(8) :: xxxs(nsite,nom),yyys(nsite,nom),zzzs(nsite,nom)
    real(8), allocatable :: xxx(:), yyy(:), zzz(:)
    real(8), allocatable :: xxxs(:, :), yyys(:, :), zzzs(:, :)
    real(8) :: boxii
    real(8) :: pomij,ppomij,fdmpij,dfdmpij,r_av
    real(8) :: pomik,ppomik,fdmpik,dfdmpik
    real(8) :: pomjk,ppomjk,fdmpjk,dfdmpjk
    real(8) :: dfdmpijv,dfdmpikv,dfdmpjkv
    real(8) :: fdmp,dfdmp,dfdmpv

    integer :: is,js,ks,i,j,k,ii,idx
    real(8) :: qis,qjs,rc2_3b
    
    real(8) :: dxij,dyij,dzij,dxik,dyik,dzik,dxjk,dyjk,dzjk
    real(8) :: xij,yij,zij,xik,yik,zik,xjk,yjk,zjk

    real(8) :: dij2,dik2,djk2

    real(8) :: xsij,ysij,zsij,dsij2,dsij,dsij1i
    real(8) :: xsik,ysik,zsik,dsik2,dsik,dsik1i
    real(8) :: xsjk,ysjk,zsjk,dsjk2,dsjk,dsjk1i

    real(8) :: dampsij,gsij_s,r0sij_s,bsij_s
    real(8) :: dampsik,gsik_s,r0sik_s,bsik_s
    real(8) :: dampsjk,gsjk_s,r0sjk_s,bsjk_s

    real(8) :: gsij_l,r0sij_l,bsij_l
    real(8) :: gsik_l,r0sik_l,bsik_l
    real(8) :: gsjk_l,r0sjk_l,bsjk_l

    real(8) :: eksij,eksik,eksjk,eks
    
    real(8) :: qqfs3, lqqfs3,fs3_vir

    integer, parameter :: ns3=17, lmax=1
    integer, parameter :: nzero = ns3*ns3*ns3*(lmax+1)*(lmax+1)*(lmax+1)
!    integer :: ind_lin (0:lmax,0:lmax,0:lmax,ns3,ns3,ns3)
!    integer :: ind_lin_2 (nzero)
    integer, allocatable :: ind_lin(:,:,:,:,:,:)
!    integer, pointer :: ind_lin_2(:) => NULL()
    integer, allocatable :: ind_lin_2(:)
!    equivalence (ind_lin,ind_lin_2)
    logical :: init
    integer :: iat,jat,kat,ip,jp,kp,new,nlin,itab,iind,nl
!    real(8) :: aj(5000),term1,term2,term3
    real(8) :: term1,term2,term3 
    real(8), allocatable :: aj(:) ! (5000)

    real(8) :: fsij,fsik,fsjk

    real(8) :: erf_ij_ik_jk, der_erfij,der_erfjk,der_erfik
    real(8) :: alp3,r0dmp 
    real(8), parameter :: epsprim = 1.d-10
    real(8), parameter :: delta33 = 10.d0 !log(1.d0/epsprim)/log(2.d0)

    qqfs3 = 0.0D0
    
    alp3 = 6.51292d0*bohr2a

    r0dmp = 7.40446d0/bohr2a

    rmin = 4.5d0/bohr2a    
    rmax = 5.5d0/bohr2a

  
    ! *** converting positions in AU ***
    ! *** ends here ***
   
!    allocate(ind_lin(0:lmax,0:lmax,0:lmax,ns3,ns3,ns3))
    allocate(ind_lin_2(nzero))
    allocate(aj(5000))

    allocate(xxx(nom));
    allocate(yyy(nom));
    allocate(zzz(nom));
    allocate(xxxs(nsite,nom));
    allocate(yyys(nsite,nom));
    allocate(zzzs(nsite,nom));

    lqqfs3 = 0.0D0;

    do i=1,nm-2

       do j=i+1,nm-1
 
          do k=j+1,nm

             dxij = xx(i)-xx(j)   ! com sep of  i-j
             dyij = yy(i)-yy(j)
             dzij = zz(i)-zz(j)

             dxik = xx(i)-xx(k)   ! com sep of  i-k
             dyik = yy(i)-yy(k)
             dzik = zz(i)-zz(k)

             dxjk = xx(j)-xx(k)   ! com sep of  j-k
             dyjk = yy(j)-yy(k)
             dzjk = zz(j)-zz(k)

             xij = dxij 
             yij = dyij 
             zij = dzij 

             xik = dxik 
             yik = dyik 
             zik = dzik 

             xjk = dxjk 
             yjk = dyjk 
             zjk = dzjk 

             dij2 = xij**2 + yij**2 + zij**2   ! rij2
             dik2 = xik**2 + yik**2 + zik**2   ! rik2
             djk2 = xjk**2 + yjk**2 + zjk**2   ! rjk2

             RAB = sqrt(dij2)
             RBC = sqrt(djk2)
             RAC = sqrt(dik2)

             Ravg = (RAB+RBC+RAC)/3.d0

             Ravg3b =  Ravg*bohr2a

             pomij = exp(alp3*(Ravg-r0dmp))
             ppomij = 1.d0/(1.d0+pomij)

             if (lrd_3b .eqv. .true.) then 
                fdmp = ppomij**delta33 
             elseif  (lrd_3b .eqv. .false.) then  
                fdmp = 1.d0 
             end if 

             !write(*,*)fdmpij,dfdmpij,fdmpik,dfdmpik,fdmpjk,dfdmpjk

                ! **************
                init = .true. ! GRU:This variable is useless.
                
                if(init) then
                  init = .false.
!                  ind_lin = 0
                ind_lin_2 = 0
                  new = 0  
 
                  do iat=1,ns3
                  do jat=1,ns3
                  do kat=1,ns3
                  do ip=0,lmax
                  do jp=0,lmax
                  do kp=0,lmax

                     idx = (kp + 1)
                     idx = idx + jp * (lmax + 1)
                     idx = idx + ip * (lmax + 1)**2
                     idx = idx + (kat - 1) * (lmax + 1)**3
                     idx = idx + (jat - 1) * (lmax + 1)**3 * ns3
                     idx = idx + (iat - 1) * (lmax + 1)**3 * ns3**2

!                     if (ind_lin(kp,jp,ip,kat,jat,iat).eq.0) then
                     if (ind_lin_2(idx) .EQ. 0) then
                        new=new+1
!                        call fill(ind_lin, iat,jat,kat, ip,jp,kp, new)
                        call fill(ind_lin_2, iat,jat,kat, ip,jp,kp, new)
!			ind_lin_2 = PACK(ind_lin, .TRUE.)
                     end if

                  end do
                  end do
                  end do
                  end do
                  end do
                  end do              

                end if 

!! ****************************  Short range !************************************** 
                if (Ravg .le. rmin) then     

                nlin = new
                aj = 0.d0 
                itab = 0
                !***************

                do is=1,ns3
                   do js=1,ns3
                      do ks=1,ns3

                         bsij_s  = bS3_s(is,js)          ! sites is-js 
                         gsij_s  = gS3_s(is,js)
                         r0sij_s = r0S3_s(is,js)

                         xsij = (xxs(is,i) - xxs(js,j))
                         ysij = (yys(is,i) - yys(js,j))
                         zsij = (zzs(is,i) - zzs(js,j))

                         dsij2 = xsij**2 + ysij**2 + zsij**2

                         dsij = sqrt(dsij2)
                         dsij1i = 1.d0/dsij ! GRU: this variable isn't used

                         dampsij = 1.d0 /(1.d0 +exp(-gsij_s*(dsij-r0sij_s)))

                         bsik_s  = bS3_s(is,ks)          ! sites is-ks 
                         gsik_s  = gS3_s(is,ks)
                         r0sik_s = r0S3_s(is,ks)

                         xsik = (xxs(is,i) - xxs(ks,k))
                         ysik = (yys(is,i) - yys(ks,k))
                         zsik = (zzs(is,i) - zzs(ks,k))

                         dsik2 = xsik**2 + ysik**2 + zsik**2

                         dsik = sqrt(dsik2)
                         dsik1i = 1.d0/dsik

                         dampsik = 1.d0 /(1.d0 +exp(-gsik_s*(dsik-r0sik_s)))

                         bsjk_s  = bS3_s(js,ks)          ! sites js-ks
                         gsjk_s  = gS3_s(js,ks)
                         r0sjk_s = r0S3_s(js,ks)

                         xsjk = (xxs(js,j) - xxs(ks,k))
                         ysjk = (yys(js,j) - yys(ks,k))
                         zsjk = (zzs(js,j) - zzs(ks,k))

                         dsjk2 = xsjk**2 + ysjk**2 + zsjk**2

                         dsjk = sqrt(dsjk2)
                         dsjk1i = 1.d0/dsjk

                         dampsjk = 1.d0 /(1.d0 +exp(-gsjk_s*(dsjk-r0sjk_s)))
 
                         eksij = exp(-bsij_s*dsij)*dampsij 
                         eksik = exp(-bsik_s*dsik)*dampsik
                         eksjk = exp(-bsjk_s*dsjk)*dampsjk  

                         eks = eksij*eksik*eksjk

                         
                         do ip=0,lmax

                            term1 = eks*dsjk**ip

                         do jp=0,lmax

                            term2 = term1*dsik**jp
                            
                         do kp=0,lmax

                            itab = itab + 1 

                            iind = ind_lin_2(itab) 

                            term3 = term2*dsij**kp

                            lqqfs3 = lqqfs3 + cs3_s(iind)*term3*fdmp 


                         end do 
                         end do
                         end do 

                         !end if ! check whether r_avg is with in rcut3 !! ends

                      end do  ! ks
                   end do     ! js
                end do        ! is 

!                do nl=1,nlin
!                qqfs3 = qqfs3 + cs3(nl)*aj(nl)
!                end do

!! ****************************  Long range !**************************************
                else if (Ravg .ge. rmax) then    

                nlin = new
                aj = 0.d0 
                itab = 0
                !***************

                do is=1,ns3
                   do js=1,ns3
                      do ks=1,ns3

                         bsij_l  = bS3_l(is,js)          ! sites is-js 
                         gsij_l  = gS3_l(is,js)
                         r0sij_l = r0S3_l(is,js)

                         xsij = (xxs(is,i) - xxs(js,j))
                         ysij = (yys(is,i) - yys(js,j))
                         zsij = (zzs(is,i) - zzs(js,j))

                         dsij2 = xsij**2 + ysij**2 + zsij**2

                         dsij = sqrt(dsij2)
                         dsij1i = 1.d0/dsij ! GRU: this variable isn't used

                         dampsij = 1.d0 /(1.d0 +exp(-gsij_l*(dsij-r0sij_l)))

                         bsik_l  = bS3_l(is,ks)          ! sites is-ks 
                         gsik_l  = gS3_l(is,ks)
                         r0sik_l = r0S3_l(is,ks)

                         xsik = (xxs(is,i) - xxs(ks,k))
                         ysik = (yys(is,i) - yys(ks,k))
                         zsik = (zzs(is,i) - zzs(ks,k))

                         dsik2 = xsik**2 + ysik**2 + zsik**2

                         dsik = sqrt(dsik2)
                         dsik1i = 1.d0/dsik

                         dampsik = 1.d0 /(1.d0 +exp(-gsik_l*(dsik-r0sik_l)))

                         bsjk_l  = bS3_l(js,ks)          ! sites js-ks
                         gsjk_l  = gS3_l(js,ks)
                         r0sjk_l = r0S3_l(js,ks)

                         xsjk = (xxs(js,j) - xxs(ks,k))
                         ysjk = (yys(js,j) - yys(ks,k))
                         zsjk = (zzs(js,j) - zzs(ks,k))

                         dsjk2 = xsjk**2 + ysjk**2 + zsjk**2

                         dsjk = sqrt(dsjk2)
                         dsjk1i = 1.d0/dsjk

                         dampsjk = 1.d0 /(1.d0 +exp(-gsjk_l*(dsjk-r0sjk_l)))
 
                         eksij = exp(-bsij_l*dsij)*dampsij 
                         eksik = exp(-bsik_l*dsik)*dampsik
                         eksjk = exp(-bsjk_l*dsjk)*dampsjk  

                         eks = eksij*eksik*eksjk

                         
                         do ip=0,lmax

                            term1 = eks*dsjk**ip

                         do jp=0,lmax

                            term2 = term1*dsik**jp
                            
                         do kp=0,lmax

                            itab = itab + 1 

                            iind = ind_lin_2(itab) 

                            term3 = term2*dsij**kp

                            lqqfs3 = lqqfs3 + cs3_l(iind)*term3*fdmp 


                         end do 
                         end do
                         end do 

                         !end if ! check whether r_avg is with in rcut3 !! ends

                      end do  ! ks
                   end do     ! js
                end do        ! is 


!! ****************************  Region of overlap !**************************************
                else    

                cfact = (Ravg - rmin)/(rmax-rmin) 
                smfun = 1.d0 + 2.d0*cfact**3.d0 - 3.d0*cfact**2.d0
                smfun2 = 1.d0 - smfun
                cf = 2.d0/(rmax-rmin)*(cfact**2.d0-cfact)
      

!                write(*,*) Ravg,rmin,rmax 

                nlin = new
                aj = 0.d0 
                itab = 0
                !***************

                do is=1,ns3
                   do js=1,ns3
                      do ks=1,ns3

                         bsij_s  = bS3_s(is,js)          ! sites is-js 
                         gsij_s  = gS3_s(is,js)
                         r0sij_s = r0S3_s(is,js)

                         xsij = (xxs(is,i) - xxs(js,j))
                         ysij = (yys(is,i) - yys(js,j))
                         zsij = (zzs(is,i) - zzs(js,j))

                         dsij2 = xsij**2 + ysij**2 + zsij**2

                         dsij = sqrt(dsij2)
                         dsij1i = 1.d0/dsij ! GRU: this variable isn't used

                         dampsij = 1.d0 /(1.d0 +exp(-gsij_s*(dsij-r0sij_s)))

                         bsik_s  = bS3_s(is,ks)          ! sites is-ks 
                         gsik_s  = gS3_s(is,ks)
                         r0sik_s = r0S3_s(is,ks)

                         xsik = (xxs(is,i) - xxs(ks,k))
                         ysik = (yys(is,i) - yys(ks,k))
                         zsik = (zzs(is,i) - zzs(ks,k))

                         dsik2 = xsik**2 + ysik**2 + zsik**2

                         dsik = sqrt(dsik2)
                         dsik1i = 1.d0/dsik

                         dampsik = 1.d0 /(1.d0 +exp(-gsik_s*(dsik-r0sik_s)))

                         bsjk_s  = bS3_s(js,ks)          ! sites js-ks
                         gsjk_s  = gS3_s(js,ks)
                         r0sjk_s = r0S3_s(js,ks)

                         xsjk = (xxs(js,j) - xxs(ks,k))
                         ysjk = (yys(js,j) - yys(ks,k))
                         zsjk = (zzs(js,j) - zzs(ks,k))

                         dsjk2 = xsjk**2 + ysjk**2 + zsjk**2

                         dsjk = sqrt(dsjk2)
                         dsjk1i = 1.d0/dsjk

                         dampsjk = 1.d0 /(1.d0 +exp(-gsjk_s*(dsjk-r0sjk_s)))
 
                         eksij = exp(-bsij_s*dsij)*dampsij 
                         eksik = exp(-bsik_s*dsik)*dampsik
                         eksjk = exp(-bsjk_s*dsjk)*dampsjk  

                         eks = eksij*eksik*eksjk


                         do ip=0,lmax

                            term1 = eks*dsjk**ip

                         do jp=0,lmax

                            term2 = term1*dsik**jp
                            
                         do kp=0,lmax

                            itab = itab + 1 

                            iind = ind_lin_2(itab) 

                            term3 = term2*dsij**kp
                           
                            lqqfs3 = lqqfs3 + smfun*(cs3_s(iind)*term3*fdmp) 
  

                         end do 
                         end do
                         end do 

                         !end if ! check whether r_avg is with in rcut3 !! ends

                      end do  ! ks
                   end do     ! js
                end do        ! is 

                nlin = new
                aj = 0.d0 
                itab = 0
                !***************

                do is=1,ns3
                   do js=1,ns3
                      do ks=1,ns3

                         bsij_l  = bS3_l(is,js)          ! sites is-js 
                         gsij_l  = gS3_l(is,js)
                         r0sij_l = r0S3_l(is,js)

                         xsij = (xxs(is,i) - xxs(js,j))
                         ysij = (yys(is,i) - yys(js,j))
                         zsij = (zzs(is,i) - zzs(js,j))

                         dsij2 = xsij**2 + ysij**2 + zsij**2

                         dsij = sqrt(dsij2)
                         dsij1i = 1.d0/dsij ! GRU: this variable isn't used

                         dampsij = 1.d0 /(1.d0 +exp(-gsij_l*(dsij-r0sij_l)))

                         bsik_l  = bS3_l(is,ks)          ! sites is-ks 
                         gsik_l  = gS3_l(is,ks)
                         r0sik_l = r0S3_l(is,ks)

                         xsik = (xxs(is,i) - xxs(ks,k))
                         ysik = (yys(is,i) - yys(ks,k))
                         zsik = (zzs(is,i) - zzs(ks,k))

                         dsik2 = xsik**2 + ysik**2 + zsik**2

                         dsik = sqrt(dsik2)
                         dsik1i = 1.d0/dsik

                         dampsik = 1.d0 /(1.d0 +exp(-gsik_l*(dsik-r0sik_l)))

                         bsjk_l  = bS3_l(js,ks)          ! sites js-ks
                         gsjk_l  = gS3_l(js,ks)
                         r0sjk_l = r0S3_l(js,ks)

                         xsjk = (xxs(js,j) - xxs(ks,k))
                         ysjk = (yys(js,j) - yys(ks,k))
                         zsjk = (zzs(js,j) - zzs(ks,k))

                         dsjk2 = xsjk**2 + ysjk**2 + zsjk**2

                         dsjk = sqrt(dsjk2)
                         dsjk1i = 1.d0/dsjk

                         dampsjk = 1.d0 /(1.d0 +exp(-gsjk_l*(dsjk-r0sjk_l)))
 
                         eksij = exp(-bsij_l*dsij)*dampsij 
                         eksik = exp(-bsik_l*dsik)*dampsik
                         eksjk = exp(-bsjk_l*dsjk)*dampsjk  

                         eks = eksij*eksik*eksjk


                         do ip=0,lmax

                            term1 = eks*dsjk**ip

                         do jp=0,lmax

                            term2 = term1*dsik**jp
                           
                         do kp=0,lmax

                            itab = itab + 1 

                            iind = ind_lin_2(itab) 

                            term3 = term2*dsij**kp
                            
                            lqqfs3 = lqqfs3 + (smfun2)*(cs3_l(iind)*term3*fdmp)    
   
 
                         end do 
                         end do
                         end do 

                         !end if ! check whether r_avg is with in rcut3 !! ends

                      end do  ! ks
                   end do     ! js
                end do        ! is 


             end if

          end do             ! k

       end do                ! j

    enddo                   ! i



!    deallocate(ind_lin)
    deallocate(ind_lin_2)
    deallocate(aj)

    deallocate(xxx)
    deallocate(yyy)
    deallocate(zzz)
    deallocate(xxxs)
    deallocate(yyys)
    deallocate(zzzs)

    qqfs3 = qqfs3 + lqqfs3 

    qfs3 = h2kcal*qqfs3
 
  end subroutine FS3_term

  subroutine fill (ntab,iat,jat,kat,ip,jp,kp,nr)

    implicit none

    integer :: iat,jat,kat,ip,jp,kp,nr,is
    integer :: ifrst,ilast,jfrst,jlast,kfrst,klast 
    integer :: ii,jj,kk
    logical :: ilook,jlook,klook 
    integer, parameter :: nsite3=17, lmax=1
    integer :: ntab(0:lmax,0:lmax,0:lmax,nsite3,nsite3,nsite3)   
    integer :: ntype(nsite3)

    data ntype /1,2,2,3,3,4,4,4,4,5,5,5,5,6,6,6,6/  ! must be blocks with consecutive numbers

    ilook=.true.
    jlook=.true.
    klook=.true.

    do is=1,nsite3
         if ((ntype(iat).eq.ntype(is)).and.ilook) then
            ifrst=is
            ilook=.false.
         end if
         if ((ntype(jat).eq.ntype(is)).and.jlook) then
            jfrst=is
            jlook=.false.
         end if
         if ((ntype(kat).eq.ntype(is)).and.klook) then
            kfrst=is
            klook=.false.
         end if
         if (ntype(iat).eq.ntype(is)) then
            if (is.eq.nsite3) then
               ilast=is
            else
               if (ntype(iat).ne.ntype(is+1)) ilast=is
            end if
         end if
         if (ntype(jat).eq.ntype(is)) then
            if (is.eq.nsite3) then
               jlast=is
            else
               if (ntype(jat).ne.ntype(is+1)) jlast=is
            end if
         end if
         if (ntype(kat).eq.ntype(is)) then
            if (is.eq.nsite3) then
               klast=is
            else
               if (ntype(kat).ne.ntype(is+1)) klast=is
            end if
         end if
    end do

    do ii=ifrst,ilast
    do jj=jfrst,jlast
    do kk=kfrst,klast
         ntab(kp,jp,ip, kk,jj,ii)=nr
         ntab(jp,kp,ip, jj,kk,ii)=nr
         ntab(ip,kp,jp, ii,kk,jj)=nr
         ntab(kp,ip,jp, kk,ii,jj)=nr
         ntab(ip,jp,kp, ii,jj,kk)=nr
         ntab(jp,ip,kp, jj,ii,kk)=nr
    end do
    end do
    end do  

  end subroutine fill

end module threebody_potential
