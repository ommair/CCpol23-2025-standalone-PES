module common_arrays

   implicit none;

!******************************************************************************
  ! unit of time      (to)    = 1.000000000 x 10**(-12) seconds   (pico-seconds)
  ! unit of length    (lo)    = 1.000000000 x 10**(-10) metres    (angstroms)
  ! unit of mass      (mo)    = 1.660540200 x 10**(-27) kilograms (daltons)
  ! unit of charge    (qo)    = 1.602177330 x 10**(-19) coulombs  (electrons)
  ! unit of energy    (eo)    = 1.660540200 x 10**(-23) joules    (10 j mol^-1)
  ! unit of pressure  (po)    = 1.660540200 x 10**(  7) pascals  (163.882576
  ! atm).
  !                           = 1.638825760 x 10**(  2) atmospheres
!*****************************************************************************

!************************************************************************
!*
!*          index of parameters and common variables
!*          ========================================
!*
  integer, parameter :: nom=5000     ! max number of molecules (32,108,256,500,864,1372,2048,2916)
  integer, parameter :: nsite=50   ! max number of sites 
  integer, parameter :: npimax=50  ! number of pair interactions
  integer, parameter :: ncsmax=3   ! max coreshell pairs 
  integer, parameter :: npolc=3  ! number of polarization sites per molecule
  integer, parameter :: nhist = 1000 !dimension of rdf histogram 
  integer, parameter :: ksqmax = 200 ! max squred vlues of k vectors indices in FOurier sum 
  real(8), parameter :: pi=dacos(-1.d0)  ! pi
  real(8), parameter :: boltz=8.31451115d-1  ! Boltzman Const eo/(k*mol) (dlpoly)
  real(8), parameter :: avsno=6.0222e+23    ! Avogadro number
  real(8), parameter :: r4pie0=138935.4835d0 ! 1/4pi*epsilon in internal units
  real(8), parameter :: ewp = 1.0d-6   ! precision for Ewald sum
  real(8), parameter :: prsunt=0.163882576d0  ! internal units to katm
  real(8), parameter :: bohr2a = 0.529177249d0
  real(8), parameter :: h2kcal = 627.510d0   ! Hartree to kcal/mol conversion
  real(8), parameter :: twosqpi= 2.0d0/sqrt(pi)  

  integer  :: nm                   ! number of molecule in the system
  integer  :: ns                   ! sites per molecules
  integer  :: natom                ! number  of massive sites per molecules  
  integer :: npi                   ! number of pair interactions in field file
  character(len=15) :: system_name ! name of the system in the field file
  character(len=3) :: an(nsite)    ! name of atoms in a molecule
  character(len=3) :: ann(nsite,nom)    ! name of atoms in a molecule

  real(8) :: alpha,alpha3, r0_tilda, delta33  ! for ewald summation (intraction of permanent charges and point induced dipoles)
  !real(8), parameter :: alphad = 0.35818465241608965 !0.28d0  ! for induction model with GSF 0.28
  real(8) :: alphad
  integer :: kmax1,kmax2,kmax3
  character(len=15) :: initial_config
  real(8) :: hsize  ! size of step for force check
  logical :: check_forces ! (T or F)
  logical :: rand_vcom  ! random velocities
  logical :: rand_quat  ! random quaternios
  logical :: restart

!  RDF varaibles **** starts here
  integer ::  hist_freq ! frequncy of history file data
  integer :: ityp(nsite) ! types of sites (for RDF)
  integer :: numrdf ! number steps used fpr rdf statistics 
  real(8), parameter :: drdf = 0.05d0
  integer :: mxrdf 
  real(8) :: delrdf ! default bin width
! RDF varaiables **** ends here  

  real(8) :: diameter ! molecular diameter
  character(len=3) :: an1(npimax),an2(npimax) ! site pair names in vdw interactions  
  real(8) :: massites(nsite) ! mass of sites 
  real(8) :: coords(3,nsite) ! coordinates of given molecular geometry 
  real(8) :: pcoord(3,nsite) ! principal corrdiates of each sites in molecule
  real(8) :: pmi(3)          ! principal moment of inertia
  real(8) :: charge(nsite)   ! charge on sites
  real(8) :: sig(npimax),eps(npimax) ! lj parameters 
  real(8) :: ssig(nsite,nsite),eeps(nsite,nsite)  ! test
  real(8) :: temp0  ! desired temperature
  real(8) :: pres0  ! desired pressure
  real(8) :: pres   
  real(8) :: ttemp  ! temp related to translation KE
  real(8) :: rtemp  ! temp rlated to rotational KE
  real(8) :: deltmp ! tolerence in temp
  real(8) :: volm,volm0,volmt   ! volume of simulation box
  real(8) :: totm   ! mass of single molecule
  real(8) :: dt     ! time step size
  real(8) :: rcut,rcutsd ! real space cutoff radius / squared
  real(8) :: rcut3,rcutsd3
  real(8) :: box,box0,boxt,hbox    ! box length in each direction / half of box 
  integer :: nsteps,nequil ! total steps, no of equilibration steps
  integer :: step,stepr,navge,laststp ! step number , no of steps over which avg is taken

 !  arrays for molecular centers  
  real(8) :: x(nom),y(nom),z(nom) ! position of com of molecule
  real(8) :: q0(nom),q1(nom),q2(nom),q3(nom) ! quaternions   
  real(8) :: vx(nom),vy(nom),vz(nom)  !  velocity of com of molecule 
  real(8) :: jx(nom),jy(nom),jz(nom)  ! angular momentum of com of molecule
  real(8) :: fx(nom),fy(nom),fz(nom)  ! force on com of molecule
  real(8) :: tx(nom),ty(nom),tz(nom)  ! torque on com of molecule

  real(8) :: xs(nsite,nom), ys(nsite,nom), zs(nsite,nom) ! sites corrds
  real(8) :: xst(nsite,nom),yst(nsite,nom),zst(nsite,nom)
  real(8) :: xstmp(nsite,nom),ystmp(nsite,nom),zstmp(nsite,nom)
  real(8) :: vxs(nsite,nom), vys(nsite,nom), vzs(nsite,nom) ! sites velocites
  real(8) :: fxs(nsite,nom),fys(nsite,nom),fzs(nsite,nom) ! forces on sites
  real(8) :: chgs(nsite,nom)                              ! charges on sites

  real(8) :: xcom,ycom,zcom,vxcom,vycom,vzcom,sysmass

  integer :: integrator ! integer to select integrator
  character(len=15) :: pot_name   ! water potential
  character(len=10) :: ensemble   ! ensemble for simulation

  real(8) :: vir_test,vircom
  real(8) :: vdwpe,vir,tke,rke ! PE, virial, trans. KE, rot. KE
  real(8) :: rcpe,kcpe,scpe ! real, kspace, self coulumb pe
  real(8) :: vir_vdw,vir_rc,vir_kc,vir_self,vir_rind,vir_kind,vir_selfind,vir_shtind
  real(8) :: indpe          ! induction energy (point dipole polarization model)
  real(8) :: fs2pe,fs3pe          ! 3b energy  
  real(8) :: vdwlrc,virlrc  ! long corrections for VDW and virial
  real(8) :: rho,enth
  real(8) :: sqf,sqt
  real(8) :: stphv,stprho,stpenth        
  real(8) :: stpvdwpe ! step vdw pe
  real(8) :: stpcpe   ! sum of real, kspace, self coulumb pe
  real(8) :: stpindpe ! step induction energy
  real(8) :: stp3bpe ! step 3b energy
  real(8) :: stppe          ! step total pe
  real(8) :: stptke ! step inst trans. ke
  real(8) :: stpttp ! step inst trans. temp
  real(8) :: stprke ! inst rot. ke
  real(8) :: stprtp ! inst rot. temp
  real(8) :: stpvir ! step vir
  real(8) :: stpte  ! step total energy
  real(8) :: stpprs ! step presurre 
  real(8) :: stptotke ! step tot kinetic energy
  real(8) :: stptmp ! step total temperature
  real(8) :: stpvolm ! step volume


  ! sapt5s,  ccpol5 or ccpol8s parameters variables

  real(8) :: a0(npimax),a1(npimax),a2(npimax),a3(npimax)
  real(8) :: expalp(npimax), beta(npimax),aalp(npimax)
  real(8) :: d1(npimax), d6(npimax),d8(npimax),d10(npimax)
  real(8) :: c0(npimax),c1(npimax),c2(npimax),c3(npimax)
  real(8) :: c6(npimax),c8(npimax),c10(npimax)
  real(8) :: domo(npimax),qa(npimax),qb(npimax)

  ! induction model parameter

  logical :: induction      ! true or false (wants to include induction model) 
  character(len=10) :: induction_type         ! (damped or undamped)
  character(len=18) :: ind_model              ! physical dipole or point dipole  
  character(len=10) :: method       ! method for long range correctin to avoid total energy drift 

  real(8) :: E0x(nsite,nom),E0y(nsite,nom),E0z(nsite,nom)    ! efield of static charges
  real(8) :: Eidmx(nsite,nom),Eidmy(nsite,nom),Eidmz(nsite,nom) ! efld of ind dipoles

  real(8) :: rE0x(nsite,nom),rE0y(nsite,nom),rE0z(nsite,nom)    ! efield of static charges (real space)
  real(8) :: kE0x(nsite,nom),kE0y(nsite,nom),kE0z(nsite,nom)    ! efield of static charges (recp space)
  real(8) :: slfE0x(nsite,nom),slfE0y(nsite,nom),slfE0z(nsite,nom)    ! efield of static charges (self correction)
  real(8) :: srfE0x(nsite,nom),srfE0y(nsite,nom),srfE0z(nsite,nom) 
  real(8) :: shtE0x(nsite,nom),shtE0y(nsite,nom),shtE0z(nsite,nom)

  real(8) :: rEidmx(nsite,nom),rEidmy(nsite,nom),rEidmz(nsite,nom) ! efld of ind dipoles (real space)
  real(8) :: kEidmx(nsite,nom),kEidmy(nsite,nom),kEidmz(nsite,nom) ! efld of ind dipoles (recp space)
  real(8) :: slfEidmx(nsite,nom),slfEidmy(nsite,nom),slfEidmz(nsite,nom) ! efld of ind dipoles (self correction)
  real(8) :: srfEidmx(nsite,nom),srfEidmy(nsite,nom),srfEidmz(nsite,nom)
  real(8) :: shtEidmx(nsite,nom),shtEidmy(nsite,nom),shtEidmz(nsite,nom)
  real(8) :: Etotx(nsite,nom),Etoty(nsite,nom),Etotz(nsite,nom)

  real(8) :: idmx(nsite,nom),idmy(nsite,nom),idmz(nsite,nom) ! induced dipoles

  real(8) :: polarizability(nsite)
  real(8) :: apol(nsite,nom)
  integer :: npol

  real(8), parameter :: delta2 = 1.65d0  ! 1/Bohr to 1/Ang
  real(8), parameter :: delta3 = 1.55d0

  ! 3B interaction arrays
  character(len=3) :: a3b(npimax),b3b(npimax) ! site pair names in 3b interactions
  logical :: nonadditive_3B      ! true or false (wants to include 3b potential)
  integer :: ipar
  real(8) :: params(10)
  real(8) :: eta3, beta3, g3 ! FS^2 parameters

  real(8) :: bS3(nsite,nsite),gS3(nsite,nsite),r0S3(nsite,nsite)
  real(8) :: bS3_s(nsite,nsite),gS3_s(nsite,nsite),r0S3_s(nsite,nsite)
  real(8) :: bS3_l(nsite,nsite),gS3_l(nsite,nsite),r0S3_l(nsite,nsite)
  character(len=3) ::  buf1(nsite), buf2(nsite)
  real(8) :: cs3(5000),cs3_s(5000),cs3_l(5000)
  real(8) :: ind_beta(nsite,nsite),ind_gamma(nsite,nsite),ind_r0(nsite,nsite) 
  
  real(8) :: Ravg3b,coff(25,3) ! cofficients for offatomic sites
  character(len=3) :: atom(nsite) 

  logical, parameter :: lrd_3b = .true.   

end module common_arrays
