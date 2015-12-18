! Name : GBLUP
! Written by Ignacy Mistzal
! Modified by : Joonho Lee - PhD working in Hankyong National Univ.(Korea) zoonolee@gmail.com, +82-10-3408-2895
! Version : Linux 1.0 (2012-04-16)
! Contents : 		1. Estimation of breeding value using pedigree information and genomic relationship matrix

program GBLUP
! As lsqmt but with support for random effects
integer,parameter::effcross=0,& !effects can be cross-classified
                   effcov=1     !or covariables
!Types of random effects
integer, parameter :: g_fixed=1,& ! fixed effect
                       g_diag=2,& ! diagonal
                          g_A=3,& ! additive animal
                         g_As=4,& ! additive sire
                          g_grm=5 ! additive inverse genomic

integer::neff,ntrait,ncoreff,coreff,RGR,miss,neq,io,i,j,k,l
integer,allocatable::address(:,:),nlev(:),indata(:),effecttype(:),nestedcov(:),randomtype(:),randomnumb(:)
real, allocatable:: xx(:,:),xy(:),sol(:),r(:,:),rinv(:,:),weight_cov(:,:),y(:),g(:,:,:),igrm(:,:)  !storage for the equations
 character(len=3),allocatable::NAOTR(:),NAOEFF(:),NAOET(:),ETYPE(:),RTYPE(:),RCTYPE(:)
 character(len=8),allocatable::NOAF(:)
 character(len=50),allocatable::MN(:)
 character(len=8)::OFT,NOMF,NOFEF,NOTRF
 character(len=16)::WF1,WF2,WF3,WF4,WF5
 character(len=220)::MDD,ODD,GDD,PDD
 character(len=30)::GDFN,PDFN,MDFN
 character(len=250)::MDFNIP,GDFNIP,ODFNIP,PDFNIP,PF

!reading parameters
print *, "***********************************************************************************************************************" 
print *, "GBLUP_GRM_Lv1 : GBLUP using GRM version 1.0"
print *, "Developed by Joonho Lee : zoonolee@gmail.com, +82-10-3408-2895"
print *, "Genotype data form is the results of Post_fastPHASE and of converted genotype 0,1,2"
print *, "Maximum number of trait is 20"
print *, "Phenotype data form is animalID, Fixedeffect...,Trait..."
print *, "Fixed effect must be numbered as integer"
print *, "Trait name must be consisted of 3 character ex) EMA,BFT,..."
print *, "Please enter the name of parameter file including path(ex) ./parGRDH"
print *, "***********************************************************************************************************************"
read *, PF
PF=adjustl(PF)
open (99, file=PF) 
read (99,*) OFT
print *, "file tag ", OFT
read (99,*) neff
print *, "Number of effects in model ", neff
allocate(nlev(neff),effecttype(neff),nestedcov(neff),NAOEFF(neff),randomtype(neff),randomnumb(neff),g(neff,20,20))
allocate(NAOEFF(i),nlev(i),ETYPE(i),RTYPE(i),RCTYPE(i),RGR
do i=1,neff
  read (99,*) NAOEFF(i),nlev(i),ETYPE(i),RTYPE(i),RCTYPE(i),RGR
  if(ETYPE(i).eq."cro") effecttype(i)=effcross
  if(ETYPE(i).eq."cov") effecttype(i)=effcov
  ! Assign random effects and G and R matrices
  !types of random effect
  if(RTYPE(i).eq."fix") randomtype(i)=g_fixed
  if(RTYPE(i).eq."dia") randomtype(i)=g_diag
  if(RTYPE(i).eq."ani") then 
    randomtype(i)=g_A
    read(99,*) PDD
    read(99,*) PDFN
    print *, "directory of pedigree data file ", PDD
    print *, "file name of pedigree ", PDFN
    PDD=adjustr(PDD)
    PDFN=adjustl(PDFN)
    PDFNIP=PDD//PDFN
    PDFNIP=adjustl(PDFNIP)
    open(97,file=PDFNIP,status='old')
  if(RTYPE(i).eq."a_s") randomtype(i)=g_As
  if(RTYPE(i).eq."grm") then 
    randomtype(i)=g_grm
    read(99,*) GDD
    read(99,*) GDFN
    print *, "directory of genomic relationship matrix data file ", GDD
    print *, "file name of genomic relationship matrix ", GDFN
    allocate(igrm(nlev(neff),nlev(neff)))
    GDD=adjustr(GDD)
    GDFN=adjustl(GDFN)
    GDFNIP=GDD//GDFN
    GDFNIP=adjustl(GDFNIP)
    open(96,file=GDFNIP,status='old')
    do j=1,nlev(neff)
      read(96,*) igrm(nlev(neff),:)
    enddo
    print *, "reading genomic relationship matrix is finished"
  endif
  !number of correlated effects per effect
  if(RCTYPE(i).eq."not") randomtype(i)=0
  if(RCTYPE(i).eq."cor") randomtype(i)=RGR
  print *, "Effect ( "//NAOEFF(i)//" ) : level = ", nlev(i),", effecttype = "//ETYPE(i)//", randomtype = "//RTYPE(i)//&
           " ,randomgroup = ",RGR
enddo
read (99,*) ntrait
print *, "Number of traits ", ntrait
allocate(r(ntrait,ntrait),rinv(ntrait,ntrait),y(ntrait),weight_cov(neff,ntrait),NAOTR(ntrait))
read (99,*) NAOTR(:)
print *, "Name of traits ", NAOTR(:)
g=0
ncoreff=sum(randomnumb)
do i=1,ncoreff
  do j=1,ntraits
    read (99,*) coreff
    read (99,*) g(coreff,j,:)
  enddo
enddo
read (99,*) ODD
print *, "directory of output file ", ODD
close(99)

ODFNIP=ODD//"Solutions"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(95,file=ODFNIP,status='unknown')
ODFNIP=ODD//"Incident_Matrix"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(94,file=ODFNIP,status='unknown')
ODFNIP=ODD//"XX"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(93,file=ODFNIP,status='unknown')
ODFNIP=ODD//"XY"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(92,file=ODFNIP,status='unknown')

neq=ntrait*sum(nlev)
allocate (xx(neq,neq), xy(neq), sol(neq),indata(neff*ntrait),&
         address(neff,ntrait))
xx=0; xy=0

call setup_g                  ! invert G matrices
open(2,file='pedname')        !pedigree file
open(1,file='data_pr3')       !data file

do
  read(1,*,iostat=io)indata,y
  if (io.ne.0) exit
  call find_addresses
  call find_rinv
  do i=1,neff
    do j=1,neff
      do k=1,ntrait
        do l=1,ntrait
          xx(address(i,k),address(j,l))=xx(address(i,k),address(j,l))+&
          weight_cov(i,k)*weight_cov(j,l)*rinv(k,l)
        enddo
      enddo
    enddo
    do k=1,ntrait
      do l=1,ntrait
        xy(address(i,k))=xy(address(i,k))+rinv(k,l)*y(l)*weight_cov(i,k)
      enddo
    enddo
  enddo
enddo
print*,'left hand side'

do i=1,neq
  print '(100f5.1)',xx(i,:)
enddo
!
print '( '' right hand side:'' ,100f6.1)',xy

! Random effects' contributions
do i=1,neff
  select case (randomtype(i))
    case (g_fixed)
      continue  ! fixed effect, do nothing
    case (g_diag)
      call add_g_diag(i)
    case (g_A, g_As)
      call add_g_add(randomtype(i),i)
    case (g_grm)
      call add_g_grm(i)
    case default
      print*,'unimplemented random type',randomtype(i)
  endselect
enddo

call solve_dense_gs(neq,xx,xy,sol)
!solution by Gauss-Seidel

do i=1,neff
  do j=1,nlev(i)
    do k=1,ntrait
      l=l+1
      print *, "Effect ( "//NAOEFF(i)//" ) : level = ",j,", trait = "//NAOTR(k)//", solution = ",sol(l)
    enddo
  enddo
enddo

contains

subroutine setup_g
! inverts g matrices
real :: w(20)
integer :: rank
do i=1,neff
  if (randomnumb(i).ne.0) then
    call printmat('original G',g(i,:,:),20,randomnumb(i)*ntrait)
    call ginv(g(i,:,:),20,1e-5,rank)
    call printmat('inverted G',g(i,:,:),20,randomnumb(i)*ntrait)
  endif
enddo
end subroutine

subroutine ginv(a,n,tol,rank)
! returns generalized inverse of x(n,n). tol is working zero
! and irank returns the rank of the matrix. rework of rohan fernando's
! f77 subroutine by i. misztal 05/05/87-05/23/00
implicit none
integer n, rank
real a(n,n),w(n),tol
integer i,ii,j
rank=n
do i=1,n
  do j=1,i-1
    a(i:n,i)=a(i:n,i)-a(i,j)*a(i:n,j)
  enddo
  if (a(i,i).lt.tol) then
    a(i:n,i)=0.0
    rank=rank-1
  else
    a(i,i)=sqrt(a(i,i))
    a(i+1:n,i)=a(i+1:n,i)/a(i,i)
  endif
enddo
do i=1,n
  if (a(i,i).eq.0.) then
    a(i+1:n,i)=0
  else
    a(i,i)=1.0/ a(i,i)
    w(i+1:n)=0
    do ii=i+1,n
      w(ii:n)=w(ii:n)-a(ii:n,ii-1)*a(ii-1,i)
      if (a(ii,ii).eq.0.) then
        a(ii,i)=0.
      else
        a(ii,i)=w(ii)/a(ii,ii)
      endif
    enddo
  endif
enddo
do j=1,n
  do i=j,n
    a(i,j)=dot_product(a(i:n,j),a(i:n,i))
  enddo
enddo
do i=1,n
  a(i,i+1:n)=a(i+1:n,i)
enddo
end subroutine

function address1(e,l,t)
! returns address for given level l of effect e and trait t
integer :: e,l,t, address1
address1= sum(nlev(1:e-1))*ntrait+(l-1)*ntrait+t
end function

subroutine add_g_diag(eff)
! adds diagonal (IID) contributions to MME
integer :: eff, i,j,k,l,m,t1,t2
do i=1,nlev(eff)
  do j=0,randomnumb(eff)-1
    do k=0,randomnumb(eff)-1
      do t1=1,ntrait
        do t2=1,ntrait
          m=address1(eff+j,i,t1); l=address1(eff+k,i,t2)
          xx(m,l)=xx(m,l)+g(eff,t1+j*ntrait,t2+k*ntrait)
        enddo
      enddo
    enddo
  enddo
enddo
end subroutine

subroutine add_g_add(type,eff)
! generates contributions for additive sire or animal effect
integer :: type,eff,i,j,t1,t2,k,l,m,n,io,animal,sire,dam,par_stat,p(3)
real ::w(3),res(4)
!
select case (type)
  case (g_A)
    w=(/1., -.5, -.5/)
    res=(/2., 4/3., 1., 0./)
  case (g_As)
    w=(/1., -.5, -.25/)
    res=(/16/11., 4/3., 15/15., 1./)
  end select
do
  read(2,*,iostat=io) animal, sire, dam,par_stat !status of parents
  if (io /= 0) exit
  p(1)=animal
  p(2)=sire
  p(3)=dam
  print*,p
  do i=0,randomnumb(eff) - 1
    do j=0,randomnumb(eff) - 1
      do t1=1,ntrait
        do t2=1,ntrait
          do k=1,3
            do l=1,3
              if (p(k)/=0 .and.p(l)/=0) then
                m=address1(eff+i,p(k),t1)
                n=address1(eff+j,p(l),t2)
                xx(m,n)=xx(m,n)+g(eff,t1+i*ntrait,t2+j*ntrait)*&
                w(k)*w(l)*res(par_stat)
              endif
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
enddo
end subroutine

subroutine add_g_grm(eff)
! genomic relationship matrix contributions to MME
!real :: igrm(nlev(eff),nlev(eff))
integer :: eff,h,i,j,k,l,m,t1,t2
do h=1,nlev(eff)
  do i=1,nlev(eff)
    do j=0,randomnumb(eff)-1
      do k=0,randomnumb(eff)-1
        do t1=1,ntrait
          do t2=1,ntrait
            m=address1(eff+j,h,t1); l=address1(eff+k,i,t2)
            xx(m,l)=xx(m,l)+g(eff,t1+j*ntrait,t2+k*ntrait)*igrm(h,i)
          enddo
        enddo
      enddo
    enddo
  enddo
enddo
end subroutine

!subroutine find_addresses
!integer :: i,j,baseaddr
!do i=1,neff
!  do j=1,ntrait
!    if (datum(i,j) == miss) then  !missing effect
!      address(i,j)=0 !dummy address
!      weight_cov(i,j)=0.0
!      cycle
!    endif
!    select case (effecttype(i))
!      case (effcross)
!        address(i,j)=address1(i,datum(i,j),j)
!        weight_cov(i,j)=1.0
!      case (effcov)
!        weight_cov(i,j)=datum(i,j)
!        if (nestedcov(i) == 0) then
!          address(i,j)=address1(i,1,j)
!        elseif (nestedcov(i)>0 .and. nestedcov(i).lt.neff) then
!          address(i,j)=address1(i,nestedcov(i),j)
!        else
!          print*,'wrong description of nested covariable'
!          stop
!        endif
!      case default
!        print*,'unimplemented effect ',i
!        stop
!    end select
!  enddo
!enddo
!end subroutine

subroutine find_addresses 
integer :: i,j,baseaddr
do i=1,neff
  do j=1,ntrait
    if (datum(i,j) == miss) then            !missing effect
      address(i,j)=0    !dummy address
      weight_cov(i,j)=0.0
      cycle
    endif
    baseaddr=sum(nlev(1:i-1))*ntrait+j       !base address (start)
    select case (effecttype(i))
      case (effcross)
 !       address(i,j)=baseaddr+(datum(i,j)-1)*ntrait
        address(i,j)=address1(i,int(datum(i,j)),j)
        weight_cov(i,j)=1.0
      case  (effcov)
        weight_cov(i,j)=datum(i,j)
        if (nestedcov(i) == 0) then
          address(i,j)=baseaddr
        elseif (nestedcov(i)>0 .and. nestedcov(i).lt.neff) then
          address(i,j)=baseaddr+(datum(nestedcov(i),j)-1)*ntrait
        else
          print*,'wrong description of nested covariable'
          stop
        endif
      case default
        print*,'unimplemented effect ',i
        stop
    end select
  enddo 
enddo
  end subroutine

function datum(ef,tr)
real :: datum
integer :: ef,tr
! calculates the value effect ef and trait tr
datum=indata(ef +(tr-1)*neff)
end function

subroutine find_rinv
! calculates inv(Q R Q), where Q is an identity matrix zeroed for
! elements corresponding to y(i)=miss
integer :: i,irank
real:: w(10)
rinv=r
do i=1,neff
  if (y(i) == miss) then
    rinv(i,:)=0; rinv(:,i)=0
  endif
enddo
call ginv(rinv,ntrait,1e-5,irank)
end subroutine

subroutine solve_dense_gs(n,lhs,rhs,sol)
! finds sol in the system of linear equations: lhs*sol=rhs
! the solution is iterative by Gauss-Seidel
integer :: n
real :: lhs(n,n),rhs(n),sol(n),eps
integer :: round
!
round=0
do
  eps=0; round=round+1
  do i=1,n
    if (lhs(i,i).eq.0) cycle
    solnew=sol(i)+(rhs(i)-sum(lhs(i,:)*sol))/lhs(i,i)
    eps=eps+ (sol(i)-solnew)**2
    sol(i)=solnew
  end do
  if (eps.lt. 1e-10) exit
end do
print*,'solutions computed in ',round,' rounds of iteration'
end subroutine

subroutine printmat(text,matrix,n,m)
! prints text and matrix
character text*(*)
integer :: n,m,i
real :: matrix(n,n)
print*,text
do i=1,m
   print '(11F7.2)',matrix(i,1:m)
enddo
end subroutine

end program GBLUP
