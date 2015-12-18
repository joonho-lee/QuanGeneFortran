! Name : NRM(Numerator Relationship Matrix)
! Developed by : Joonho Lee - PhD working in Hankyong National Univ.(Korea) zoonolee@gmail.com, +82-10-3408-2895
! Version : Windows 2.0 (2012-05-24)
! Contents : 		1. Making NRM using renumbered & ordered pedigree information(OUTPUT41 of renumf90)
!		 	2. Printing NRM of selected animals out with same order
program MakingNRM
implicit none
INTEGER :: NOA,NOA_S,ANO1,ANO2,ANO_S1,ANO_S2,CNOA,io
INTEGER, allocatable :: LOSA(:)
INTEGER :: AID, SID, DID
INTEGER :: sn, dn, k
REAL, allocatable :: NRM(:,:),NRM_S(:,:)
real :: max1
character(len=8)::OFT
character(len=30)::PDFN,SADFN
character(len=220)::PDD,SADD,ODD
character(len=250)::ODFNIP,PDFNIP,SADFNIP,PF
character(len=16),allocatable::ID(:),DCV(:),ID_S(:)
!reading parameters
print *, "***********************************************************************************************************************" 
print *, "NRM_Wv1 : Making NRM & print NRM of selected animals Windows version 1.0"
print *, "Developed by Joonho Lee : zoonolee@gmail.com, +82-10-3408-2895"
print *, "Pedigree information must be renumbered as integer and ordered by generation"
print *, "1st, 2nd and 3rd columns of pedigree file must be renumbered(as integer) animalID, sireID and DamID"
print *, "Every element of data must be less then 16 chracters"
print *, "All of data file(pedigree, selected animals) must have header for 1 first row"
print *, "Please enter the name of parameter file including path(ex) ./parGRDH"
print *, "***********************************************************************************************************************"
read *, PF
PF=adjustl(PF)
open (99, file=PF) 
read (99,*) OFT
print *, "file tag ", OFT
read (99,*) CNOA
print *, "column number of animal ID is ", CNOA
read (99,*) PDD
print *, "directory of pedigree data file ", PDD
read (99,*) PDFN
print *, "file name of pedigree data ", PDFN
read (99,*) SADD
print *, "directory of Selected_Animals file ", SADD
read (99,*) SADFN
print *, "file name of Selected_Animals ", SADFN
read (99,*) ODD
print *, "directory of output file ", ODD
close(99)
if(CNOA.le.3) print *,"parameter error : 1-3 column of pedigree file must be renumbered(as integer) animalID, sireID and DamID"
if(CNOA.gt.4) then 
  allocate(DCV(CNOA-4))
endif
ODD=adjustr(ODD)
PDD=adjustr(PDD)
SADD=adjustr(SADD)
PDFN=adjustl(PDFN)
SADFN=adjustl(SADFN)
PDFNIP=PDD//PDFN
PDFNIP=adjustl(PDFNIP)
SADFNIP=SADD//SADFN
SADFNIP=adjustl(SADFNIP)
open(98,file=PDFNIP,status='old')
open(97,file=SADFNIP,status='old')
!ODFNIP=ODD//"NRM"//OFT//".out"
!ODFNIP=adjustl(ODFNIP)
!open(96,file=ODFNIP,status='unknown')
ODFNIP=ODD//"NRM_selected"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(95,file=ODFNIP,status='unknown')
ODFNIP=ODD//"IBC"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(94,file=ODFNIP,status='unknown')
ODFNIP=ODD//"RC_selected"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(93,file=ODFNIP,status='unknown')
ODFNIP=ODD//"RENUM_selected_animals"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(92,file=ODFNIP,status='unknown')
NOA=0
do 
  read(98,*,iostat=io) 
  if(io/=0) exit
  NOA=NOA+1 !number of sire
enddo
print *, "number of animal is ", NOA
rewind(98)
allocate(NRM(NOA,NOA),ID(NOA))
NRM=0
do ANO1=1, NOA
  if(CNOA.gt.4) then
    read(98,*) AID,SID,DID,DCV(:),ID(ANO1)
  elseif(CNOA.eq.4) then
    read(98,*) AID,SID,DID,ID(ANO1)
  else
    print *, "reading error - pedigree file"
    stop
  endif
  ! both known parents
  if (SID.gt.0.and.DID.gt.0) then
    do ANO2=1, AID-1
      NRM(ANO2,AID) = 0.5 * ( NRM(ANO2,SID) + NRM(ANO2,DID) )
      NRM(AID,ANO2) = NRM(ANO2,AID)
    end do    
      NRM(AID,AID) = 1.0 + 0.5*NRM(SID,DID)  
  else
    ! one known parent
    if (SID+DID.gt.0.and.SID*DID.eq.0) then
      do ANO2=1,AID-1
        NRM(ANO2,AID) = 0.5 * (NRM(ANO2,SID+DID))
        NRM(AID,ANO2) = NRM(ANO2,AID)
      end do
    end if
      NRM(AID,AID) = 1.0
  end if
end do
close(98)
do ANO1=1,NOA
!  write(96,*) ID(ANO1),NRM(ANO1,:)
  do ANO2=ANO1,NOA
  if(ANO1.eq.ANO2) write(94,*) ID(ANO1),NRM(ANO1,ANO2)-1
  enddo
enddo
NOA_S=0
do 
  read(97,*,iostat=io) 
  if(io/=0) exit
  NOA_S=NOA_S+1 !number of sire
enddo
NOA_S=NOA_S-1
print *, "number of selected animal is ", NOA_S
rewind(97)
read(97,*,iostat=io) 
allocate(NRM_S(NOA_S,NOA_S),ID_S(NOA_S),LOSA(NOA_S))
NRM_S=0
LOSA=0
do ANO_S1=1, NOA_S
  read(97,*) ID_S(ANO_S1)
  do ANO1=1,NOA
    if(ID_S(ANO_S1).eq.ID(ANO1).and.LOSA(ANO_S1).ne.0) print *, "Animal ",ID_S(ANO_S1)," is dupliccated" 
    if(ID_S(ANO_S1).eq.ID(ANO1)) then 
      LOSA(ANO_S1)=ANO1
      write(92,*) LOSA(ANO_S1)
    endif 
  enddo
  if(LOSA(ANO_S1).eq.0) print *, "Animal ",ID_S(ANO_S1)," is not in the pedigreefile"
enddo
do ANO_S1=1, NOA_S
  do ANO_S2=1, NOA_S
    NRM_S(ANO_S1,ANO_S2)=NRM(LOSA(ANO_S1),LOSA(ANO_S2))
    if(ANO_S2.ge.ANO_S1) write(93,fmt='(2i10,f15.10,2a22)') LOSA(ANO_S1),&
                         &LOSA(ANO_S2),NRM_S(ANO_S1,ANO_S2),ID_S(ANO_S1),ID_S(ANO_S2)
  enddo
  write(95,*) ID_S(ANO_S1),NRM_S(ANO_S1,:)
enddo
end program MakingNRM