! Name : Merging 50k Genomic Data and 777k Genomic Data
! Developed by : Joonho Lee - PhD working in Hankyong National Univ.(Korea) zoonolee@gmail.com, +82-10-3408-2895
! Version : Linux 1.0 (2012-06-12)
! Contents : 		1. Merging 50k Genomic Data and 777k Genomic Data

program M57
 integer::NOA5,NOA7,NOA57,NOM5,NOM7,NOM57,MNO5,MNO7,MNO57,NOOM  !Numbers
 integer,allocatable::CHNO5(:),MKPO5(:),CHNO7(:),MKPO7(:)
 character(len=2),allocatable::GT5(:,:),GT7(:,:)
 character(len=8)::OFT="_DEF0000"
 character(len=10)::SNOAFWF
 character(len=20),allocatable::ID5(:),ID7(:)
 character(len=30)::GDFN5,GDFN7,WF_ID,WF_GT
 character(len=50)::Dummy
 character(len=50),allocatable::MKNA5(:),MKNA7(:)
 character(len=200)::GDD5,GDD7,ODD
 character(len=230)::GDFNIP5,GDFNIP7,ODFNIP,PFIP

!reading parameters
print *, "***********************************************************************************************************************" 
print *, "M57_Lv1 : Merging 50k Genomic Data and 777k Genomic Data Linux version 1.0"
print *, "Developed by Joonho Lee : zoonolee@gmail.com, +82-10-3408-2895"
print *, "Data format : genotype data file - row is markers column is animal"
print *, "Data format : marker data file - (2nd row) ALGA0000009	1	52297	AG	GG	GG	... (Genotypes)"
print *, "Input data must be sorted by chromosome No.(1st) and Marker position(2nd)"
print *, "Sex chromosome and markers of unknown chromosome must be numbered more then the biggest chromosome number"
print *, "Length of Output file tag must be 8 characters"
print *, "Please enter the name of parameter file including directory(ex) ./parGRDH"
print *, "***********************************************************************************************************************"
read *, PFIP
PFIP=adjustl(PFIP)
open (99, file=PFIP) 
read (99,*) OFT ! output file tag
read (99,*) NOA5 ! number of animals on 50k genomic data file
read (99,*) NOA7 ! number of animals on 777k genomic data file
read (99,*) NOM5 ! number of markers on 50k genomic data file
read (99,*) NOM7 ! number of markers on 777k genomic data file
read (99,*) GDD5 ! directory of 50k genomic data file
read (99,*) GDFN5 ! 50k genomic data file name
read (99,*) GDD7 ! directory of 777k genomic data file
read (99,*) GDFN7 ! 777k genomic data file name
read (99,*) ODD ! output data directory
close(99)

print *, "output file tag : ", OFT 
print *, "number of animals on 50k genomic data file : ", NOA5 
print *, "number of animals on 777k genomic data file : ", NOA7
print *, "number of markers on 50k genomic data file : ", NOM5
print *, "number of markers on 777k genomic data file : ", NOM7
print *, "directory of 50k genomic data file : ", GDD5
print *, "50k genomic data file name : ", GDFN5
print *, "directory of 777k genomic data file : ", GDD7 
print *, "777k genomic data file name : ", GDFN7
print *, "output data directory : ", ODD

GDD5=adjustr(GDD5)
GDD7=adjustr(GDD7)
GDFNIP5=GDD5//GDFN5
GDFNIP7=GDD7//GDFN7
GDFNIP5=adjustl(GDFNIP5)
GDFNIP7=adjustl(GDFNIP7)
open(98,file=GDFNIP5,status='old')
print *, "50k genomic file - ", GDFN5," is opened"
open(97,file=GDFNIP7,status='old')
print *, "777k genomic file - ", GDFN7," is opened"
ODD=adjustr(ODD)
ODFNIP=ODD//"GT"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(96,file=ODFNIP,status='unknown')
ODFNIP=ODD//"odd-matching"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(95,file=ODFNIP,status='unknown')

allocate(GT5(NOM5,NOA5),GT7(NOM7,NOA7))
allocate(ID5(NOA5),ID7(NOA7))
allocate(MKNA5(NOM5),MKNA7(NOM7),CHNO5(NOM5),CHNO7(NOM7),MKPO5(NOM5),MKPO7(NOM7))

NOA57=NOA5+NOA7
write(SNOAFWF,fmt='(i10)') NOA57
SNOAFWF=adjustr(SNOAFWF)
WF_ID='(a50, a4, a15, '//SNOAFWF//'a20)'
WF_GT='(a50, i4, i15, '//SNOAFWF//'a3 )'

read(98,*) Dummy,Dummy,Dummy,ID5(:)
read(97,*) Dummy,Dummy,Dummy,ID7(:)
write(96, fmt=WF_ID) "Name","Chr","Position ",ID5(:),ID7(:)

do MNO5=1, NOM5
  read(98,*) MKNA5(MNO5),CHNO5(MNO5),MKPO5(MNO5),GT5(MNO5,:)
enddo
print *, "reading 50k genomic file is finished "
do MNO7=1, NOM7
  read(97,*) MKNA7(MNO7),CHNO7(MNO7),MKPO7(MNO7),GT7(MNO7,:)
enddo
print *, "reading 777k genomic file is finished "

NOM57=0
NOOM=0
do MNO5=1, NOM5
  do MNO7=1, NOM7
    if(MKNA5(MNO5).eq.MKNA7(MNO7)) then
       NOM57=NOM57+1
       write(96, fmt=WF_GT) MKNA7(MNO7),CHNO7(MNO7),MKPO7(MNO7),GT5(MNO5,:),GT7(MNO7,:)
       if(CHNO5(MNO5).ne.CHNO7(MNO7).or.MKPO5(MNO5).ne.MKPO7(MNO7)) then
         NOOM=NOOM+1
         write(95, *) MKNA5(MNO5),CHNO5(MNO5),MKPO5(MNO5),GT5(MNO5,1:5),"  ",MKNA7(MNO7),CHNO7(MNO7),MKPO7(MNO7),GT7(MNO7,1:5)
       endif
    endif
  enddo
enddo
print *, "number of matched markers are ",NOM57
print *, "number of odd_matched markers are ",NOOM
end program M57




