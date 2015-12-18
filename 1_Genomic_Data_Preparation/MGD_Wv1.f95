! Name : Merging Genomic Data
! Developed by : Joonho Lee - PhD working in Hankyong National Univ.(Korea) zoonolee@gmail.com, +82-10-3408-2895
! Version : Windows 1.0 (2013-12-16)
! Contents : 		1. Merging Genomic Data

program MGD
 integer::NOF,FNO,TNOA=0,SANO,EANO,NOM,MNO,CNOFGT,CNOMKN,CNOCHN,CNOMKP,NOCNC,CCNO  !Numbers
 integer,allocatable::NOABF(:),CHNO(:),MKP(:)
 character(len=2),allocatable::GT(:)
 character(len=8)::OFT="_DEF0000"
 character(len=10)::SNOAFWF
 character(len=30),allocatable::ID(:)
 character(len=30)::WF_ID,WF_GT
 character(len=50)::GDFN
 character(len=50),allocatable::IOM(:),MKNA(:),CCN(:,:)
 character(len=200)::GDD,MDD,ODD
 character(len=230)::MDFNIP,ODFNIP,PFIP
 character(len=250)::GDFNIP

!reading parameters
print *, "***********************************************************************************************************************" 
print *, "MGD_Wv1 : Merging Genomic Data Windows version 4.0"
print *, "Developed by Joonho Lee : zoonolee@gmail.com, +82-10-3408-2895"
print *, "Data format : genotype data file - row is markers column is animal"
print *, "Data format : marker data file - (2nd row) ALGA0000009	1	52297	AG	GG	GG	... (Genotypes)"
print *, "Input data must be sorted by chromosome No.(1st) and Marker position(2nd)"
print *, "Sex chromosome and markers of unknown chromosome must be numbered more then the biggest chromosome number"
print *, "Length of Output file tag must be 8 characters"
print *, "Animal ID should be less than 30 characters"
print *, "Please enter the name of parameter file including directory(ex) ./parGRDH"
print *, "***********************************************************************************************************************"
read *, PFIP
PFIP=adjustl(PFIP)
open (99, file=PFIP) 
read (99,*) OFT ! output file tag
read (99,*) NOF ! number of files
read (99,*) CNOFGT ! column number of first genotype 
read (99,*) CNOMKN ! column number of marker name 
read (99,*) CNOCHN ! column number of chromosome number 
read (99,*) CNOMKP ! column number of marker position
read (99,*) NOCNC  ! number of chromosome number changes
allocate(CCN(NOCNC,2))
do CCNO=1,NOCNC
  read (99,*) CCN(CCNO,1),CCN(CCNO,2) 
  print *, CCN(CCNO,1),"chromosome will be numbered as ", CCN(CCNO,2)
enddo
allocate(NOABF(NOF))
read (99,*) GDD ! genomic data directory
GDD=adjustr(GDD)
do FNO=1,NOF
  read (99,*) GDFN, NOABF(FNO) !reading genomic data file name and number of animals by file
  GDFNIP=GDD//GDFN
  GDFNIP=adjustl(GDFNIP)
  open(FNO+20,file=GDFNIP,status='old')
  print *, GDFN,"file is opened"
  TNOA=TNOA+NOABF(FNO)
enddo
read (99,*) ODD ! genomic data directory
close(99)

ODD=adjustr(ODD)
ODFNIP=ODD//"GT"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(98,file=ODFNIP,status='unknown')

NOM=0
do 
  read(21,*,iostat=io)
  if(io/=0) exit
  NOM=NOM+1
enddo
rewind(unit=21)
NOM=NOM-1
print *, "number of markers : ", NOM
allocate(IOM(CNOFGT-1))
allocate(ID(TNOA),MKNA(NOF),CHNO(NOF),MKP(NOF),GT(TNOA))
write(SNOAFWF,fmt='(i10)') TNOA
SNOAFWF=adjustr(SNOAFWF)
WF_ID='(a50, a4, a15, '//SNOAFWF//'a30)'
WF_GT='(a50, a4, a15, '//SNOAFWF//'a3 )'

SANO=1
EANO=0
do FNO=1,NOF
  EANO=EANO+NOABF(FNO)
  read(FNO+20,*) IOM(:), ID(SANO:EANO)
  SANO=SANO+NOABF(FNO)
enddo
write(98, fmt=WF_ID) "Name","Chr","Position ",ID(:)

do MNO=1,NOM
  SANO=1
  EANO=0
  do FNO=1,NOF
    EANO=EANO+NOABF(FNO)
    read(FNO+20,*) IOM(:), GT(SANO:EANO)
    do CCNO=1,NOCNC
      if(IOM(CNOCHN).eq.CCN(CCNO,1)) IOM(CNOCHN)=CCN(CCNO,2)
    enddo
    read(IOM(CNOMKN),fmt='(a50)') MKNA(FNO)
    read(IOM(CNOCHN),fmt='(i15)') CHNO(FNO)
    read(IOM(CNOMKP),fmt='(i15)') MKP(FNO)
    if(FNO.ne.1.and.MKNA(FNO).ne.MKNA(FNO-1)) print *,MNO,"th marker name is unmatched between file",FNO-1,"and",FNO
    if(FNO.ne.1.and.CHNO(FNO).ne.CHNO(FNO-1)) print *,MNO,"th marker chromosome is unmatched between file",FNO-1,"and",FNO
    if(FNO.ne.1.and.MKP(FNO).ne.MKP(FNO-1)) print *, MNO,"th marker position is unmatched between file",FNO-1,"and",FNO
    SANO=SANO+NOABF(FNO)
  enddo
  write(98, fmt=WF_GT) IOM(CNOMKN),IOM(CNOCHN),IOM(CNOMKP),GT(:)
enddo
print *, "sorting command is"
print *, "awk '{if(NR>1) print | ",'"sort -n -k 2 -k 3"',"; else print}' ","GT"//OFT//".out"," > ","GT"//OFT//"_awksorted.out"
end program MGD



