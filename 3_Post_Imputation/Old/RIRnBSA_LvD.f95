! Name : Reformation of Imputation Results & Basic Statistics Analysis 
! Developed by : Joonho Lee - PhD working in Hankyong National Univ.(Korea) zoonolee@gmail.com, +82-10-3408-2895
! Version : Linux 1.0 (2011-10-25)
! Contents : 		1. Reformation of Imputation Results phenotype(name, records))
!		 	2. Calculation of gene & genotype frequency
!			4. Statistics of data & population genetics 

program RIRnBSA 
integer::i,j,k,l,m,n,o,p,q,r,s,iz,i1,i2,AMNO=0,MNO=0,SMNBC=1,EMNBC,CMN,io 
integer::NOCEXC=0,ANO=0 !number of chromosomes & chromosome number for X, Y, Unknown
integer::NOMWMAF1=0,NOMWMAF2=0,NOMWMAF3=0,NOMWMAF4=0,NOMWMAF5=0,CMP=0,EXMP=0,CN=0,EXCN=999, SODBAM=0,NOMBCN=0 ! for statistics
integer::NODMBOFOMGT=0,NODMBOAH=0,NODMBOFOMIG=0,NODMBOCSOHWT=0,NODMBONP=0,NODMBOAHE=0,EXCNO=0 
integer,allocatable::GN(:,:),NOMBC(:),COGT(:,:),COG(:,:) !data input & missing data output
real::FOMAH,FOMIH,FOHE,FOMAG,FOMIG,SOPIC=0,AOPIC=0,PV !frequencies
real::HWD=0,ENO0=0,ENO1=0,ENO2=0,CSOHWT=0,EH=0,PIC=0,TMGR=0,ADOAMBCN=0,AHWD=0,AHE=0,AEH=0,SOHWD=0,SOHE=0,SOEH=0 !statistics for population genetics
real::SOIBC=0,AOIBC=0,IBC=0 !statistics for population genetics
real,allocatable::FOGT(:,:),FOG(:,:)
character(len=1)::SNP
character(len=2)::CHNO
character(len=300)::GDD,ODD,DR,Dummy
character(len=30)::GDFN,WF,ODFN
character(len=15)::MWF
character(len=330)::GDFNIP,ODFNIP,PF
character(len=1),allocatable::GTDAEOM(:,:),JMIMA(:,:)

!reading parameters
print *, "***********************************************************************************************************************" 
print *, "RIRnBSA_Lv1 : Reformation of Imputation Results & Basic Statistics Analysis Linux version 3.0"
print *, "Developed by Joonho Lee : zoonolee@gmail.com, +82-10-3408-2895"
print *, "Data form is the results of fastPHASE"
print *, "Please enter the name of parameter file including path(ex) ./parGRDH"
print *, "***********************************************************************************************************************"
read *, PF
PF=adjustl(PF)
open (99, file=PF) 
read (99,*) ANO
read (99,*) GDD
read (99,*) GDFN
read (99,*) ODD
read (99,*) PV
read (99,*) NOCEXC
allocate(NOMBC(NOCEXC))
do i=1,NOCEXC
  read (99,*) NOMBC(i)
  MNO=MNO+NOMBC(i)
enddo
close(99)

print *, "number of markers", MNO
print *, "number of animals", ANO
print *, "directory of genomic data file", GDD
print *, "name of genomic data file?", GDFN
print *, "directory of output file", ODD
print *, "number of autosomal chromosomes : ", NOCEXC
do j=1,NOCEXC
  print *, "Number of markers on chromosome No.",j," : ", NOMBC(j)
enddo

GDD=adjustr(GDD)
allocate(GTDAEOM(ANO*2,MNO))
allocate(GN(ANO,MNO))
allocate(FOGT(MNO,3),FOG(MNO,2))
allocate(COGT(MNO,3),COG(MNO,2))
allocate(JMIMA(MNO,2))

do k=1, NOCEXC
  iz=ichar('0')
  i1=floor(k/10.)+iz
  i2=mod(k,10)+iz
  CHNO(1:1)=achar(i1)
  CHNO(2:2)=achar(i2)
  GDFN(1:2)=CHNO(1:2)
  GDFNIP=GDD//GDFN
  GDFNIP=adjustl(GDFNIP)
  open(66,file=GDFNIP,status='old')
  print*,"reading",GDFNIP
  do l=1,21
    read(66,fmt='(a300)') Dummy
  enddo
  EMNBC=SMNBC+NOMBC(k)-1
  do m=1, ANO
    read(66,fmt='(a300)') Dummy
    read(66,*,iostat=io) GTDAEOM(m*2-1,SMNBC:EMNBC)
    if(m.eq.1) JMIMA(SMNBC:EMNBC,1)=GTDAEOM(m*2-1,SMNBC:EMNBC)
    read(66,*,iostat=io) GTDAEOM(m*2,SMNBC:EMNBC)
  enddo
  SMNBC=SMNBC+NOMBC(k)
  print*,"Reading file",k,"finished"
  close(66)
enddo
  print*,"Reading file end"
!  print*,JMIMA
COGT=0
COG=0

SMNBC=1
do i=1, NOCEXC
  do j=1, ANO
    do k=1, NOMBC(i)
!  print*,i,j,k
      CMN=SMNBC+k-1
      if(j.eq.1) JMIMA(CMN,1)=GTDAEOM(j*2-1,CMN)
!      if(GTDAEOM(j*2-1,CMN).ne.JMIMA(CMN,1)) GTDAEOM(j*2-1,CMN) =JMIMA(CMN,2)
!      if(GTDAEOM(j*2-1,CMN).ne.GTDAEOM(j*2,CMN)) then
!        if(GTDAEOM(j*2-1,CMN).eq.JMIMA(CMN,1)) then
!          JMIMA(CMN,2)=GTDAEOM(j*2,CMN)
!        else
!          JMIMA(CMN,2)=GTDAEOM(j*2-1,CMN)
!        endif
!      endif
      if(GTDAEOM(j*2-1,CMN).eq.GTDAEOM(j*2,CMN).and.GTDAEOM(j*2-1,CMN).eq.JMIMA(CMN,1)) then
        GN(j,CMN)=2
        COGT(CMN,3)=COGT(CMN,3)+1
        COG(CMN,2)=COG(CMN,2)+2
      elseif(GTDAEOM(j*2-1,CMN).ne.GTDAEOM(j*2,CMN)) then
        GN(j,CMN)=1
        COGT(CMN,2)=COGT(CMN,2)+1
        COG(CMN,1)=COG(CMN,1)+1
        COG(CMN,2)=COG(CMN,2)+1
      else
        GN(j,CMN)=0
        COGT(CMN,1)=COGT(CMN,1)+1
        COG(CMN,1)=COG(CMN,1)+2
      endif
    enddo
  enddo
  SMNBC=SMNBC+NOMBC(i)
enddo

ODD=adjustr(ODD)
ODFNIP=ODD//"Genotype012.out"
ODFNIP=adjustl(ODFNIP)
open(71,file=ODFNIP,status='unknown')
ODFNIP=ODD//"Genefrequency.out"
ODFNIP=adjustl(ODFNIP)
open(72,file=ODFNIP,status='unknown')
write(72,fmt='(2a15)') "G0","G1"
ODFNIP=ODD//"Genotypefrequency.out"
ODFNIP=adjustl(ODFNIP)
open(73,file=ODFNIP,status='unknown')
write(73,fmt='(3a15)') "GT0","GT1","GT2"
ODFNIP=ODD//"GT_maching_012.out"
ODFNIP=adjustl(ODFNIP)
open(74,file=ODFNIP,status='unknown')
write(74,fmt='(3a10)') "GT0","GT1","GT2"
ODFNIP=ODD//"G_maching_01.out"
ODFNIP=adjustl(ODFNIP)
open(75,file=ODFNIP,status='unknown')
write(75,fmt='(2a10)') "G0","G1"
ODFNIP=ODD//"Pop_Gen_Stat.out"
ODFNIP=adjustl(ODFNIP)
open(76,file=ODFNIP,status='unknown')
write(76,fmt='(a5,5a15)')  "CNO","ObservedH","ExpectedH","HWT_Chi2","PIC","IBC"
ODFNIP=ODD//"Pop_Gen_Stat_CH.out"
ODFNIP=adjustl(ODFNIP)
open(77,file=ODFNIP,status='unknown')
write(77,fmt='(a17,3a15)') "NO","AOEH","AOPIC","AOIBC"

write(MWF,fmt='(a1,i10,a3)') '(',CMN,'a2)'

do i=1,ANO
  write(71,fmt='(100000a2)') GN(i,1:NOM)
enddo
AMNO=0
do i=1,NOCEXC
  do j=1,NOMBC(i)
    AMNO=AMNO+1
    FOGT(AMNO,1)=real(COGT(AMNO,1))/real(ANO)
    FOGT(AMNO,2)=real(COGT(AMNO,2))/real(ANO)
    FOGT(AMNO,3)=real(COGT(AMNO,3))/real(ANO)
    FOG(AMNO,1)=real(COG(AMNO,1))/real(ANO*2)
    FOG(AMNO,2)=real(COG(AMNO,2))/real(ANO*2)
    write(72,fmt='(2f15.10)') FOG(AMNO,1), FOG(AMNO,2)
    write(73,fmt='(3f15.10)') FOGT(AMNO,1), FOGT(AMNO,2), FOGT(AMNO,3)
    write(74,fmt='(3a10)') JMIMA(AMNO,2)//JMIMA(AMNO,2), JMIMA(AMNO,2)//JMIMA(AMNO,1), JMIMA(AMNO,1)//JMIMA(AMNO,1)
    write(75,fmt='(2a10)') JMIMA(AMNO,2), JMIMA(AMNO,1)
    ENO0=(FOG(AMNO,1)**2)*ANO
    ENO1=(2*FOG(AMNO,1)*FOG(AMNO,2))*ANO
    ENO2=(FOG(AMNO,2)**2)*ANO
    CSOHWT=((COGT(AMNO,1)-ENO0)**2/ENO0)+((COGT(AMNO,2)-ENO1)**2/ENO1)+((COGT(AMNO,3)-ENO2)**2/ENO2)
    EH=1-(FOG(AMNO,1)**2+FOG(AMNO,2)**2)
    PIC=1-(FOG(AMNO,1)**2+FOG(AMNO,2)**2)-2*((FOG(AMNO,1)**2)*(FOG(AMNO,2)**2))
    IBC=1-(real(COGT(AMNO,2))/ENO1)
    SOHE=SOHE+EH
    SOPIC=SOPIC+PIC
    SOIBC=SOIBC+IBC
    if(CSOHWT>PV) NODMBOCSOHWTBC=NODMBOCSOHWTBC+1
    write(76,fmt='(i5,2f15.8,3f15.5)') i,FOGT(AMNO,2),EH,CSOHWT,PIC,IBC 
  enddo
  AEH=real(SOHE)/real(NOMBC(i))
  AOPIC=real(SOPIC)/real(NOMBC(i))
  AOIBC=real(SOIBC)/real(NOMBC(i))
  write(77,fmt='(a12,i5,2f15.8,3f15.5)') "chromosome", i,AEH,AOPIC,AOIBC
  NODMBOCSOHWTBC=0
  SOHE=0
  SOPIC=0
  SOIBC=0
enddo
end program RIRnBSA

