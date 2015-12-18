! Name : Reformation of Imputation Results & Basic Statistics Analysis 
! Developed by : Joonho Lee - PhD working in Hankyong National Univ.(Korea) zoonolee@gmail.com, +82-10-3408-2895
! Version : Linux 1.0 (2011-10-25)
! Contents : 		1. Reformation of Imputation Results phenotype(name, records))
!		 	2. Calculation of gene & genotype frequency
!			4. Statistics of data & population genetics 

program RIRnBSA 
integer::i,j,k,iz,i1,i2,AMNO=0,MAXNOMBC=0,MNO=0,SMNBC=1,p,q,r,s !iostate & do_loop
integer::NOCEXC=0,ANO=0 !number of chromosomes & chromosome number for X, Y, Unknown
integer::NOMWMAF1=0,NOMWMAF2=0,NOMWMAF3=0,NOMWMAF4=0,NOMWMAF5=0,CMP=0,EXMP=0,CN=0,EXCN=999, SODBAM=0,NOMBCN=0 ! for statistics
integer::NODMBOFOMGT=0,NODMBOAH=0,NODMBOFOMIG=0,NODMBOCSOHWT=0,NODMBONP=0,NODMBOAHE=0,EXCNO=0 
integer,allocatable::GN(:,:),NOMBC(:),COGT(:,:),COG(:,:) !data input & missing data output
real::FOMAH,FOMIH,FOHE,FOMAG,FOMIG,SOPIC=0,AOPIC=0,PV !frequencies
real::HWD=0,ENO0=0,ENO1=0,ENO2=0,CSOHWT=0,EH=0,PIC=0,TMGR=0,ADOAMBCN=0,AHWD=0,AHE=0,AEH=0,SOHWD=0,SOHE=0,SOEH=0 !statistics for population genetics
real::SOIBC=0,AOIBC=0,IBC=0 !statistics for population genetics
real,allocatable::FOGT(:,:),FOG(:,:)
character(len=2)::CHNO
character(len=300)::GDD,ODD,DR,Dummy
character(len=30)::GDFN,WF,ODFN,DR1,DR2,DR3
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

print *, "number of animals", ANO
print *, "directory of genomic data file", GDD
print *, "name of genomic data file?", GDFN
print *, "directory of output file", ODD
print *, "number of autosomal chromosomes : ", NOCEXC
do p=1,NOCEXC
  print *, "Number of markers on chromosome No.",p," : ", NOMBC(p)
enddo
GDD=adjustr(GDD)
allocate(GTDAEOM(ANO*2,MNO))
allocate(GN(ANO,MNO))
allocate(FOGT(MNO,3),FOG(MNO,2))
allocate(COGT(MNO,3),COG(MNO,2))
allocate(JMIMA(MNO,2))
ODD=adjustr(ODD)
ODFNIP=ODD//"Genotype012.out"
ODFNIP=adjustl(ODFNIP)
open(98,file=ODFNIP,status='unknown')
ODFNIP=ODD//"Genefrequency.out"
ODFNIP=adjustl(ODFNIP)
open(97,file=ODFNIP,status='unknown')
write(97,fmt='(2a15)') "0","1"
ODFNIP=ODD//"Genotypefrequency.out"
ODFNIP=adjustl(ODFNIP)
open(96,file=ODFNIP,status='unknown')
write(96,fmt='(3a15)') "0","1","2"
ODFNIP=ODD//"GT_maching_012.out"
ODFNIP=adjustl(ODFNIP)
open(95,file=ODFNIP,status='unknown')
write(95,fmt='(3a3)') "0","1","2"
ODFNIP=ODD//"G_maching_01.out"
ODFNIP=adjustl(ODFNIP)
open(94,file=ODFNIP,status='unknown')
write(94,fmt='(2a3)') "0","1"
ODFNIP=ODD//"Pop_Gen_Stat.out"
ODFNIP=adjustl(ODFNIP)
open(93,file=ODFNIP,status='unknown')
write(93,fmt='(a5,5a10)')  "CHNO.","ObservedH","ExpectedH","HWT_Chi2","PIC","IBC"
ODFNIP=ODD//"Pop_Gen_Stat_CH.out"
ODFNIP=adjustl(ODFNIP)
open(92,file=ODFNIP,status='unknown')
write(92,fmt='(a17,3a10)') "NO","AOEH","AOPIC","AOIBC"

do q=1, NOCEXC
  iz=ichar('0')
  i1=floor(q/10.)+iz
  i2=mod(q,10)+iz
  CHNO(1:1)=achar(i1)
  CHNO(2:2)=achar(i2)
  GDFN(1:2)=CHNO(1:2)
  GDFNIP=GDD//GDFN
  GDFNIP=adjustl(GDFNIP)
  open(q,file=GDFNIP,status='old')
  print*,"reading",GDFNIP
    read(q,fmt='(a300)',iostat=io) Dummy
    read(q,fmt='(a300)',iostat=io) Dummy
    read(q,fmt='(a300)',iostat=io) Dummy
    read(q,fmt='(a300)',iostat=io) Dummy
    read(q,fmt='(a300)',iostat=io) Dummy
    read(q,fmt='(a300)',iostat=io) Dummy
    read(q,fmt='(a300)',iostat=io) Dummy
    read(q,fmt='(a300)',iostat=io) Dummy
    read(q,fmt='(a300)',iostat=io) Dummy
    read(q,fmt='(a300)',iostat=io) Dummy
    read(q,fmt='(a300)',iostat=io) Dummy
    read(q,fmt='(a300)',iostat=io) Dummy
    read(q,fmt='(a300)',iostat=io) Dummy
    read(q,fmt='(a300)',iostat=io) Dummy
    read(q,fmt='(a300)',iostat=io) Dummy
    read(q,fmt='(a300)',iostat=io) Dummy
    read(q,fmt='(a300)',iostat=io) Dummy
    read(q,fmt='(a300)',iostat=io) Dummy
    read(q,fmt='(a300)',iostat=io) Dummy
    read(q,fmt='(a300)',iostat=io) Dummy
    read(q,fmt='(a300)',iostat=io) Dummy
  do k=1, NOA
      
    read(q,fmt='(a300)',iostat=io) Dummy
    print*,Dummy
    read(q,*,iostat=io) GTDAEOM(k*2-1,SMNBC:(SMNBC+NOMBC(q)-1))
    if(io/=0) exit
    read(q,*,iostat=io) GTDAEOM(k*2,SMNBC:(SMNBC+NOMBC(q)-1))
  enddo
  SMNBC=SMNBC+NOMBC(q)
  print*,"Reading file",q,"finished"

enddo

SMNBC=1
do i=1, NOCEXC
  do j=1, ANO
    do k=1, NOMBC(i)
      AMNO=AMNO+1
      if(j.eq.1) JMIMA(AMNO,1)=GTDAEOM(j*2-1,SMNBC+k-1)
      if(GTDAEOM(j*2-1,SMNBC+k-1).ne.GTDAEOM(j*2,SMNBC+k-1)) then
        if(GTDAEOM(j*2-1,SMNBC+k-1).eq.JMIMA(AMNO,1)) then
          JMIMA(AMNO,2)=GTDAEOM(j*2,SMNBC+k-1)
        else
          JMIMA(AMNO,2)=GTDAEOM(j*2-1,SMNBC+k-1)
        endif
      endif
      if(GTDAEOM(j*2-1,SMNBC+k-1).eq.GTDAEOM(j*2,SMNBC+k-1).and.GTDAEOM(j*2-1,SMNBC+k-1).ne.JMIMA(AMNO,1)) then
        GN(j,AMNO)=0
        COGT(AMNO,1)=COGT(AMNO,1)+1
        COG(AMNO,1)=COG(AMNO,1)+2
      elseif(GTDAEOM(j*2-1,SMNBC+k-1).ne.GTDAEOM(j*2,SMNBC+k-1)) then
        GN(j,AMNO)=1
        COGT(AMNO,2)=COGT(AMNO,2)+1
        COG(AMNO,1)=COG(AMNO,1)+1
        COG(AMNO,2)=COG(AMNO,2)+1
      elseif(GTDAEOM(j*2-1,SMNBC+k-1).eq.GTDAEOM(j*2,SMNBC+k-1).and.GTDAEOM(j*2-1,SMNBC+k-1).eq.JMIMA(AMNO,1)) then
        GN(j,AMNO)=2
        COGT(k,3)=COGT(k,3)+1
        COG(k,2)=COG(k,2)+2
      endif
    enddo
  enddo
  SMNBC=SMNBC+NOMBC(i)
enddo
FOGT=real(COGT)/real(ANO)
FOG=real(COG)/real(ANO)
write(MWF,fmt='(a1,i10,a3)') '(',AMNO,'a2)'

do i=1,ANO
  write(98,fmt=MWF) GN(i,:)
enddo
AMNO=0
do i=1,NOCEXC
  do j=1,NOMBC(i)
    AMNO=AMNO+1
    write(97,fmt='(2f15.10)') FOG(AMNO,1), FOG(AMNO,2)
    write(96,fmt='(3f15.10)') FOGT(AMNO,1), FOGT(AMNO,2), FOGT(AMNO,3)
    write(95,fmt='(3a3)') JMIMA(AMNO,2)//JMIMA(AMNO,2), JMIMA(AMNO,2)//JMIMA(AMNO,1), JMIMA(AMNO,1)//JMIMA(AMNO,1)
    write(94,fmt='(2a3)') JMIMA(AMNO,2), JMIMA(AMNO,1)
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
    write(93,fmt='(i5,2f10.8,3f10.5)') i,FOGT(AMNO,2),EH,CSOHWT,PIC,IBC 
  enddo
  AEH=real(SOHE)/real(NOMBC(i))
  AOPIC=real(SOPIC)/real(NOMBC(i))
  AOIBC=real(SOIBC)/real(NOMBC(i))
  write(92,fmt='(a12,i5,2f10.8,3f10.5)') "chromosome", i,AEH,AOPIC,AOIBC
  NODMBOCSOHWTBC=0
  SOHE=0
  SOPIC=0
  SOIBC=0
enddo
end program RIRnBSA

