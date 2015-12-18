! Name : Reformation of Imputation Results & Basic Statistics Analysis 
! Developed by : Joonho Lee - PhD working in Hankyong National Univ.(Korea) zoonolee@gmail.com, +82-10-3408-2895
! Version : Linux 2.0 (2012-10-25)
! Contents : 		1. Reformation of Imputation Results phenotype(name, records))
!		 	2. Calculation of gene & genotype frequency
!			4. Statistics of data & population genetics 

program RIRnBSA 
integer::i,j,k,iz,i1,i2,p,q,r,s !iostate & do_loop
integer::NOAC=0,SMNBC=1,NOA=0,ANO=0,NOM=0,MNO=0 !number of chromosomes, animals and markers & chromosome, animal and marekr number
integer::BC(4)
integer::NOMWMAF1=0,NOMWMAF2=0,NOMWMAF3=0,NOMWMAF4=0,NOMWMAF5=0,CMP=0,EXMP=0,CN=0,EXCN=999, SODBAM=0,NOMBCN=0 ! for statistics
integer::NODMBOFOMGT=0,NODMBOAH=0,NODMBOFOMIG=0,NODMBOCSOHWT=0,NODMBONP=0,NODMBOAHE=0,EXCNO=0 
integer,allocatable::GN(:,:),NOSMBC(:),COGT(:,:),COG(:,:) !data input & missing data output
real::FOMAH,FOMIH,FOHE,FOMAG,FOMIG,SOPIC=0,AOPIC=0,PV !frequencies
real::HWD=0,ENO0=0,ENO1=0,ENO2=0,CSOHWT=0,EH=0,PIC=0,TMGR=0,ADOAMBCN=0,AHWD=0,AHE=0,AEH=0,SOHWD=0,SOHE=0,SOEH=0 !statistics for population genetics
real::SOIBC=0,AOIBC=0,IBC=0 !statistics for population genetics
real,allocatable::FOGT(:,:),FOG(:,:)
character(len=1)::BA(4)=(/'A','T','G','C'/)
character(len=2)::CHNO
character(len=8)::OFT
character(len=24)::GDFN='_SSL049301_genotypes.out'
character(len=300)::Dummy
character(len=200)::GDD,ODD
character(len=30)::WF,ODFN
character(len=15)::MWF
character(len=230)::GDFNIP,ODFNIP,PF
character(len=1),allocatable::GDS(:,:),JMIMA(:,:)

!reading parameters
print *, "***********************************************************************************************************************" 
print *, "RIRnBSA_Lv2 : Reformation of Imputation Results & Basic Statistics Analysis Linux version 2.0"
print *, "Developed by Joonho Lee : zoonolee@gmail.com, +82-10-3408-2895"
print *, "Data form is the results of fastPHASE"
print *, "Please enter the name of parameter file including path(ex) ./parGRDH"
print *, "***********************************************************************************************************************"
read *, PF
PF=adjustl(PF)
open (99, file=PF) 
read (99,*) OFT
read (99,*) NOA
read (99,*) GDD
read (99,*) ODD
read (99,*) PV
read (99,*) NOAC
allocate(NOSMBC(NOAC))
do i=1,NOAC
  read (99,*) NOSMBC(i)
  NOM=NOM+NOSMBC(i)
enddo
close(99)

print *, "file tag", OFT
print *, "number of animals", NOA
print *, "number of markers", NOM
print *, "directory of imputed genomic data file", GDD
print *, "directory of output file", ODD
print *, "number of autosomal chromosomes : ", NOAC
do p=1,NOAC
  print *, "Number of markers on chromosome No.",p," : ", NOSMBC(p)
enddo
GDD=adjustr(GDD)
allocate(GDS(NOA*2,NOM))
allocate(GN(NOA,NOM))
allocate(FOGT(NOM,3),FOG(NOM,2))
allocate(COGT(NOM,3),COG(NOM,2))
allocate(JMIMA(NOM,2))
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

do q=1, NOAC
  iz=ichar('0')
  i1=floor(q/10.)+iz
  i2=mod(q,10)+iz
  CHNO(1:1)=achar(i1)
  CHNO(2:2)=achar(i2)
  GDFN(1:8)=OFT(1:8)
  GDFN(9:10)=CHNO(1:2)
  GDFNIP=GDD//GDFN
  GDFNIP=adjustl(GDFNIP)
  open(q,file=GDFNIP,status='old')
  print*,"reading",GDFNIP
  do r=1, 21
    read(q,fmt='(a300)',iostat=io) Dummy
    if(io/=0) exit
  enddo
  do k=1, NOA
    read(q,fmt='(a300)',iostat=io) Dummy
    print*,Dummy
    read(q,*,iostat=io) GDS(k*2-1,SMNBC:(SMNBC+NOSMBC(q)-1))
    if(io/=0) exit
    read(q,*,iostat=io) GDS(k*2,SMNBC:(SMNBC+NOSMBC(q)-1))
  enddo
  SMNBC=SMNBC+NOSMBC(q)
  print*,"Reading imputed genotype file",q,"finished"
enddo

SMNBC=1
do i=1, NOAC
  do j=1, NOA
    do k=1, NOSMBC(i)

  do j=1,NOA
  if(GTS(j).EQ.GT(11)) GTS(j)="??"
    do k=1,2
      select case (GTS(j)(k:k))
        case("A")
          BC(1)=BC(1)+1
        case("T")
          BC(2)=BC(2)+1
        case("G")
          BC(3)=BC(3)+1
        case("C")
          BC(4)=BC(4)+1
        case("?")
          BC(5)=BC(5)+1
      end select
    enddo

      MNO=MNO+1
      if(j.eq.1) JMIMA(MNO,1)=GDS(j*2-1,SMNBC+k-1)
      if(GDS(j*2-1,SMNBC+k-1).ne.GDS(j*2,SMNBC+k-1)) then
        if(GDS(j*2-1,SMNBC+k-1).eq.JMIMA(MNO,1)) then
          JMIMA(MNO,2)=GDS(j*2,SMNBC+k-1)
        else
          JMIMA(MNO,2)=GDS(j*2-1,SMNBC+k-1)
        endif
      endif
      if(GDS(j*2-1,SMNBC+k-1).eq.GDS(j*2,SMNBC+k-1).and.GDS(j*2-1,SMNBC+k-1).ne.JMIMA(MNO,1)) then
        GN(j,MNO)=0
        COGT(MNO,1)=COGT(MNO,1)+1
        COG(MNO,1)=COG(MNO,1)+2
      elseif(GDS(j*2-1,SMNBC+k-1).ne.GDS(j*2,SMNBC+k-1)) then
        GN(j,MNO)=1
        COGT(MNO,2)=COGT(MNO,2)+1
        COG(MNO,1)=COG(MNO,1)+1
        COG(MNO,2)=COG(MNO,2)+1
      elseif(GDS(j*2-1,SMNBC+k-1).eq.GDS(j*2,SMNBC+k-1).and.GDS(j*2-1,SMNBC+k-1).eq.JMIMA(MNO,1)) then
        GN(j,MNO)=2
        COGT(k,3)=COGT(k,3)+1
        COG(k,2)=COG(k,2)+2
      endif
    enddo
  enddo
  SMNBC=SMNBC+NOSMBC(i)
enddo
FOGT=real(COGT)/real(NOA)
FOG=real(COG)/real(NOA)
write(MWF,fmt='(a1,i10,a3)') '(',MNO,'a2)'

do i=1,NOA
  write(98,fmt=MWF) GN(i,:)
enddo
MNO=0
do i=1,NOAC
  do j=1,NOSMBC(i)
    MNO=MNO+1
    write(97,fmt='(2f15.10)') FOG(MNO,1), FOG(MNO,2)
    write(96,fmt='(3f15.10)') FOGT(MNO,1), FOGT(MNO,2), FOGT(MNO,3)
    write(95,fmt='(3a3)') JMIMA(MNO,2)//JMIMA(MNO,2), JMIMA(MNO,2)//JMIMA(MNO,1), JMIMA(MNO,1)//JMIMA(MNO,1)
    write(94,fmt='(2a3)') JMIMA(MNO,2), JMIMA(MNO,1)
    ENO0=(FOG(MNO,1)**2)*NOA
    ENO1=(2*FOG(MNO,1)*FOG(MNO,2))*NOA
    ENO2=(FOG(MNO,2)**2)*NOA
    CSOHWT=((COGT(MNO,1)-ENO0)**2/ENO0)+((COGT(MNO,2)-ENO1)**2/ENO1)+((COGT(MNO,3)-ENO2)**2/ENO2)
    EH=1-(FOG(MNO,1)**2+FOG(MNO,2)**2)
    PIC=1-(FOG(MNO,1)**2+FOG(MNO,2)**2)-2*((FOG(MNO,1)**2)*(FOG(MNO,2)**2))
    IBC=1-(real(COGT(MNO,2))/ENO1)
    SOHE=SOHE+EH
    SOPIC=SOPIC+PIC
    SOIBC=SOIBC+IBC
    if(CSOHWT>PV) NODMBOCSOHWTBC=NODMBOCSOHWTBC+1
    write(93,fmt='(i5,2f10.8,3f10.5)') i,FOGT(MNO,2),EH,CSOHWT,PIC,IBC 
  enddo
  AEH=real(SOHE)/real(NOSMBC(i))
  AOPIC=real(SOPIC)/real(NOSMBC(i))
  AOIBC=real(SOIBC)/real(NOSMBC(i))
  write(92,fmt='(a12,i5,2f10.8,3f10.5)') "chromosome", i,AEH,AOPIC,AOIBC
  NODMBOCSOHWTBC=0
  SOHE=0
  SOPIC=0
  SOIBC=0
enddo
end program RIRnBSA

