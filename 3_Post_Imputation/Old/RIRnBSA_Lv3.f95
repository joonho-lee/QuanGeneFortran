! Name : Reformation of Imputation Results & Basic Statistics Analysis 
! Developed by : Joonho Lee - PhD working in Hankyong National Univ.(Korea) zoonolee@gmail.com, +82-10-3408-2895
! Version : Linux 3.0 (2012-03-18)
! Contents : 		1. Reformation of Imputation Results phenotype(name, records))
!		 	2. Calculation of gene & genotype frequency
!			4. Statistics of data & population genetics 

program RIRnBSA 
integer::i,j,k,iz,i1,i2,p,q,r,s !iostate & do_loop
integer::NOAC=0,SMNBC=1,NOA=0,ANO=0,NOM=0,MNO=0 !number of chromosomes, animals and markers & chromosome, animal and marekr number
integer::BC(4)
integer::GC(10)
integer::NOMWMAF1=0,NOMWMAF2=0,NOMWMAF3=0,NOMWMAF4=0,NOMWMAF5=0,CMP=0,EXMP=0,CN=0,EXCN=999, SODBAM=0,NOMBCN=0 ! for statistics
integer::NODMBOFOMGT=0,NODMBOAH=0,NODMBOFOMIG=0,NODMBOCSOHWT=0,NODMBONP=0,NODMBOAHE=0,EXCNO=0 
integer,allocatable::G012DS(:,:),NOSMBC(:) !data input & missing data output
real::GF(5),GTF(11),FOMAH,FOMIH,FOHE,FOMAG,FOMIG,SOGF,SOGTF,FOMG,FOMGT,IOL(4),SOPIC=0,AOPIC=0,PV !frequencies
real::HWD=0,ENO0=0,ENO1=0,ENO2=0,CSOHWT=0,EH=0,PIC=0,TMGR=0,ADOAMBCN=0,AHWD=0,AHE=0,AEH=0,SOHWD=0,SOHE=0,SOEH=0 !statistics for population genetics
real::SOIBC=0,AOIBC=0,IBC=0 !statistics for population genetics
real::SOMAFBC=0,SOIBCBC=0,SOCSOHWTBC=0,SOPICBC=0,SOHEBC=0,SOEHBC=0,SOTMAF=0,SOTIBC=0,SOTCSOHWT=0,SOTPIC=0,SOTHE=0,SOTEH=0 !statistics for population genetics by chromosome
real::AOMAFBC=0,AOIBCBC=0,AOCSOHWTBC=0,AOPICBC=0,AOHEBC=0,AOEHBC=0,AOTMAF=0,AOTIBC=0,AOTCSOHWT=0,AOTPIC=0,AOTHE=0,AOTEH=0 !statistics for population genetics by chromosome
real::SDOMAFBC=0,SDOIBCBC=0,SDOCSOHWTBC=0,SDOPICBC=0,SDOHEBC=0,SDOEHBC=0
real::SDOTMAF=0,SDOTIBC=0,SDOTCSOHWT=0,SDOTPIC=0,SDOTHE=0,SDOTEH=0 !statistics for population genetics by chromosome
real::MINOMAFBC=0,MINOIBCBC=0,MINOCSOHWTBC=0,MINOPICBC=0,MINOHEBC=0,MINOEHBC=0
real::MINOTMAF=0,MINOTIBC=0,MINOTCSOHWT=0,MINOTPIC=0,MINOTHE=0,MINOTEH=0 !statistics for population genetics by chromosome
real::MAXOMAFBC=0,MAXOIBCBC=0,MAXOCSOHWTBC=0,MAXOPICBC=0,MAXOHEBC=0,MAXOEHBC=0
real::MAXOTMAF=0,MAXOTIBC=0,MAXOTCSOHWT=0,MAXOTPIC=0,MAXOTHE=0,MAXOTEH=0 !statistics for population genetics by chromosome


character(len=1)::BA(4)=(/'A','T','G','C'/)
character(len=2)::GT(10)=(/'AA','AT','AG','AC','TT','TG','TC','GG','GC','CC'/),MAH,MIH,HE,CHNO
character(len=8)::OFT
character(len=24)::GDFN='01_DEF0000_genotypes.out'
character(len=300)::Dummy
character(len=200)::GDD,ODD
character(len=30)::WF,ODFN
character(len=15)::MWF
character(len=230)::GDFNIP,ODFNIP,PF
character(len=1),allocatable::GDS(:,:)
character(len=2),allocatable::GTDS(:,:)

!reading parameters
print *, "***********************************************************************************************************************" 
print *, "RIRnBSA_Lv3 : Reformation of Imputation Results & Basic Statistics Analysis Linux version 3.0"
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
allocate(GTDS(NOM,NOA))
allocate(G012DS(NOM,NOA))
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
write(93,fmt='(a5,6a10)')  "CHNO.","NOmarkers","ObservedH","ExpectedH","HWT_Chi2","PIC","IBC"
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
  GDFN(3:10)=OFT(1:8)
  GDFN(1:2)=CHNO(1:2)
  GDFNIP=GDD//GDFN
  GDFNIP=adjustl(GDFNIP)
  open(q,file=GDFNIP,status='old')
  print*,"reading",GDFNIP
  do r=1, 21
    read(q,fmt='(a300)',iostat=io) Dummy
    if(io/=0) exit
  enddo
  do s=1, NOA
    read(q,fmt='(a300)',iostat=io) Dummy
    print*,Dummy
    read(q,*,iostat=io) GDS(s*2-1,SMNBC:(SMNBC+NOSMBC(q)-1))
    if(io/=0) exit
    read(q,*,iostat=io) GDS(s*2,SMNBC:(SMNBC+NOSMBC(q)-1))
  enddo
  SMNBC=SMNBC+NOSMBC(q)
  print*,"Reading imputed genotype file",q,"finished"
close(q)
enddo

do i=1, NOAC
  do j=1, NOSMBC(i)
    MNO=MNO+1
    do k=1, NOA
      GTDS(MNO,k)=GDS(k*2-1,MNO)//GDS(k*2,MNO)
      select case (GDS(k*2-1,MNO))
        case("A")
          BC(1)=BC(1)+1
        case("T")
          BC(2)=BC(2)+1
        case("G")
          BC(3)=BC(3)+1
        case("C")
          BC(4)=BC(4)+1
      end select
      select case (GDS(k*2,MNO))
        case("A")
          BC(1)=BC(1)+1
        case("T")
          BC(2)=BC(2)+1
        case("G")
          BC(3)=BC(3)+1
        case("C")
          BC(4)=BC(4)+1
      end select
      select case (GTDS(MNO,k))
        case("AA")
          GC(1)=GC(1)+1
        case("AT")
          GC(2)=GC(2)+1
        case("TA")
          GC(2)=GC(2)+1
        case("AG")
          GC(3)=GC(3)+1
        case("GA")
          GC(3)=GC(3)+1
        case("AC")
          GC(4)=GC(4)+1
        case("CA")
          GC(4)=GC(4)+1
        case("TT")
          GC(5)=GC(5)+1
        case("TG")
          GC(6)=GC(6)+1
        case("GT")
          GC(6)=GC(6)+1
        case("TC")
          GC(7)=GC(7)+1
        case("CT")
          GC(7)=GC(7)+1
        case("GG")
          GC(8)=GC(8)+1
        case("GC")
          GC(9)=GC(9)+1
        case("CG")
          GC(9)=GC(9)+1
        case("CC")
          GC(10)=GC(10)+1
      end select
    enddo
    GF=real(BC)/real(sum(BC(:4)))
    GTF=real(GC)/real(sum(GC(:10)))
    FOMIH=0.0; FOHE=0.0; FOMAH=0.0; FOMIG=0.0; FOMAG=0.0; 
    MIH='??';MAH='??';HE='??';MIBA='?';MABA='?'
    do l=1,4
      if(GF(l).NE.0 .AND. GF(l)<0.5) then
        MIBA=BA(l)
        FOMIG=GF(l)
        NOMIG=BC(l)
      elseif(GF(l).NE.0 .AND. GF(l)>0.5) then
        MABA=BA(l)
        FOMAG=GF(l)
        NOMAG=BC(l)
      elseif(GF(l)==0.5) then
        if(sum(GF(:l)).eq.0.5) then
          NOMIG=BC(l)
          MIBA=BA(l)
          FOMIG=GF(l)
        else 
          MABA=BA(l)
          FOMAG=GF(l)
          NOMAG=BC(l)
        endif
      endif
    enddo
    MIH=MIBA//MIBA
    MAH=MABA//MABA
    HE=MABA//MIBA
    if (HE=='TA') then 
      HE="AT"
    elseif(HE=='GA') then 
      HE="AG"
    elseif(HE=='CA') then 
      HE="AC"
    elseif(HE=='GT') then 
      HE="TG"
    elseif(HE=='CT') then 
      HE="TC"
    elseif(HE=='CG') then 
      HE="GC"
    endif
    do l=1,10
      if(GT(l)==MIH) then
        FOMIH=GTF(l)
        NOMIH=GC(l) 
      elseif(GT(l)==MAH) then
        FOMAH=GTF(l)
        NOMAH=GC(l)
      elseif(GT(l)==HE) then
        FOHE=GTF(l)
        NOHE=GC(l)
      endif
    enddo
    do l=1,NOA
      if(GTDS(MNO,l)==MIH) then
        G012DS(MNO,l)=0
      elseif(GTDS(MNO,l)==HE) then
        G012DS(MNO,l)=1
      elseif(GTDS(MNO,l)==MAH) then
        G012DS(MNO,l)=2
      endif
    enddo
    ENOMIH=(FOMIG**2)*real((NOMIH+NOHE+NOMAH))
    ENOHE=(2*FOMIG*FOMAG)*real((NOMIH+NOHE+NOMAH))
    ENOMAH=(FOMAG**2)*real((NOMIH+NOHE+NOMAH))
    CSOHWT=((real(NOMIH)-ENOMIH)**2/ENOMIH)+((real(NOHE)-ENOHE)**2/ENOHE)+((real(NOMAH)-ENOMAH)**2/ENOMAH)
    EH=1-(FOMIG**2+FOMAG**2)
    PIC=1-(FOMIG**2+FOMAG**2)-2*((FOMIG**2)*(FOMAG**2))
    IBC=1-(real(NOHE)/ENOHE)
    SOMAFBC=SOMAFBC+FOMIG;SOIBCBC=SOIBCBC+IBC;SOCSOHWTBC=SOCSOHWTBC+CSOHWT;SOPICBC=SOPICBC+PIC;SOHEBC=SOHEBC+FOHE;SOEHBC=SOEHBC+EH
    SOTMAF=SOTMAF+FOMIG;SOTIBC=SOTIBC+IBC;SOTCSOHWT=SOTCSOHWT+CSOHWT;SOTPIC=SOTPIC+PIC;SOTHE=SOTHE+FOHE;SOTEH=SOTEH+EH
    if(FOMIG .LT. 0.1) then 
      NOSMWMAF1=NOSMWMAF1+1
    elseif(FOMIG .GE. 0.1  .and. FOMIG .LT. 0.2) then 
      NOSMWMAF2=NOSMWMAF2+1
    elseif(FOMIG .GE. 0.2 .and. FOMIG .LT. 0.3) then 
      NOSMWMAF3=NOSMWMAF3+1    
    elseif(FOMIG .GE. 0.3 .and. FOMIG .LT. 0.4) then 
      NOSMWMAF4=NOSMWMAF4+1
    elseif(FOMIG .GE. 0.4 ) then 
      NOSMWMAF5=NOSMWMAF5+1
    endif
  enddo
  AOMAFBC=SOMAFBC/real(NOSMBC(i))
  AOIBCBC=SOIBCBC/real(NOSMBC(i))
  AOCSOHWTBC=SOCSOHWTBC/real(NOSMBC(i))
  AOPICBC=SOPICBC/real(NOSMBC(i))
  AOHEBC=SOHEBC/real(NOSMBC(i))
  AOEHBC=SOEHBC/real(NOSMBC(i))
  AODBAMBC=SODBAMBC/real(NOSMBC(i)-1)
enddo
AOTMAF=SOTMAF/real(NOM)
AOTIBC=SOTIBC/real(NOM)
AOTCSOHWT=SOTCSOHWT/real(NOM)
AOTPIC=SOTPIC/real(NOM)
AOTHE=SOTHE/real(NOM)
AOTEH=SOTEH/real(NOM)
AOTDBAM=SOTDBAM/real(NOM-NOAC)






end program RIRnBSA

