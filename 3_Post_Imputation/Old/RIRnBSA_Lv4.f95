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
integer,allocatable::G012DS(:,:),NOSMBC(:) !data input
real,allocatable::MAFS(:),IBCS(:),CSOHWTS(:),PICS(:),HES(:),EHS(:) !data input
real,allocatable::MAFSS(:),IBCSS(:),CSOHWTSS(:),PICSS(:),HESS(:),EHSS(:) !data input
real::GF(4),GTF(10),FOMAH,FOMIH,FOHE,FOMAG,FOMIG,PV !frequencies
real::HWD=0,CSOHWT=0,EH=0,PIC=0,IBC=0 !statistics for population genetics
real::SOMAFBC=0,SOIBCBC=0,SOCSOHWTBC=0,SOPICBC=0,SOHEBC=0,SOEHBC=0,SOTMAF=0,SOTIBC=0,SOTCSOHWT=0,SOTPIC=0,SOTHE=0,SOTEH=0 !statistics for population genetics by chromosome
real::AOMAFBC=0,AOIBCBC=0,AOCSOHWTBC=0,AOPICBC=0,AOHEBC=0,AOEHBC=0,AOTMAF=0,AOTIBC=0,AOTCSOHWT=0,AOTPIC=0,AOTHE=0,AOTEH=0 !statistics for population genetics by chromosome
real::SDOMAFBC=0,SDOIBCBC=0,SDOCSOHWTBC=0,SDOPICBC=0,SDOHEBC=0,SDOEHBC=0
real::SDOTMAF=0,SDOTIBC=0,SDOTCSOHWT=0,SDOTPIC=0,SDOTHE=0,SDOTEH=0 !statistics for population genetics by chromosome
real::MINOMAFBC=0,MINOIBCBC=0,MINOCSOHWTBC=0,MINOPICBC=0,MINOHEBC=0,MINOEHBC=0
real::MINOTMAF=0,MINOTIBC=0,MINOTCSOHWT=0,MINOTPIC=0,MINOTHE=0,MINOTEH=0 !statistics for population genetics by chromosome
real::MAXOMAFBC=0,MAXOIBCBC=0,MAXOCSOHWTBC=0,MAXOPICBC=0,MAXOHEBC=0,MAXOEHBC=0
real::MAXOTMAF=0,MAXOTIBC=0,MAXOTCSOHWT=0,MAXOTPIC=0,MAXOTHE=0,MAXOTEH=0 !statistics for population genetics by chromosome


character(len=1)::BA(4)=(/'A','T','G','C'/),MABA,MIBA
character(len=2)::GT(10)=(/'AA','AT','AG','AC','TT','TG','TC','GG','GC','CC'/),MAH,MIH,HE,CHNO
character(len=8)::OFT,NOAF
character(len=12)::WF1,WF2
character(len=24)::GDFN='01_DEF0000_genotypes.out'
character(len=300)::Dummy
character(len=200)::GDD,ODD
character(len=30)::ODFN
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
rewind(unit=99)
read (99,*) OFT
read (99,*) NOAF
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
allocate(MAFS(NOM),IBCS(NOM),CSOHWTS(NOM),PICS(NOM),HES(NOM),EHS(NOM))
allocate(MAFSS(NOM),IBCSS(NOM),CSOHWTSS(NOM),PICSS(NOM),HESS(NOM),EHSS(NOM))
WF1='            ';WF2='            ';
write(WF1,fmt='(a12)') '('//adjustr(NOAF)//'i1)'
write(WF2,fmt='(a12)') '('//adjustr(NOAF)//'a3)'
print *, WF1
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
write(93,fmt='(a5,6a15)')  "MAF","IBC","CSOHWT","PIC","HE","EH"
ODFNIP=ODD//"Pop_Gen_Stat_BC.out"
ODFNIP=adjustl(ODFNIP)
open(92,file=ODFNIP,status='unknown')
write(92,fmt='(a17,3a10)') "NO","AOEH","AOPIC","AOIBC"
ODFNIP=ODD//"Genotype.out"
ODFNIP=adjustl(ODFNIP)
open(91,file=ODFNIP,status='unknown')
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
  open(q+20,file=GDFNIP,status='old')
  print*,"reading",GDFNIP
  do r=1, 21
    read(q+20,fmt='(a300)',iostat=io) Dummy
    if(io/=0) exit
  enddo
  do s=1, NOA
    read(q+20,fmt='(a300)',iostat=io) Dummy
    read(q+20,*,iostat=io) GDS(s*2-1,SMNBC:(SMNBC+NOSMBC(q)-1))
    read(q+20,*,iostat=io) GDS(s*2,SMNBC:(SMNBC+NOSMBC(q)-1))
  enddo
  print *,"Reading imputed genotype from",SMNBC,"th marker to ",(SMNBC+NOSMBC(q)-1),"th marker"
  SMNBC=SMNBC+NOSMBC(q)
  print *,"Reading imputed genotype file",q,"finished"
close(q+20)
enddo
SMNBC=1
do i=1, NOAC
SOMAFBC=0;SOIBCBC=0;SOCSOHWTBC=0;SOPICBC=0;SOHEBC=0;SOEHBC=0
AOMAFBC=0;AOIBCBC=0;AOCSOHWTBC=0;AOPICBC=0;AOHEBC=0;AOEHBC=0
  do j=1, NOSMBC(i)
    MNO=MNO+1
    GF=0.0;GTF=0.0;BC=0;GC=0;
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
    GF=real(BC)/real(sum(BC(1:4)))
    GTF=real(GC)/real(sum(GC(1:10)))
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
      elseif(GTDS(MNO,l)==MAH) then
        G012DS(MNO,l)=2
      else
        G012DS(MNO,l)=1
      endif
    enddo
    write(98,fmt=WF1) G012DS(MNO,1:NOA) 
    write(97,fmt='(2F15.10)') FOMIG,FOMAG
    write(96,fmt='(3F15.10)') FOMIH,FOHE,FOMAH
    write(95,fmt='(3a3)') MIH,HE,MAH
    write(94,fmt='(2a3)') MIBA,MABA
    write(91,fmt=WF2) GTDS(MNO,1:NOA)
    ENOMIH=(FOMIG**2)*real((NOMIH+NOHE+NOMAH))
    ENOHE=(2*FOMIG*FOMAG)*real((NOMIH+NOHE+NOMAH))
    ENOMAH=(FOMAG**2)*real((NOMIH+NOHE+NOMAH))
    CSOHWT=((real(NOMIH)-ENOMIH)**2/ENOMIH)+((real(NOHE)-ENOHE)**2/ENOHE)+((real(NOMAH)-ENOMAH)**2/ENOMAH)
    EH=1-(FOMIG**2+FOMAG**2)
    PIC=1-(FOMIG**2+FOMAG**2)-2*((FOMIG**2)*(FOMAG**2))
    IBC=1-(real(NOHE)/ENOHE)
    MAFS(MNO)=FOMIG;IBCS(MNO)=IBC;CSOHWTS(MNO)=CSOHWT;PICS(MNO)=PIC;HES(MNO)=FOHE;EHS(MNO)=EH
    MAFSS(MNO)=FOMIG**2;IBCSS(MNO)=IBC**2;CSOHWTSS(MNO)=CSOHWT**2;PICSS(MNO)=PIC**2;HESS(MNO)=FOHE**2;EHSS(MNO)=EH**2
    write(93,fmt='(6F15.7)') MAFS(MNO),IBCS(MNO),CSOHWTS(MNO),PICS(MNO),HES(MNO),EHS(MNO)
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
  SOMAFBC=SOMAFBC+MAFS(MNO)
  SOIBCBC=SOIBCBC+IBCS(MNO)
  SOCSOHWTBC=SOCSOHWTBC+CSOHWTS(MNO)
  SOPICBC=SOPICBC+PICS(MNO)
  SOHEBC=SOHEBC+HES(MNO)
  SOEHBC=SOEHBC+EHS(MNO)
  enddo
  AOMAFBC=SOMAFBC/real(NOSMBC(i))
  AOIBCBC=SOIBCBC/real(NOSMBC(i))
  AOCSOHWTBC=SOCSOHWTBC/real(NOSMBC(i))
  AOPICBC=SOPICBC/real(NOSMBC(i))
  AOHEBC=SOHEBC/real(NOSMBC(i))
  AOEHBC=SOEHBC/real(NOSMBC(i))

!  SDOMAFBC=(sum(MAFSS(SMNBC:NOSMBC(i)-1))-sum(MAFS(SMNBC:NOSMBC(i)-1))**2/real(NOSMBC(i)))/real(NOSMBC(i)-1) 
!  SDOIBCBC=(sum(IBCSS(SMNBC:NOSMBC(i)-1))-sum(IBCS(SMNBC:NOSMBC(i)-1))**2/real(NOSMBC(i)))/real(NOSMBC(i)-1) 
!  SDOCSOHWTBC=(sum(CSOHWTSS(SMNBC:NOSMBC(i)-1))-sum(CSOHWTS(SMNBC:NOSMBC(i)-1))**2/real(NOSMBC(i)))/real(NOSMBC(i)-1) 
!  SDOPICBC=(sum(PICSS(SMNBC:NOSMBC(i)-1))-sum(PICS(SMNBC:NOSMBC(i)-1))**2/real(NOSMBC(i)))/real(NOSMBC(i)-1) 
!  SDOHEBC=(sum(HESS(SMNBC:NOSMBC(i)-1))-sum(HES(SMNBC:NOSMBC(i)-1))**2/real(NOSMBC(i)))/real(NOSMBC(i)-1) 
!  SDOEHBC=(sum(EHSS(SMNBC:NOSMBC(i)-1))-sum(EHS(SMNBC:NOSMBC(i)-1))**2/real(NOSMBC(i)))/real(NOSMBC(i)-1) 
!  MINOMAFBC=min(MAFS(SMNBC:NOSMBC(i)-1,1))
!  MINOIBCBC=min(IBCS(SMNBC:NOSMBC(i)-1,1))
!  MINOCSOHWTBC=min(CSOHWTS(SMNBC:NOSMBC(i)-1,1))
!  MINOPICBC=min(PICS(SMNBC:NOSMBC(i)-1,1))
!  MINOHEBC=min(HES(SMNBC:NOSMBC(i)-1,1))
!  MINOEHBC=min(EHS(SMNBC:NOSMBC(i)-1,1))
!  MAXOMAFBC=MAX(MAFS(SMNBC:NOSMBC(i)-1,1))
!  MAXOIBCBC=MAX(IBCS(SMNBC:NOSMBC(i)-1,1))
!  MAXOCSOHWTBC=MAX(CSOHWTS(SMNBC:NOSMBC(i)-1,1))
!  MAXOPICBC=MAX(PICS(SMNBC:NOSMBC(i)-1,1))
!  MAXOHEBC=MAX(HES(SMNBC:NOSMBC(i)-1,1))
!  MAXOEHBC=MAX(EHS(SMNBC:NOSMBC(i)-1,1))
!  write(93,fmt='(24F15.7)') AOMAFBC,SDOMAFBC,MINOMAFBC,MAXOMAFBC,AOIBCBC,SDOIBCBC,MINOIBCBC,MAXOIBCBC,AOCSOHWTBC,&
!  &SDOCSOHWTBC,MINOCSOHWTBC,MAXOCSOHWTBC,AOPICBC,SDOPICBC,MINOPICBC,MAXOPICBC,AOHEBC,SDOHEBC,MINOHEBC,MAXOHEBC,&
!  &AOEHBC,SDOEHBC,MINOEHBC,MAXOEHBC
  write(92,fmt='(12F15.5)') AOMAFBC,SDOMAFBC,AOIBCBC,SDOIBCBC,AOCSOHWTBC,SDOCSOHWTBC,AOPICBC,&
  &SDOPICBC,AOHEBC,SDOHEBC,AOEHBC,SDOEHBC
  SMNBC=SMNBC+NOSMBC(i)
  print *,"Converting imputed genotype file and analyzing statiscics of chromosome",i,"finished"
enddo

end program RIRnBSA

