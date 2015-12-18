! Name : Post Imputation & Basic Statistics Analysis 
! Developed by : Joonho Lee - PhD working in Hankyong National Univ.(Korea) zoonolee@gmail.com, +82-10-3408-2895
! Version : Linux 5.0 (2012-06-25)
! Contents : 1. Converting Imputation Results -> 0,1,2
!	2. Calculation of gene & genotype frequency
!	4. Statistics of data & population genetics 

program RIRnBSA 
integer::i,j,k,iz,i1,i2,p,q,r,s,ss=0,checkID=0 !iostate & do_loop
integer::NOAC=0,SMNBC=1,NOA=0,ANO=0,NOM=0,MNO=0 !number of chromosomes, animals and markers & chromosome, animal and marekr number
integer::BC(4)
integer::GC(10)
integer::CMPBC=0,EXMPBC=0,MDG(4)
integer::NODBAM1=0,NODBAM2=0,NODBAM3=0,NODBAM4=0,NODBAM5=0
integer::NOMWMAF1=0,NOMWMAF2=0,NOMWMAF3=0,NOMWMAF4=0,NOMWMAF5=0,CMP=0,EXMP=0,CN=0,EXCN=999, SODBAM=0,NOMBCN=0 ! for statistics
integer::NODMBOFOMGT=0,NODMBOAH=0,NODMBOFOMIG=0,NODMBOCSOHWT=0,NODMBONP=0,NODMBOAHE=0,EXCNO=0 
integer,allocatable::G012DS(:,:),CG101DS(:,:),NOSMBC(:),IOM(:,:),ON(:) !data input

real::GF(4),GTF(10),FOMAH,FOMIH,FOHE,FOMAG,FOMIG,IOL(4) !frequencies
real::HWD=0,CSOHWT=0,EH=0,PIC=0,IBC=0 !statistics for population genetics
real::SOMAFBC=0,SOIBCBC=0,SOCSOHWTBC=0,SOPICBC=0,SOHEBC=0,SOEHBC=0,SOTMAF=0,SOTIBC=0,SOTCSOHWT=0,SOTPIC=0,SOTHE=0,SOTEH=0 !statistics for population genetics by chromosome
real::SSOMAFBC=0,SSOIBCBC=0,SSOCSOHWTBC=0,SSOPICBC=0,SSOHEBC=0,SSOEHBC=0,SSOTMAF=0,SSOTIBC=0,SSOTCSOHWT=0,SSOTPIC=0,SSOTHE=0 !statistics for population genetics by chromosome
real::SSOTEH=0 !statistics for population genetics by chromosome
real::AOMAFBC=0,AOIBCBC=0,AOCSOHWTBC=0,AOPICBC=0,AOHEBC=0,AOEHBC=0,AOTMAF=0,AOTIBC=0,AOTCSOHWT=0,AOTPIC=0,AOTHE=0,AOTEH=0 !statistics for population genetics by chromosome
real::SDOMAFBC=0,SDOIBCBC=0,SDOCSOHWTBC=0,SDOPICBC=0,SDOHEBC=0,SDOEHBC=0
real::SDOTMAF=0,SDOTIBC=0,SDOTCSOHWT=0,SDOTPIC=0,SDOTHE=0,SDOTEH=0 !statistics for population genetics by chromosome
real::MINOMAFBC=0,MINOIBCBC=0,MINOCSOHWTBC=0,MINOPICBC=0,MINOHEBC=0,MINOEHBC=0
real::MINOTMAF=0,MINOTIBC=0,MINOTCSOHWT=0,MINOTPIC=0,MINOTHE=0,MINOTEH=0 !statistics for population genetics by chromosome
real::MAXOMAFBC=0,MAXOIBCBC=0,MAXOCSOHWTBC=0,MAXOPICBC=0,MAXOHEBC=0,MAXOEHBC=0
real::MAXOTMAF=0,MAXOTIBC=0,MAXOTCSOHWT=0,MAXOTPIC=0,MAXOTHE=0,MAXOTEH=0 !statistics for population genetics by chromosome
real(10)::SODBAMBC=0,SOTDBAM=0,SSODBAMBC=0,SSOTDBAM=0,MINODBAMBC=0,MAXODBAMBC=0,MINOTDBAM=0,MAXOTDBAM=0
real(10)::DBAMBC=0,AODBAMBC=0,SDODBAMBC=0,AOTDBAM=0,SDOTDBAM=0

character(len=1)::BA(4)=(/'A','T','G','C'/),MABA,MIBA
character(len=2)::GT(10)=(/'AA','AT','AG','AC','TT','TG','TC','GG','GC','CC'/),MAH,MIH,HE,CHNO
character(len=8)::OFT,NOFT,NOAF
character(len=12)::WF1,WF2,WF3
character(len=300)::Dummy
character(len=220)::IDD,ODD,GDD,DDD,ADD
character(len=30)::GDFN,DDFN,ADFN
character(len=15)::MWF
character(len=280)::IDFNIP
character(len=250)::GDFNIP,ODFNIP,DDFNIP,ADFNIP,PF
character(len=1),allocatable::GDS(:,:),GDS1(:)
character(len=2),allocatable::GTDS(:,:)
character(len=20),allocatable::ID(:),DID(:)
character(len=38),allocatable::IDFNS(:)
character(len=50),allocatable::MN(:)

!reading parameters
print *, "***********************************************************************************************************************" 
print *, "RIRnBSA_Lv5 : Reformation of Imputation Results & Basic Statistics Analysis Linux version 5.0"
print *, "Developed by Joonho Lee : zoonolee@gmail.com, +82-10-3408-2895"
print *, "Data form is the results of BEAGLE"
print *, "Please enter the name of parameter file including path(ex) ./parGRDH"
print *, "***********************************************************************************************************************"
read *, PF
PF=adjustl(PF)
open (99, file=PF) 
read (99,*) OFT
read (99,*) NOTA
read (99,*) NOA
rewind(unit=99)
read (99,*) OFT
print *, "file tag", OFT
read (99,*) NOTA
print *, "number of total animals on genomic data file", NOTA
read (99,*) NOAF
print *, "number of selected animals", NOA
read (99,*) IDD
print *, "directory of imputed genomic data file ", IDD
read (99,*) GDD
print *, "directory of previous data file ", GDD
read (99,*) GDFN
print *, "name of selected marker information file ", GDFN
read (99,*) ODD
print *, "directory of output file", ODD
if(NOTA.ne.NOA) then
  read (99,*) NOFT
  read (99,*) ADD
  print *, "directory of selected animal file ", ADD
  read (99,*) ADFN
  print *, "name of selected animal file ", ADFN
  read (99,*) DDD
  print *, "directory of being deleted animal list file ", DDD
  read (99,*) DDFN
  print *, "name of being deleted animal list file ", DDFN
  allocate(ID(NOTA),DID(NOTA-NOA))
endif  
read (99,fmt='(4f10.5)') IOL(:)
print *, "Outliers : MAF<",IOL(2),"HW-chisquare>",IOL(3)  
read (99,fmt='(4i15)') MDG(:)
read (99,*) NOAC
print *, "number of autosomal chromosomes : ", NOAC
allocate(NOSMBC(NOAC),IDFNS(NOAC),GDS1(NOTA*2))
do i=1,NOAC
  read (99,*) NOSMBC(i),IDFNS(i)
  NOM=NOM+NOSMBC(i)
  print*,IDFNS(i)
enddo
close(99)
print *, "number of markers", NOM
do p=1,NOAC
  print *, "Number of markers on chromosome No.",p," : ", NOSMBC(p)
enddo
IDD=adjustr(IDD)
allocate(GDS(NOM,NOA*2))
allocate(GTDS(NOM,NOA))
allocate(G012DS(NOM,NOA))
allocate(CG101DS(NOM,NOA))
allocate(IOM(NOM,3))
allocate(MN(NOM))
allocate(ON(NOTA))
WF1='            ';WF2='            ';WF3='            ';
write(WF1,fmt='(a12)') '('//adjustr(NOAF)//'i1)'
write(WF2,fmt='(a12)') '('//adjustr(NOAF)//'a3)'
write(WF3,fmt='(a12)') '('//adjustr(NOAF)//'i3)'
print *, WF1
ODD=adjustr(ODD)
GDD=adjustr(GDD)
ODFNIP=ODD//"Genotype012"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(98,file=ODFNIP,status='unknown')
ODFNIP=ODD//"Genefrequency"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(97,file=ODFNIP,status='unknown')
write(97,fmt='(2a15)') "0","1"
ODFNIP=ODD//"Genotypefrequency"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(96,file=ODFNIP,status='unknown')
write(96,fmt='(3a15)') "0","1","2"
ODFNIP=ODD//"GT_maching_012"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(95,file=ODFNIP,status='unknown')
write(95,fmt='(3a3)') "0","1","2"
ODFNIP=ODD//"G_maching_01"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(94,file=ODFNIP,status='unknown')
write(94,fmt='(2a3)') "0","1"
ODFNIP=ODD//"Pop_Gen_Stat"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(93,file=ODFNIP,status='unknown')
write(93,fmt='(8a15)') "MARKER NO.","MAF","IBC","CSOHWT","PIC","O_HE","E_HE","DBAM"
ODFNIP=ODD//"Pop_Gen_Stat_BC"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(92,file=ODFNIP,status='unknown')
write(92,fmt='(a5,28a15)') "CHNO","AMAF","SDMAF","MINMAF","MAXMAF","AIBC","SDIBC","MINIBC","MAXIBC","ACSOHWT","SDCSOHWT"&
&,"MINCSOHWT","MAXCSOHWT","APIC","SDPIC","MINPIC","MAXPIC","AHE","SDHE","MINHE","MAXHE","AEH","SDEH","MINEH","MAXEH"&
&,"ADBAM","SDDBAM","MINDBAM","MAXDBAM"
ODFNIP=ODD//"Genotype"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(91,file=ODFNIP,status='unknown')
ODFNIP=ODD//"CenteredGT_101"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(90,file=ODFNIP,status='unknown')
GDFNIP=GDD//GDFN
GDFNIP=adjustl(GDFNIP)
open(89,file=GDFNIP,status='old')
ODFNIP=ODD//"Statistics"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(88,file=ODFNIP,status='unknown')
if(NOTA.ne.NOA) then
  DDD=adjustr(DDD)
  ADD=adjustr(ADD)
  ADFNIP=ADD//ADFN
  ADFNIP=adjustl(ADFNIP)
  open(87,file=ADFNIP,status='old')
  DDFNIP=DDD//DDFN
  DDFNIP=adjustl(DDFNIP)
  open(86,file=DDFNIP,status='old')
  ODFNIP=ODD//"Selected_Animals"//OFT//".out"
  ODFNIP=adjustl(ODFNIP)
  open(85,file=ODFNIP,status='unknown')
  write(85,fmt='(a20,2a10)') "Animal ID    ","Animal NO.","Animal NNO."
  read(87,*) 
  do q=1,NOTA
    read(87,*) ID(q),ON(q)
  enddo
  do q=1,NOTA-NOA
    read(86,*) DID(q)
  enddo
endif  

do q=1, NOAC
  if(NOTA.ne.NOA) OFT=NOFT
!  iz=ichar('0')
!  i1=floor(q/10.)+iz
!  i2=mod(q,10)+iz
!  CHNO(1:1)=achar(i1)
!  CHNO(2:2)=achar(i2)
!  IDFN(3:10)=OFT(1:8)
!  IDFN(1:2)=CHNO(1:2)
  IDFNIP=IDD//IDFNS(q)
  IDFNIP=adjustl(IDFNIP)
  open(q+20,file=IDFNIP,status='old')
  print*,"reading",IDFNIP
  read(q+20,fmt='(a300)') Dummy
  print *,"Reading imputed genotype from",SMNBC,"th marker to ",(SMNBC+NOSMBC(q)-1),"th marker"
  do r=1,NOSMBC(q)
    MNO=MNO+1
    read(q+20,*) Dummy,MN(MNO),GDS1(:)
    if(NOTA.ne.NOA) then
      do s=1, NOTA
        do i=1, NOTA-NOA
          if(ID(s).eq.DID(i)) checkID=checkID+1
        enddo  
        if(checkID.ge.2) print *, s,"th animal was found ",checkID,"times in delete animal list please check and run again" 
        if(checkID.gt.0) then
          if(q.eq.1.and.r.eq.1) print *, "animal ",ID(s),"was deleted"
        else
          ss=ss+1
          GDS(MNO,ss*2-1)=GDS1(s*2-1)
          GDS(MNO,ss*2)=GDS1(s*2)
          if(q.eq.1.and.r.eq.1) write(85,fmt='(a20,2i10)') ID(s),ON(s),ss
        endif
        checkID=0
      enddo
      ss=0
    else
      do s=1, NOTA
        GDS(MNO,s*2-1)=GDS1(s*2-1)
        GDS(MNO,s*2)=GDS1(s*2)
      enddo
    endif  
  enddo
  SMNBC=SMNBC+NOSMBC(q)
  print *,"Reading imputed genotype file",q,"finished"
close(q+20)
enddo
MNO=0
SOTMAF=0;SOTIBC=0;SOTCSOHWT=0;SOTPIC=0;SOTHE=0;SOTEH=0
SSOTMAF=0;SSOTIBC=0;SSOTCSOHWT=0;SSOTPIC=0;SSOTHE=0;SSOTEH=0
AOTMAF=0;AOTIBC=0;AOTCSOHWT=0;AOTPIC=0;AOTHE=0;AOTEH=0
SDOTMAF=0;SDOTIBC=0;SDOTCSOHWT=0;SDOTPIC=0;SDOTHE=0;SDOTEH=0
MINOTMAF=0;MINOTIBC=0;MINOTCSOHWT=0;MINOTPIC=0;MINOTHE=0;MINOTEH=0
MAXOTMAF=0;MAXOTIBC=0;MAXOCSOTHWT=0;MAXOTPIC=0;MAXTOHE=0;MAXOTEH=0
read(89,fmt='(a300)',iostat=io) Dummy
print *, Dummy
do i=1, NOAC
SOMAFBC=0;SOIBCBC=0;SOCSOHWTBC=0;SOPICBC=0;SOHEBC=0;SOEHBC=0
SSOMAFBC=0;SSOIBCBC=0;SSOCSOHWTBC=0;SSOPICBC=0;SSOHEBC=0;SSOEHBC=0
AOMAFBC=0;AOIBCBC=0;AOCSOHWTBC=0;AOPICBC=0;AOHEBC=0;AOEHBC=0
SDOMAFBC=0;SDOIBCBC=0;SDOCSOHWTBC=0;SDOPICBC=0;SDOHEBC=0;SDOEHBC=0
MINOMAFBC=0;MINOIBCBC=0;MINOCSOHWTBC=0;MINOPICBC=0;MINOHEBC=0;MINOEHBC=0
MAXOMAFBC=0;MAXOIBCBC=0;MAXOCSOHWTBC=0;MAXOPICBC=0;MAXOHEBC=0;MAXOEHBC=0
SODBAMBC=0;SSODBAMBC=0;AODBAMBC=0;SDODBAMBC=0;MINODBAMBC=0;MAXODBAMBC=0;
DBAMBC=0;CMPBC=0;EXMPBC=0
  do j=1, NOSMBC(i)
    MNO=MNO+1
    read(89,fmt='(a50,i10,2i15)') MN(MNO),IOM(MNO,1),IOM(MNO,2),IOM(MNO,3)
    GF=0.0;GTF=0.0;BC=0;GC=0;
    do k=1, NOA
      GTDS(MNO,k)=GDS(MNO,k*2-1)//GDS(MNO,k*2)
      select case (GDS(MNO,k*2-1))
        case("A")
          BC(1)=BC(1)+1
        case("T")
          BC(2)=BC(2)+1
        case("G")
          BC(3)=BC(3)+1
        case("C")
          BC(4)=BC(4)+1
      end select
      select case (GDS(MNO,k*2))
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
        CG101DS(MNO,l)=-1
      elseif(GTDS(MNO,l)==MAH) then
        G012DS(MNO,l)=2
        CG101DS(MNO,l)=1
      else
        G012DS(MNO,l)=1
        CG101DS(MNO,l)=0
      endif
    enddo
    write(98,fmt=WF1) G012DS(MNO,1:NOA)
    write(90,fmt=WF3) CG101DS(MNO,1:NOA) 
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
    if(FOMIG<IOL(2).or.CSOHWT>IOL(3).or.FOHE.eq.1) print *, "Marker",MNO,"looks unusual",FOMIG,FOHE,CSOHWT
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
    CMPBC=IOM(MNO,3)
    if(j.ne.1) then
      DBAMBC=real(CMPBC-EXMPBC)/1000
      SODBAMBC=SODBAMBC+DBAMBC
      SSODBAMBC=SSODBAMBC+DBAMBC**2
      if(DBAMBC .LT. MDG(1)) then 
        NODBAM1=NODBAM1+1
      elseif(DBAMBC .GE. MDG(1) .and. DBAMBC .LT. MDG(2)) then 
        NODBAM2=NODBAM2+1
      elseif(DBAMBC .GE. MDG(2) .and. DBAMBC .LT. MDG(3)) then 
        NODBAM3=NODBAM3+1    
      elseif(DBAMBC .GE. MDG(3) .and. DBAMBC .LT. MDG(4)) then 
        NODBAM4=NODBAM4+1
      elseif(DBAMBC .GE. MDG(4) ) then 
        NODBAM5=NODBAM5+1
      endif
      if(j.eq.2) then; 
        MINODBAMBC=DBAMBC
        MAXODBAMBC=DBAMBC
      endif
      MINODBAMBC=MIN(MINODBAMBC,DBAMBC)
      MAXODBAMBC=MAX(MAXODBAMBC,DBAMBC)
    endif
    EXMPBC=CMPBC
    if(j.eq.1) then; 
      MINOMAFBC=FOMIG;MINOIBCBC=IBC;MINOCSOHWTBC=CSOHWT;MINOPICBC=PIC;MINOHEBC=FOHE;MINOEHBC=EH
      MAXOMAFBC=FOMIG;MAXOIBCBC=IBC;MAXOCSOHWTBC=CSOHWT;MAXOPICBC=PIC;MAXOHEBC=FOHE;MAXOEHBC=EH
    endif
    MINOMAFBC=MIN(MINOMAFBC,FOMIG);MINOIBCBC=MIN(MINOIBCBC,IBC);MINOCSOHWTBC=MIN(MINOCSOHWTBC,CSOHWT);
    MINOPICBC=MIN(MINOPICBC,PIC);MINOHEBC=MIN(MINOHEBC,FOHE);MINOEHBC=MIN(MINOEHBC,EH)
    MAXOMAFBC=MAX(MAXOMAFBC,FOMIG);MAXOIBCBC=MAX(MAXOIBCBC,IBC);MAXOCSOHWTBC=MAX(MAXOCSOHWTBC,CSOHWT);
    MAXOPICBC=MAX(MAXOPICBC,PIC);MAXOHEBC=MAX(MAXOHEBC,FOHE);MAXOEHBC=MAX(MAXOEHBC,EH)
    SOMAFBC=SOMAFBC+FOMIG
    SOIBCBC=SOIBCBC+IBC
    SOCSOHWTBC=SOCSOHWTBC+CSOHWT
    SOPICBC=SOPICBC+PIC
    SOHEBC=SOHEBC+FOHE
    SOEHBC=SOEHBC+EH
    SSOMAFBC=SSOMAFBC+FOMIG**2
    SSOIBCBC=SSOIBCBC+IBC**2
    SSOCSOHWTBC=SSOCSOHWTBC+CSOHWT**2
    SSOPICBC=SSOPICBC+PIC**2
    SSOHEBC=SSOHEBC+FOHE**2
    SSOEHBC=SSOEHBC+EH**2
    write(93,fmt='(i15,7F15.7)') MNO,FOMIG,IBC,CSOHWT,PIC,FOHE,EH,DBAMBC
  enddo
  AOMAFBC=SOMAFBC/real(NOSMBC(i))
  AOIBCBC=SOIBCBC/real(NOSMBC(i))
  AOCSOHWTBC=SOCSOHWTBC/real(NOSMBC(i))
  AOPICBC=SOPICBC/real(NOSMBC(i))
  AOHEBC=SOHEBC/real(NOSMBC(i))
  AOEHBC=SOEHBC/real(NOSMBC(i))
  AODBAMBC=SODBAMBC/real(NOSMBC(i)-1)
  SDOMAFBC=(SSOMAFBC-SOMAFBC**2/real(NOSMBC(i)))/real(NOSMBC(i))
  SDOIBCBC=(SSOIBCBC-SOIBCBC**2/real(NOSMBC(i)))/real(NOSMBC(i))
  SDOCSOHWTBC=(SSOCSOHWTBC-SOCSOHWTBC**2/real(NOSMBC(i)))/real(NOSMBC(i))
  SDOPICBC=(SSOPICBC-SOPICBC**2/real(NOSMBC(i)))/real(NOSMBC(i))
  SDOHEBC=(SSOHEBC-SOHEBC**2/real(NOSMBC(i)))/real(NOSMBC(i))
  SDOEHBC=(SSOEHBC-SOEHBC**2/real(NOSMBC(i)))/real(NOSMBC(i))
  SDODBAMBC=(SSODBAMBC-SODBAMBC**2/real(NOSMBC(i)-1))/real(NOSMBC(i)-1)
  write(92,fmt='(i5,28f15.5)') i,AOMAFBC,SDOMAFBC,MINOMAFBC,MAXOMAFBC,AOIBCBC,SDOIBCBC,MINOIBCBC,MAXOIBCBC,&
  AOCSOHWTBC,SDOCSOHWTBC,MINOCSOHWTBC,MAXOCSOHWTBC,AOPICBC,SDOPICBC,MINOPICBC,MAXOPICBC,AOHEBC,SDOHEBC,&
  MINOHEBC,MAXOHEBC,AOEHBC,SDOEHBC,MINOEHBC,MAXOEHBC,AODBAMBC,SDODBAMBC,MINODBAMBC,MAXODBAMBC
  SOTMAF=SOTMAF+SOMAFBC
  SOTIBC=SOTIBC+SOIBCBC
  SOTCSOHWT=SOTCSOHWT+SOCSOHWTBC
  SOTPIC=SOTPIC+SOPICBC
  SOTHE=SOTHE+SOHEBC
  SOTEH=SOTEH+SOEHBC
  SOTDBAM=SOTDBAM+SODBAMBC
  SSOTMAF=SSOTMAF+SSOMAFBC
  SSOTIBC=SSOTIBC+SSOIBCBC
  SSOTCSOHWT=SSOTCSOHWT+SSOCSOHWTBC
  SSOTPIC=SSOTPIC+SSOPICBC
  SSOTHE=SSOTHE+SSOHEBC
  SSOTEH=SSOTEH+SSOEHBC
  SSOTDBAM=SSOTDBAM+SSODBAMBC
  if(i.eq.1) then; 
    MINOTMAF=MINOMAFBC;MINOTIBC=MINOIBCBC;MINOTCSOHWT=MINOCSOHWTBC;MINOTPIC=MINOPICBC;MINOTHE=MINOHEBC;MINOTEH=MINOEHBC
    MAXOTMAF=MAXOMAFBC;MAXOTIBC=MAXOIBCBC;MAXOTCSOHWT=MAXOCSOHWTBC;MAXOTPIC=MAXOPICBC;MAXOTHE=MAXOHEBC;MAXOTEH=MAXOEHBC
    MINOTDBAM=MINODBAMBC;MAXOTDBAM=MAXODBAMBC;
  endif
  MINOTMAF=MIN(MINOMAFBC,MINOTMAF);MINOTIBC=MIN(MINOIBCBC,MINOTIBC);MINOTCSOHWT=MIN(MINOCSOHWTBC,MINOTCSOHWT);
  MINOTPIC=MIN(MINOPICBC,MINOTPIC);MINOTHE=MIN(MINOHEBC,MINOTHE);MINOTEH=MIN(MINOEHBC,MINOTEH)
  MAXOTMAF=MAX(MAXOMAFBC,MAXOTMAF);MAXOTIBC=MAX(MAXOIBCBC,MAXOTIBC);MAXOTCSOHWT=MAX(MAXOCSOHWTBC,MAXOTCSOHWT);
  MAXOTPIC=MAX(MAXOPICBC,MAXOTPIC);MAXOTHE=MAX(MAXOHEBC,MAXOTHE);MAXOTEH=MAX(MAXOEHBC,MAXOTEH)
  MINOTDBAM=MIN(MINODBAMBC,MINOTDBAM);MAXOTDBAM=MAX(MAXODBAMBC,MAXOTDBAM);
  print *,"Reading marker information on",IOM(MNO,2),"th chromosome"
  print *,"Converting imputed genotype file and analyzing statiscics of chromosome",i,"finished"
enddo
AOTMAF=SOTMAF/real(NOM)
AOTIBC=SOTIBC/real(NOM)
AOTCSOHWT=SOTCSOHWT/real(NOM)
AOTPIC=SOTPIC/real(NOM)
AOTHE=SOTHE/real(NOM)
AOTEH=SOTEH/real(NOM)
AOTDBAM=SOTDBAM/real(NOM-NOAC)
SDOTMAF=(SSOTMAF-SOTMAF**2/real(NOM))/real(NOM)
SDOTIBC=(SSOTIBC-SOTIBC**2/real(NOM))/real(NOM)
SDOTCSOHWT=(SSOTCSOHWT-SOTCSOHWT**2/real(NOM))/real(NOM)
SDOTPIC=(SSOTPIC-SOTPIC**2/real(NOM))/real(NOM)
SDOTHE=(SSOTHE-SOTHE**2/real(NOM))/real(NOM)
SDOTEH=(SSOTEH-SOTEH**2/real(NOM))/real(NOM)
SDOTDBAM=(SSOTDBAM-SOTDBAM**2/real(NOM-NOAC))/real(NOM-NOAC)
write(92,fmt='(a5,28f15.5)') "Tot",AOTMAF,SDOTMAF,MINOTMAF,MAXOTMAF,AOTIBC,SDOTIBC,MINOTIBC,MAXOTIBC,&
AOTCSOHWT,SDOTCSOHWT,MINOTCSOHWT,MAXOTCSOHWT,AOTPIC,SDOTPIC,MINOTPIC,MAXOTPIC,AOTHE,SDOTHE,MINOTHE,MAXOTHE,&
&AOTEH,SDOTEH,MINOTEH,MAXOTEH,AOTDBAM,SDOTDBAM,MINOTDBAM,MAXOTDBAM
    write(88,fmt='(a120)') "***** frequency of MAF *****"
    write(88,fmt='(a120)') "number of markers with Minor allel frequency"
    write(88,fmt='(a120,i15)') "       minor allele frequency < 0.1  :  ", NOSMWMAF1
    write(88,fmt='(a120,i15)') "0.1 <= minor allele frequency < 0.2  :  ", NOSMWMAF2
    write(88,fmt='(a120,i15)') "0.2 <= minor allele frequency < 0.4  :  ", NOSMWMAF3
    write(88,fmt='(a120,i15)') "0.3 <= minor allele frequency < 0.4  :  ", NOSMWMAF4
    write(88,fmt='(a120,i15)') "0.4 <= minor allele frequency        :  ", NOSMWMAF5
    write(88,fmt='(a120)') " "
    write(88,fmt='(a120)') "***** frequency of distance *****"
    write(88,fmt='(a120)') "number of total markers with distance between adjacent markers"
    write(88,fmt='(50x,a55,2i15)') "               distance between adjacent markers < ",MDG(1), NODBAM1
    write(88,fmt='(50x,i15,a40,2i15)') MDG(1),"<= distance between adjacent markers <",MDG(2), NODBAM2
    write(88,fmt='(50x,i15,a40,2i15)') MDG(2),"<= distance between adjacent markers <",MDG(3), NODBAM3
    write(88,fmt='(50x,i15,a40,2i15)') MDG(3),"<= distance between adjacent markers <",MDG(4), NODBAM4
    write(88,fmt='(50x,i15,a55,i15)') MDG(4),"<= distance between adjacent markers                 ", NODBAM5
end program RIRnBSA
