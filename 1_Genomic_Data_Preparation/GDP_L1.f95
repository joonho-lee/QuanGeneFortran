! Name : Outlier Elimination, Basic Statistics of Genomic Raw Data & Preparation for BEAGLE 
! Developed by : Joonho Lee - PhD working in Hankyong National Univ.(Korea) zoonolee@gmail.com
! Version : Linux version 1.0 (2014-07-26)
! Contents : 1. Data reading(animal ID, SNP marker name, chromosome number, marker position, phenotype(name, records))
!	2. Calculation of gene & genotype frequency
!	3. Unusual data elimination
!	4. Statistics of data & population genetics 
!	5. Transformation of genotype data for BEAGLE 

program GDP_L1
integer::io,i,j,k,l,m,n,o,p=0,q=0,r,iz,i1,i2,SMNBC=1 !iostate & do_loop
integer::NOAC=0,XCN=19,YCN=20,UCN=99,ANO=0,EXCN=999 !number of autosoal chromosomes & chromosome number for X, Y, Unknown
integer::BC(5)
integer::GC(11)
integer::NOA,NOZ=0,NOZBC=0,NANO,NOMOAC=0,NOM,NOSM,NOMAG=0,NOMIG=0,NOMIH=0,NOHE=0,NOMAH=0,NOMG=0,NOMOXC=0,NOMOYC=0,NOMOUC=0,TNOMG=0 !counting
integer::NOAMM=0,NOAMMBC=0,NONMM=0,NODM=0,NODMBC=0,NODA=0,NODMBOFOMGT=0,NODMBOAH=0,NODMBOFOMIG=0,NODMBOCSOHWT=0,NODMBOAHE=0 !counting
integer::ACNONMM=0,ACNOAMM=0,ACNODMBOFOMGT=0,ACNODMBOAH=0,ACNODMBOFOMIG=0,ACNODMBOCSOHWT=0,ACNODMBOAHE=0 !counting
integer::UCNONMM=0,UCNOAMM=0,UCNODMBOFOMGT=0,UCNODMBOAH=0,UCNODMBOFOMIG=0,UCNODMBOCSOHWT=0,UCNODMBOAHE=0 !counting
integer::SCNONMM=0,SCNOAMM=0,SCNODMBOFOMGT=0,SCNODMBOAH=0,SCNODMBOFOMIG=0,SCNODMBOCSOHWT=0,SCNODMBOAHE=0 !counting
integer::NODMBO_M_H=0,NODMBO_F_H=0,NODMBO_M_F=0,NODMBO_MHF=0,NODMBO_HOM=0,NODMBO_HEM=0 !counting
integer::ACNODMBO_M_H=0,ACNODMBO_F_H=0,ACNODMBO_M_F=0,ACNODMBO_MHF=0,ACNODMBO_HOM=0,ACNODMBO_HEM=0 !counting
integer::UCNODMBO_M_H=0,UCNODMBO_F_H=0,UCNODMBO_M_F=0,UCNODMBO_MHF=0,UCNODMBO_HOM=0,UCNODMBO_HEM=0 !counting
integer::SCNODMBO_M_H=0,SCNODMBO_F_H=0,SCNODMBO_M_F=0,SCNODMBO_MHF=0,SCNODMBO_HOM=0,SCNODMBO_HEM=0 !counting
integer::NODMBOFOMGTBC=0,NODMBOAHBC=0,NODMBOFOMIGBC=0,NODMBOCSOHWTBC=0,NODMBOAHEBC=0 !counting
integer::NODMBO_M_HBC=0,NODMBO_F_HBC=0,NODMBO_M_FBC=0,NODMBO_MHFBC=0,NODMBO_HOMBC=0,NODMBO_HEMBC=0 !counting
integer::NOSMWMAF1=0,NOSMWMAF2=0,NOSMWMAF3=0,NOSMWMAF4=0,NOSMWMAF5=0,NOTMWMAF1=0,NOTMWMAF2=0,NOTMWMAF3=0,NOTMWMAF4=0,NOTMWMAF5=0
integer::NODBAM1=0,NODBAM2=0,NODBAM3=0,NODBAM4=0,NODBAM5=0,NODBASM1=0,NODBASM2=0,NODBASM3=0,NODBASM4=0,NODBASM5=0
integer::CMPBC=0,EXMPBC=0
integer::CMPBC_S=0,EXMPBC_S=0
integer*8::SODBAMBC=0,SOTDBAM=0,SODBAMBC_S=0,SOTDBAM_S=0
integer::NOMBC=0,MDG(4) ! for statistics
integer,allocatable::IOM(:,:),GN(:,:),NOMGBA(:),NOSMBC(:) !data input & missing data output
real::GF(5),GTF(11),FOMAH,FOMIH,FOHE,FOMAG,FOMIG,SOGF,SOGTF,FOMG,FOMGT,IOL(4),SOPIC=0,AOPIC=0 !frequencies
real::DBAMBC=0,DBAMBC_S=0
real::ENOMIH=0,ENOHE=0,ENOMAH=0,CSOHWT=0,EH=0,PIC=0,IBC=0,TMGR=0,AODBAMBC_S=0,AODBAMBC=0,AOTDBAM_S=0,AOTDBAM=0 !statistics for population genetics
real::SOMAFBC_S=0,SOIBCBC_S=0,SOCSOHWTBC_S=0,SOPICBC_S=0,SOHEBC_S=0,SOEHBC_S=0
real::AOMAFBC_S=0,AOIBCBC_S=0,AOCSOHWTBC_S=0,AOPICBC_S=0,AOHEBC_S=0,AOEHBC_S=0
real::SOTMAF_S=0,SOTIBC_S=0,SOTCSOHWT_S=0,SOTPIC_S=0,SOTHE_S=0,SOTEH_S=0 !statistics for population genetics by chromosome
real::AOTMAF_S=0,AOTIBC_S=0,AOTCSOHWT_S=0,AOTPIC_S=0,AOTHE_S=0,AOTEH_S=0 !statistics for population genetics by chromosome
real::SOMAFBC=0,SOIBCBC=0,SOCSOHWTBC=0,SOPICBC=0,SOHEBC=0,SOEHBC=0,SOTMAF=0,SOTIBC=0,SOTCSOHWT=0,SOTPIC=0,SOTHE=0,SOTEH=0 !statistics for population genetics by chromosome
real::AOMAFBC=0,AOIBCBC=0,AOCSOHWTBC=0,AOPICBC=0,AOHEBC=0,AOEHBC=0,AOTMAF=0,AOTIBC=0,AOTCSOHWT=0,AOTPIC=0,AOTHE=0,AOTEH=0 !statistics for population genetics by chromosome
character(len=50)::MN, D1, D2, D3
character(len=50),allocatable::MKN(:)
character(len=1)::BA(5)=(/'A','T','G','C','Z'/),MABA,MIBA,MG
character(len=2)::GT(11)=(/'AA','AT','AG','AC','TT','TG','TC','GG','GC','CC','ZZ'/),MAH,MIH,HE,CHNO
character(len=200)::GDD,ODD,IDD,PDD
character(len=30)::GDFN,WF
character(len=80)::r_BEAGLE
character(len=27)::MWF
character(len=6)::SOSM="NORMAL"
character(len=8)::OFT="_DEF0000"
character(len=230)::GDFNIP,ODFNIP,IDFNIP,PDFNIP,PFIP
character(len=30),allocatable::ID(:)
character(len=2),allocatable::GTS(:)
character(len=1),allocatable::GTDFFP(:,:),GTDFFP2(:,:)

!reading parameters
print *, "***********************************************************************************************************************" 
print *, "GDP_L1 : Outlier Elimination of Genomic Raw Data & Preparation for BEAGLE Linux version 1.0"
print *, "Developed by Joonho Lee : zoonolee@gmail.com"
print *, "Data format : (1st row) Markername Chromosome Position LL09066048	LL08131059	LL09147068 ... (Animal IDs)"
print *, "Data format : (2nd row) ALGA0000009	1	52297	AG	GG	GG	... (Genotypes)"
print *, "Input data must be sorted by chromosome No.(1st) and Marker position(2nd)"
print *, "Sex chromosome and markers of unknown chromosome must be numbered more then the biggest chromosome number"
print *, "Length of Output file tag must be 8 characters"
print *, "Please enter the name of parameter file including directory(ex) ./parGRDH"
print *, "***********************************************************************************************************************"
read *, PFIP
PFIP=adjustl(PFIP)
open (99, file=PFIP) 
read (99,*) OFT
read (99,*) NOA
read (99,*) GDD
read (99,*) GDFN
read (99,*) ODD
read (99,*) IDD
read (99,*) PDD
read (99,*) GT(11)
read (99,*) IOL(:)
read (99,*) NOAC
read (99,*) XCN, YCN, UCN
read (99,*) DXCN, DYCN, DUCN
read (99,*) MDG(:)

print *, "number of animals", NOA
print *, "directory of genomic data file", GDD
print *, "name of genomic data file?", GDFN
print *, "directory of output file", ODD
print *, "missing genotype", GT(11)
print *, "range of outliers"
print *, "   marker-missing% : ", IOL(1)
print *, "   minor allele frequency : ", IOL(2)
print *, "   HW-chisquare(df=1) : ", IOL(3)
print *, "   animal-missing% : ", IOL(4)
print *, "number of autosomal chromosomes : ", NOAC
print *, "X chromosome No. , Y chromosome No. Unknown chromosome No. ", XCN, YCN, UCN

GDFN=adjustl(GDFN)
GDD=adjustr(GDD)
GDFNIP=GDD//GDFN
GDFNIP=adjustl(GDFNIP)

print *, "input genomic data file is ",GDFNIP 
print *, "Output file for Statistics will be created in",ODD
print *, "Output file for Imputation will be created in",IDD
print *, "Parameter file for Reformation will be created in",PDD  
ODD=adjustr(ODD)
IDD=adjustr(IDD)
PDD=adjustr(PDD)
close(99)

open(1,file=GDFNIP,status='unknown')
ODFNIP=ODD//"Deleted_Markers"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(11,file=ODFNIP,status='unknown')
write(11,fmt='(a50,6a10,2a15)') "Marker Name","MarkerNO","SOSM","FOMGT","FOMIG","FOHE","CSOHWT","chromosome","position"
ODFNIP=ODD//"Selected_Markers"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(12,file=ODFNIP,status='unknown')
write(12,fmt='(a50,a10,2a15)') "Marker Name","MarkerNO","chromosome","position"
ODFNIP=ODD//"Deleted_Animals"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(13,file=ODFNIP,status='unknown')
write(13,fmt='(a30,a10,2a15)') "Animal ID    ","Animal NO.","NO.Missing","Missing rate"
ODFNIP=ODD//"Selected_Animals"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(14,file=ODFNIP,status='unknown')
write(14,fmt='(a30,2a10,2a15)') "Animal ID    ","Animal NO.","Animal NNO.","NO.Missing","Missing rate"
IDFNIP=IDD//"run_BEAGLE"//OFT//".sh"
IDFNIP=adjustl(IDFNIP)
open(10,file=IDFNIP,status='unknown')
IDFNIP=IDD//"run_gzip"//OFT//".sh"
IDFNIP=adjustl(IDFNIP)
open(15,file=IDFNIP,status='unknown')
write(10,*) "cd "//adjustl(IDD)
do i=1,NOAC
  iz=ichar('0')
  i1=floor(i/10.)+iz
  i2=mod(i,10)+iz
  CHNO(1:1)=achar(i1)
  CHNO(2:2)=achar(i2)
  IDFNIP=IDD//"To_BEAGLE"//CHNO(1:2)//OFT//".bgl"
  IDFNIP=adjustl(IDFNIP)
  write(10,*) "java -Xmx1000m -jar beagle.jar unphased="//"To_BEAGLE"//CHNO(1:2)//OFT//".bgl missing=? out=imputed"
  open(i+20,file=IDFNIP,status='unknown')
enddo
ODFNIP=ODD//"Basic_Statistics"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(16,file=ODFNIP,status='unknown')
ODFNIP=ODD//"By_Chromosome_S"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(17,file=ODFNIP,status='unknown')
write(17,fmt='(a14,19a10)') "Chromosome NO.","NOM","NOSM","NODM","NODM_M","NODM_F","NODM_H","NODM_M_H","NODM_F_H",&
&"NODM_M_F","NODM_MHF","NODM_AHO","NODM_AHE","AODAMB","AOMAF","AOIBC","AOCSOHWT","AOPIC","AOHE","AOEH"
ODFNIP=ODD//"By_Chromosome_A"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(20,file=ODFNIP,status='unknown')
write(20,fmt='(a14,8a10)') "Chromosome NO.","NOM","AODAMB","AOMAF","AOIBC","AOCSOHWT","AOPIC","AOHE","AOEH"
PDFNIP=PDD//"PIP"//OFT//".par"
PDFNIP=adjustl(PDFNIP)
open(18,file=PDFNIP,status='unknown')
ODFNIP=ODD//"Whole_Results"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(19,file=ODFNIP,status='unknown')
write(19,fmt='(a50,a5,a15,3a10,3a5,3a6,3a10,2a5,2a6,6a10)') "Marker Name","CHno","position","SOSM","FOMGT","NOMGT","MAH","HE",&
"MIH","NOMAH","NOHE","NOMIH","FOMAH","FOHE","FOMIH","MABA","MIBA","NOMAG","NOMIG","FOMAG","FOMIG","EH","CSOHWT","PIC","IBC"
NOM=0
do 
 read(1,*,iostat=io)
  if(io/=0) exit
 NOM=NOM+1
enddo
rewind(unit=1)
NOM=NOM-1
allocate(IOM(NOM,2),MKN(NOM),GN(NOM,NOA),GTDFFP(NOM,NOA*2))
allocate(ID(NOA),GTS(NOA),NOMGBA(NOA),NOSMBC(NOAC))
NOSMBC=0
IOM=0
GN=0
NOMGBA=0
NODA=0
NOSM=0
NANO=0
GTDFFP=' '
MKN='                                                  '
read(1,*)D1,D2,D3,ID
MG=GT(11)(1:1)
NODM=0
do i=1,NOM
  GTS='  '
  GF=0.0
  GTF=0.0
  FOMGT=0.0; FOMG=0.0
  read(1,*,iostat=io) MN, IOM(i,:), GTS
  if(io/=0) exit
  BC=0; GC=0; nMI=0; nMA=0; NOMIH=0; NOHE=0; NOMAH=0; NOMG=0
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
    select case (GTS(j))
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
      case("??")
        GC(11)=GC(11)+1
        TNOMG=TNOMG+1
    end select
  enddo
  GF=real(BC)/real(sum(BC(:4)))
  GTF=real(GC)/real(sum(GC(:10)))
  FOMG=real(BC(5))/real(NOA*2)
  FOMGT=real(GC(11))/(real(NOA))
  FOMIH=0.0; FOHE=0.0; FOMAH=0.0; FOMIG=0.0; FOMAG=0.0; 
  MIH='??';MAH='??';HE='??';MIBA='?';MABA='?';SOSM='NORMAL'
  do k=1,4
    if(GF(k).NE.0 .AND. GF(k)<0.5) then
      MIBA=BA(k)
      FOMIG=GF(k)
      NOMIG=BC(k)
    elseif(GF(k).NE.0 .AND. GF(k)>0.5) then
      MABA=BA(k)
      FOMAG=GF(k)
      NOMAG=BC(k)
    elseif(GF(k)==0.5) then
      if(sum(GF(:k)).eq.0.5) then
        NOMIG=BC(k)
        MIBA=BA(k)
        FOMIG=GF(k)
      else 
        MABA=BA(k)
        FOMAG=GF(k)
        NOMAG=BC(k)
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
  do k=1,10
    if(GT(k)==MIH) then
      FOMIH=GTF(k)
      NOMIH=GC(k) 
    elseif(GT(k)==MAH) then
      FOMAH=GTF(k)
      NOMAH=GC(k)
    elseif(GT(k)==HE) then
      FOHE=GTF(k)
      NOHE=GC(k)
    endif
  enddo
  if(EXCN.ne.XCN.and.EXCN.ne.YCN.and.EXCN.ne.UCN.and.EXCN.ne.999.and.IOM(i,1).ne.EXCN) then 
    AOMAFBC_S=SOMAFBC_S/real(NOSMBC(IOM(i-1,1)))
    AOIBCBC_S=SOIBCBC_S/real(NOSMBC(IOM(i-1,1)))
    AOCSOHWTBC_S=SOCSOHWTBC_S/real(NOSMBC(IOM(i-1,1)))
    AOPICBC_S=SOPICBC_S/real(NOSMBC(IOM(i-1,1)))
    AOHEBC_S=SOHEBC_S/real(NOSMBC(IOM(i-1,1)))
    AOEHBC_S=SOEHBC_S/real(NOSMBC(IOM(i-1,1)))
    AODBAMBC_S=SODBAMBC_S/real(NOSMBC(IOM(i-1,1))-1)
    AOMAFBC=SOMAFBC/real(NOMBC)
    AOIBCBC=SOIBCBC/real(NOMBC-NOZBC)
    AOCSOHWTBC=SOCSOHWTBC/real(NOMBC-NOZBC)
    AOPICBC=SOPICBC/real(NOMBC)
    AOHEBC=SOHEBC/real(NOMBC)
    AOEHBC=SOEHBC/real(NOMBC)
    AODBAMBC=SODBAMBC/real(NOMBC-1)
    write(19,fmt='(a1)') " "
    write(17,fmt='(a10,i4,12i10,f10.1,6f10.5)') "chromosome",EXCN,NOMBC,NOSMBC(IOM(i-1,1)),NODMBC,NODMBOFOMGTBC,NODMBOFOMIGBC,&
    NODMBOCSOHWTBC,NODMBO_M_HBC,NODMBO_F_HBC,NODMBO_M_FBC,NODMBO_MHFBC,NODMBO_HOM,NODMBO_HEM,AODBAMBC_S,AOMAFBC_S,AOIBCBC_S,&
    AOCSOHWTBC_S,AOPICBC_S,AOHEBC_S,AOEHBC_S
    write(20,fmt='(a10,i4,i10,f10.1,6f10.5)') "chromosome",EXCN,NOMBC,AODBAMBC,AOMAFBC,AOIBCBC,AOCSOHWTBC,AOPICBC,AOHEBC,AOEHBC
    NODMBOCSOHWTBC=0;NODMBOFOMGTBC=0;NODMBOFOMIGBC=0;NODMBO_M_HBC=0;NODMBO_F_HBC=0;NODMBO_M_FBC=0;NODMBO_MHFBC=0;NOZBC=0
    NODMBO_HOM=0;NODMBO_HEM=0
    SODBAMBC=0;SOMAFBC=0;SOIBCBC=0;SOCSOHWTBC=0;SOPICBC=0;SOHEBC=0;SOEHBC=0;NOMBC=0;NODMBC=0
    SODBAMBC_S=0;SOMAFBC_S=0;SOIBCBC_S=0;SOCSOHWTBC_S=0;SOPICBC_S=0;SOHEBC_S=0;SOEHBC_S=0
  end if
  ENOMIH=(FOMIG**2)*real((NOMIH+NOHE+NOMAH))
  ENOHE=(2*FOMIG*FOMAG)*real((NOMIH+NOHE+NOMAH))
  ENOMAH=(FOMAG**2)*real((NOMIH+NOHE+NOMAH))
  CSOHWT=((real(NOMIH)-ENOMIH)**2/ENOMIH)+((real(NOHE)-ENOHE)**2/ENOHE)+((real(NOMAH)-ENOMAH)**2/ENOMAH)
  EH=1-(FOMIG**2+FOMAG**2)
  PIC=1-(FOMIG**2+FOMAG**2)-2*((FOMIG**2)*(FOMAG**2))
  IBC=1-(real(NOHE)/ENOHE)
  if(ENOMIH.eq.0.or.ENOMAH.eq.0.or.ENOHE.eq.0) then; NOZ=NOZ+1;NOZBC=NOZBC+1;CSOHWT=0;IBC=0;endif;   
  if(CSOHWT>IOL(3)) then; NODMBOCSOHWT=NODMBOCSOHWT+1;SOSM="OL_HWT";NODMBOCSOHWTBC=NODMBOCSOHWTBC+1;endif;
  if(FOMIG<IOL(2)) then; NODMBOFOMIG=NODMBOFOMIG+1;SOSM="OL_MAF";NODMBOFOMIGBC=NODMBOFOMIGBC+1;endif;
  if(FOMIG<IOL(2).and.CSOHWT>IOL(3)) then; NODMBO_F_H=NODMBO_F_H+1;SOSM="OL_F_H";NODMBO_F_HBC=NODMBO_F_HBC+1;endif;
  if(FOMGT>IOL(1)) then; NODMBOFOMGT=NODMBOFOMGT+1;SOSM="OL_MGR";NODMBOFOMGTBC=NODMBOFOMGTBC+1;endif;
  if(FOMGT>IOL(1).and.CSOHWT>IOL(3)) then; NODMBO_M_H=NODMBO_M_H+1;SOSM="OL_M_H";NODMBO_M_HBC=NODMBO_M_HBC+1;endif;
  if(FOMGT>IOL(1).and.FOMIG<IOL(2))  then; NODMBO_M_F=NODMBO_M_F+1;SOSM="OL_M_F";NODMBO_M_FBC=NODMBO_M_FBC+1;endif;
  if(FOMGT>IOL(1).and.FOMIG<IOL(2).and.CSOHWT>IOL(3)) then; NODMBO_MHF=NODMBO_MHF+1;SOSM="OL_MHF";NODMBO_MHFBC=NODMBO_MHFBC+1;endif;
  if(FOMGT.eq.1) then; NOAMM=NOAMM+1;SOSM="OL_AMG";NOAMMBC=NOAMMBC+1;endif;
  if(FOMAH.eq.1) then; NODMBOAH=NODMBOAH+1;SOSM="OL_AHO";NODMBOAHBC=NODMBOAHBC+1;endif;
  if(FOHE.eq.1) then; NODMBOAHE=NODMBOAHE+1;SOSM="OL_AHE";NODMBOAHEBC=NODMBOAHEBC+1;endif;
  if(FOMAH.eq.1 .and. FOMGT>IOL(1)) then; NODMBO_HOM=NODMBO_HOM+1;SOSM="OL_HOM";NODMBO_HOMBC=NODMBO_HOMBC+1;endif;
  if(FOHE.eq.1 .and. FOMGT>IOL(1)) then; NODMBO_HEM=NODMBO_HEM+1;SOSM="OL_HEM";NODMBO_HEMBC=NODMBOA_HEMBC+1;endif;
  if(FOMGT.eq.0) NONMM=NONMM+1
  if(IOM(i,1).ne.XCN.and.IOM(i,1).ne.YCN.and.IOM(i,1).ne.UCN) then
    if(CSOHWT>IOL(3)) ACNODMBOCSOHWT=ACNODMBOCSOHWT+1
    if(FOMIG<IOL(2)) ACNODMBOFOMIG=ACNODMBOFOMIG+1
    if(FOMGT>IOL(1)) ACNODMBOFOMGT=ACNODMBOFOMGT+1
    if(FOMIG<IOL(2).and.CSOHWT>IOL(3)) ACNODMBO_F_H=ACNODMBO_F_H+1
    if(FOMGT>IOL(1).and.CSOHWT>IOL(3)) ACNODMBO_M_H=ACNODMBO_M_H+1
    if(FOMGT>IOL(1).and.FOMIG<IOL(2))  ACNODMBO_M_F=ACNODMBO_M_F+1
    if(FOMGT>IOL(1).and.FOMIG<IOL(2).and.CSOHWT>IOL(3)) ACNODMBO_MHF=ACNODMBO_MHF+1
    if(FOMGT.eq.1) ACNOAMM=ACNOAMM+1
    if(FOMGT.eq.0) ACNONMM=ACNONMM+1
    if(FOMAH.eq.1) ACNODMBOAH=ACNODMBOAH+1
    if(FOHE.eq.1) ACNODMBOAHE=ACNODMBOAHE+1
    if(FOMAH.eq.1.and.FOMGT>IOL(1)) ACNODMBO_HOM=ACNODMBO_HOM+1
    if(FOHE.eq.1.and.FOMGT>IOL(1)) ACNODMBO_HEM=ACNODMBO_HEM+1
  elseif(IOM(i,1).eq.UCN) then 
    if(CSOHWT>IOL(3)) UCNODMBOCSOHWT=UCNODMBOCSOHWT+1
    if(FOMIG<IOL(2)) UCNODMBOFOMIG=UCNODMBOFOMIG+1
    if(FOMGT>IOL(1)) UCNODMBOFOMGT=UCNODMBOFOMGT+1
    if(FOMIG<IOL(2).and.CSOHWT>IOL(3)) UCNODMBO_F_H=UCNODMBO_F_H+1
    if(FOMGT>IOL(1).and.CSOHWT>IOL(3)) UCNODMBO_M_H=UCNODMBO_M_H+1
    if(FOMGT>IOL(1).and.FOMIG<IOL(2))  UCNODMBO_M_F=UCNODMBO_M_F+1
    if(FOMGT>IOL(1).and.FOMIG<IOL(2).and.CSOHWT>IOL(3)) UCNODMBO_MHF=UCNODMBO_MHF+1
    if(FOMGT.eq.1) UCNOAMM=UCNOAMM+1
    if(FOMGT.eq.0) UCNONMM=UCNONMM+1
    if(FOMAH.eq.1) UCNODMBOAH=UCNODMBOAH+1
    if(FOHE.eq.1) UCNODMBOAHE=UCNODMBOAHE+1
    if(FOMAH.eq.1.and.FOMGT>IOL(1)) UCNODMBO_HOM=UCNODMBO_HOM+1
    if(FOHE.eq.1.and.FOMGT>IOL(1)) UCNODMBO_HEM=UCNODMBO_HEM+1
  elseif(IOM(i,1).eq.XCN.or.IOM(i,1).eq.YCN) then 
    if(CSOHWT>IOL(3)) SCNODMBOCSOHWT=SCNODMBOCSOHWT+1
    if(FOMIG<IOL(2)) SCNODMBOFOMIG=SCNODMBOFOMIG+1
    if(FOMGT>IOL(1)) SCNODMBOFOMGT=SCNODMBOFOMGT+1
    if(FOMIG<IOL(2).and.CSOHWT>IOL(3)) SCNODMBO_F_H=SCNODMBO_F_H+1
    if(FOMGT>IOL(1).and.CSOHWT>IOL(3)) SCNODMBO_M_H=SCNODMBO_M_H+1
    if(FOMGT>IOL(1).and.FOMIG<IOL(2))  SCNODMBO_M_F=SCNODMBO_M_F+1
    if(FOMGT>IOL(1).and.FOMIG<IOL(2).and.CSOHWT>IOL(3)) SCNODMBO_MHF=SCNODMBO_MHF+1
    if(FOMGT.eq.1) SCNOAMM=SCNOAMM+1
    if(FOMGT.eq.0) SCNONMM=SCNONMM+1
    if(FOMAH.eq.1) SCNODMBOAH=SCNODMBOAH+1
    if(FOHE.eq.1) SCNODMBOAHE=SCNODMBOAHE+1
    if(FOMAH.eq.1.and.FOMGT>IOL(1)) SCNODMBO_HOM=SCNODMBO_HOM+1
    if(FOHE.eq.1.and.FOMGT>IOL(1)) SCNODMBO_HEM=SCNODMBO_HEM+1
  endif
  if(IOM(i,1) .eq. XCN) then; NOMOXC=NOMOXC+1;SOSM="OL_XCH";endif; 
  if(IOM(i,1) .eq. YCN) then; NOMOYC=NOMOYC+1;SOSM="OL_YCH";endif;
  if(IOM(i,1) .eq. UCN) then; NOMOUC=NOMOUC+1;SOSM="OL_UCH";endif;
  write(19,fmt='(a50,i5,i15,a10,f10.6,i10,3a5,3i6,3f10.6,2a5,2i6,3f10.6,f10.3,2f10.6)') MN,IOM(i,1),IOM(i,2),SOSM,FOMGT,GC(11),&
  &MAH,HE,MIH,NOMAH,NOHE,NOMIH,FOMAH,FOHE,FOMIH,MABA,MIBA,NOMAG,NOMIG,FOMAG,FOMIG,EH,CSOHWT,PIC,IBC
  if(IOM(i,1).ne.DXCN.and.IOM(i,1).ne.DYCN.and.IOM(i,1).ne.DUCN) then
    NOMOAC=NOMOAC+1 
    NOMBC=NOMBC+1
    if(FOMIG .LT. 0.1) then 
      NOTMWMAF1=NOTMWMAF1+1
    elseif(FOMIG .GE. 0.1  .and. FOMIG .LT. 0.2) then 
      NOTMWMAF2=NOTMWMAF2+1
    elseif(FOMIG .GE. 0.2 .and. FOMIG .LT. 0.3) then 
      NOTMWMAF3=NOTMWMAF3+1    
    elseif(FOMIG .GE. 0.3 .and. FOMIG .LT. 0.4) then 
      NOTMWMAF4=NOTMWMAF4+1
    elseif(FOMIG .GE. 0.4 ) then 
      NOTMWMAF5=NOTMWMAF5+1
    endif
    CMPBC=IOM(i,2)
    if(NOMBC.ne.1) then
      DBAMBC=real(CMPBC-EXMPBC)/1000.0
      SODBAMBC=SODBAMBC+DBAMBC
      SOTDBAM=SOTDBAM+DBAMBC
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
    endif
    EXMPBC=CMPBC
    SOMAFBC=SOMAFBC+FOMIG;SOIBCBC=SOIBCBC+IBC;SOCSOHWTBC=SOCSOHWTBC+CSOHWT;SOPICBC=SOPICBC+PIC;SOHEBC=SOHEBC+FOHE;SOEHBC=SOEHBC+EH
    SOTMAF=SOTMAF+FOMIG;SOTIBC=SOTIBC+IBC;SOTCSOHWT=SOTCSOHWT+CSOHWT;SOTPIC=SOTPIC+PIC;SOTHE=SOTHE+FOHE;SOTEH=SOTEH+EH
  endif
  if(FOMGT>IOL(1).or.FOMIG<IOL(2).or.CSOHWT>IOL(3).or.FOHE.EQ.1.or.IOM(i,1).eq.DXCN.or.IOM(i,1).eq.DYCN&
  &.or.IOM(i,1).eq.DUCN) then 
    NODM=NODM+1
    NODMBC=NODMBC+1
    write(11,fmt='(a50,i10,a10,3f10.7,f10.3,2i15)') MN,i,SOSM,FOMGT,FOMIG,FOHE,CSOHWT,IOM(i,1),IOM(i,2) 
  else
    SOMAFBC_S=SOMAFBC_S+FOMIG;SOIBCBC_S=SOIBCBC_S+IBC;SOCSOHWTBC_S=SOCSOHWTBC_S+CSOHWT;
    SOPICBC_S=SOPICBC_S+PIC;SOHEBC_S=SOHEBC_S+FOHE;SOEHBC_S=SOEHBC_S+EH
    SOTMAF_S=SOTMAF_S+FOMIG;SOTIBC_S=SOTIBC_S+IBC;SOTCSOHWT_S=SOTCSOHWT_S+CSOHWT;
    SOTPIC_S=SOTPIC_S+PIC;SOTHE_S=SOTHE_S+FOHE;SOTEH_S=SOTEH_S+EH
    NOSMBC(IOM(i,1))=NOSMBC(IOM(i,1))+1
    CMPBC_S=IOM(i,2)
    if(NOSMBC(IOM(i,1)).ne.1) then;
      DBAMBC_S=real(CMPBC_S-EXMPBC_S)/1000.0
      SODBAMBC_S=SODBAMBC_S+DBAMBC_S
      SOTDBAM_S=SOTDBAM_S+DBAMBC_S
      if(DBAMBC_S .LT. MDG(1)) then 
        NODBASM1=NODBASM1+1
      elseif(DBAMBC_S .GE. MDG(1) .and. DBAMBC_S .LT. MDG(2)) then 
        NODBASM2=NODBASM2+1
      elseif(DBAMBC_S .GE. MDG(2) .and. DBAMBC_S .LT. MDG(3)) then 
        NODBASM3=NODBASM3+1    
      elseif(DBAMBC_S .GE. MDG(3) .and. DBAMBC_S .LT. MDG(4)) then 
        NODBASM4=NODBASM4+1
      elseif(DBAMBC_S .GE. MDG(4) ) then 
        NODBASM5=NODBASM5+1
      endif
    endif
    EXMPBC_S=CMPBC_S
    NOSM=NOSM+1
    MKN(NOSM)=MN
    write(12,fmt='(a50,i10,2i15)') MN,i,IOM(i,1),IOM(i,2)
    do l=1,NOA
      GTDFFP(NOSM,l*2-1)=GTS(l)(1:1)
      GTDFFP(NOSM,l*2)=GTS(l)(2:2)
      if(GTS(l).eq.'??') NOMGBA(l)=NOMGBA(l)+1 
    enddo
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
  endif
  EXCN=IOM(i,1)
enddo
AOTMAF_S=SOTMAF_S/real(NOM-NODM)
AOTIBC_S=SOTIBC_S/real(NOM-NODM)
AOTCSOHWT_S=SOTCSOHWT_S/real(NOM-NODM)
AOTPIC_S=SOTPIC_S/real(NOM-NODM)
AOTHE_S=SOTHE_S/real(NOM-NODM)
AOTEH_S=SOTEH_S/real(NOM-NODM)
AOTDBAM_S=SOTDBAM_S/real(NOM-NODM-NOAC)
AOTMAF=SOTMAF/real(NOMOAC)
AOTIBC=SOTIBC/real(NOMOAC-NOZ)
AOTCSOHWT=SOTCSOHWT/real(NOMOAC-NOZ)
AOTPIC=SOTPIC/real(NOMOAC)
AOTHE=SOTHE/real(NOMOAC)
AOTEH=SOTEH/real(NOMOAC)
AOTDBAM=SOTDBAM/real(NOMOAC-NOAC)
write(17,fmt='(a14,12i10,f10.1,6f10.5)') "total  ",NOMOAC,NOM-NODM,NOMOAC-(NOM-NODM),ACNODMBOFOMGT,ACNODMBOFOMIG,&
ACNODMBOCSOHWT,ACNODMBO_M_H,ACNODMBO_F_H,ACNODMBO_M_F,ACNODMBO_MHF,ACNODMBO_HOM,ACNODMBO_HEM,AOTDBAM_S,AOTMAF_S,&
&AOTIBC_S,AOTCSOHWT_S,AOTPIC_S,AOTHE_S,AOTEH_S
write(20,fmt='(a14,i10,f10.1,6f10.5)') "total  ",NOMOAC,AOTDBAM,AOTMAF,AOTIBC,AOTCSOHWT,AOTPIC,AOTHE,AOTEH

print *, "Done with the outlier marker elimination : ",NODM,"markers were deleted" 

allocate(GTDFFP2(NOSM,NOA*2))
GTDFFP2=' '

do m=1,NOA
  if((real(NOMGBA(m))/real((NOM-NODM))).GT.IOL(4)) then
    NODA=NODA+1
    write(13,fmt='(a30,i10,i15,f15.7)') ID(m),m,NOMGBA(m),real(NOMGBA(m))/real((NOM-NODM))
  else
    NANO=NANO+1
    write(14,fmt='(a30,2i10,i15,f15.7)') ID(m),m,NANO,NOMGBA(m),real(NOMGBA(m))/real((NOM-NODM))
    do p=1, NOSM
      GTDFFP2(p,NANO*2-1)=GTDFFP(p,m*2-1)
      GTDFFP2(p,NANO*2)=GTDFFP(p,m*2)
    enddo
  endif
enddo

write(MWF,fmt='(a14,i10,a3)') '(a1,1x,a50,1x,',NANO*2,'a2)'
MWF=trim(MWF) 
MWF=adjustl(MWF) 

do o=1,NOAC
  print *, "Number of markers on chromosome No.",o,NOSMBC(o)
  do m=SMNBC,SMNBC+NOSMBC(o)-1
    write(o+20,fmt=MWF) "M",MKN(m),GTDFFP2(m,:)
  enddo
  SMNBC=SMNBC+NOSMBC(o)
close(o+20)
enddo  
print *, "Done with the outlier animal elimination : ",NODA,"animals were deleted"
    write(16,fmt='(a120)') "*****Raw data & parameter description*****"
    write(16,fmt='(a120,i15)') "number of animals", NOA
    write(16,fmt='(a120,a230)') "directory of genomic data file", GDD
    write(16,fmt='(a120,a230)') "name of genomic data file", GDFN
    write(16,fmt='(a120,a230)') "directory of output file", ODD
    write(16,fmt='(a120,a15)') "missing genotype", GT(11)
    write(16,fmt='(a120,i15)') "range of outliers"
    write(16,fmt='(a120,f15.3)') "   marker-missing% :", IOL(1)
    write(16,fmt='(a120,f15.3)') "   minor allele frequency :", IOL(2)
    write(16,fmt='(a120,f15.3)') "   HW-chisquare(df=1) :", IOL(3)
    write(16,fmt='(a120,f15.3)') "   animal-missing% :", IOL(4)
    write(16,fmt='(a120,i15)') "number of autosomal chromosomes :", NOAC
    write(16,fmt='(a120,3i5)') "X chromosome No. , Y chromosome No. Unknown chromosome No. ", XCN, YCN, UCN
    write(16,fmt='(a120)') " "
    write(16,fmt='(a120)') "*****Basic statistics*****"
    write(16,fmt='(a120,i15)') "Number of markers on chromosome   X",NOMOXC
    write(16,fmt='(a120,i15)') "Number of markers on chromosome   Y",NOMOYC
    write(16,fmt='(a120,i15)') "Number of markers with Unknown marker location(position)",NOMOUC    
    write(16,fmt='(a120,i15)') "Total number of markers", NOM
    write(16,fmt='(a120,i15)') "Number of markers on atuosomal chromosomes", NOMOAC
    write(16,fmt='(a120,i15)') "Number of selected(useful) markers", NOM-NODM
    write(16,fmt='(a120,i15)') "Number of outlier markers", NODM
    write(16,fmt='(a109,f10.2,a1,i15)') "Number of animals with missing over ",IOL(4),"%", NODA
    TMGR=real(TNOMG)/real(NOM*NOA)
    write(16,fmt='(a120,f15.10)') "Total missing rate(Number of Missing genotypes/Number of genotypes(NO.animal*NO.markers))" , TMGR
    write(16,fmt='(a120)') " "
    write(16,fmt='(a120)') "*** total ***"
    write(16,fmt='(a120,i15)') "Number of all missing markers", NOAMM
    write(16,fmt='(a120,i15)') "Number of all homo_genotype markers", NODMBOAH
    write(16,fmt='(a120,i15)') "Number of all hetero_genotype markers", NODMBOAHE
    write(16,fmt='(a109,f10.2,a1,i15)') "Number of all homo_genotype markers and missing over ", IOL(1)*100,"%", NODMBO_HOM
    write(16,fmt='(a109,f10.2,a1,i15)') "Number of all hetero_genotype markers and missing over ", IOL(1)*100,"%", NODMBO_HEM
    write(16,fmt='(a120,i15)') "Number of markers with no missing", NONMM         
    write(16,fmt='(a109,f10.2,a1,i15)') "Number of markers with missing over",IOL(1)*100,"%", NODMBOFOMGT
    write(16,fmt='(a110,f10.2,i15)') "Number of markers with minor allele frequency less than",IOL(2), NODMBOFOMIG
    write(16,fmt='(a110,f10.2,i15)') "Number of markers with HWchisquare over",IOL(3), NODMBOCSOHWT
    write(16,fmt='(a78,f10.2,a22,f10.2,i15)') "Number of markers with missing over",IOL(1)*100,"% and HWchisquare over",IOL(3),&
    NODMBO_M_H
    write(16,fmt='(a79,f10.2,a21,f10.2,i15)') "Number of markers with minor allele frequency less than",IOL(2),&
    " and HWchisquare over",IOL(3),NODMBO_F_H
    write(16,fmt='(a62,f10.2,a38,f10.2,i15)') "Number of markers with missing over ",IOL(1)*100,&
    "% and minor allele frequency less than" ,IOL(2),NODMBO_M_F
    write(16,fmt='(a35,f9.2,a35,f10.2,a21,f10.2,i15)') "Number of markers with missing over",IOL(1)*100,&
    "%, minor allele frequency less than" ,IOL(2)," and HWchisquare over",IOL(3), NODMBO_MHF
    write(16,fmt='(a120)') " "
    write(16,fmt='(a120)') "*** on autosomal chromosome ***"
    write(16,fmt='(a120,i15)') "Number of all missing markers", ACNOAMM
    write(16,fmt='(a120,i15)') "Number of all homo_genotype markers", ACNODMBOAH
    write(16,fmt='(a120,i15)') "Number of all hetero_genotype markers", ACNODMBOAHE
    write(16,fmt='(a109,f10.2,a1,i15)') "Number of all homo_genotype markers and missing over ", IOL(1)*100,"%", ACNODMBO_HOM
    write(16,fmt='(a109,f10.2,a1,i15)') "Number of all hetero_genotype markers and missing over ", IOL(1)*100,"%", ACNODMBO_HEM
    write(16,fmt='(a120,i15)') "Number of markers with no missing", ACNONMM         
    write(16,fmt='(a109,f10.2,a1,i15)') "Number of markers with missing over",IOL(1)*100,"%", ACNODMBOFOMGT
    write(16,fmt='(a110,f10.2,i15)') "Number of markers with minor allele frequency less than",IOL(2), ACNODMBOFOMIG
    write(16,fmt='(a110,f10.2,i15)') "Number of markers with HWchisquare over",IOL(3), ACNODMBOCSOHWT
    write(16,fmt='(a78,f10.2,a22,f10.2,i15)') "Number of markers with missing over",IOL(1)*100,"% and HWchisquare over",IOL(3),&
    ACNODMBO_M_H
    write(16,fmt='(a79,f10.2,a21,f10.2,i15)') "Number of markers with minor allele frequency less than",IOL(2),&
    " and HWchisquare over",IOL(3),ACNODMBO_F_H
    write(16,fmt='(a62,f10.2,a38,f10.2,i15)') "Number of markers with missing over ",IOL(1)*100,&
    "% and minor allele frequency less than" ,IOL(2),ACNODMBO_M_F
    write(16,fmt='(a35,f9.2,a35,f10.2,a21,f10.2,i15)') "Number of markers with missing over",IOL(1)*100,&
    "%, minor allele frequency less than" ,IOL(2)," and HWchisquare over",IOL(3), ACNODMBO_MHF
    write(16,fmt='(a120)') " "
    write(16,fmt='(a120)') "*** on unknown chromosome ***"
    write(16,fmt='(a120,i15)') "Number of all missing markers", UCNOAMM
    write(16,fmt='(a120,i15)') "Number of all homo_genotype markers", UCNODMBOAH
    write(16,fmt='(a120,i15)') "Number of all hetero_genotype markers", UCNODMBOAHE
    write(16,fmt='(a109,f10.2,a1,i15)') "Number of all homo_genotype markers and missing over ", IOL(1)*100,"%", UCNODMBO_HOM
    write(16,fmt='(a109,f10.2,a1,i15)') "Number of all hetero_genotype markers and missing over ", IOL(1)*100,"%", UCNODMBO_HEM
    write(16,fmt='(a120,i15)') "Number of markers with no missing", UCNONMM         
    write(16,fmt='(a109,f10.2,a1,i15)') "Number of markers with missing over",IOL(1)*100,"%", UCNODMBOFOMGT
    write(16,fmt='(a110,f10.2,i15)') "Number of markers with minor allele frequency less than",IOL(2), UCNODMBOFOMIG
    write(16,fmt='(a110,f10.2,i15)') "Number of markers with HWchisquare over",IOL(3), UCNODMBOCSOHWT
    write(16,fmt='(a78,f10.2,a22,f10.2,i15)') "Number of markers with missing over",IOL(1)*100,"% and HWchisquare over",IOL(3),&
    UCNODMBO_M_H
    write(16,fmt='(a79,f10.2,a21,f10.2,i15)') "Number of markers with minor allele frequency less than",IOL(2),&
    " and HWchisquare over",IOL(3),UCNODMBO_F_H
    write(16,fmt='(a62,f10.2,a38,f10.2,i15)') "Number of markers with missing over ",IOL(1)*100,&
    "% and minor allele frequency less than" ,IOL(2),UCNODMBO_M_F
    write(16,fmt='(a35,f9.2,a35,f10.2,a21,f10.2,i15)') "Number of markers with missing over",IOL(1)*100,&
    "%, minor allele frequency less than" ,IOL(2)," and HWchisquare over",IOL(3), UCNODMBO_MHF
    write(16,fmt='(a120)') " "
    write(16,fmt='(a120)') "*** on sex chromosome ***"
    write(16,fmt='(a120,i15)') "Number of all missing markers", SCNOAMM
    write(16,fmt='(a120,i15)') "Number of all homo_genotype markers", SCNODMBOAH
    write(16,fmt='(a120,i15)') "Number of all hetero_genotype markers", SCNODMBOAHE
    write(16,fmt='(a109,f10.2,a1,i15)') "Number of all homo_genotype markers and missing over ", IOL(1)*100,"%", SCNODMBO_HOM
    write(16,fmt='(a109,f10.2,a1,i15)') "Number of all hetero_genotype markers and missing over ", IOL(1)*100,"%", SCNODMBO_HEM
    write(16,fmt='(a120,i15)') "Number of markers with no missing", SCNONMM         
    write(16,fmt='(a109,f10.2,a1,i15)') "Number of markers with missing over",IOL(1)*100,"%", SCNODMBOFOMGT
    write(16,fmt='(a110,f10.2,i15)') "Number of markers with minor allele frequency less than",IOL(2), SCNODMBOFOMIG
    write(16,fmt='(a110,f10.2,i15)') "Number of markers with HWchisquare over",IOL(3), SCNODMBOCSOHWT
    write(16,fmt='(a78,f10.2,a22,f10.2,i15)') "Number of markers with missing over",IOL(1)*100,"% and HWchisquare over",IOL(3),&
    SCNODMBO_M_H
    write(16,fmt='(a79,f10.2,a21,f10.2,i15)') "Number of markers with minor allele frequency less than",IOL(2),&
    " and HWchisquare over",IOL(3),SCNODMBO_F_H
    write(16,fmt='(a62,f10.2,a38,f10.2,i15)') "Number of markers with missing over ",IOL(1)*100,&
    "% and minor allele frequency less than" ,IOL(2),SCNODMBO_M_F
    write(16,fmt='(a35,f9.2,a35,f10.2,a21,f10.2,i15)') "Number of markers with missing over",IOL(1)*100,&
    "%, minor allele frequency less than" ,IOL(2)," and HWchisquare over",IOL(3), SCNODMBO_MHF
    write(16,fmt='(a120)') " "
    write(16,fmt='(a120)') "***** frequency of MAF *****"
    write(16,fmt='(a120)') "number of total markers with Minor allel frequency"
    write(16,fmt='(a120,i15)') "       minor allele frequency < 0.1  :  ", NOTMWMAF1
    write(16,fmt='(a120,i15)') "0.1 <= minor allele frequency < 0.2  :  ", NOTMWMAF2
    write(16,fmt='(a120,i15)') "0.2 <= minor allele frequency < 0.3  :  ", NOTMWMAF3
    write(16,fmt='(a120,i15)') "0.3 <= minor allele frequency < 0.4  :  ", NOTMWMAF4
    write(16,fmt='(a120,i15)') "0.4 <= minor allele frequency        :  ", NOTMWMAF5
    write(16,fmt='(a120)') " "
    write(16,fmt='(a120)') "number of selected markers with Minor allel frequency(except for unusual data)"
    write(16,fmt='(a120,i15)') "       minor allele frequency < 0.1  :  ", NOSMWMAF1
    write(16,fmt='(a120,i15)') "0.1 <= minor allele frequency < 0.2  :  ", NOSMWMAF2
    write(16,fmt='(a120,i15)') "0.2 <= minor allele frequency < 0.3  :  ", NOSMWMAF3
    write(16,fmt='(a120,i15)') "0.3 <= minor allele frequency < 0.4  :  ", NOSMWMAF4
    write(16,fmt='(a120,i15)') "0.4 <= minor allele frequency        :  ", NOSMWMAF5
    write(16,fmt='(a120)') " "
    write(16,fmt='(a120)') "***** frequency of distance *****"
    write(16,fmt='(a120)') "number of total markers with distance between adjacent markers"
    write(16,fmt='(50x,a55,2i15)') "               distance between adjacent markers < ",MDG(1), NODBAM1
    write(16,fmt='(50x,i15,a40,2i15)') MDG(1),"<= distance between adjacent markers <",MDG(2), NODBAM2
    write(16,fmt='(50x,i15,a40,2i15)') MDG(2),"<= distance between adjacent markers <",MDG(3), NODBAM3
    write(16,fmt='(50x,i15,a40,2i15)') MDG(3),"<= distance between adjacent markers <",MDG(4), NODBAM4
    write(16,fmt='(50x,i15,a55,i15)') MDG(4),"<= distance between adjacent markers                 ", NODBAM5
    write(16,fmt='(a120)') " "
    write(16,fmt='(a120)') "number of selected markers with distance between adjacent markers(except for unusual data)"
    write(16,fmt='(50x,a55,2i15)') "          distance between adjacent markers <",MDG(1), NODBASM1
    write(16,fmt='(50x,i15,a40,2i15)') MDG(1),"<= distance between adjacent markers <",MDG(2), NODBASM2
    write(16,fmt='(50x,i15,a40,2i15)') MDG(2),"<= distance between adjacent markers <",MDG(3), NODBASM3
    write(16,fmt='(50x,i15,a40,2i15)') MDG(3),"<= distance between adjacent markers <",MDG(4), NODBASM4
    write(16,fmt='(50x,i15,a55,i15)') MDG(4),"<= distance between adjacent markers                 ", NODBASM5
    write(18,fmt='(a10)') "'"//OFT//"'"
    write(18,*) NANO
    write(18,*) NANO
    IDD="'"//adjustl(IDD)
    write(18,*) IDD//"'"
    ODD="'"//adjustl(ODD)
    write(18,*) ODD//"'"
    write(18,fmt='(a30)') "'"//"Selected_Markers"//OFT//".out"//"'"
    PDD="'"//adjustl(PDD)
    write(18,*) PDD//"'"
    write(18,fmt='(4f10.5)') IOL(:)
    write(18,fmt='(4i15)') MDG(:)
    write(18,*) NOAC
do r=1,NOAC
  iz=ichar('0')
  i1=floor(r/10.)+iz
  i2=mod(r,10)+iz
  CHNO(1:1)=achar(i1)
  CHNO(2:2)=achar(i2)
  IDFNIP=IDD//"To_BEAGLE"//CHNO(1:2)//OFT//".bgl"
  r_BEAGLE="cd ../2_Imputation/|nohup java -Xmx1000m -jar beagle.jar unphased=To_BEAGLE"//CHNO(1:2)&
  &//OFT//".bgl missing=? out=imputed &"
  write(18,*) NOSMBC(r),"imputed.To_BEAGLE"//CHNO(1:2)//OFT//".bgl.phased"
!  call system(r_BEAGLE)
  write(15,*) "gzip -d imputed.To_BEAGLE"//CHNO(1:2)//OFT//".bgl.phased.gz"
enddo
print*, "Done with Outlier Elimination"
print*, "Run BEAGLE ex) java -Xmx1000m -jar beagle.jar unphased=To_BEAGLE01"//OFT//".bgl missing=? out=imputed "
print*, "Parameter file for 3_Post_Imputation is ",PDD
print*, "Parameter file name for 3_Post_Imputation is ","PIP"//OFT//".par"

end program GDP_L1


