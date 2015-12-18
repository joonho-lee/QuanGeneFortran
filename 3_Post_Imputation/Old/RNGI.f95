! Name : Reformation of Imputation Results & Basic Statistics Analysis 
! Developed by : Joonho Lee - PhD working in Hankyong National Univ.(Korea) zoonolee@gmail.com, +82-10-3408-2895
! Version : Linux 1.0 (2011-10-25)
! Contents : 		1. Data reading(animal ID, SNP marker name, chromosome number, marker position, phenotype(name, records))
!		 	2. Calculation of gene & genotype frequency
!			3. Unusual data elimination
!			4. Statistics of data & population genetics 
!			5. Transformation of genotype data for fastPHASE 

program OEGRDPF 
integer::io,i,j,k,l,m,n,o,p=0,q=0,iz,i1,i2 !iostate & do_loop
integer::NOCEXC=0,XCN=19,YCN=20,UCN=99,ANO=0 !number of chromosomes & chromosome number for X, Y, Unknown
integer::BC(5)
integer::GC(11)
integer::NOA,NOM,NOMIH=0,NOHE=0,NOMAH=0,NOMG=0,NOMOXC=0,NOMOYC=0,NOMOUC=0,TNOMM=0 !counting
integer::NOAMM=0,NONMM=0,NODM=0,NODA=0,NOMWMRO=0,NOMWMAFO=0,NOMWHCO=0,NOAWMRO=0 !counting
integer::NOMWMAF1=0,NOMWMAF2=0,NOMWMAF3=0,NOMWMAF4=0,NOMWMAF5=0,CMP=0,EXMP=0,CN=0,EXCN=999, SODBAM=0,NOMBCN=0 ! for statistics
integer::NODMBOFOMGT=0,NODMBOAH=0,NODMBOFOMIG=0,NODMBOCSOHWT=0,NODMBONP=0,NODMBOAHE=0,EXCNO=0 
integer,allocatable::IOM(:,:),GN(:,:),NOMGBA(:),NOMBC(:) !data input & missing data output
real::GF(5),GTF(11),FOMAH,FOMIH,FOHE,FOMAG,FOMIG,SOGF,SOGTF,FOMG,FOMGT,IOL(5),SOPIC=0,AOPIC=0 !frequencies
real::HWD=0,ENOMIH=0,ENOHE=0,ENOMAH=0,CSOHWT=0,EH=0,PIC=0,TMGR=0,ADOAMBCN=0,AHWD=0,AHE=0,AEH=0,SOHWD=0,SOHE=0,SOEH=0 !statistics for population genetics
real::SOIBC=0,AIBC=0,IBC=0 !statistics for population genetics
character(len=16)::MN, D1, D2, D3
character(len=1)::BA(5)=(/'A','T','G','C','Z'/),MABA,MIBA,MG
character(len=2)::GT(11)=(/'AA','AT','AG','AC','TT','TG','TC','GG','GC','CC','ZZ'/),MAH,MIH,HE,CHNO
character(len=300)::GDD,ODD
character(len=30)::GDFN,WF,ODFN
character(len=15)::MWF(5)
character(len=330)::GDFNIP,ODFNIP,PF
character(len=20),allocatable::ID(:)
character(len=2),allocatable::GTS(:)
character(len=1),allocatable::GTDAEOM(:,:,:)

print *, "***********************************************************************************************************************" 
print *, "OEGRDPF_Lv1 : Outlier Elimination of Genomic Raw Data & Preparation for fastPHASE Linux version 3.0"
print *, "Developed by Joonho Lee : zoonolee@gmail.com, +82-10-3408-2895"
print *, "Data format : (1st row) Markername Chromosome Position LL09066048	LL08131059	LL09147068 ... (Animal IDs)"
print *, "Data format : (2nd row) ALGA0000009	1	52297	AG	GG	GG	... (Genotypes)"
print *, "Input data must be sorted by chromosome No.(1st) and Marker position(2nd)"
print *, "Sex chromosome and markers of unknown chromosome must be numbered more then the biggest chromosome number"
print *, "Please enter the name of parameter file including path(ex) ./parGRDH"
print *, "***********************************************************************************************************************"
read *, PF
PF=adjustl(PF)
open (99, file=PF) 
read (99,*) NOA
read (99,*) GDD
read (99,*) GDFN
read (99,*) ODD
read (99,*) GT(11)
read (99,*) IOL(:)
read (99,*) NOCEXC
read (99,*) NOMOTBC
read (99,*) XCN, YCN, UCN
read (99,*) DXCN, DYCN, DUCN

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
print *, "   markers with no position(if you want to delete unknown position marker;enter 0) : ", IOL(5)
print *, "X chromosome No. , Y chromosome No. Unknown chromosome No. ", XCN, YCN, UCN

GDFN=adjustl(GDFN)
GDD=adjustr(GDD)
GDFNIP=GDD//GDFN
GDFNIP=adjustl(GDFNIP)

print *, "input genomic data file is ",GDFNIP 
print *, "Output file will be created in",ODD 
ODD=adjustr(ODD)
close(99)

do i=1,NOCEXC
  open(i,file=GDFNIP,status='unknown')
ODFN="Delected_Markers.dat"
ODFN=adjustl(ODFN)
ODFNIP=ODD//ODFN
ODFNIP=adjustl(ODFNIP)
open(11,file=ODFNIP,status='replace')
ODFNIP=ODD//"Selected_Markers.dat"
ODFNIP=adjustl(ODFNIP)
write(11,fmt='(a20,6a10)') "Marker_Name","FOMGT","FOMIG","CSOHWT","position","FOHE","chno" 
open(12,file=ODFNIP,status='unknown')
ODFNIP=ODD//"Delected_Animals.dat"
ODFNIP=adjustl(ODFNIP)
open(13,file=ODFNIP,status='unknown')
ODFNIP=ODD//"Selected_Animals.dat"
ODFNIP=adjustl(ODFNIP)
open(14,file=ODFNIP,status='unknown')
do i=1,NOCEXC
  iz=ichar('0')
  i1=floor(i/10.)+iz
  i2=mod(i,10)+iz
  CHNO(1:1)=achar(i1)
  CHNO(2:2)=achar(i2)
  ODFNIP=ODD//"To_fastPHASE"//CHNO(1:2)//".inp"
  ODFNIP=adjustl(ODFNIP)
  open(i+20,file=ODFNIP,status='unknown')
enddo
ODFNIP=ODD//"Result.dat"
ODFNIP=adjustl(ODFNIP)
open(16,file=ODFNIP,status='unknown')
ODFNIP=ODD//"Statistics.dat"
ODFNIP=adjustl(ODFNIP)
open(17,file=ODFNIP,status='unknown')
write(17,fmt='(a12,4a10)') "          ","EXCN","NOMBCN","ADOAMBCN","AEH","AOPIC"
NOM=0
do 
 read(1,*,iostat=io)
  if(io/=0) exit
 NOM=NOM+1
enddo
rewind(unit=1)
NOM=NOM-1
allocate(IOM(NOM,2),GN(NOM,NOA),GTDAEOM(NOCEXC,NOA*2,NOMOTBC))
allocate(ID(NOA),GTS(NOA),NOMGBA(NOA),NOMBC(NOCEXC))
NOMBC=0
IOM=0
GN=0
NOMGBA=0
NODA=0
GTDAEOM=' '
read(1,*)D1,D2,D3,ID
MG=GT(11)(1:1)
NODM=0
do i=1,NOM
    GTS='  '
    GF=0.0
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
              TNOMM=TNOMM+1
    end select
    enddo
    GF=real(BC)/real(sum(BC(:4)))
    GTF=real(GC)/real(sum(GC(:10)))
    FOMG=real(BC(5))/real(NOA*2)
    FOMGT=real(GC(11))/(real(NOA))
    FOMIH=0.0; FOHE=0.0; FOMAH=0.0; FOMIG=0.0; FOMAG=0.0; 
    MIH='  ';MAH='  ';HE='  ';MIBA=' ';MABA=' '
    do k=1,4
        if(GF(k).NE.0 .AND. GF(k)<0.5) then
          MIBA=BA(k)
          FOMIG=GF(k)
        elseif(GF(k).NE.0 .AND. GF(k)>0.5) then
          MABA=BA(k)
          FOMAG=GF(k)
          nMA=BC(k)
        elseif(GF(k)==0.5) then
          if(sum(GF(:k)).eq.0.5) then
            MIBA=BA(k)
            FOMIG=GF(k)
          else 
            MABA=BA(k)
            FOMAG=GF(k)
          endif
        endif
    enddo
    MIH=MIBA//MIBA
    MAH=MABA//MABA
    HE=MIBA//MABA
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
  ENOMIH=(FOMIG**2)*(NOMIH+NOHE+NOMAH)
  ENOHE=(2*FOMIG*FOMAG)*(NOMIH+NOHE+NOMAH)
  ENOMAH=(FOMAG**2)*(NOMIH+NOHE+NOMAH)
  CSOHWT=((NOMIH-ENOMIH)**2/ENOMIH)+((NOHE-ENOHE)**2/ENOHE)+((NOMAH-ENOMAH)**2/ENOMAH)
  EH=1-(FOMIG**2+FOMAG**2)
  PIC=1-(FOMIG**2+FOMAG**2)-2*((FOMIG**2)*(FOMAG**2))
  IBC=1-(real(NOHE)/ENOHE)
  if(FOMGT>IOL(1)) NODMBOFOMGT=NODMBOFOMGT+1
  if(FOMAH.eq.1) NODMBOAH=NODMBOAH+1
  if(FOMIG<IOL(2)) NODMBOFOMIG=NODMBOFOMIG+1
  if(CSOHWT>IOL(3)) NODMBOCSOHWT=NODMBOCSOHWT+1
  if(IOM(i,2).le.IOL(5)) NODMBONP=NODMBONP+1
  if(FOHE.EQ.1) NODMBOAHE=NODMBOAHE+1
  if(IOM(i,1) .eq. XCN) NOMOXC=NOMOXC+1
  if(IOM(i,1) .eq. YCN) NOMOYC=NOMOYC+1
  if(IOM(i,1) .eq. UCN) NOMOUC=NOMOUC+1
  if(FOMGT>IOL(1).or.FOMIG<IOL(2).or.CSOHWT>IOL(3).or.IOM(i,2).le.IOL(5).or.FOHE.EQ.1.or.IOM(i,1).eq.DXCN.or.IOM(i,1).eq.DYCN&
  &.or.IOM(i,1).eq.DUCN) then 
    NODM=NODM+1
    write(11,fmt='(a20,3f10.7,i10,f10.7,i10)') MN,FOMGT,FOMIG,CSOHWT,IOM(i,2),FOHE,IOM(i,1) 
  else
    if(IOM(i,1).ne.EXCNO) then
      p=0
    endif
      EXCNO=IOM(i,1)
    NOMBC(IOM(i,1))=NOMBC(IOM(i,1))+1
    p=p+1
    write(12,fmt='(a20)') MN
    do l=1,NOA
      GTDAEOM(IOM(i,1),l*2-1,p)=GTS(l)(1:1)
      GTDAEOM(IOM(i,1),l*2,p)=GTS(l)(2:2)
      if(GTS(l).eq.'??') NOMGBA(l)=NOMGBA(l)+1 
    enddo
    if(FOMIG .LT. 0.1) then 
      NOMWMAF1=NOMWMAF1+1
    elseif(0.1 .GE. FOMIG .and. FOMIG .LT. 0.2) then 
      NOMWMAF2=NOMWMAF2+1
    elseif(0.2 .GE. FOMIG .and. FOMIG .LT. 0.3) then 
      NOMWMAF3=NOMWMAF3+1    
    elseif(0.3 .GE. FOMIG .and. FOMIG .LT. 0.4) then 
      NOMWMAF4=NOMWMAF4+1
    elseif(0.4 .GE. FOMIG) then 
      NOMWMAF5=NOMWMAF5+1
    endif
  endif
  if(FOMGT.eq.1) NOAMM=NOAMM+1
  if(FOMGT.eq.0) NONMM=NONMM+1
  if(EXCN.ne.XCN.and.EXCN.ne.YCN.and.EXCN.ne.UCN.and.EXCN.ne.999.and.IOM(i,1).ne.EXCN) then 
    ADOAMBCN=real(SODBAM)/real(NOMBCN)
    AEH=real(SOHE)/real(NOMBCN)
    AOPIC=real(SOPIC)/real(NOMBCN)
    write(17,fmt='(a10,i2,i10,f10.2,2f10.7)') "chromosome",EXCN,NOMBCN,ADOAMBCN,AEH,AOPIC
  end if
  if(EXCN.ne.DXCN.and.EXCN.ne.DYCN.and.EXCN.ne.DUCN.and.IOM(i,1).ne.EXCN) then 
    CN=IOM(i,1)
    NOMBCN=1
    SODBAM=0
    ADOAMBCN=0
    EXMP=0
    SOHWD=HWD
    SOHE=EH
    SOPIC=PIC
    SOIBC=IBC
    CMP=0
  else
    SOHWD=SOHWD+HWD
    SOHE=SOHE+EH
    SOPIC=SOPIC+PIC
    if(EXMP.ne.0) CMP=IOM(i,2)-EXMP
    NOMBCN=NOMBCN+1
    SODBAM=SODBAM+CMP
  end if
  EXMP=IOM(i,2)
  EXCN=CN
enddo
print *, "Done with the outlier marker elimination : ",NODM,"markers were deleted" 
  do m=1,NOA
    if(real(NOMGBA(m)/(NOM-NODM)).GT.IOL(4)) then 
      NODA=NODA+1
      write(13,fmt='(a20,i10,f10.7)') ID(m),NOMGBA(m),real(NOMGBA(m)/(NOM-NODM))
    else
      write(14,fmt='(a20)') ID(m)
    endif
  enddo
do o=1,NOCEXC
  if(NOA-NODA.lt.10) then
    write(MWF(3),fmt='(a4)') '(a1)'
  elseif(NOA-NODA.ge.10.and.NOA-NODA.lt.100) then
    write(MWF(3),fmt='(a4)') '(a2)'
  elseif(NOA-NODA.ge.100.and.NOA-NODA.lt.1000) then
    write(MWF(3),fmt='(a4)') '(a3)'
  elseif(NOA-NODA.ge.1000.and.NOA-NODA.lt.10000) then
    write(MWF(3),fmt='(a4)') '(a4)'
  elseif(NOA-NODA.ge.10000.and.NOA-NODA.lt.100000) then
    write(MWF(3),fmt='(a4)') '(a5)'
  endif
  write(MWF(1),fmt='(i10)') NOA-NODA
  MWF(1)=adjustl(MWF(1)) 
  write(o+20,fmt=MWF(3)) MWF(1)
  if(NOMBC(o).lt.10) then
    write(MWF(3),fmt='(a4)') '(a1)'
  elseif(NOMBC(o).ge.10.and.NOMBC(o).lt.100) then
    write(MWF(3),fmt='(a4)') '(a2)'
  elseif(NOMBC(o).ge.100.and.NOMBC(o).lt.1000) then
    write(MWF(3),fmt='(a4)') '(a3)'
  elseif(NOMBC(o).ge.1000.and.NOMBC(o).lt.10000) then
    write(MWF(3),fmt='(a4)') '(a4)'
  elseif(NOMBC(o).ge.10000.and.NOMBC(o).lt.100000) then
    write(MWF(3),fmt='(a4)') '(a5)'
  elseif(NOMBC(o).ge.100000.and.NOMBC(o).lt.1000000) then
    write(MWF(3),fmt='(a4)') '(a6)'
  endif
  write(MWF(2),fmt='(i10)') NOMBC(o)
  MWF(2)=adjustl(MWF(2)) 
  write(o+20,fmt=MWF(3)) MWF(2)
  write(MWF(5),fmt='(a1,i10,a3)') '(',NOMBC(o),'a1)'
  MWF(5)=trim(MWF(5)) 
  MWF(5)=adjustl(MWF(5)) 
  print *, "Number of markers on chromosome No.",o,NOMBC(o),MWF(5)
  ANO=0
  do m=1,NOA
    if(real(NOMGBA(m)/NOM-NODM).LE.IOL(4)) then
      ANO=ANO+1
      write(MWF(4),fmt='(i10)') ANO
      MWF(4)=adjustl(MWF(4))
      if(ANO.lt.10) then
        write(MWF(3),fmt='(a10)') '(a4,1x,a1)'
      elseif(ANO.ge.10.and.ANO.lt.100) then
        write(MWF(3),fmt='(a10)') '(a4,1x,a2)'
      elseif(ANO.ge.100.and.ANO.lt.1000) then
        write(MWF(3),fmt='(a10)') '(a4,1x,a3)'
      elseif(ANO.ge.1000.and.ANO.lt.10000) then
        write(MWF(3),fmt='(a10)') '(a4,1x,a4)'
      elseif(ANO.ge.10000.and.ANO.lt.100000) then
        write(MWF(3),fmt='(a10)') '(a4,1x,a5)'
      endif
      write(o+20,fmt=MWF(3)) "# id",MWF(4)
      write(o+20,fmt=MWF(5)) GTDAEOM(o,m*2-1,:NOMBC(o))
      write(o+20,fmt=MWF(5)) GTDAEOM(o,m*2,:NOMBC(o))
    endif
  enddo
enddo  
print *, "Done with the outlier animal elimination : ",NODA,"animals were deleted"
    write(16,*) "number of markers on chromosome   X  ",NOMOXC
    write(16,*) "number of markers on chromosome   Y  ",NOMOYC
    write(16,*) "number of markers on chromosome   Unknown  ",NOMOUC    
    write(16,*) "Total number of markers   ", NOM
    write(16,*) "All missing markers   ", NOAMM
    write(16,*) "All homo_genotype markers   ", NODMBOAH
    write(16,*) "All hetero_genotype markers   ", NODMBOAHE
    write(16,*) "Number of markers with no missing  ", NONMM         
    write(16,*) "Number of useful markers  ", NOM-NODM
    write(16,*) "Number of outlier markers  ", NODM
    write(16,*) "Number of markers with missing over ",IOL(1),"%   ", NODMBOFOMGT
    write(16,*) "Number of markers with minor allele frequency less than ",IOL(2),"   ", NODMBOFOMIG
    write(16,*) "Number of markers with HWchisquare over ",IOL(3),"   ", NODMBOCSOHWT
    write(16,*) "Number of animals with missing over ",IOL(4),"%   ", NODA
    write(16,*) "Number of animals with no position information ", NODMBONP
    TMGR=real(TNOMM)/real(NOM*NOA)
    write(16,*) "Total missing rate   " , TMGR
    write(16,*) "number of markers with Minor allel frequency(except for unusual data)"
    write(16,*) "       minor allele frequency < 0.1  :  ", NOMWMAF1
    write(16,*) "0.1 <= minor allele frequency < 0.2  :  ", NOMWMAF2
    write(16,*) "0.2 <= minor allele frequency < 0.4  :  ", NOMWMAF3
    write(16,*) "0.3 <= minor allele frequency < 0.4  :  ", NOMWMAF4
    write(16,*) "0.4 <= minor allele frequency        :  ", NOMWMAF5
end program OEGRDPF





