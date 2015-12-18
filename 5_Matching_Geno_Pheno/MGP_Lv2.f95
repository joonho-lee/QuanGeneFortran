! Name : Matching Genotype data with Phenotype data 
! Developed by : Joonho Lee - PhD working in Hankyong National Univ.(Korea) zoonolee@gmail.com, +82-10-3408-2895
! Version : Linux 1.0 (2012-04-02)
! Contents : 		1. Matching Genotype data with Phenotype data
!		 	2. Elimination of genotype data for animals with no phenotype 

program MGP 
integer::CNO,MNO1,MNO2,NEJ,iter,DummyI,iz,i1,i2,NONEG !iostate & do_loop
integer::NOA=0,ANO=0,NOM=0,MNO=0,NOR=0,RNO=0,NOTR=0,TRNO=0,ENO=0,NOE=0,MNODPR=0,TNODPR=0 !numbers and counting variables
integer,allocatable::G012DS(:,:),TG012DS(:,:),EDS(:,:),MNOMPR(:),MNOOPR(:),MNOPRBT(:),TNOMPR(:),TNOOPR(:),TNOPRBT(:) !data storage
real::MV
real,allocatable::PRDS(:,:),IOL(:,:)
real,allocatable::MSOPRBT(:),MAOPRBT(:),MSSOPRBT(:),MSDOPRBT(:),MMINOPRBT(:),MMAXOPRBT(:) !output data
real,allocatable::TSOPRBT(:),TAOPRBT(:),TSSOPRBT(:),TSDOPRBT(:),TMINOPRBT(:),TMAXOPRBT(:) !output data
character(len=2)::CHNO
character(len=3),allocatable::NAOTR(:),NAOE(:)
character(len=50)::MN
character(len=8)::OFT,NOAF,NOMF,NOEF,NOTRF
character(len=16)::WF1,WF2,WF3,WF4
character(len=27)::WF5
character(len=300)::Dummy
character(len=220)::MDD,ODD,GDD,PDD
character(len=30)::GDFN,PDFN,MDFN
character(len=250)::MDFNIP,GDFNIP,ODFNIP,PDFNIP,PF
character(len=20)::PID
character(len=20),allocatable::IDG(:),IDP(:)

!reading parameters
print *, "***********************************************************************************************************************" 
print *, "MGP_Lv1 : Matching Genotype data with Phenotype data Linux version 1.0"
print *, "Developed by Joonho Lee : zoonolee@gmail.com, +82-10-3408-2895"
print *, "Genotype data form is the results of Post_fastPHASE and of converted genotype 0,1,2"
print *, "Phenotype data form is animalID, Fixedeffect...,Trait..."
print *, "Fixed effect must be numbered as integer"
print *, "Trait name must be consisted of 3 character ex) EMA,BFT,..."
print *, "Please enter the name of parameter file including path(ex) ./parGRDH"
print *, "***********************************************************************************************************************"
read *, PF
PF=adjustl(PF)
open (99, file=PF) 
read (99,*) OFT
print *, "file tag ", OFT
read (99,*) NOA
read (99,*) NOM
read (99,*) NOE
if(NOE .NE. 0) then 
  read (99,*)
endif
read (99,*) NOTR
rewind(unit=99)
read (99,*) OFT
read (99,*) NOAF
print *, "number of animals", NOA
read (99,*) NOMF
print *, "number of markers", NOM
allocate(G012DS(NOM,NOA),TG012DS(NOA,NOM))
allocate(IDG(NOA))
read (99,*) NOEF
print *, "number of effects", NOE
if(NOE .NE. 0) then 
  print *, "name of fixed effects "
  allocate(NAOE(NOE))
  read (99,*) NAOE(:)
  do ENO=1,NOE
    print *, "Effect",ENO,":",NAOE(ENO)
  enddo
endif
read (99,*) NOTRF
print *, "number of traits", NOTR
allocate(NAOTR(NOTR))
allocate(IOL(NOTR,2))
allocate(MNOMPR(NOTR),MNOOPR(NOTR),MNOPRBT(NOTR))
allocate(TNOMPR(NOTR),TNOOPR(NOTR),TNOPRBT(NOTR))
allocate(MSOPRBT(NOTR),MAOPRBT(NOTR),MSSOPRBT(NOTR),MSDOPRBT(NOTR),MMINOPRBT(NOTR),MMAXOPRBT(NOTR))
allocate(TSOPRBT(NOTR),TAOPRBT(NOTR),TSSOPRBT(NOTR),TSDOPRBT(NOTR),TMINOPRBT(NOTR),TMAXOPRBT(NOTR))
read (99,*) NAOTR(:)
print *, "name of traits and outlier ranges"
do TRNO=1,NOTR
  read (99,*) IOL(TRNO,:)
  print *, NAOTR(TRNO), "(",IOL(TRNO,1)," ~ ",IOL(TRNO,2)," )" 
enddo
read (99,*) MV
print *, "missing value ", MV
read (99,*) GDD
print *, "directory of Post_Imputation output file ", GDD
read (99,*) GDFN
print *, "name of allelic coding data file(example : Genotype012_OUT_TAG.out) ", GDFN
read (99,*) MDD
print *, "directory of Selected animals file ", MDD
read (99,*) MDFN
print *, "name of Selected animals file(example : Selected_Animals_OUT_TAG.out) ", MDFN
read (99,*) PDD
print *, "directory of Phenotype file ", PDD
read (99,*) PDFN
print *, "name of Phenotype file(example : Phenotype_OUT_TAG.out) ", PDFN
read (99,*) ODD
print *, "directory of output file ", ODD
close(99)

PDD=adjustr(PDD)
ODD=adjustr(ODD)
GDD=adjustr(GDD)
MDD=adjustr(MDD)

WF1='            '
write(WF1,fmt='(a16)') '('//adjustr(NOAF)//'i1    )'
print *, "reading format for Selected markers file is : ",WF1
WF2='            '
write(WF2,fmt='(a16)') '('//adjustr(NOMF)//'i1    )'
print *, "writing format for Selected markers file is : ",WF2
WF3='            '
if(NOE.NE.0)  write(WF3,fmt='(a16)') '('//adjustr(NOEF)//'i10   )'
print *, "writing format for effect file is : ",WF3
WF4='            '
write(WF4,fmt='(a16)') '('//adjustr(NOTRF)//'f10.3 )'
print *, "writing format for trait observations is : ",WF4
WF5='                           '
write(WF5,fmt='(a27)') '('//adjustr(NOMF)//'i1,'//adjustr(NOTRF)//'f10.3 )'
print *, "writing format for whole data is : ",WF5
GDFNIP=GDD//GDFN
GDFNIP=adjustl(GDFNIP)
open(98,file=GDFNIP,status='old')
PDFNIP=PDD//PDFN
PDFNIP=adjustl(PDFNIP)
open(97,file=PDFNIP,status='old')
MDFNIP=MDD//MDFN
MDFNIP=adjustl(MDFNIP)
open(96,file=MDFNIP,status='old')
ODFNIP=ODD//"Whole_Data"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(95,file=ODFNIP,status='unknown')
ODFNIP=ODD//"MGP_Stat"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(94,file=ODFNIP,status='unknown')
do TRNO=1, NOTR
  ODFNIP=ODD//"GT"//NAOTR(TRNO)//OFT//".out"
  ODFNIP=adjustl(ODFNIP)
  open(TRNO+10,file=ODFNIP,status='unknown')
  ODFNIP=ODD//"PR"//NAOTR(TRNO)//OFT//".out"
  ODFNIP=adjustl(ODFNIP)
  open(TRNO+30,file=ODFNIP,status='unknown')
  if(NOE.NE.0) then
    ODFNIP=ODD//"FE"//NAOTR(TRNO)//OFT//".out"
    ODFNIP=adjustl(ODFNIP)
    open(TRNO+50,file=ODFNIP,status='unknown')
  endif
  ODFNIP=ODD//"ID"//NAOTR(TRNO)//OFT//".out"
  ODFNIP=adjustl(ODFNIP)
  open(TRNO+70,file=ODFNIP,status='unknown')
enddo

do MNO=1, NOM
  read(98,fmt=WF1) G012DS(MNO,:)
enddo
TG012DS=transpose(G012DS)
print *, "reading genotype data is finished"
read(96,*)
do ANO=1, NOA
  read(96,*) IDG(ANO), DummyI, DummyI
enddo
do 
  read(97,*,iostat=io)
  if(io/=0) exit
  NOR=NOR+1
enddo
rewind(unit=97)
print *, "number of records", NOR
write(94,*) "number of genotype records  : ", NOA
write(94,*) "number of phenotype records : ", NOR
write(94,*) " "
write(94,*) "*** statistics for total phenotype records by traits ***"
write(94,fmt='(8a15)') "Trait","reco_P","miss_P","outl_P","mean_P","vari_P","minV_P","maxV_P"
if(NOE.NE.0) allocate(EDS(NOR,NOE))
allocate(PRDS(NOR,NOTR))
allocate(IDP(NOR))
if(NOE.EQ.0) then
  do RNO=1,NOR
    read(97,*) IDP(RNO),PRDS(RNO,:)
  enddo
else
  do RNO=1,NOR
    read(97,*) IDP(RNO),EDS(RNO,:),PRDS(RNO,:)
  enddo
endif
print *, "reading phenotype data is finished"

MNOMPR=0;MNOOPR=0;MNOPRBT=0;
TNOMPR=0;TNOOPR=0;TNOPRBT=0;
MSOPRBT=0;MAOPRBT=0;MSSOPRBT=0;MSDOPRBT=0;MMINOPRBT=0;MMAXOPRBT=0
TSOPRBT=0;TAOPRBT=0;TSSOPRBT=0;TSDOPRBT=0;TMINOPRBT=0;TMAXOPRBT=0
do TRNO=1,NOTR
  do ANO=1,NOA
    do RNO=1,NOR
      if(ANO.eq.1.and.PRDS(RNO,TRNO).ne.MV.and.PRDS(RNO,TRNO).gt.IOL(TRNO,1).and.PRDS(RNO,TRNO).lt.IOL(TRNO,2)) then
        TNOPRBT(TRNO)=TNOPRBT(TRNO)+1
        if(TNOPRBT(TRNO).eq.1) then
          TMINOPRBT(TRNO)=PRDS(RNO,TRNO)
          TMAXOPRBT(TRNO)=PRDS(RNO,TRNO)
        endif
        TSOPRBT(TRNO)=TSOPRBT(TRNO)+PRDS(RNO,TRNO)
        TSSOPRBT(TRNO)=TSSOPRBT(TRNO)+(PRDS(RNO,TRNO))**2
        TMINOPRBT(TRNO)=min(TMINOPRBT(TRNO),PRDS(RNO,TRNO))
        TMAXOPRBT(TRNO)=max(TMAXOPRBT(TRNO),PRDS(RNO,TRNO))
      endif
      if(ANO.eq.1.and.PRDS(RNO,TRNO).eq.MV) TNOMPR(TRNO)=TNOMPR(TRNO)+1
      if(ANO.eq.1.and.PRDS(RNO,TRNO).ne.MV.and.PRDS(RNO,TRNO).le.IOL(TRNO,1)) TNOOPR(TRNO)=TNOOPR(TRNO)+1
      if(ANO.eq.1.and.PRDS(RNO,TRNO).ne.MV.and.PRDS(RNO,TRNO).ge.IOL(TRNO,2)) TNOOPR(TRNO)=TNOOPR(TRNO)+1
      if(TRNO.eq.1.and.IDG(ANO).eq.IDP(RNO)) write(95,WF5) TG012DS(ANO,:),PRDS(RNO,:)
      if(IDG(ANO).eq.IDP(RNO).and.PRDS(RNO,TRNO).ne.MV.and.PRDS(RNO,TRNO).gt.IOL(TRNO,1).and.PRDS(RNO,TRNO).lt.IOL(TRNO,2)) then
        if(TRNO.eq.1.and.IDG(ANO).eq.PID) then 
          MNODPR=MNODPR+1
          print *, "Animal ",IDP(RNO),"has duplicated phenotypes"  
        endif
        PID=IDG(ANO)
        MNOPRBT(TRNO)=MNOPRBT(TRNO)+1
        if(MNOPRBT(TRNO).eq.1) then
          MMINOPRBT(TRNO)=PRDS(RNO,TRNO)
          MMAXOPRBT(TRNO)=PRDS(RNO,TRNO)
        endif
        MSOPRBT(TRNO)=MSOPRBT(TRNO)+PRDS(RNO,TRNO)
        MSSOPRBT(TRNO)=MSSOPRBT(TRNO)+(PRDS(RNO,TRNO))**2
        MMINOPRBT(TRNO)=min(MMINOPRBT(TRNO),PRDS(RNO,TRNO))
        MMAXOPRBT(TRNO)=max(MMAXOPRBT(TRNO),PRDS(RNO,TRNO))
        write(TRNO+10,WF2) TG012DS(ANO,:)
        write(TRNO+30,WF4) PRDS(RNO,TRNO)
        if(NOE.NE.0) write(TRNO+50,WF3) EDS(RNO,:)
        write(TRNO+70,fmt='(2a20)') IDG(ANO),IDP(RNO)
      endif
      if(IDG(ANO).eq.IDP(RNO).and.PRDS(RNO,TRNO).eq.MV) MNOMPR(TRNO)=MNOMPR(TRNO)+1
      if(IDG(ANO).eq.IDP(RNO).and.PRDS(RNO,TRNO).ne.MV.and.PRDS(RNO,TRNO).le.IOL(TRNO,1)) MNOOPR(TRNO)=MNOOPR(TRNO)+1
      if(IDG(ANO).eq.IDP(RNO).and.PRDS(RNO,TRNO).ne.MV.and.PRDS(RNO,TRNO).ge.IOL(TRNO,2)) MNOOPR(TRNO)=MNOOPR(TRNO)+1
    enddo
  enddo
  MAOPRBT(TRNO)=MSOPRBT(TRNO)/real(MNOPRBT(TRNO))
  MSDOPRBT(TRNO)=(MSSOPRBT(TRNO)-MSOPRBT(TRNO)**2/real(MNOPRBT(TRNO)))/real(MNOPRBT(TRNO))
  TAOPRBT(TRNO)=TSOPRBT(TRNO)/real(TNOPRBT(TRNO))
  TSDOPRBT(TRNO)=(TSSOPRBT(TRNO)-TSOPRBT(TRNO)**2/real(TNOPRBT(TRNO)))/real(TNOPRBT(TRNO))
  write(94,fmt='(a15,3i15,4f15.3)') NAOTR(TRNO),TNOPRBT(TRNO),TNOMPR(TRNO),TNOOPR(TRNO),&
  &TAOPRBT(TRNO),TSDOPRBT(TRNO),TMINOPRBT(TRNO),TMAXOPRBT(TRNO)
enddo
write(94,*) " "
write(94,*) "*** statistics for matched phenotype records with genotype data by traits ***"
write(94,fmt='(8a15)') "Trait","reco_M","miss_M","outl_M","mean_M","vari_M","minV_M","maxV_M"
do TRNO=1,NOTR
  write(94,fmt='(a15,3i15,4f15.3)') NAOTR(TRNO),MNOPRBT(TRNO),MNOMPR(TRNO),MNOOPR(TRNO),&
  &MAOPRBT(TRNO),MSDOPRBT(TRNO),MMINOPRBT(TRNO),MMAXOPRBT(TRNO)
enddo
write(94,*) "number of duplicated phenotype records(in matched phenotype records with genotype data", MNODPR
end program MGP

