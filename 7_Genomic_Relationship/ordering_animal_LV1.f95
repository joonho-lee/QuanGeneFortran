! Name : Making Genomic Relationship Matrix 
! Developed by : Joonho Lee - PhD working in Hankyong National Univ.(Korea) zoonolee@gmail.com, +82-10-3408-2895
! Version : Linux 1.0 (2012-04-30)
! Contents : 		1. Ordering Animal information for Genomic Relationship Matrix

program ordering_animal
integer::i,j,l,mk,an,ErrorFlag !iostate & do_loop
integer::NOA=0,ANO=0,NOM=0,MNO=0 !number of animals and markers & animal and marekr number
integer,allocatable::MT(:,:),M(:,:),IDOD(:) !data storage

real(10)::GF1=0, GF2=0, SEH=0, S=0  !gnen frequencies & expected heterozygosity
real(10)::SOIBC=0,SSOIBC=0,AOIBC=0,SDOIBC=0,MINOIBC=0,MAXOIBC=0 !statistics for population genetics
real,allocatable::GRM(:,:),IGRM(:,:),W(:,:),IBC(:),k(:) !output data

character(len=8)::OFT,NOAF
character(len=16)::WF1,WF2,WF3
character(len=50),allocatable::dummyV(:)
character(len=300)::Dummy
character(len=220)::ADD,ODD,GDD,PDD
character(len=30)::GDFN1,GDFN2,ADFN,PDFN1,PDFN2
character(len=250)::ADFNIP,GDFNIP,ODFNIP,PDFNIP,PF
character(len=20),allocatable::ID(:)

!reading parameters
print *, "***********************************************************************************************************************" 
print *, "MGRM_Lv2 : Making Genomic Relationship Matrix & Basic Statistics Analysis Linux version 2.0"
print *, "Developed by Joonho Lee : zoonolee@gmail.com, +82-10-3408-2895"
print *, "Data form is the results of Post_fastPHASE and of converted genotype 0,1,2"
print *, "Please enter the name of parameter file including path(ex) ./parGRDH"
print *, "***********************************************************************************************************************"
read *, PF
PF=adjustl(PF)
open (99, file=PF) 
read (99,*) OFT
print *, "file tag ", OFT
read (99,*) NOA
rewind(unit=99)
read (99,*) OFT
read (99,*) NOAF
print *, "number of animals", NOA
read (99,*) NOM
print *, "number of markers", NOM
read (99,*) GDD
print *, "directory of Post_Imputation output file ", GDD
read (99,*) GDFN1
print *, "name of allelic coding data file(example : Genotype012_OUT_TAG.out) ", GDFN1
read (99,*) GDFN2
print *, "name of allele frequency data file(example : Genefrequency_OUT_TAG.out) ", GDFN2
read (99,*) PDD
print *, "directory of renum output file ", PDD
read (99,*) PDFN1
print *, "name of renumbered phenotype data file(example : OUTPUT31) ", PDFN1
read (99,*) PDFN2
print *, "name of renumbered animal data file(example : OUTPUT41) ", PDFN2
read (99,*) ACNO
print *, "column number of animal ", ACNO
allocate(dummyV(ACNO-1))
read (99,*) ADD
print *, "directory of Selected animals file ", ADD
read (99,*) ADFN
print *, "name of Selected animals file(example : Selected_Animals_OUT_TAG.out) ", ADFN
read (99,*) ODD
print *, "directory of output file ", ODD
close(99)

allocate(MT(NOM,NOA))
allocate(IDOD(NOA))
allocate(M(NOA,NOM))

WF1='            ';WF2='            ';WF3='            ';
write(WF1,fmt='(a16)') '('//adjustr(NOAF)//'f15.11)'
write(WF2,fmt='(a16)') '('//adjustr(NOAF)//'f15.8 )'
write(WF3,fmt='(a16)') '('//adjustr(NOAF)//'i1    )'
print *, WF1,WF2,WF3
ODD=adjustr(ODD)
GDD=adjustr(GDD)
ADD=adjustr(ADD)
PDD=adjustr(PDD)
GDFNIP=GDD//GDFN1
GDFNIP=adjustl(GDFNIP)
open(98,file=GDFNIP,status='old')
GDFNIP=GDD//GDFN2
GDFNIP=adjustl(GDFNIP)
open(97,file=GDFNIP,status='old')
ADFNIP=ADD//ADFN
ADFNIP=adjustl(ADFNIP)
open(96,file=ADFNIP,status='old')
ODFNIP=ODD//"GRM"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(95,file=ODFNIP,status='unknown')
ODFNIP=ODD//"IGRM"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(94,file=ODFNIP,status='unknown')
ODFNIP=ODD//"inbreeding"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(93,file=ODFNIP,status='unknown')
write(93,fmt='(2a20)') "Animal ID","Inbreeding coeff"
PDFNIP=PDD//PDFN1
PDFNIP=adjustl(PDFNIP)
open(92,file=PDFNIP,status='old')
PDFNIP=PDD//PDFN2
PDFNIP=adjustl(PDFNIP)
open(91,file=PDFNIP,status='old')

read(96,*)
do mk=1, NOM
  read(98,fmt=WF3) MT(mk,:)
enddo
M=transpose(MT)

do j=1, NOA
  read(96,*) ID(j),Dummy,Dummy
enddo

read(92,*)
do ANO=1, NOA
  read(92,*) DummyV, IDNO(:)
enddo

 



print *, "reading genotype is finished"

do i=1,NOA
  do j=1,NOA
    S=0
    do l=1,NOM
      S=S+W(MT(l,i)+1,i)*W(MT(l,j)+1,l)
    enddo
    GRM(i,j)=S/sqrt(k(i)*k(j))
    GRM(j,i)=GRM(i,j)
  enddo
enddo

print *, "making genomic relationship matrix is finished "
CALL FindInv(GRM, IGRM, NOA, ErrorFlag)
print *, "making inverse genomic relationship matrix is finished "

do j=1, NOA
  read(96,*) ID(j),Dummy,Dummy
  IBC(j)=GRM(j,j)-1.0
  SSOIBC=SSOIBC+IBC(j)**2
  write(93,fmt='(a20,f20.10)') ID(j),IBC(j)
  write(95,fmt=WF1) GRM(j,:)
  write(94,fmt=WF2) IGRM(j,:)
enddo
AOIBC=sum(IBC)/NOA;SDIBC=(SSOIBC-(sum(IBC)**2/NOA))/NOA
MINOIBC=minval(IBC);MAXOIBC=maxval(IBC)
write(93,fmt='(4a20)') "AOIBC","SDIBC","MINOIBC","MAXOIBC"
write(93,fmt='(4f20.10)') AOIBC,SDIBC,MINOIBC,MAXOIBC

print *, "estimating inbreeding coefficient is finished "

end program ordering_animal


