! Name : Extent of Linkage Disequilibrium 
! Developed by : Joonho Lee - PhD working in Hankyong National Univ.(Korea) zoonolee@gmail.com, +82-10-3408-2895
! Version : Linux 1.0 (2012-03-26)
! Contents : 		1. Calculation of Haplotype frequency
!		 	2. Estimation of Haplotype frequency using EM algorithm 
!		 	3. Calculation of Linkage Disequribrium
!			4. Estimation of Effective Population Size

program EOLD 
integer::CNO,MNO1,MNO2,NEJ,iter,DummyI,iz,i1,i2,NONEG !iostate & do_loop
integer::NOA=0,ANO=0,NOM=0,MNO=0,NOHTBC=0,CBPTCM=1000000 !number of animals and markers & animal and marekr number
integer::nMIHMIH,nMIHHE,nMIHMAH,nHEMIH,nHEHE,nHEMAH,nMAHMIH,nMAHHE,nMAHMAH 
integer,allocatable::G012DS(:,:),SEMNBC(:,:),MP(:),NOSMBC(:),NEGA(:,:),CNEG(:) !data storage
real::D=0,rsquare=0,CONV=0,DLFLDD=0,HG=0,DFNEG=0,EPS=0
real::nMIMI,nMIMA,nMAMI,nMAMA
real::MIMIf,MIMAf,MAMIf,MAMAf,newMIMIf,newMIMAf,newMAMIf,newMAMAf,prMIMIMAMA,prMIMAMAMI,diff=0.0
real(10)::SLDBC=0,ALDBC=0,SALDBC=0,AALDBC=0,STLD=0,ATLD=0,STALD=0,ATALD=0
real,allocatable::GF(:,:),NEG(:,:),SNEG(:),ANEG(:) !output data
character(len=2)::CHNO
character(len=50)::MN
character(len=8)::OFT,NOAF
character(len=16)::WF1
character(len=300)::Dummy
character(len=220)::MDD,ODD,GDD
character(len=30)::GDFN1,GDFN2,MDFN
character(len=250)::MDFNIP,GDFNIP,ODFNIP,PF
character(len=20),allocatable::ID(:)

!reading parameters
print *, "***********************************************************************************************************************" 
print *, "EOLD_Lv1 : Extent of Linkage Disequilibrium Linux version 1.0"
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
read (99,*) MDD
print *, "directory of Selected markers file ", MDD
read (99,*) MDFN
print *, "name of Selected markers file(example : Selected_markers_OUT_TAG.out) ", MDFN
read (99,*) ODD
print *, "directory of output file ", ODD
read (99,*) CONV
print *, "convergence of haplotype frequency estimation ", CONV
read (99,*) CBPTCM
print *, "1 centi morgan(cM) to base pair(bp) (example 1000000)) ", CBPTCM
read (99,*) DLFLDD
print *, "distance limit for LD decay (unit : cM) ", DLFLDD
read (99,*) NOAC
print *, "number of autosomal chromosomes : ", NOAC

allocate(SEMNBC(NOAC,2))
allocate(NOSMBC(NOAC))
allocate(G012DS(NOM,NOA))
allocate(MP(NOM))
allocate(GF(NOM,2))
ODD=adjustr(ODD)
GDD=adjustr(GDD)
MDD=adjustr(MDD)
MNO=0

do CNO=1,NOAC
  read (99,*) NOSMBC(CNO)
  SEMNBC(CNO,1)=MNO+1
  MNO=MNO+NOSMBC(CNO)
  SEMNBC(CNO,2)=MNO
  print *, "Marker nomber on chromosome No.",CNO," from ", SEMNBC(CNO,1)," to ", SEMNBC(CNO,2)  
  print *, "Number of markers on chromosome No.",CNO," : ", NOSMBC(CNO)
  iz=ichar('0')
  i1=floor(CNO/10.)+iz
  i2=mod(CNO,10)+iz
  CHNO(1:1)=achar(i1)
  CHNO(2:2)=achar(i2)
  print *, "Name of LD analysis result of chromosome No.",CNO," : ","LD_CH_"//CHNO(1:2)//OFT//".out"
  ODFNIP=ODD//"LD_CH_"//CHNO(1:2)//OFT//".out"
  ODFNIP=adjustl(ODFNIP)
  open(CNO+20,file=ODFNIP,status='unknown')
enddo
read (99,*) NONEG
print *, "number of groups for NE estimation : ", NONEG
allocate(NEG(NONEG,2))
allocate(NEGA(NONEG,2))
allocate(CNEG(NONEG))
allocate(SNEG(NONEG))
allocate(ANEG(NONEG))
do i=1,NONEG 
  read(99,*) NEGA(i,:)
  NEG(i,1)=1/real(NEGA(i,1)*2)*100;NEG(i,2)=1/real(NEGA(i,2)*2)*100
  print *, "effective population size from ",NEGA(i,1)," to ",NEGA(i,2)," generations ago "
  print *, "using average LD of distance between ",NEG(i,1)," and ",NEG(i,2)
enddo
close(99)
WF1='            '
write(WF1,fmt='(a16)') '('//adjustr(NOAF)//'i1    )'
print *, "reading format for Selected markers file is : ",WF1
GDFNIP=GDD//GDFN1
GDFNIP=adjustl(GDFNIP)
open(98,file=GDFNIP,status='old')
GDFNIP=GDD//GDFN2
GDFNIP=adjustl(GDFNIP)
open(97,file=GDFNIP,status='old')
MDFNIP=MDD//MDFN
MDFNIP=adjustl(MDFNIP)
open(96,file=MDFNIP,status='old')
ODFNIP=ODD//"LD_Decay"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(95,file=ODFNIP,status='unknown')
ODFNIP=ODD//"LD_STAT_BC"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(94,file=ODFNIP,status='unknown')
ODFNIP=ODD//"NE_estimation"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(93,file=ODFNIP,status='unknown')
ODFNIP=ODD//"LD_Adjacent"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(92,file=ODFNIP,status='unknown')

read(97,*) DummyI,DummyI
read(96,*) 

do MNO=1, NOM
  read(98,fmt=WF1) G012DS(MNO,:)
  read(97,*) GF(MNO,:)
  read(96,*) MN,DummyI,DummyI,MP(MNO)
enddo

do CNO=1, NOAC
NOHTBC=0
write(92,*) "chromosome No. ", CNO
SLDBC=0;ALDBC=0;SALDBC=0;AALDBC=0
  do MNO1=SEMNBC(CNO,1), SEMNBC(CNO,2)
    do MNO2=MNO1+1,SEMNBC(CNO,2)
    HG=real(MP(MNO2)-MP(MNO1))/CBPTCM
    NOHTBC=NOHTBC+1
    nMIHMIH=0;nMIHHE=0;nMIHMAH=0;nHEMIH=0;nHEHE=0;nHEMAH=0;nMAHMIH=0;nMAHHE=0;nMAHMAH=0
      do ANO=1,NOA
        select case (G012DS(MNO1,ANO))
          case(0)
          select case (G012DS(MNO2,ANO))
            case(0)
              nMIHMIH=nMIHMIH+1
            case(1)
              nMIHHE=nMIHHE+1
            case(2)
              nMIHMAH=nMIHMAH+1
          end select
          case(1)
          select case (G012DS(MNO2,ANO))
            case(0)
              nHEMIH=nHEMIH+1
            case(1)
              nHEHE=nHEHE+1
            case(2)
              nHEMAH=nHEMAH+1
          end select
          case(2)
          select case (G012DS(MNO2,ANO))
            case(0)
              nMAHMIH=nMAHMIH+1
            case(1)
              nMAHHE=nMAHHE+1
            case(2)
              nMAHMAH=nMAHMAH+1
          end select
        end select
      enddo
      nMIMI=2*nMIHMIH+nMIHHE+nHEMIH+real(nHEHE)/2
      nMIMA=nMIHHE+2*nMIHMAH+real(nHEHE)/2+nHEMAH
      nMAMI=nMAHHE+2*nMAHMIH+real(nHEHE)/2+nHEMIH
      nMAMA=2*nMAHMAH+nMAHHE+nHEMAH+real(nHEHE)/2
      do iter=1,50
        MIMIf=nMIMI/(nMIMI+nMIMA+nMAMI+nMAMA)
        MIMAf=nMIMA/(nMIMI+nMIMA+nMAMI+nMAMA)
        MAMIf=nMAMI/(nMIMI+nMIMA+nMAMI+nMAMA)
        MAMAf=nMAMA/(nMIMI+nMIMA+nMAMI+nMAMA)
        prMIMIMAMA=(MIMIf*MAMAf)/((MIMIf*MAMAf)+(MIMAf*MAMIf))
        prMIMAMAMI=(MIMAf*MAMIf)/((MIMIf*MAMAf)+(MIMAf*MAMIf))
        nMIMI=2*nMIHMIH+nMIHHE+nHEMIH+(nHEHE*prMIMIMAMA)
        nMIMA=nMIHHE+2*nMIHMAH+(nHEHE*prMIMAMAMI)+nHEMAH
        nMAMI=nMAHHE+2*nMAHMIH+(nHEHE*prMIMAMAMI)+nHEMIH
        nMAMA=2*nMAHMAH+nMAHHE+nHEMAH+(nHEHE*prMIMIMAMA)
        newMIMIf=nMIMI/(nMIMI+nMIMA+nMAMI+nMAMA)
        newMIMAf=nMIMA/(nMIMI+nMIMA+nMAMI+nMAMA)
        newMAMIf=nMAMI/(nMIMI+nMIMA+nMAMI+nMAMA)
        newMAMAf=nMAMA/(nMIMI+nMIMA+nMAMI+nMAMA)
        DIFF=2*abs(newMIMIf-MIMIf)
        if( DIFF < CONV ) exit
      enddo
    D=(newMIMIf*newMAMAf)-(newMIMAf*newMAMIf)
    rsquare=(D**2)/(GF(MNO1,1)*GF(MNO1,2)*GF(MNO2,1)*GF(MNO2,2))
    SLDBC=SLDBC+rsquare
    write(CNO+20,fmt='(2i15,f15.10,f15.12)') MNO1,MNO2,HG,rsquare
    if(HG.le.DLFLDD) write(95,fmt='(2f15.10)') HG, rsquare 
    if(MNO2-MNO1.eq.1) then 
      SALDBC=SALDBC+rsquare
      write(92,fmt='(2f15.10)') HG, rsquare 
    endif
      do NEJ=1, NONEG
        if(HG.gt.NEG(NEJ,2).and.HG.gt.NEG(NEJ,2)) then
          CNEG(NEJ)=CNEG(NEJ)+1
          SNEG(NEJ)=SNEG(NEJ) + rsquare
        endif
      enddo
    enddo
  enddo
  print *,"haplotype frequency calculation of chromosome", CNO, "is finished" 
  print *,"haplotype number on chromosome", CNO, "is", NOHTBC
  close(CNO+20)
  ALDBC=SLDBC/real(NOHTBC);AALDBC=SALDBC/real(NOSMBC(CNO)+1)
  STLD=STLD+SLDBC;STALD=STALD+SALDBC
  NOTHT=NOTHT+NOHTBC
  write(94,fmt='(a12,i5,i10,2f15.10)') 'chromosome',CNO,NOHTBC,ALDBC,AALDBC
enddo
ATLD=STLD/real(NOTHT);ATALD=STALD/real(NOM-NOAC)
write(94,fmt='(a27, 2f15.10)') 'total', ATLD, ATALD
write(93,*) 'Estimated effective population size'
do NEJ=1, NONEG
  ANEG(NEJ)=SNEG(NEJ)/CNEG(NEJ)
  EPS=(1/((NEG(NEJ,1)+NEG(NEJ,2))/200))*((1/ANEG(NEJ))-1.0)
  write(93,fmt='(i4,a3,i4,10x,f20.10)') NEGA(NEJ,1)," ~ ",NEGA(NEJ,2),EPS
enddo
end program EOLD

