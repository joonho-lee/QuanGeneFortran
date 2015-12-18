! Name : Making Genomic Relationship Matrix 
! Developed by : Joonho Lee - PhD working in Hankyong National Univ.(Korea) zoonolee@gmail.com, +82-10-3408-2895
! Version : Linux 4.0 (2012-06-15)
! Contents : 		1. Ordering Animal information for Genomic Relationship Matrix
!		 	2. Making Genomic Relationship Matrix
!		 	3. Inversion of Genomic Relationship Matrix
!			4. Statistics of Relationship

program MGRM 
integer::i,j,l,mk,an,ErrorF !iostate & do_loop
integer::NOA=0,ANO=0,NOM=0,MNO=0 !number of animals and markers & animal and marekr number
integer,allocatable::MT(:,:) !data storage

real(10)::GF1=0,GF2=0,SEH=0,SG05=0,SGMF=0,AOMAF=0,SGF1,SSGF1,AGF1,VGF1,ALPHA,BETA,P0,Q0,KGOFS,TRG05  !gnen frequencies & expected heterozygosity
real(10)::SOIBCG05=0,SSOIBCG05=0,AOIBCG05=0,SDOIBCG05=0,MINOIBCG05=0,MAXOIBCG05=0 !statistics for population genetics
real(10)::SOIBCGMF=0,SSOIBCGMF=0,AOIBCGMF=0,SDOIBCGMF=0,MINOIBCGMF=0,MAXOIBCGMF=0 !statistics for population genetics
real(10)::SOIBCGOF=0,SSOIBCGOF=0,AOIBCGOF=0,SDOIBCGOF=0,MINOIBCGOF=0,MAXOIBCGOF=0 !statistics for population genetics
real(10)::SOIBCGOFS=0,SSOIBCGOFS=0,AOIBCGOFS=0,SDOIBCGOFS=0,MINOIBCGOFS=0,MAXOIBCGOFS=0 !statistics for population genetics
real(10)::SOIBCGN=0,SSOIBCGN=0,AOIBCGN=0,SDOIBCGN=0,MINOIBCGN=0,MAXOIBCGN=0 !statistics for population genetics
real,allocatable::G05(:,:),IG05(:,:),IBCG05(:),W05(:,:),WMF(:,:) !output data
real,allocatable::GMF(:,:),IGMF(:,:),IBCGMF(:)
real,allocatable::GOF(:,:),IGOF(:,:),IBCGOF(:)
real,allocatable::GOFS(:,:),IGOFS(:,:),IBCGOFS(:)
real,allocatable::GN(:,:),IGN(:,:),IBCGN(:)
real,allocatable::k(:),DII(:)
character(len=8)::OFT,NOAF
character(len=16)::WF1,WF2,WF3
character(len=300)::Dummy
character(len=220)::ADD,ODD,GDD
character(len=30)::GDFN1,GDFN2,ADFN
character(len=250)::ADFNIP,GDFNIP,ODFNIP,PF
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
read (99,*) ADD
print *, "directory of Selected animals file ", ADD
read (99,*) ADFN
print *, "name of Selected animals file(example : Selected_Animals_OUT_TAG.out) ", ADFN
read (99,*) ODD
print *, "directory of output file ", ODD
close(99)

allocate(ID(NOA))
allocate(MT(NOM,NOA))
allocate(G05(NOA,NOA))
allocate(IG05(NOA,NOA))
allocate(IBCG05(NOA))
allocate(W05(3,NOM))
allocate(k(NOA))

allocate(GMF(NOA,NOA))
allocate(IGMF(NOA,NOA))
allocate(IBCGMF(NOA))
allocate(WMF(3,NOM))

allocate(GOF(NOA,NOA))
allocate(IGOF(NOA,NOA))
allocate(IBCGOF(NOA))
allocate(DII(NOM))

allocate(GOFS(NOA,NOA))
allocate(IGOFS(NOA,NOA))
allocate(IBCGOFS(NOA))

allocate(GN(NOA,NOA))
allocate(IGN(NOA,NOA))
allocate(IBCGN(NOA))

G05=0;IG05=0
GMF=0;IGMF=0
GOF=0;IGOF=0;
GOFS=0;IGOFS=0;
GN=0;IGN=0

WF1='            ';WF2='            ';WF3='            ';
write(WF1,fmt='(a16)') '('//adjustr(NOAF)//'f15.11)'
write(WF2,fmt='(a16)') '('//adjustr(NOAF)//'f22.15)'
write(WF3,fmt='(a16)') '('//adjustr(NOAF)//'i1    )'
print *, WF1,WF2,WF3
ODD=adjustr(ODD)
GDD=adjustr(GDD)
ADD=adjustr(ADD)
GDFNIP=GDD//GDFN1
GDFNIP=adjustl(GDFNIP)
open(98,file=GDFNIP,status='old')
GDFNIP=GDD//GDFN2
GDFNIP=adjustl(GDFNIP)
open(97,file=GDFNIP,status='old')
ADFNIP=ADD//ADFN
ADFNIP=adjustl(ADFNIP)
open(96,file=ADFNIP,status='old')
ODFNIP=ODD//"G05"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(95,file=ODFNIP,status='unknown')
ODFNIP=ODD//"IG05"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(94,file=ODFNIP,status='unknown')
ODFNIP=ODD//"IBCG05"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(93,file=ODFNIP,status='unknown')
write(93,fmt='(2a20)') "Animal ID","Inbreeding coeff"
ODFNIP=ODD//"G05RC_selected"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(92,file=ODFNIP,status='unknown')

ODFNIP=ODD//"GMF"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(91,file=ODFNIP,status='unknown')
ODFNIP=ODD//"IGMF"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(90,file=ODFNIP,status='unknown')
ODFNIP=ODD//"IBCGMF"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(89,file=ODFNIP,status='unknown')
write(89,fmt='(2a20)') "Animal ID","Inbreeding coeff"
ODFNIP=ODD//"GMFRC_selected"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(88,file=ODFNIP,status='unknown')

ODFNIP=ODD//"GOF"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(87,file=ODFNIP,status='unknown')
ODFNIP=ODD//"IGOF"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(86,file=ODFNIP,status='unknown')
ODFNIP=ODD//"IBCGOF"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(85,file=ODFNIP,status='unknown')
write(85,fmt='(2a20)') "Animal ID","Inbreeding coeff"
ODFNIP=ODD//"GOFRC_selected"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(84,file=ODFNIP,status='unknown')

ODFNIP=ODD//"GOFS"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(83,file=ODFNIP,status='unknown')
ODFNIP=ODD//"IGOFS"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(82,file=ODFNIP,status='unknown')
ODFNIP=ODD//"IBCGOFS"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(81,file=ODFNIP,status='unknown')
write(81,fmt='(2a20)') "Animal ID","Inbreeding coeff"
ODFNIP=ODD//"GOFSRC_selected"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(80,file=ODFNIP,status='unknown')

ODFNIP=ODD//"GN"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(79,file=ODFNIP,status='unknown')
ODFNIP=ODD//"IGN"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(78,file=ODFNIP,status='unknown')
ODFNIP=ODD//"IBCGN"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(77,file=ODFNIP,status='unknown')
write(77,fmt='(2a20)') "Animal ID","Inbreeding coeff"
ODFNIP=ODD//"GNRC_selected"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(76,file=ODFNIP,status='unknown')

read(97,*) GF1,GF2

SGF1=0;SSGF1=0
do mk=1, NOM
  read(98,fmt=WF3) MT(mk,:)
  read(97,*) GF1,GF2
  AOMAF=AOMAF+GF1
  SEH=SEH+(GF1*GF2*2)
  DII(mk)=1.0/(NOM*GF1*GF2*2)
  W05(1,mk)=GF2*(-2)
  W05(2,mk)=1.0-GF2*2
  W05(3,mk)=2.0-GF2*2
  SGF1=SGF1+GF1
!  SGF2=SGF2+GF2
  SSGF1=SSGF1+GF1**2
!  SSGF2=SSGF2+GF2**2
enddo
AGF1=SGF1/NOM
!AGF2=SGF2/NOM
VGF1=(SSGF1-(SGF1**2/NOM))/NOM
!VGF2=(SSGF2-(SGF2**2/NOM))/NOM
ALPHA=AGF1*(AGF1*(1-AGF1)/VGF1-1)
BETA=(1-AGF1)*(AGF1*(1-AGF1)/VGF1-1)
AOMAF=AOMAF/NOM
Q0=ALPHA/(ALPHA+BETA)
P0=BETA/(ALPHA+BETA)
KGOFS=NOM*((P0-Q0)**2+(SEH/NOM)*((ALPHA+BETA+2)/(ALPHA+BETA)))

close(98)
close(97)
k=SEH

do mk=1, NOM
  WMF(1,mk)=AOMAF*2-GF1*2-1.0
  WMF(2,mk)=AOMAF*2-GF1*2
  WMF(3,mk)=AOMAF*2-GF1*2+1.0  
enddo

print *, "reading genotype is finished"

TRG05=0
read(96,*)
do i=1,NOA
  read(96,*) ID(i)
  do j=1,NOA
    SG05=0
    SGMF=0
    SGOF=0
    do l=1,NOM
      SG05=SG05+W05(MT(l,i)+1,l)*W05(MT(l,j)+1,l)
      SGMF=SGMF+WMF(MT(l,i)+1,l)*WMF(MT(l,j)+1,l)
      SGOF=SGOF+W05(MT(l,i)+1,l)*W05(MT(l,j)+1,l)*DII(l)
    enddo
    IF(i.EQ.j) TRG05=TRG05+SG05
    G05(i,j)=SG05/sqrt(k(i)*k(j))
    G05(j,i)=G05(i,j)
    GMF(i,j)=SGMF/sqrt(k(i)*k(j))
    GMF(j,i)=GMF(i,j)
    GOF(i,j)=SGOF
    GOF(j,i)=GOF(i,j)
    GOFS(i,j)=SG05/KGOFS
    GOFS(j,i)=GOFS(i,j)
    GN(i,j)=SG05
    GN(j,i)=GN(i,j)     
  enddo
enddo

GN=GN/(TRG05/NOA)
close(96)

print *, "making genomic relationship matrix is finished "
CALL FindInv(G05, IG05, NOA, ErrorF)
if(ErrorF.lt.0) print *, "error status during inversion"
print *, "making inverse genomic relationship matrix is finished -G05"
CALL FindInv(GMF, IGMF, NOA, ErrorF)
if(ErrorF.lt.0) print *, "error status during inversion"
print *, "making inverse genomic relationship matrix is finished -GMF"
CALL FindInv(GOF, IGOF, NOA, ErrorF)
if(ErrorF.lt.0) print *, "error status during inversion"
print *, "making inverse genomic relationship matrix is finished -GOF"
CALL FindInv(GOFS, IGOFS, NOA, ErrorF)
if(ErrorF.lt.0) print *, "error status during inversion"
print *, "making inverse genomic relationship matrix is finished -GOFS"
CALL FindInv(GN, IGN, NOA, ErrorF)
if(ErrorF.lt.0) print *, "error status during inversion"
print *, "making inverse genomic relationship matrix is finished -GN"

do j=1, NOA
  IBCG05(j)=G05(j,j)-1.0
  SSOIBCG05=SSOIBCG05+IBCG05(j)**2
  write(93,fmt='(a20,f20.10)') ID(j),IBCG05(j)
  write(95,fmt=WF1) G05(j,:)
  write(94,fmt=WF2) IG05(j,:)
  do i=1,NOA
    if(i.ge.j) write(92,fmt='(2i10,f15.10,2a22)') j,i,G05(j,i),ID(j),ID(i)
  enddo
enddo
AOIBCG05=sum(IBCG05)/NOA;SDIBCG05=(SSOIBCG05-(sum(IBCG05)**2/NOA))/NOA
MINOIBCG05=minval(IBCG05);MAXOIBCG05=maxval(IBCG05)
write(93,fmt='(4a20)') "AOIBC","SDIBC","MINOIBC","MAXOIBC"
write(93,fmt='(4f20.10)') AOIBCG05,SDIBCG05,MINOIBCG05,MAXOIBCG05

do j=1, NOA
  IBCGMF(j)=GMF(j,j)-1.0
  SSOIBCGMF=SSOIBCGMF+IBCGMF(j)**2
  write(89,fmt='(a20,f20.10)') ID(j),IBCGMF(j)
  write(91,fmt=WF1) GMF(j,:)
  write(90,fmt=WF2) IGMF(j,:)
  do i=1,NOA
    if(i.ge.j) write(88,fmt='(2i10,f15.10,2a22)') j,i,GMF(j,i),ID(j),ID(i)
  enddo
enddo
AOIBCGMF=sum(IBCGMF)/NOA;SDIBCGMF=(SSOIBCGMF-(sum(IBCGMF)**2/NOA))/NOA
MINOIBCGMF=minval(IBCGMF);MAXOIBCGMF=maxval(IBCGMF)
write(89,fmt='(4a20)') "AOIBC","SDIBC","MINOIBC","MAXOIBC"
write(89,fmt='(4f20.10)') AOIBCGMF,SDIBCGMF,MINOIBCGMF,MAXOIBCGMF

do j=1, NOA
  IBCGOF(j)=GOF(j,j)-1.0
  SSOIBCGOF=SSOIBCGOF+IBCGOF(j)**2
  write(85,fmt='(a20,f20.10)') ID(j),IBCGOF(j)
  write(87,fmt=WF1) GOF(j,:)
  write(86,fmt=WF2) IGOF(j,:)
  do i=1,NOA
    if(i.ge.j) write(84,fmt='(2i10,f15.10,2a22)') j,i,GOF(j,i),ID(j),ID(i)
  enddo
enddo
AOIBCGOF=sum(IBCGOF)/NOA;SDIBCGOF=(SSOIBCGOF-(sum(IBCGOF)**2/NOA))/NOA
MINOIBCGOF=minval(IBCGOF);MAXOIBC=maxval(IBCGOF)
write(85,fmt='(4a20)') "AOIBC","SDIBC","MINOIBC","MAXOIBC"
write(85,fmt='(4f20.10)') AOIBCGOF,SDIBCGOF,MINOIBCGOF,MAXOIBCGOF

do j=1, NOA
  IBCGOFS(j)=GOFS(j,j)-1.0
  SSOIBCGOFS=SSOIBCGOFS+IBCGOFS(j)**2
  write(81,fmt='(a20,f20.10)') ID(j),IBCGOFS(j)
  write(83,fmt=WF1) GOFS(j,:)
  write(82,fmt=WF2) IGOFS(j,:)
  do i=1,NOA
    if(i.ge.j) write(80,fmt='(2i10,f15.10,2a22)') j,i,GOFS(j,i),ID(j),ID(i)
  enddo
enddo
AOIBCGOFS=sum(IBCGOFS)/NOA;SDIBCGOFS=(SSOIBCGOFS-(sum(IBCGOFS)**2/NOA))/NOA
MINOIBCGOFS=minval(IBCGOFS);MAXOIBCGOFS=maxval(IBCGOFS)
write(81,fmt='(4a20)') "AOIBC","SDIBC","MINOIBC","MAXOIBC"
write(81,fmt='(4f20.10)') AOIBCGOFS,SDIBCGOFS,MINOIBCGOFS,MAXOIBCGOFS

do j=1, NOA
  IBCGN(j)=GN(j,j)-1.0
  SSOIBCGN=SSOIBCGN+IBCGN(j)**2
  write(77,fmt='(a20,f20.10)') ID(j),IBCGN(j)
  write(79,fmt=WF1) GN(j,:)
  write(78,fmt=WF2) IGN(j,:)
  do i=1,NOA
    if(i.ge.j) write(76,fmt='(2i10,f15.10,2a22)') j,i,GN(j,i),ID(j),ID(i)
  enddo
enddo
AOIBCGN=sum(IBCGN)/NOA;SDIBCGN=(SSOIBCGN-(sum(IBCGN)**2/NOA))/NOA
MINOIBCGN=minval(IBCGN);MAXOIBCGN=maxval(IBCGN)
write(77,fmt='(4a20)') "AOIBC","SDIBC","MINOIBC","MAXOIBC"
write(77,fmt='(4f20.10)') AOIBCGN,SDIBCGN,MINOIBCGN,MAXOIBCGN


print *, "estimating inbreeding coefficient is finished "

end program MGRM

!Subroutine comes here 

!Subroutine to find the inverse of a square matrix
!Author : Louisda16th a.k.a Ashwith J. Rego
!Reference : Algorithm has been well explained in:
!http://math.uww.edu/~mcfarlat/inverse.htm           
!http://www.tutor.ms.unimelb.edu.au/matrix/matrix_inverse.html
SUBROUTINE FINDInv(matrix, inverse, n, errorflag)
	IMPLICIT NONE
	!Declarations
	INTEGER, INTENT(IN) :: n
	INTEGER, INTENT(OUT) :: errorflag  !Return error status. -1 for error, 0 for normal
	REAL, INTENT(IN), DIMENSION(n,n) :: matrix  !Input matrix
	REAL, INTENT(OUT), DIMENSION(n,n) :: inverse !Inverted matrix
	
	LOGICAL :: FLAG = .TRUE.
	INTEGER :: i, j, k, l
	REAL :: m
	REAL, DIMENSION(n,2*n) :: augmatrix !augmented matrix
	
	!Augment input matrix with an identity matrix
	DO i = 1, n
		DO j = 1, 2*n
			IF (j <= n ) THEN
				augmatrix(i,j) = matrix(i,j)
			ELSE IF ((i+n) == j) THEN
				augmatrix(i,j) = 1
			Else
				augmatrix(i,j) = 0
			ENDIF
		END DO
	END DO
	
	!Reduce augmented matrix to upper traingular form
	DO k =1, n-1
		IF (augmatrix(k,k) == 0) THEN
			FLAG = .FALSE.
			DO i = k+1, n
				IF (augmatrix(i,k) /= 0) THEN
					DO j = 1,2*n
						augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
					END DO
					FLAG = .TRUE.
					EXIT
				ENDIF
				IF (FLAG .EQV. .FALSE.) THEN
					PRINT*, "Matrix is non - invertible"
					inverse = 0
					errorflag = -1
					return
				ENDIF
			END DO
		ENDIF
		DO j = k+1, n			
			m = augmatrix(j,k)/augmatrix(k,k)
			DO i = k, 2*n
				augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
			END DO
		END DO
	END DO
	
	!Test for invertibility
	DO i = 1, n
		IF (augmatrix(i,i) == 0) THEN
			PRINT*, "Matrix is non - invertible"
			inverse = 0
			errorflag = -1
			return
		ENDIF
	END DO
	
	!Make diagonal elements as 1
	DO i = 1 , n
		m = augmatrix(i,i)
		DO j = i , (2 * n)				
			   augmatrix(i,j) = (augmatrix(i,j) / m)
		END DO
	END DO
	
	!Reduced right side half of augmented matrix to identity matrix
	DO k = n-1, 1, -1
		DO i =1, k
		m = augmatrix(i,k+1)
			DO j = k, (2*n)
				augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m
			END DO
		END DO
	END DO				
	
	!store answer
	DO i =1, n
		DO j = 1, n
			inverse(i,j) = augmatrix(i,j+n)
		END DO
	END DO
	errorflag = 0
END SUBROUTINE FINDinv



