! Name : Making Genomic Relationship Matrix 
! Developed by : Joonho Lee - PhD working in Hankyong National Univ.(Korea) zoonolee@gmail.com, +82-10-3408-2895
! Version : Linux 3.0 (2012-05-25)
! Contents : 		1. Ordering Animal information for Genomic Relationship Matrix
!		 	2. Making Genomic Relationship Matrix
!		 	3. Inversion of Genomic Relationship Matrix
!			4. Statistics of Relationship

program MGRM 
integer::i,j,l,mk,an,ErrorF !iostate & do_loop
integer::NOA=0,ANO=0,NOM=0,MNO=0 !number of animals and markers & animal and marekr number
integer,allocatable::MT(:,:) !data storage

real(10)::GF1=0, GF2=0, SEH=0, S=0  !gnen frequencies & expected heterozygosity
real(10)::SOIBC=0,SSOIBC=0,AOIBC=0,SDOIBC=0,MINOIBC=0,MAXOIBC=0 !statistics for population genetics
real,allocatable::GRM(:,:),IGRM(:,:),W(:,:),IBC(:),k(:) !output data

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

allocate(GRM(NOA,NOA))
allocate(IGRM(NOA,NOA))
allocate(ID(NOA))
allocate(IBC(NOA))
allocate(MT(NOM,NOA))
allocate(W(3,NOM))
allocate(k(NOA))
GRM=0;IGRM=0
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
ODFNIP=ODD//"GRC_selected"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(92,file=ODFNIP,status='unknown')

read(97,*) GF1,GF2

do mk=1, NOM
  read(98,fmt=WF3) MT(mk,:)
  read(97,*) GF1,GF2
  SEH=SEH+(GF1*GF2*2)
  W(1,mk)=GF2*(-2)
  W(2,mk)=1.0-GF2*2
  W(3,mk)=2.0-GF2*2  
enddo
close(98)
close(97)
k=SEH

print *, "reading genotype is finished"

read(96,*)
do i=1,NOA
  read(96,*) ID(i)
  do j=1,NOA
    S=0
    do l=1,NOM
      S=S+W(MT(l,i)+1,l)*W(MT(l,j)+1,l)
    enddo
    GRM(i,j)=S/sqrt(k(i)*k(j))
    GRM(j,i)=GRM(i,j)
  enddo
enddo
close(96)

print *, "making genomic relationship matrix is finished "
CALL FindInv(GRM, IGRM, NOA, ErrorF)
if(ErrorF.lt.0) print *, "error status during inversion"
print *, "making inverse genomic relationship matrix is finished "

do j=1, NOA
  IBC(j)=GRM(j,j)-1.0
  SSOIBC=SSOIBC+IBC(j)**2
  write(93,fmt='(a20,f20.10)') ID(j),IBC(j)
  write(95,fmt=WF1) GRM(j,:)
  write(94,fmt=WF2) IGRM(j,:)
  do i=1,NOA
    if(i.ge.j) write(92,fmt='(2i10,f15.10,2a22)') j,i,GRM(j,i),ID(j),ID(i)
  enddo
enddo
AOIBC=sum(IBC)/NOA;SDIBC=(SSOIBC-(sum(IBC)**2/NOA))/NOA
MINOIBC=minval(IBC);MAXOIBC=maxval(IBC)
write(93,fmt='(4a20)') "AOIBC","SDIBC","MINOIBC","MAXOIBC"
write(93,fmt='(4f20.10)') AOIBC,SDIBC,MINOIBC,MAXOIBC

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



