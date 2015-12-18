program HFC !Haplotype Frequency Calculation

integer,parameter::NID=525
integer,parameter::nchromosome=19
character(len=10)::Variables(8)
character,allocatable::Mi(:),Ma(:)
character(len=60)::file_name='C:\My Research\Swine\fromKSU\Output\HF\HFCR00.dat'
integer,allocatable::GN(:,:),noMK(:),CHno(:),MKpo(:)
integer::nMIHMIH,nMIHHE,nMIHMAH,nHEMIH,nHEHE,nHEMAH,nMAHMIH,nMAHHE,nMAHMAH
real::nMIMI,nMIMA,nMAMI,nMAMA
real,allocatable::MIAf(:),MAAf(:)
integer::io,nnMK,c,i,k,l,chinfo(nchromosome,2)=0,Dist=0
open(99,file='C:\My Research\Swine\fromKSU\Output\genotype012UDE.dat',status='old')
nnMK=0
do 
 read(99,*,iostat=io) 
  if(io/=0) exit
 nnMK=nnMK+1
enddo
rewind(unit=99)
nnMK=nnMK-1
read(99,*) Variables
allocate(GN(nnMK,NID),noMK(nnMK),CHno(nnMK),MKpo(nnMK),Mi(nnMK),Ma(nnMK),MAAf(nnMK),MIAf(nnMK))
chinfo(1,1)=1;chinfo(nchromosome,2)=nnMK
do i=1, nnMK
  read(99,fmt='(3i10,2F10.7,2A4,1x,600i1)')noMK(i),CHno(i),MKpo(i),MIAf(i),MAAf(i),Mi(i),Ma(i),GN(i,:)
  if(CHno(i)>1.AND.CHno(i-1).NE.CHno(i)) then
    chinfo(CHno(i-1),2)=noMK(i-1)
    chinfo(CHno(i),1)=noMK(i)
  endif
enddo
close(99)
nMIHMIH=0;nMIHHE=0;nMIHMAH=0;nHEMIH=0;nHEHE=0;nHEMAH=0;nMAHMIH=0;nMAHHE=0;nMAHMAH=0;nMIMI=0;nMIMA=0;nMAMI=0;nMAMA=0
do c=1, nchromosome
  write(file_name(44:45),'(I2.2)') c
  open(c,file=file_name,status='unknown')
  write(c,fmt='(24A10)') "Dist","i","MIAf_i","MAAf_i","Mi_i","Ma_i","j","MIAf_j","MAAf_j","Mi_j","Ma_j","nMIHMIH",&
  &"nMIHHE","nMIHMAH","nHEMIH","nHEHE","nHEMAH","nMAHMIH","nMAHHE","nMAHMAH","nMIMI","nMIMA","nMAMI","nMAMA"
  do i=chinfo(c,1), chinfo(c,2)
    do j=i+1,chinfo(c,2)
    dist=abs(MKpo(j)-MKpo(i))
      do k=1,NID
        select case (GN(i,k))
          case(0)
          select case (GN(j,k))
            case(0)
              nMIHMIH=nMIHMIH+1
            case(1)
              nMIHHE=nMIHHE+1
            case(2)
              nMIHMAH=nMIHMAH+1
          end select
          case(1)
          select case (GN(j,k))
            case(0)
              nHEMIH=nHEMIH+1
            case(1)
              nHEHE=nHEHE+1
            case(2)
              nHEMAH=nHEMAH+1
          end select
          case(2)
          select case (GN(j,k))
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
  write(c,fmt='(2i10,2F10.7,2a10,i10,2F10.7,2a10,9i10,4f10.1)') Dist,i,MIAf(i),MAAf(i),Mi(i),Ma(i),j,MIAf(j),MAAf(j),Mi(j),Ma(j),&
      &nMIHMIH,nMIHHE,nMIHMAH,nHEMIH,nHEHE,nHEMAH,nMAHMIH,nMAHHE,nMAHMAH,nMIMI,nMIMA,nMAMI,nMAMA
      nMIHMIH=0;nMIHHE=0;nMIHMAH=0;nHEMIH=0;nHEHE=0;nHEMAH=0;nMAHMIH=0;nMAHHE=0;nMAHMAH=0;nMIMI=0;nMIMA=0;nMAMI=0;nMAMA=0;dist=0
    enddo
  enddo
!  print *,"finishing chromosome no. ", c, chinfo(c,1), chinfo(c,2)
  close(c)
enddo
end program HFC
