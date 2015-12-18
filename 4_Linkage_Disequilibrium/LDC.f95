program LDC !Linkage Disequilibrium Calculation among markers

integer,parameter::nchromosome=19
character(len=10)::Variables(21)
character::Mii,Mai,Mij,Maj
character(len=64)::file_namer='C:\My Data\Research\Genomic data\Swine\DATA\HF\EHF\EHFCR00.dat'
character(len=59)::file_nameW='C:\My Data\Research\Genomic data\Swine\DATA\LD\LDCR00.dat'
integer::Ai,Aj,Dist
real::nMIMI,nMIMA,nMAMI,nMAMA,newMIMIf,newMIMAf,newMAMIf,newMAMAf,diff,MIAfi,MAAfi,MIAfj,MAAfj,D,rsquare
integer::io,nnnMK,c,i,iter,k
do c=1, nchromosome
  k=c+50
  write(file_namer(57:58),'(I2.2)') c
  write(file_namew(52:53),'(I2.2)') c
  open(c,file=file_namer,status='old')
  open(k,file=file_namew,status='unknown')
  write(k,fmt='(12A10)') "Dist","Ai","MIAfi","MAAfi","Mii","Mai","Aj","MIAfj","MAAfj","Mij","Maj","r2"
  nnnMK=0
  do 
    read(c,*,iostat=io) 
      if(io/=0) exit
    nnnMK=nnnMK+1
  enddo
  rewind(unit=c)
  nnnMK=nnnMK-1
  read(c,*) Variables
  do i=1,nnnMK
    read(c,fmt='(2i10,2f10.7,2a10,i10,2f10.7,2a10,9f10.5,i10)') Dist,Ai,MIAfi,MAAfi,Mii,Mai,Aj,MIAfj,MAAfj,Mij,Maj,&
    &nMIMI,nMIMA,nMAMI,nMAMA,newMIMIf,newMIMAf,newMAMIf,newMAMAf,diff,iter
    D=(newMIMIf*newMAMAf)-(newMIMAf*newMAMIf)
    rsquare=(D**2)/(MIAfi*MAAfi*MIAfj*MAAfj)
    write(k,fmt='(2i10,2f10.7,2a10,i10,2f10.7,2a10,f10.7)') Dist,Ai,MIAfi,MAAfi,Mii,Mai,Aj,MIAfj,MAAfj,Mij,Maj,rsquare
  enddo
  close(c)
  close(k)
!  print*, "file NO.",c,"finished"
enddo
end program LDC
