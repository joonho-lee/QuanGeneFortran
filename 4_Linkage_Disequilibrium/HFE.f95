program HFE !Haplotype Frequency Estimation using EM Algorithm

integer,parameter::NID=525
integer,parameter::nchromosome=19
real,parameter::conv=0.00001
character(len=10)::Variables(20)
character::Mii,Mai,Mij,Maj
character(len=59)::file_namer='C:\My Research\Swine\fromKSU\Output\HF\HFCR00.dat'
character(len=64)::file_namew='C:\My Research\Swine\fromKSU\Output\HF\EHF\EHFCR00.dat'
integer::Ai,Aj,Dist
real::nMIHMIH,nMIHHE,nMIHMAH,nHEMIH,nHEHE,nHEMAH,nMAHMIH,nMAHHE,nMAHMAH,nMIMI,nMIMA,nMAMI,nMAMA
real::MIMIf,MIMAf,MAMIf,MAMAf,newMIMIf,newMIMAf,newMAMIf,newMAMAf,prMIMIMAMA,prMIMAMAMI,diff=0.0,MIAfi,MAAfi,MIAfj,MAAfj
integer::io,nnnMK,c,i,iter,k
do c=1, nchromosome
  k=c+50
  write(file_namer(44:45),'(I2.2)') c
  write(file_namew(49:50),'(I2.2)') c
  open(c,file=file_namer,status='old')
  open(k,file=file_namew,status='unknown')
  write(k,fmt='(21A10)') "Dist","Ai","MIAfi","MAAfi","Mii","Mai","Aj","MIAfj","MAAfj","Mij","Maj","nMIMI","nMIMA","nMAMI","nMAMA",&
  &"newMIMIf","newMIMAf","newMAMIf","newMAMAf","diff","iter"
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
    read(c,fmt='(2i10,2f10.7,2a10,i10,2f10.7,2a10,9F10.0,4F10.1)') Dist,Ai,MIAfi,MAAfi,Mii,Mai,Aj,MIAfj,MAAfj,Mij,Maj,&
    &nMIHMIH,nMIHHE,nMIHMAH,nHEMIH,nHEHE,nHEMAH,nMAHMIH,nMAHHE,nMAHMAH,nMIMI,nMIMA,nMAMI,nMAMA
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
      diff=2*abs(newMIMIf-MIMIf)
      if( diff < conv ) goto 10
    enddo
10  write(k,fmt='(2i10,2f10.7,2a10,i10,2f10.7,2a10,9f10.5,i10)') Dist,Ai,MIAfi,MAAfi,Mii,Mai,Aj,MIAfj,MAAfj,Mij,Maj,&
    &nMIMI,nMIMA,nMAMI,nMAMA,newMIMIf,newMIMAf,newMAMIf,newMAMAf,diff,iter
  enddo
  close(c)
  close(k)
!  print*, "file NO.",c,"finished"
enddo
end program HFE
