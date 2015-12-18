! Name : EM_REML using GRM
! Developed by : Joonho Lee - PhD working in Hankyong National Univ.(Korea) zoonolee@gmail.com, +82-10-3408-2895
! Version : Linux 1.0 (2012-04-16)
! Contents : 		1. Estimation of marker effect
!		 	2. Estimation of breeding value using marker effect

program RRB 
integer::ANO=0,NOM=0,MNO=0,NOR=0,RNO=0,NOTR=0,TRNO=0,FENO=0,NOFE=0,Dummy,PCHNO=1,MAXP=0,SG !numbers and counting variables
integer,allocatable::TG012DS(:,:),G012DS(:,:),NOA(:),FEDS(:,:),IOM(:,:),NOMWSG(:,:),OMN(:),NOSM(:) !data storage
real::JtJ,JtX,XtJ,XtX,JtY,XtY,YtY,det,mu_hat,g_hat,F_val,tmp,P_val,log10P_val,CLOP
real,allocatable::PRDS(:)
character(len=3),allocatable::NAOTR(:)
character(len=8),allocatable::NOAF(:)
character(len=50),allocatable::MN(:)
character(len=8)::OFT,NOMF,NOFEF,NOTRF
character(len=4)::MAXPF
character(len=16)::WF1,WF2,WF3,WF4,WF5
character(len=220)::MDD,ODD,GDD1,GDD2,PDD
character(len=30)::GDFN1,GDFN2,PDFN,MDFN
character(len=250)::MDFNIP,GDFNIP1,GDFNIP2,ODFNIP,PDFNIP,PF

!reading parameters
print *, "***********************************************************************************************************************" 
print *, "MGP_Lv1 : Extent of Linkage Disequilibrium Linux version 1.0"
print *, "Developed by Joonho Lee : zoonolee@gmail.com, +82-10-3408-2895"
print *, "Genotype data form is the results of Post_fastPHASE and of converted genotype 0,1,2"
print *, "Maximum number of trait is 20"
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
read (99,*) NOTR
allocate(NAOTR(NOTR))
allocate(NOA(NOTR))
allocate(NOAF(NOTR))
allocate(NOSM(NOTR))

read (99,*) NOM
read (99,*) MAXP
read (99,*) CLOP
read (99,*) NOTR
allocate(MN(NOM))
allocate(OMN(NOM))
allocate(IOM(NOM,2))
rewind(unit=99)
read (99,*) OFT
read (99,*) NOMF
read (99,*) MAXPF
read (99,*) CLOP
read (99,*) NOTRF
print *, "number of markers", NOM
print *, "number of traits", NOTR
allocate(NAOTR(NOTR))
allocate(NOA(NOTR))
allocate(NOAF(NOTR))
allocate(NOSM(NOTR))
allocate(NOMWSG(MAXP+1,NOTR))
read (99,*) NAOTR(:)
print *, "name of traits and number of records"
do TRNO=1,NOTR
  read (99,*) NOA(TRNO)
  backspace(unit=99)
  read (99,*) NOAF(TRNO)
  print *, NAOTR(TRNO),NOA(TRNO),NOAF(TRNO) 
enddo
read (99,*) GDD1
print *, "directory of matched genotype and phenotype data file ", GDD1
read (99,*) GDD2
print *, "directory of Selected_Markers file ", GDD2
read (99,*) GDFN2
print *, "file name of Selected_Markers ", GDFN2
read (99,*) ODD
print *, "directory of output file ", ODD
close(99)

ODD=adjustr(ODD)
GDD1=adjustr(GDD1)
GDD2=adjustr(GDD2)
GDFN2=adjustl(GDFN2)
GDFNIP2=GDD2//GDFN2
GDFNIP2=adjustl(GDFNIP2)
open(98,file=GDFNIP2,status='old')
read(98,*)
do MNO=1,NOM
  read(98,*) MN(MNO),OMN(MNO),IOM(MNO,1),IOM(MNO,2)
enddo

WF1='            ';WF2='            ';WF3='            ';WF4='            '
write(WF1,fmt='(a16)') '('//adjustr(NOMF)//'i1    )'
write(WF2,fmt='(a16)') '(a9,'//adjustr(NOTRF)//'a15)'
write(WF3,fmt='(a16)') '(i9,'//adjustr(NOTRF)//'i15)'
write(WF4,fmt='(a16)') '(a9,'//adjustr(NOTRF)//'i15)'
print *, "reading format for matched genotype file is : ",WF1
print *, "writing format for number of signigicant markers file is : ",WF2,WF3,WF4

ODFNIP=ODD//"SM_frequency"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(97,file=ODFNIP,status='unknown')
write(97,fmt=WF2) "-log10_P", NAOTR(:)
ODFNIP=ODD//"SM_numbers"//OFT//".out"
ODFNIP=adjustl(ODFNIP)
open(96,file=ODFNIP,status='unknown')

NOMWSG=0
NOSM=0
do TRNO=1, NOTR
  WF5='            '
  write(WF5,fmt='(a16)') '('//adjustr(NOAF(TRNO))//'i1    )'
  GDFNIP1=GDD1//"GT"//NAOTR(TRNO)//OFT//".out"
  GDFNIP1=adjustl(GDFNIP1)
  open(TRNO+10,file=GDFNIP1,status='old')
  GDFNIP1=GDD1//"PR"//NAOTR(TRNO)//OFT//".out"
  GDFNIP1=adjustl(GDFNIP1)
  open(TRNO+30,file=GDFNIP1,status='old')
  ODFNIP=ODD//"ME"//NAOTR(TRNO)//OFT//".out"
  ODFNIP=adjustl(ODFNIP)
  open(TRNO+70,file=ODFNIP,status='unknown')
  write(TRNO+70,fmt='(a50,7a15)')"Marker Name","Chromosome","Position","Mu_hat","Me_hat","F_val","P_val","-log10P_val"
  allocate(TG012DS(NOA(TRNO),NOM))
  allocate(G012DS(NOM,NOA(TRNO)))
  allocate(PRDS(NOA(TRNO)))
  do ANO=1,NOA(TRNO)
    read(TRNO+10,fmt=WF1) TG012DS(ANO,:)
    read(TRNO+30,*) PRDS(ANO)
  enddo
  G012DS=transpose(TG012DS)
  close(TRNO+10)
  close(TRNO+30)
  ODFNIP=ODD//"RG"//NAOTR(TRNO)//OFT//".out"
  ODFNIP=adjustl(ODFNIP)
  open(TRNO+10,file=ODFNIP,status='unknown')
  ODFNIP=ODD//"RM"//NAOTR(TRNO)//OFT//".out"
  ODFNIP=adjustl(ODFNIP)
  open(TRNO+30,file=ODFNIP,status='unknown')
  JtJ=NOA(TRNO)
!  print *, "JtJ",JtJ
  JtY=sum(PRDS)
!  print *, "JtY",JtY
  YtY=sum(PRDS**2)
!  print *, "YtY",YtY
  do MNO=1,NOM
    JtX=sum(TG012DS(:,MNO))
    XtJ=JtX
!    print *, "XtJ and JtX",XtJ, JtX
    XtX=sum(TG012DS(:,MNO)**2)
!    print *, "XtX",XtX
    XtY=dot_product(TG012DS(:,MNO),PRDS)
!    print *, "XtY",XtY
    det=1/(JtJ*XtX-XtJ**2)
!    print *, "det",det
    mu_hat=(XtX*JtY-JtX*XtY)*det
    g_hat=(JtJ*XtY-XtJ*JtY)*det
    F_val=(NOA(TRNO)-1)*abs(g_hat)*XtY/(YtY-abs(g_hat)*XtY)
    tmp=(NOA(TRNO)-1)/((NOA(TRNO)-1)+1*F_val)
    P_val=BETAI((NOA(TRNO)-1)/2.,1/2.,tmp)
    log10P_val=log10(P_val)*(-1.0)
    if(P_val.lt.CLOP) then
      NOSM(TRNO)=NOSM(TRNO)+1
      write(TRNO+30,fmt='(a50,i10,2i15,f15.10)') MN(MNO),OMN(MNO),IOM(MNO,1),IOM(MNO,2),P_val
      write(TRNO+10,fmt=WF5) G012DS(MNO,:)
    endif
    do SG=1, MAXP
      if(log10P_val.ge.(SG-1).and.log10P_val.lt.SG) NOMWSG(SG,TRNO)=NOMWSG(SG,TRNO)+1
    enddo
    if(log10P_val.ge.MAXP) NOMWSG(MAXP+1,TRNO)=NOMWSG(MAXP+1,TRNO)+1
    if(PCHNO.ne.IOM(MNO,1)) then
      print *, "Association test of trait-",NAOTR(TRNO)," using markers on chromosome",PCHNO,"is finished" 
    endif
    PCHNO=IOM(MNO,1)
    write(TRNO+70,fmt='(a50,2i15,5f15.8)') MN(MNO),IOM(MNO,1),IOM(MNO,2),mu_hat, g_hat, F_val, P_val, log10P_val
!    print '(5a15)',"mu_hat","g_hat","F_val","P_val","log10P_val"
!    print *, mu_hat, g_hat, F_val, P_val, log10P_val
  enddo
  deallocate(TG012DS,G012DS,PRDS)
enddo
print *, "Association test of trait-",NAOTR(NOTR)," using markers on chromosome",PCHNO,"is finished"
do SG=1, MAXP
  write(97,fmt=WF3) SG, NOMWSG(SG,:)
enddo
write(97,fmt=WF4) "over"//MAXPF,NOMWSG(MAXP+1,:)
write(97,fmt=WF4) "recordsNo",NOA(:)
write(97,fmt=WF4) "No.SM",NOSM(:)
write(96,fmt='(2a20)') "Trait","Significant markers"
do TRNO=1, NOTR
  write(96,fmt='(a20,i20)') NAOTR(TRNO),NOSM(TRNO)
enddo

end program RRB

!*******************************************************************
!Returns the value ln(Gamma(xx)) for xx>0. Full accuracy is obtained
!for xx > 1. For 0<xx<1, the reflection formula can be used first:
!
!    Gamma(1-z) = pi/Gamma(z)/sin(pi*z) = pi*z/Gamma(1+z)/sin(pi*z)
!******************************************************************* 
Function GAMMLN(xx)
real*8 cof(6),stp,half,one,fpf,x,tmp,ser !internal arithmetic in double precision
data cof,stp/76.18009173D0,-86.50532033D0,24.01409822D0,         &
     -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
data half,one,fpf/0.5d0,1.d0,5.5d0/
  x=xx-one
  tmp=x+fpf
  tmp=(x+half)*log(tmp)-tmp
  ser=one
  do j=1,6
    x=x+one
	ser=ser+cof(j)/x
  end do
  GAMMLN=tmp+log(stp*ser)
  return
End

!********************************************
!Returns the incomplete beta function Ix(a,b)
!********************************************
Function BETAI(a,b,x)
  if(x.lt.0..or.x.gt.1.) print *, ' Bad argument x in function BetaI.'
  if(x.eq.0..or.x.eq.1.) then
    bt=0.
  else
    tmp=gammln(a+b)-gammln(a)-gammln(b)+a*alog(x)+b*alog(1.-x)
    bt=exp(gammln(a+b)-gammln(a)-gammln(b) &
	   +a*alog(x)+b*alog(1.-x))
  end if
  if (x.lt.(a+1.)/(a+b+2.)) then !use continued fraction directly
    BETAI=bt*BETACF(a,b,x)/a
	return
  else
    BETAI=1.-bt*BETACF(b,a,1.-x)/b !use continued fraction after
	return						   !making the symmetry transformation
  end if
End

!******************************************************************
!  Continued fraction for incomplete beta function (used by BETAI)
!******************************************************************
Function BETACF(a,b,x)
Parameter(MAXIT=100,EPS=3e-7,FPMIN=1e-30)
  qab=a+b; qap=a+1.0; qam=a-1.0
  c=1.0
  d=1.0-qab*x/qap
  if (abs(d) < FPMIN) d=FPMIN
  d=1.0/d; h=d
  do m=1, MAXIT
    m2=2*m
	aa=m*(b-m)*x/((qam+m2)*(a+m2))
    d=1.0+aa*d
    if (abs(d) < FPMIN) d=FPMIN
    c=1.0+aa/c
	if (abs(c) < FPMIN) c=FPMIN
	d=1.0/d
	h=h*d*c
    aa=-(a+m)*(qab+m)*x/((qap+m2)*(a+m2))
    d=1.0+aa*d
    if (abs(d) < FPMIN) d=FPMIN
	c=1.0+aa/c
	if (abs(c) < FPMIN) c=FPMIN
	d=1.0/d
	del=d*c
	h=h*del
    if (abs(del-1.0) < EPS) goto 10
  end do !m loop
10 if (m>MAXIT) then 
	print *,' a or b too big, or MAXIT too small.'
  else 
!    print *,' Number of iterations: ', m
  end if
  BETACF=h
  return
END

