program ordering_animal
integer::neff,ntrait,nlev,nsol,io,trait,eff,lev,anieffno,nani,n,m
integer,allocatable::eff_lev(:,:)
real::sol
real,allocatable::solutions(:,:)
character(len=16),allocatable::aniname(:)
character(len=16)::dummy
print *,"number of traits?"
read *, ntrait
print *,"which effect number is animal?"
read *, anieffno

open(1,file="solutions",status='old')
open(3,file="OUTPUT41",status='old')
open(2,file="nsolutions",status='unknown')
open(4,file="rsolutions",status='unknown')
write(2,*) "effect level  sol/trait"
nsol=0;nani=0
do 
  read(1,*,iostat=io)
  if(io.ne.0) exit
  nsol=nsol+1
enddo
rewind(unit=1)
read (1,*)
nsol=nsol-1
nsol=nsol/ntrait
do 
  read(3,*,iostat=io)
  if(io.ne.0) exit
  nani=nani+1
enddo
nani=nani-1
allocate(aniname(nani))
rewind(unit=3)
read (3,*)
m=1
do
  read (3,*,iostat=io) dummy,dummy,dummy,dummy,dummy,aniname(m)
  if(io.ne.0) exit
  m=m+1
enddo
allocate(solutions(nsol,ntrait),eff_lev(nsol,2))
n=1;m=1
do 
  read (1,*,iostat=io) trait, eff, lev, sol
  if(io.ne.0) exit
  solutions(n,trait)=sol
  if(trait.eq.ntrait) then
    write(2,*) eff,lev,solutions(n,:)
    if(eff.eq.anieffno) then
      write(4,*) aniname(m),solutions(n,:)
      m=m+1
    endif
    n=n+1
  endif
enddo
end program ordering_animal

