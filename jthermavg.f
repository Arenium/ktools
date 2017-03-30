c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c copyright (c) 2017 john r. barker, jason a. sonk
c
c john r. barker
c jrbarker@umich.edu
c university of michigan
c ann arbor, mi 48109-2143
c (734) 763 6239
c
c this program is free software; you can redistribute it and/or
c modify it under the terms of the gnu general public license (version 2)
c as published by the free software foundation.
c
c this program is distributed in the hope that it will be useful,
c but without any warranty; without even the implied warranty of
c merchantability or fitness for a particular purpose. see the
c gnu general public license for more details.
c
c see the 'readme' file for a copy of the gnu general public license,
c or contact:
c
c free software foundation, inc.
c 59 temple place - suite 330
c boston, ma 02111-1307, usa.
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine jthermavg
      include 'declare.inc'
      real*8  nusum,desum,term
      real*8  jpmat(binmax,nt),jgmat(binmax,nt)
      real*8  jkmat(binmax,nt)
 100  format(i5,2x,f10.2,2x,100(ES15.6,2x))
 200  format(19x,100(F15.2,2x))

      do k=1,nt
         do i=1,binmax
            desum=0.0d0
            do j=1, maxj+1
               term=Barr(1)*j*(j+1)
               term=-1.0d0*term
               term=term/(kb*Temps(k))
               term=Dexp(term)
               term=minpmat(j,i)*term
               desum=desum+term
            end do
            jpmat(i,k)=desum
         end do

         do i=1,binmax
            nusum=0.0d0
            do j=1, maxj+1
               term=Barr(1)*j*(j+1)
               term=-1.0d0*term
               term=term/(kb*Temps(k))
               term=Dexp(term)
               term=mingmat(j,i)*term*c
               nusum=nusum+term
            end do
            jgmat(i,k)=nusum
         end do

         do i=1,binmax
            if(i.LT.(binmax-binmax))then
               jkmat(i,k)=0.0d0
            else
               jkmat(i,k)=(jgmat(i,k)/jpmat(i,k))
            end if
         end do
      end do

c j-averaged k(E) for each temperature
c      open(unit=junit,file="javge.txt")
c      write(junit,200)(temps(k),k=1,nt)
c      do i=1,binmax
c         write(junit,100)i,de*(i-1),(jkmat(i,k),k=1,nt)
c      end do
c      close(junit)

      end subroutine
      
         
