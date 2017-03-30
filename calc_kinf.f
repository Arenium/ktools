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


      subroutine calc_kinf
      include 'declare.inc'
      real*8  qr(nt),qt(nt),qu(nt)
      real*8 t(maxj+1,maxval(nbinl))
      integer x,y


 100  format(a8,2x,a8  ,2x,2(a12,2x))
 200  format(i8,2x,f8.2,2x,2(es12.5,2x))
 201  format(a8,2x,a8,2x,2(a12,2x))
 300  format(f8.2,2x,3(es12.5,2x))
 301  format(a8,2x,3(a12,2x))

      write(lunit,*)'convergence for reactant'
      x=maxj+1
      y=maxval(nbinl)
      do j=1,x
         do k=1,y
            t(j,k)=rpmat(j,k)
         end do
      end do
      call nkinfint(temps,nt,t,x,y,v1l,sym(1),de,lunit,qr)

      write(lunit,*)      
      write(lunit,*)'convergence for transition state'
      y=binmax
      do j=1,x
         do k=1,y
            t(j,k)=minpmat(j,k)*kejmat(j,k)
         end do
      end do
      call nkinfint(temps,nt,t,x,y,vefmax,sym(2),de,lunit,qt)

      if(manymins)then
         write(lunit,*)      
         write(lunit,*)'convergence for unified transition state'
         do j=1,x
            do k=1,y
               t(j,k)=minpmat(j,k)*ukejmat(j,k)
            end do
         end do
         call nkinfint(temps,nt,t,x,y,vefmax,sym(2),de,lunit,qu)
      end if

      write(lunit,*)
      write(lunit,*)'thermally averaged microcanonical partition functio
     $ns'
      write(lunit,*)
      if(manymins)then
      write(lunit,301)'temp','q(r)','q(ts)','q(uts)'
      write(lunit,*)('-',i=1,51)
      do i=1,nt
         write(lunit,300)temps(i),qr(i),qt(i)*(h/(kbj*temps(i))),qu(i)*
     $(h/(kbj*temps(i)))
      end do

      else

      write(lunit,301)'temp','q(r)','q(ts)'
      write(lunit,*)('-',i=1,51)
      do i=1,nt
         write(lunit,300)temps(i),qr(i),qt(i)*(h/(kbj*temps(i)))
      end do
      end if

      do i=1,nt
        ratekinf(1,i)=qt(i)/qr(i)
        mrate(i)=ratekinf(1,i)
        ratekinf(2,i)=qu(i)/qr(i)
       umrate(i)=ratekinf(2,i)
      end do

      if(manymins)then
         write(lunit,*)
         write(lunit,*)'thermal averaged microcanonical rate constants'
         write(lunit,*)
         write(lunit,201)'#','temp','<kv(e)>t','<ukv(e)>t'
         write(lunit,*)('-',i=1,47)
         do i=1,nt
            write(lunit,200)i,temps(i),ratekinf(1,i),ratekinf(2,i)
         end do
      else
         write(lunit,*)
         write(lunit,*)'thermal averaged microcanonical rate constants'
         write(lunit,*)
         write(lunit,201)'#','temp','<kv(e)>t'
         write(lunit,*)('-',i=1,47)
         do i=1,nt
            write(lunit,200)i,temps(i),ratekinf(1,i)
         end do
      end if

      end subroutine 
