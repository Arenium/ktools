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

      subroutine micro_writeout(viblo)
      include 'declare.inc'
      real*8  viblo
      
      call prewrite('pmat',viblo)
      call prewrite('ming',viblo)
      if(rcnt+ntts.gt.1)then
c         call prewrite('minp')
         call prewrite('kej',viblo)
      end if

      if(manymins)then
         call prewrite('uming',viblo)
         call prewrite('ukej',viblo)
      end if


      end subroutine


      subroutine prewrite(mat,viblo)
      include 'declare.inc'
      character(len=*) mat
      character fnme*15
      real*8   t(maxj+1,maxval(nbinl)),viblo
      integer  x,y


      select case(mat)
         case ('pmat')
            x=maxj+1
            y=maxval(nbinl)
            do i=1,x
               do j=1,y
                  t(i,j)=rpmat(i,j)
               end do
            end do
            fnme='2dens'
         case ('minp')
            x=maxj+1
            y=binmax
            do i=1,x
               do j=1,y
                  t(i,j)=minpmat(i,j)
               end do
            end do
            fnme='crit2dens'
         case ('ming')
            x=maxj+1
            y=binmax
            do i=1,x
               do j=1,y
                  t(i,j)=mingmat(i,j)
               end do
            end do
            fnme='2sums'
         case ('uming')
            x=maxj+1
            y=binmax
            do i=1,x
               do j=1,y
                  t(i,j)=umingmat(i,j)
               end do
            end do
            fnme='u2sums'
         case ('kej')
            x=maxj+1
            y=binmax
            do i=1,x
               do j=1,y
                  t(i,j)=kejmat(i,j)
               end do
            end do
            fnme='kej'
         case ('ukej')
            x=maxj+1
            y=binmax
            do i=1,x
               do j=1,y
                  t(i,j)=ukejmat(i,j)
               end do
            end do
            fnme='ukej'
      end select

      fnme=trim(fnme)
      numnames=numnames+1
      call write_mat(t,x,y,de,dj,v1l,vefmax,inputfile,fnme,rcnt+ntts,
     $               viblo,tstmp,names(numnames))

      end subroutine



