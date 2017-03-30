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

      subroutine micro_rates
      include 'declare.inc'
      integer jval
      real*8  viblo
      real*8  tmat(rcnt+ntts,etop)
      real*8  gmat(rcnt+ntts,etop)
      character fname*80,namefile*80
 100  format(i7,2x,100(es15.4))
 200  format(3x,'veff',2x,100(es15.4))

      call sestamp('micro_rates',1)

      call find_jmax(lunit,etopuser,jtopuser,jtop,rcnt,ntts,pcnt,maxtemp
     $,vmax,delh,epsil,barr,delhf,maxj)

c do an initial run to generate the t(r,i) matrix
      jval=0
      call veff(jval,vef)
      do i=1,rcnt+ntts
         call sterab(i,jval,tmat)
      end do

      manymins=.false.

c run over all j values

      write(*,*)'Calculating Microcanonical Rate Contants'
      do jval=0,maxj,dj
         write(*,FMT="(A1,A,t21,F6.2,A)",ADVANCE="NO") achar(13),
     $ "  Percent Complete: ",
     $ ((real(jval)+1.0E-10)/(real(maxj)+1.0E-10))*100.0,
     $ "%"
         call veff(jval,vef)
         do i=1,rcnt+ntts
            call sterabj(i,jval,tmat,gmat)
         end do
      

      if(.false.)then
         if(
     $       (jval.eq.0)    .or.
     $       (jval.eq.maxj) .or.
     $       (mod(jval,jtopuser)).eq.0
     $     ) then
            if(fileroot.eq.'')then
               write(fname,"(a,i3.3)")'gmat-j-',jval
            else
               namefile=trim(fileroot)//"-j-"
               write(fname,"(a,i3.3)")trim(namefile),jval
            end if
            namefile=trim(fname)//".gmat"
            open(unit=12, file=namefile)
            call stamp(12,2)
            call dnt(12)
            write(12,*)
            write(12,*)
            write(12,200)(vef(j),j=rcnt+1,rcnt+ntts)
            do j=1, binmax
               write(12,100)j-1,(gmat(i,j),i=rcnt+1,rcnt+ntts)
            end do
            close(12)
 
         end if
      end if
         call find_ming(gmat,binmax,jval)
      end do
      write(*,*)

      call javerage(viblo)
      call micro_writeout(viblo)
      call calc_kinf

      call sestamp('micro_rates',2)

      return

      end subroutine
