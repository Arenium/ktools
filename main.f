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
      program ktools
      include 'declare.inc'
      real start,finish
      integer values(9),sortcnt,reason
      character dname*80
 100  format('elapsed time:',F10.2,' seconds.')
 101  format('elapsed time:',F10.2,' minutes.')
 102  format('elapsed time:',F10.2,' hours.')
 103  format(i3,2x,F9.2,2x,4(ES15.5,2x))
 104  format(a3,2x,a6,3x,3(a15,2x),a14,2x)
 107  format(i1,i3.3,3(i2.2))

      start = secnds(0.)
      sortcnt=0
      call ltime(time(),values)
     
c  read input, check dofs, calculated mw's, convert energies to wavenumbers      

      call read_input
      write(*,*)'sorting data'
      if(rcnt+ntts+pcnt.gt.1)call dsort(sortcnt)
      call check_input
      if(trim(backup).eq.'backup')call write_input
      call super_mol
      write(*,*)'converting units'
      call bsort
      call hsort
      write(*,*)'calculating zpe'
      call zpe_calc
c      call multiwell_write
c  canonical rate constant subroutine

      call canon_rates
      write(*,*)'partition function calculations done'
      if(whatdo.eq.'canon')goto 200
     

c  microcanonical rate subroutine

      call micro_rates

c  writeout canonical data to the end of the canonical file
 200  if(rcnt+ntts+pcnt.gt.1)call finalprint

c  final time stamp and logfile closed

      call sestamp('main',2)
      close(canunit)
      finish = secnds(start)
      if(finish.lt.60)then
            write(lunit,100)finish
      else if (finish.lt.3600)then
            write(lunit,101)finish/60
      else
            write(lunit,102)finish/3600
      end if
      
      close(lunit)
     
      write(*,*)'Calling Jasum'
      do i=1,numnames
         if(names(i)(len_trim(names(i))-2:len_trim(names(i)))
     $.ne."kej")then
            call jaread_input(names(i),dimax1,disize,demax2)
         end if
      end do
      end program

