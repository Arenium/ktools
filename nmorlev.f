c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c densum: a code for calculating sums and densities of states.
c copyright (c) 2001 john r. barker
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
c see the "readme" file for a copy of the gnu general public license,
c or contact:
c
c free software foundation, inc.
c 59 temple place - suite 330
c boston, ma 02111-1307, usa.
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c
c	calculates energy levels for morse oscillator to be used
c          with sterab
c
c	dele	= energy grain
c	jmax	= index of emax (top of energy space)
c	we	= harmonic vibration frequency (cm-1)
c	xe	= anharmonicity (cm-1)
c           ng          = degeneracy of vibration
c	ir	= vector of energy level indices
c	imax	= nint of array (e.g. dissociation energy for morse osc.)
c           zpe         = zero point energy (cm-1), including anharmonicity
c
c	input observed weo = we + 2*xe for 0-1 transition
c	note sign of xe:  negative for usual morse oscillator
c
c	e = we*(v+1/2) + xe*(v+1/2)^2
c
c     all energies and frequencies in wavenumbers, relative to zpe = we/2 + xe/4
c
c    2/02 modification: now, hindered rotors are counted as vibs in whitten-rabinovitch
c         parameter calculations. zero point energy is calculated based on anharmonic
c         vibrations and from the hindered rotor subroutine.
c
c    9/07 fixed small bug in calculation of rmax 
c
c--------------------------------------------------------------------------------
      subroutine morlev(at,t,nmax,dele,jmax,we,xe,ng,zpe)
      implicit none
      double precision at , t , dele
      dimension at(nmax) , t(nmax)
      integer i , ir , nmax , jmax , ng
      dimension ir(20001)
      double precision we , xe , zpe, w0 , rmax , vi , r
      integer imax , j , k , karg , l
      external sterab , gam , wrab , wrden , alngam
      save 
 
      zpe = 0.5d+00*we + 0.25d+00*xe                    ! zero point energy at v=0

      if ( xe .lt. 0.0 ) then
         w0 = we + xe
         rmax = -0.5d+00*w0/xe                          ! highest bound state, e relative to zpe
      else
         rmax = jmax*dele/we                            ! highest bound state
      endif
      imax = int(rmax) 
      if ( imax.gt.20001 ) imax = 20001

      do i = 1 , imax                               ! start at v=1
        vi = i + 0.5d+00
        r = (we + xe*vi)*vi - zpe                   ! state energy relative to zpe
c        ir(i) = 2 + int( r/dele )                   ! index of bin above zpe
        ir(i) = 1 + nint( r/dele )                ! nearest integer number of grains
        ir(i) = 1 + nint( r/dele )                ! nearest integer number of grains
      end do ! i

      do l = 1 , ng                                    ! loop for degenerate vibrations
         do j = 1 , imax                                  ! imax is the number of energy states
            do k = 1 , jmax                                  ! jmax is the number of energy grains
               karg = k - ir(j) + 1
               if ( karg.gt.0 .and. karg.le.jmax ) then
                 at(k) = at(k) + t(karg)               ! at(k) is the number of states in the kth grain
               endif
            end do  ! k
         end do  ! j
 
         do j = 1 , jmax
            t(j) = at(j) + t(j)
            at(j) = 0.0d+00
         end do
      end do   ! l
 
      return
      end
