c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c densum: a code for calculating sums and densities of states.
c copyright (c) 2013 john r. barker
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
c	calculates energy levels for 1-dimensional k-rotor to be used
c          with sterab
c
c	dele	= energy grain size
c	jmax	= index of emax
c	b	= rotational constant (cm-1)
c	jk	= quantum number for total angular momentum
c	ir	= vector of energy level indices
c	nmax	= nint of array
c
c
      subroutine krotlev(at,t,nmax,dele,jmax,b,a,jk)
 
      implicit double precision(a-h,o-z)
      double precision at(nmax) , t(nmax)  
      save 
 
      rmax = sqrt(jmax*dele/b)           ! highest rotational quantum number consistent with energy nint
      imax = int(rmax) + 1               ! integer version
      if ( imax .gt. jk ) imax = jk
      if ( imax .gt. 20001 ) imax = 20001
             
      do j = 1, imax                  ! k = quantum number of 1-dimensional free rotor
         r = b*j*j+a*jk*(jk+1)              ! eigenvalues
         r = b*j*j              ! eigenvalues
         ir = nint( r/dele ) + 1            ! nearest integer: number of grains

         if (j .eq. 0 ) then
             f = 1.0                       ! f is the multiplicity of the 1-d rotor (e.g. single-axis internal rotor) energy level (see j.l. mchale, molecular spectroscopy (prentice-hall, 1999), 216f.
         else
             f = 2.0                       ! f is the multiplicity of the 1-d rotor (e.g. single-axis internal rotor) energy level (see j.l. mchale, molecular spectroscopy (prentice-hall, 1999), 216f. 
         endif
c         write(12,*)j,r,ir,b,f

         do k = ir, jmax
            karg = k - ir + 1
            if (karg.le.jmax) at(k) = at(k) + f*t(karg)
         end do  ! k

      end do  ! j

      do j = 1 , jmax
          t(j) = at(j) + t(j)
          at(j) = 0.0d+00
      end do  ! j 

c      write(12,*)
 
      return
      end
