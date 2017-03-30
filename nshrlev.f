c
c shrlev: a code to calculate eigenvalues of symmetrical, rigid 1d-hindered internal rotation
c copyright (c) 2009 lam t. nguyen and john r. barker
c
c	date:	feb. 13, 2009
c
c john r. barker
c jrbarker@umich.edu
c university of michigan
c ann arbor, mi 48109-2143
c (734) 763 6239
c
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
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine shrlev(t,at,dele,jmax,b,v,n1,n2,imax,zpe,
     &          vhr,bhr,phav,phab)

        implicit double precision(a-h,o-z), integer(i-n)
        parameter (nn=2000, no=1)
        dimension ev(nn), cv(no), cf(no), t(jmax), at(jmax)
        dimension ir(nn)
        character vhr*5, bhr*5

        emax=(jmax-1)*dele + 1000.              ! *****raise max energy to account for later subtraction of zpe******

        ncv=1
        ncf=1
        cv(1)=v
        cf(1)=b

        call ghrlev(ev,nn,imax,emax,b,ncv,ncf,cv,cf,n1,n2,
     &          vhr,bhr,phav,phab)

        zpe=ev(1)

c*********************************************************************
c     the following lines are for output of hindered rotor eigenvalues
c*********************************************************************

        do i=1, imax+1
          tp=ev(i)-zpe
        enddo
99    format (3x, i5, 1x, f10.3 , 1x, f10.3 )
c*********************************************************************

      do i=1, imax
         tp=ev(i)-zpe
         ir(i)=nint( tp/dele )
         ir(i)=nint( tp/dele )
      enddo 

      do j=1, imax
         ll=ir(j)
         do l=1, jmax - ll
            at(l+ll) = at(l+ll) + t(l)
         enddo
      enddo
      do j=1, jmax
         t(j)=at(j)/n1
        at(j)=0.0d0
      enddo 

      return
      end


