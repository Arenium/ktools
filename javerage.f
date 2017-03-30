c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c copyright (c) 2014 john r. barker, jason a. sonk
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

      subroutine javerage(viblo)

      include 'declare.inc'
      integer jdex
      real*8    gmat(rcnt+ntts,epsmax),khsp(5*etop),totp,viblo
      real*8    khsg(5*etop),ikhsg(5*etop)
      real*8    ukhsg(5*etop),uikhsg(5*etop)
 100  format(i7,2x,30(es15.4))
 200  format("**************input data summary**************")
 201  format(a,' data as generated from ',a)
 202  format(f10.2,2x,i10,2x,f10.2,2x,i10,2x,f10.2)
 203  format(a6,2x,a10,2x,5(a15))
 204  format(i6,2x,f10.2,2x,5(es15.6,2x))



c
c   calculate rate constants from sums and densities of states

      if(rcnt+ntts.eq.1)then                                            !  if only 1 species is passed write out .dens file for that species

         do j=1,maxj+1                                                  !  klippenstein, harding, smith, & gilbert j averaging 
            do i=1,maxval(nbinl)
               jdex=nint((barr(1)*j*(j-1))/de)
               khsp(i+jdex)=khsp(i+jdex)+rpmat(j,i)
            end do
         end do

         call get_viblo(rcnt+ntts,rcnt+ntts,viblo)

      else if(rcnt+ntts.ge.2)then                                       !  if r and ts passed write out .dens file for each species

         do j=1,maxj+1                                                  !  klippenstein, harding, smith, & gilbert j averaging 
            do i=1,maxval(nbinl)
               jdex=nint((barr(1)*j*(j-1))/de)
               khsp(i+jdex)=khsp(i+jdex)+rpmat(j,i)

               jdex=nint((vefmax(j)-vefmax(1))/de)
               khsg(i+jdex)=   khsg(i+jdex) +   mingmat(j,i)
              ukhsg(i+jdex)=  ukhsg(i+jdex) +  umingmat(j,i)
                 
            end do
         end do

         call get_viblo(rcnt,rcnt,viblo)


         do j=1,maxj+1,dj
            do i=1,binmax
                 kejmat(j,i)=(   mingmat(j,i)/minpmat(j,i))*c
                ukejmat(j,i)=(  umingmat(j,i)/minpmat(j,i))*c
            end do
         end do

         call jthermavg

      end if


      end subroutine
