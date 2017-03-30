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

      subroutine sterab(z,jstep,gmat)
      include 'declare.inc'

      integer z,imax,jarg,nbin,jstep,passno
      integer r(ndof(z)),vcnt,hracnt,nrg
      real*8  summ,rmax
      real*8  emax,b,v(ndof(z)),vv
      real*8,allocatable,dimension(:)::t,at
      integer,allocatable,dimension(:)::ir
      real*8  bb,vvv,ngg,zap
      real*8  hral(ndof(z)),zzpe(ndof(z))
      real*8  gmat(rcnt+ntts,etop)
      real*8  wee,anhh
      integer nng,cnt,anh,we
      logical tt,ff
 
 100  format(i3,2x,f10.2,2x,i5) 


c      
c settinng data for initializingg t and at arrays

      emax = epsil(z)                    ! settinng emax to epsilon for point z
      b=aarr(z)                          ! 1d external rotational constant
      nbin = nint(emax/de) + 1        ! number of bins
      rmax=dsqrt((nbin*de)/b)            ! max "k" for krotor 
      imax=nint(rmax)+1               ! bins for krotor

      if(z.eq.1)nbinl(jstep+1)=nbin      ! max epsilon should always be for the reactants, this corresponds to point z=1


c      
c dynamic allocation of memory      
      allocate(t(nbin))                 
      allocate(at(nbin))
      allocate(ir(imax))
   
c       
c initalize t and at arrays
      do i=1,nbin
         t(i)=0.0d0
        at(i)=0.0d0
      end do  
      t(1)=1.0d0

c
c loop over degrees of freedom

111   continue

      do i = 1 , (ndof(z))
       if    (idofl(z,i).eq.'vib')then                   ! ********  vibrations                 ********
          wee =wel(z,i)
          anhh=anhl(z,i)
          nng =ngl(z,i)
          call morlev(at,t,nbin,de,nbin,wee,anhh,nng,zap)
       elseif(idofl(z,i).eq.'hra')then                   ! ********  quantized hindered rotors  ********
          wee =wel(z,i)
          anhh=anhl(z,i)
          nng =ngl(z,i)
          bb  = amua2nu/anhl(z,i)                        ! rotational constant (cm-1)
          vv  = ((wel(z,i)/ngl(z,i))**2)/(bb)            ! hindrance barrier (cm-1)
          ngg = ngl(z,i)                                 ! potential energy symmetry (foldedness)
          call shrlev(t,at,de,nbin,bb,vv,int(ngg),1,imax,zap,'vhrd1',
     $                                             'bhrd1',0.0d0,0.0d0)
        elseif(idofl(z,i).eq.'hrb')then
          vv=anhl(z,i)                                   ! hindrance barrier (cm-1)
          ngg = ngl(z,i)                                 ! potential energy symmetry (foldedness)
          bb  = ((wel(z,i)/ngl(z,i))**2)/(vv)            ! rotational constant (cm-1)
          call shrlev(t,at,de,nbin,bb,vv,int(ngg),1,imax,zap,'vhrd1',
     $                                             'bhrd1',0.0d0,0.0d0)
        elseif(idofl(z,i).eq.'hrc')then
          bb  = amua2nu/wel(z,i)                        ! rotational constant (cm-1)
          vv  = anhl(z,i)                               ! hindrance barrier (cm-1)
          ngg = ngl(z,i)                                ! potential energy symmetry (foldedness)
          call shrlev(t,at,de,nbin,bb,vv,int(ngg),1,imax,zap,'vhrd1',
     $                                             'bhrd1',0.0d0,0.0d0)
        elseif(idofl(z,i).eq.'hrd')then

          do k=1, ncbl(z,i)
             cbb(k)=cbl(z,i,k)
          end do

          do k=1, ncvl(z,i)
             cvv(k)=cvl(z,i,k)
          end do
          call uhrlev(t,at,de,nbin,ncbl(z,i),ncvl(z,i),cbb,cvv,
     $           nsvl(z,i),nsbl(z,i),valmax(z,i),zap,vhrl(z,i),
     $           bhrl(z,i),phavl(z,i),phabl(z,i),ngl(z,i),1,zap)

        elseif(idofl(z,i).eq.'jro')then
          continue
        elseif(idofl(z,i).eq.'kro')then                   ! ********  1-dim. free k-rotor        ********
          continue
       endif
      end do  

      do i=1,nbin
         gmat(z,i)=t(i)
      end do


c
c release dynamically allocated memory
      deallocate(ir)
      deallocate(at)
      deallocate(t)

      return

      end subroutine
