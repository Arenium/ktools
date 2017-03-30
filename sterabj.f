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

      subroutine sterabj(z,jstep,tmat,gmat)
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
      real*8  tmat(rcnt+ntts,etop)
      real*8  wee,anhh,twee
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
         t(i)=tmat(z,i)
        at(i)=0.0d0
      end do  

c
c loop over degrees of freedom

111   continue

      do i = 1 , (ndof(z))
       if(idofl(z,i).eq.'kro')then                   ! ********  1-dim. free k-rotor        ********
          wee =wel(z,i)
          anhh=anhl(z,i)
          nng =ngl(z,i)
       elseif(idofl(z,i).eq.'jro')then
          twee=wel(z,i)
       end if
       
      end do 
      wee=wee-twee
c      call krotlev(at,t,nbin,de,nbin,wee,jstep)
      call krotlev(at,t,nbin,de,nbin,wee,twee,jstep)

c
c include 2j+1 external rotation contribution
      do i=1, nbin
         t(i)=((2.0d0*jstep)+1.0d0)*t(i)
      end do

c  
c sum up densities of states in t and store in vector at
      summ=0.0d0
      do i=1,nbin
         summ=summ+t(i)
         at(i)=summ
         t(i)=t(i)/de
      end do

c
c find sums above veff max (ecrit) and store in gmat, store number of grains above vmax in smat, and densities in minpmat
      cnt=0
      if(z.eq.1)then                            ! store dos of reactant
         do i=1,nbin
            rpmat(jstep+1,i)=t(i)                              ! store full dos
            nrg=(i-1)*de
            if(nrg.ge.((vmax-vef(z))))then                    ! store info above critical energy
               cnt = cnt + 1
               if(cnt.gt.binmax)exit
               minpmat(jstep+1,cnt)=t(i)
               gmat(z,cnt)=at(i)
            end if
         end do
      else
         do i=1,nbin                                          ! store sos of trial ts
            nrg=(i-1)*de
            if(nrg.ge.((vmax-vef(z))))then
               cnt = cnt + 1
               if(cnt.gt.binmax)exit
               gmat(z,cnt)=at(i)
            end if
         end do
      end if

c
c release dynamically allocated memory
      deallocate(ir)
      deallocate(at)
      deallocate(t)

      return
      end subroutine
