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
      subroutine zpe_calc
      include 'declare.inc'
      integer nn,no,ncv,ncf,n1,n2
      integer ncb,imax
      real*8  bt(datpts),bat(datpts)
      real time
      real*8 zpel(datpts)
      real*8 b,v,ngg,zap
      parameter (nn=2000,no=1)
      real*8 cv(no),cf(no),vo,vm
      real*8 phav,phab,emax,delmin
      character vhr*5,bhr*5
      character*80 efilename
      save

 111  format(5x,22(es20.8))
 222  format(1x,a3,2x,4(a12,2x))
 333  format(1x,i3,2x,4(f12.4,2x))
 444  format(a12,2x,a12,2x,10(a12,2x))
 555  format(9x,i3,9x,f5.2,2x,10(f12.4,2x))

      call sestamp('zpe_calc',1)

c
c create array zpel containing zpe data for 
c all structures rcnt + ntts

c loop over all structures rcnt + ntts
      do i=1,rcnt+ntts+pcnt
         zpel(i)=0.0d0
         do j=1,ndof(i)
            if(idofl(i,j).eq.'vib')then
               zpel(i)=zpel(i) + 0.50d0*(wel(i,j)*ngl(i,j))    ! zpe for vibrational dof is (0.5)*vibrational frequency
            else if(idofl(i,j).eq.'hra')then                   ! zpe for hindered rotors is calculated in ghrlev
               call nrgconvert(anhl(i,j),1,keyword(i,2),'cm-1')
               b  = anhl(i,j)                                  ! rotational constant (cm-1)
               v  = ((wel(i,j)/ngl(i,j))**2)/(b)               ! hindrance barrier (cm-1)
               ngg = ngl(i,j)                                  ! potential energy symmetry (foldedness)
               vhr='vhrd1'
               bhr='bhrd1'
               ncv=1
               ncf=1
               cv(1)=v
               cf(1)=b
               emax=20.0d0*(0.6950356D+00)*maxtemp
               n1=ngg
               n2=1
               phav=0.0d0
               phab=0.0d0
               vo=0.0d0
               vm=0.0d0
c               write(*,*)
c               write(*,*)'shr in ',emax,b,v,n1,n2,valmax(i,j),zap,vhr,
c     $bhr,phav,phab,vo,vm
               call tshrlev(emax,ev,b,v,n1,n2,valmax(i,j),zap,vhr,bhr,
     $                      phav,phab,vo,vm)
c               write(*,*)'shr out',emax,b,v,n1,n2,valmax(i,j),zap,vhr,
c     $bhr,phav,phab,vo,vm
c               write(*,*)
               do k=1,valmax(i,j)
                  evhl(i,j,k)=ev(k)
               end do
               zap = ev(1)
               zpel(i)=zpel(i) + zap
            else if(idofl(i,j).eq.'hrb')then                 ! to be implimented
               v  = anhl(i,j)                                ! hindrance barrier (cm-1)
               b  = (1.0d0/v)*((wel(i,j)/ngl(i,j))**2)       ! rotational constant (cm-1)
               ngg = ngl(i,j)                                ! potential energy symmetry (foldedness)
               vhr='vhrd1'
               bhr='bhrd1'
               ncv=1
               ncf=1
               cv(1)=v
               cf(1)=b
               emax=20.0d0*(0.6950356D+00)*maxtemp
               n1=ngg
               n2=1
               phav=0.0d0
               phab=0.0d0
               vo=0.0d0
               vm=0.0d0
               call tshrlev(emax,ev,b,v,n1,n2,valmax(i,j),zap,vhr,bhr,
     $                      phav,phab,vo,vm)
               do k=1,valmax(i,j)
                  evhl(i,j,k)=ev(k)
               end do
               zap = ev(1)
               zpel(i)=zpel(i) + zap
            else if(idofl(i,j).eq.'hrc')then
               v  = anhl(i,j)                                ! hindrance barrier (cm-1)
               b  = amua2nu/wel(i,j)                         ! rotational constant (cm-1)
               ngg = ngl(i,j)                                ! potential energy symmetry (foldedness)
               vhr='vhrd1'
               bhr='bhrd1'
               ncv=1
               ncf=1
               cv(1)=v
               cf(1)=b
               emax=20.0d0*(0.6950356D+00)*maxtemp
               n1=ngg
               n2=1
               phav=0.0d0
               phab=0.0d0
               vo=0.0d0
               vm=0.0d0
               call tshrlev(emax,ev,b,v,n1,n2,valmax(i,j),zap,vhr,bhr,
     $                      phav,phab,vo,vm)
               do k=1,valmax(i,j)
                  evhl(i,j,k)=ev(k)
               end do
               zap = ev(1)
               zpel(i)=zpel(i) + zap
            else if(idofl(i,j).eq.'hrd')then

               ncb=int(anhl(i,j))
               do k=1, ncb
                  cbb(k)=cbl(i,j,k)
               enddo

               ncv=int(wel(i,j))
               do k=1, ncv
                  cvv(k)=cvl(i,j,k)
               enddo

               call uhrlev(bt,bat,de,etopuser,ncb,ncv,cbb,cvv,
     $              nsvl(i,j),nsbl(i,j),valmax(i,j),zap,vhrl(i,j),
     $              bhrl(i,j),phavl(i,j),phabl(i,j),ngl(i,j),0,ev)

               do k=1, valmax(i,j)
                  evhl(i,j,k)=ev(k)
               end do
               zpel(i)=zpel(i) + zap

            end if
         end do
      end do

      write(lunit,*)'# of structures:',rcnt+ntts+pcnt
      write(lunit,*)
      write(lunit,*)'reporting v, zpe, and electronic (cm-1)'
      write(lunit,*)
      write(lunit,222)'#','v','zpe','electronic'
      write(lunit,*)('-',i=1,45)
      do i=1, rcnt+ntts+pcnt
         write(lunit,333)i,delh(i),zpel(i),(delh(i)-zpel(i))-
     $(delh(1)-zpel(1))
      end do
      write(lunit,*)

      if(rcnt+ntts+pcnt.ge.2)then
         if((inputfile.eq.'').or.(inputfile.eq.'ktools.dat'))then
            efilename='efile.txt'
         else
            efilename=trim(fileroot)//"-efile.txt"
         endif

         efilename=trim(efilename)

         open(unit=efunit,file=efilename)
         call stamp(efunit,2)
         call dnt(efunit)
         write(efunit,*)'run id: ',tstmp
         write(efunit,*)
         write(efunit,*)
         write(efunit,444)'struct#','rxncoord(a)','dh(cm-1)','zpe(cm-1)'
     $,'delect(cm-1)','k-rot(cm-1)','j-rot(cm-1)'
         do i=1,rcnt+ntts+pcnt
            write(efunit,555)i,distl(i),delh(i),zpel(i),
     $(delh(i)-zpel(i))-(delh(1)-zpel(1)),aarr(i),barr(i)
         end do
         close(efunit)
      end if

      call sestamp('zpe_calc',2)

      return

      end subroutine zpe_calc


                                                            
