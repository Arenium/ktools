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

      subroutine find_ming(gmat,n,jval)
      include 'declare.inc'
      integer jval,fend,minct,minil(rcnt+ntts)
      real*8  gmat(rcnt+ntts,n),mingl(rcnt+ntts)
      real*8  x(3),y(3),ming
      real*8  tdat(rcnt+ntts),urate(n)
      integer tmins(rcnt+ntts),fct(n)
      integer tempx,rx,lx
      character frmt*60,namefile*80
 100  format(i3,2x,f10.2,2x,i5)
 200  format(2x,a1,6x,a6,2x,a5,2x,4(a10,2x))
 300  format(i5,2x,4(es15.5,2x))

c
c   if there are 2 or more trial ts open nmins file
c      if(ntts.ge.2)then         
c         if(inputfile.eq.'')then
c            namefile='nmins.txt'
c         else
c            namefile=trim(fileroot)//"-nmins.txt"
c         end if
c
c         open(unit=nminsunit,file=namefile)
c         if(jval.eq.0)then
c            call stamp(nminsunit,2)
c            call dnt(nminsunit)
c            write(nminsunit,*)'run id: ',tstmp
c            write(nminsunit,*)
c            write(nminsunit,*)
c            write(nminsunit,200)'j','e cm-1','# min','found@','min g'
c         end if
c      end if

c
c   check for number of trial ts
      if(ntts.le.1)then                                  ! only 1 ts structure given
         do i=1,n                                        ! loop over all energy grains
            mingmat(jval+1,i)=gmat(rcnt+ntts,i)
           umingmat(jval+1,i)=gmat(rcnt+ntts,i)
         end do

      else                                               ! multiple ts given
c         write(*,*)'jval',jval 
         do i=1,n                                        ! loop over all energy grains

            fct(i)=0                                     ! found min count

            do j=rcnt+1,rcnt+ntts                        ! only over ntts
               tdat(j-rcnt)=gmat(j,i)
            end do
           
c            write(*,*)jval,i 
            call find_tmins(tdat,tmins,ntts,fct(i),microh)


            if(fct(i).gt.1)manymins=.true.            

            urate(i)=0.0d0
            ming=gmat(rcnt+tmins(1),i)
            do j=1,fct(i)
               minil(j)=rcnt+tmins(j)
               mingl(j)=gmat(rcnt+tmins(j),i)
               urate(i)=urate(i)+(1.0d0/mingl(j))
               if(mingl(j).lt.ming)ming=mingl(j)             
            end do
            urate(i)=1.0d0/urate(i)
           
c            if(ntts.ge.2)then 
c            write(frmt,'("(i3,2x,f10.2,2x,i3,2x,",i0,"(i3,2x),",i0,
c     $           "(es10.2,2x),i3)")')fct(i),fct(i)
c            write(nminsunit,frmt)jval,(i-1)*de,fct(i),
c     $           (minil(l)-rcnt,l=1,fct(i)),(mingl(l),l=1,fct(i))
c            end if

            mingmat(jval+1,i)=ming
           umingmat(jval+1,i)=urate(i)            


         end do

      end if

      return
      end subroutine
