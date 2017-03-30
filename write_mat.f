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


      subroutine write_mat(mat,dimmi,dimmj,de,dj,vr,vt,uname,fname,sw,
     $                     vl,tstmp,nname)
      
      integer i,j,dimmi,dimmj,tunit,dj,sw,one
      real*8 de,vl
      real*8 mat(dimmi,dimmj),vr(dimmi),vt(dimmi)
      character(len=*) fname,uname,nname
      character frmt*60,namefile*60,tstmp*10
      parameter (tunit=24,one=1)
      
      if(uname.eq.'')then
         namefile='dat.'//trim(fname)
      else
         namefile=uname(1:len_trim(uname)-4)//'.'//trim(fname)
      end if

      namefile=trim(namefile)
      nname=namefile
      open(tunit,file=namefile)
c      call stamp(tunit,2)
c      call dnt(tunit)
c      write(tunit,*)
c      write(tunit,*)
      write(tunit,*)'run id: ',tstmp
      write(tunit,*)vl
      if(dimmi.eq.1)then
            write(frmt,'("(14x,i9,3x,i15,3x)")')
      else
            write(frmt,'("(14x,i9,3x,",i0,"(i15,3x))")')dimmi-1
      end if
      write(tunit,frmt)(i-1,i=1,dimmi,dj)
      write(frmt,'("(a14,",i0,"(es15.6,3x))")')dimmi

      if(sw.gt.1)then
         if(fname.eq.'kej'.or.fname.eq.'ukej')then
            write(tunit,frmt)'Erot-Re',(vr(i),i=1,dimmi,dj)
            write(tunit,frmt)'Erot-TS',(vt(i),i=1,dimmi,dj)
         elseif(fname.eq.'2dens')then
            write(tunit,frmt)'Erot-Re',(vr(i),i=1,dimmi,dj)
         elseif(fname.eq.'2sums')then
            write(tunit,frmt)'Erot-TS',(vt(i),i=1,dimmi,dj)
         elseif(fname.eq.'u2sums')then
            write(tunit,frmt)'Erot-TS',(vt(i),i=1,dimmi,dj)
         end if
      else
         if(fname.eq.'2dens')then
            write(tunit,frmt)'Erot',(vr(i),i=1,dimmi,dj)
         else if(fname.eq.'2sums')then
c            write(tunit,frmt)'Erot',(vt(i),i=1,dimmi,dj)
            write(tunit,frmt)'Erot',(vr(i),i=1,dimmi,dj)
         end if
      end if
      
      write(tunit,*)
      write(frmt,'("(f10.1,4x,",i0,"(es15.6,3x))")')dimmi
      do j=1,dimmj
         write(tunit,frmt)(j-1)*de,(mat(i,j),i=1,dimmi,dj)
      end do
      write(tunit,*)
      close(tunit)

      end subroutine


