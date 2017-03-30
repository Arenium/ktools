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
      subroutine jaread_input(inputfile,imax1,isize,temax2)
      implicit none
      character(*) inputfile
      character filetype*4,fname*80
      real*8,allocatable,dimension(:,:)::mat
      real*8,allocatable,dimension(:)::khsg,sumdens
      real*8 emax,emax1,emax2,de,rotjay(1001),e,egrain2,egrain1,temax2
      real*8 viblo
      character(60000)line
      character(15)line2,junk,arg1,arg2,arg3
      character   sumtype
      logical there
      integer i,j,jdex,nbin,tnbin
      integer iunit,ounit,ccount,imax1,isize
      integer jays(1001),maxj,eol1,eol2,reason,lnct

 100  format('*** ',A,' not found ***')
 200  format('Unsupported file. Please use a .2dens or .2sums file')
 300  format(i10,1x,f10.1,2(1x,1pe12.5))
 400  format(5(F10.2,2x)) 
 500  format(1x,A,2x,F10.2)
 600  format(1x,A,2x,I10)

      iunit=111
      ounit=222

c
c  Get name of file to open

c      call get_command_argument(1,inputfile)
c      ccount=command_argument_count()
      ccount=10
      write(*,*)
      write(*,*)"File:",trim(inputfile)
      filetype=inputfile(len_trim(inputfile)-3:len_trim(inputfile))
      if(filetype.eq.'sums')then
         sumtype=inputfile(len_trim(inputfile)-5:len_trim(inputfile)-4)
      end if

c
c  Grab command line arguments for imax and isize
c      if(ccount.gt.1)then
c         call get_command_argument(2,arg1)
c         call get_command_argument(3,arg2)
c         read(arg1,*)imax1
c         read(arg2,*)isize
         write(*,600)'imax1:',imax1
         write(*,600)'isize:',isize
c         if(ccount.gt.3)then
c            call get_command_argument(4,arg3)
c            read(arg3,*)temax2
            write(*,500)'emax2:',temax2
c         end if
c      end if
      

c
c  Check if file is present
 111  inquire(file=inputfile,exist=there)
      if(.not.there)then
          write(*,100)trim(inputfile)
          stop
      end if

c
c  Check if valid file type

      if(trim(filetype).eq.'dens')then
         fname=inputfile(1:len_trim(inputfile)-6)//"-r."//trim(filetype)
      else if(trim(filetype).eq.'sums')then
            if(sumtype.eq.'u')then
               fname=inputfile(1:len_trim(inputfile)-7)//"-uts.dens"
            else
               fname=inputfile(1:len_trim(inputfile)-6)//"-ts.dens"
            end if
      else
           write(*,200)
           stop
      end if
      
c
c  Find max J value present

      open(unit=iunit,file=inputfile)
      lnct=0
      read(iunit,*)junk
      read(iunit,*)viblo
      lnct=lnct+1

      read(iunit,"(A)",IOSTAT=eol1)line
      read(line,*,IOSTAT=eol2)jays
      maxj=maxval(jays)+1
c      write(*,*)maxj

c
c  Read in rotational energies
      read(iunit,"(A)",IOSTAT=eol1)line
      read(line,*,IOSTAT=eol2)junk,rotjay
      lnct=lnct+1


c
c  Skip blank line

      read(iunit,*)     
      lnct=lnct+1
      read(iunit,*,iostat=reason)emax1
      lnct=lnct+1
      read(iunit,*,iostat=reason)emax2
      lnct=lnct+1
      de=emax2-emax1
c      write(*,*)emax1,emax2,de

c
c  Find max E present
      do
        read(iunit,*,iostat=reason)emax
        lnct=lnct+1
        if(reason.lt.0)goto 700
      end do
 700  close(iunit)

 
c
c  Set values for nbin, egrain2
      tnbin=int(emax/de)+1
      if(temax2.ne.emax)then
         nbin=int(temax2/de)+1
         Egrain2=temax2/(isize-imax1-1)
      else
        nbin=int(emax/de)+1
        Egrain2=emax/(isize-imax1-1)
      end if

c
c  Allocate array space
      allocate(mat(maxj,tnbin))
      allocate(khsg(tnbin))
      allocate(sumdens(tnbin))

c
c  Re-openfile

      open(unit=iunit,file=inputfile)

c
c  skip 4 lines

      do i=1,4
         read(iunit,*)
      end do

c
c  read in full matrix
      do i=1,tnbin
         read(iunit,*)line,(mat(j,i),j=1,maxj)
         khsg(i)=0.0D0
      end do
      close(iunit)

c
c  do j-summing placing all mat(j,i) into khsg(i+jdex)

      do j=1,maxj
         jdex=nint((rotjay(j)-rotjay(1))/de)
         do i=1,tnbin
            if((i+jdex).le.tnbin)then
                khsg(i+jdex)=khsg(i+jdex)+mat(j,i)
            end if
         end do
      end do
       

c      if(ccount.gt.1)then
         call chkdens(khsg,tnbin,dE,imax1,tnbin,de*(tnbin-1))
c      end if

c
c  calculate sums or dens
      if(trim(filetype).eq.'dens')then
         do i=1,tnbin
            if(i.eq.1)then
                sumdens(i)=khsg(i)*de
            else
                sumdens(i)=sumdens(i-1)+khsg(i)*de
            end if
         end do
      else
         do i=1,tnbin
            if(i.gt.1)then
                sumdens(i)=(khsg(i)-khsg(i-1))/de
            else
                sumdens(i)=khsg(i)/de
            end if
         end do
      end if

c
c write out .dens file
c data order is:
c Bin Number, Energy of Bin, Density of Bin, Sum of Bins.      
      open(unit=ounit,file=fname)
      call japrewrite(ounit,fname,de,imax1,temax2,isize,viblo)
      if(ccount.lt.2)then
         if(trim(filetype).eq.'dens')then
            do i=1,tnbin
               write(ounit,300)i,de*(i-1),khsg(i),sumdens(i)
            end do
         else
            do i=1,tnbin
               write(ounit,300)i,de*(i-1),sumdens(i),khsg(i)
            end do
         end if
      else
         if(trim(filetype).eq.'dens')then
            do i=1,imax1
               E = (i-1)*de
               if(i.gt.nbin)then
                  j=nbin
                  write(ounit,300)i,e,khsg(j),sumdens(j)
               else
                  write(ounit,300)i,e,khsg(i),sumdens(i)
               end if
            end do
   
c
c  do j-summing placing all mat(j,i) into khsg(i+jdex)
         do i=imax1+1, isize
               E = (i-imax1-1)*Egrain2
               j = int(e/de) + 1
               if(j.gt.tnbin)then
                  j=tnbin
                  write(ounit,300)i,e,khsg(j),sumdens(j)
               else
                  write(ounit,300)i,e,khsg(j),sumdens(j)
               end if
            end do
         else
            do i=1,imax1
               E = (i-1)*de
               if(i.gt.nbin)then
                  j=nbin
                  write(ounit,300)i,e,sumdens(j),khsg(j)
               else
                  write(ounit,300)i,e,sumdens(i),khsg(i)
               end if
            end do
            do i=imax1+1, isize
               E = (i-imax1-1)*Egrain2
               j = int(e/de) + 1
               if(j.gt.tnbin)then
                  j=tnbin
                  write(ounit,300)i,e,sumdens(j),khsg(j)
               else
                  write(ounit,300)i,e,sumdens(j),khsg(j)
               end if
            end do
         end if
      end if
 
      deallocate(sumdens)
      deallocate(khsg)
      deallocate(mat)

      end subroutine
