
      subroutine ghrlev(ev,nn,nmax,emax,b,ncv,ncf,cv,cf,
     &            nsiv,nsif,vhr,bhr,phav,phab)
 
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c ghrlev: a code of general purposes to calculate eigenvalues of 1d-hindered internal rotation 
c copyright (c) 2009 lam t. nguyen
c
c     date: feb. 13, 2009
c
c lam t. nguyen
c nguyenlt@umich.edu
c university of michigan
c ann arbor, mi 48109-2143
c
c or contact
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
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c     input:      
c           1) emax           : maximum energy in master equation
c           2) b        : average rotational constant, defined as:
c                 b = (integral of b(x)dx from 0 to 2*pi) / (2*pi)
c           3) ncv, ncf : number of elements in vectors of cv and cf, respectively
c           4) cv       : vector of coefficients with ncv elements in torsional 
c                 potential energy function, which is expressed as = cvo/2 + sum of cvi*(1-cos(i*x*nsig))/2
c           5) cf       : vector of coefficients with ncf elements in rotational constant,
c                 which is given by a fourier series = cfo + sum of cfi*cos(i*x*nsig)
c           6) nsiv, nsif     : symmetry numbers for tortional potential e function and rotational constant, respectively.
c           7) nn       : maximum number of elements in vector ev
c           8) vhr            : model for torsional potential energy function
c           9) bhr            : model for rotational constant
c           10) phav, phab    : phases of torsional angles in vhr and bhr, respectively.         
c
c     output:     1) ev(nmax) : vector of eigenvalues with nmax elements
c           2) nmax           : number of elements of vector ev
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      implicit double precision(a-h,o-z), integer(i-n)
      parameter (nx=501)
c
c.....declare arrays.
      dimension ev(nn), cv(ncv), cf(ncf)
        character vhr*5, bhr*5
c
c.....calvo used to find vmax (vo) and vmin (vm) of torsional potential energy function
      call calvo(vo,vm,cv,ncv,nsiv,vhr,phav)

      call odqhr(ev,nn,ncv,ncf,cv,cf,nsiv,nsif,vhr,bhr,phav,phab,vm)
      nmax=nx
      call calei(ev,nn,nmax,emax,b)

        zpe = ev(1)
c        write(*,9991) vo, vm, ev(1)
9991    format(11x,'torsional potential: ',1x,'vmax (cm-1) = ',f8.1,1x,
     &  '; vmin (cm-1) = ',f8.1,1x,'; zpe (cm-1) = ',f5.1)
     
c*********************************************************************
c     the following lines are for output of hindered rotor eigenvalues
c*********************************************************************
c        write(12,*) '   '
c        write(12,9991) vo, vm, zpe

c        do i=1, imax+1
c          write(12, 99) i, ev(i), (ev(i)-zpe)
c        enddo
c99      format (3x, i5, 1x, f10.3 , 1x, f10.3 )
c*********************************************************************
     


      return 
      end


c****************************************************************************************************************
c     name: subroutine odqhr 
c     use:  to obtain a vector of quantum eigenvalues of 1d-hindered internal rotor, in which  
c           both torsional potential energy function and moment inertia are explicitly treated as 
c           functions of internal rotational angle, i.e. non-rigid rotor, based on meyer's algorithm:
c           r. meyer j. chem. phys. 52 (1970) 2053-2059.  
c
c     authors:    lam t. nguyen and john r. barker 
c
c     date: jan. 26, 2009
c
c****************************************************************************************************************
c     
c      nx  : number of grid points
c        xa  : position on a grid.
c        ar  : hamiltonian matrix, i.e. meyer's symmetrical matrix
c        er  : vector of eigenvalues 
c        zr  : eigenvectors (x,y); where
c                                      x : wavefunction.
c                                      y : energy level.       
c        npr      : 0 for eigenvalue only and 1 for both eigenvalue and eigenvector 
c        rmin     : starting point of grid
c        rmax     : end point of grid
c        zl       : grid length
c        dx       : grid spacings
c
c****************************************************************************************************************
c
c     input:      1) ncv, ncf = number of coefficients derived from fitting rotational energy profile and moment inertia
c           2) cv(ncv)  = vector of coefficients from rotational energy profile, with ncv elements
c           3) cf(ncf)  = vector of coefficients from moment inertia (rotational constant), with ncf elements
c           4) nsig           = rotational symmetry number
c     output:     1) er             = vector of eigenvalues of 1d-hindered internal rotor, with nx=501 elements 
c
c*****************************************************************************************************************

      subroutine odqhr(er,nn,ncv,ncf,cv,cf,nsiv,nsif,
     &      vhr,bhr,phav,phab,vm)
      implicit double precision(a-h,o-z), integer(i-n)
      parameter (nx=501)
c
c.....declare arrays.
      dimension ar(nx,nx),er(nn),zr(nx,nx),xa(nx) 
      dimension d(nx,nx), f(nx), v(nx), dt(nx,nx),
     & t(nx,nx), ft(nx,nx)
      dimension  cf(ncf), cv(ncv), fv1(nx), fv2(nx)
        character vhr*5, bhr*5

c.....      dimension work(nx*3)

c
c.....variable data input.

      data npr/0/

c.....now compute hamiltonian matrix:

      call cmm(ar,ncv,ncf,cv,cf,nsiv,nsif,vhr,bhr,phav,phab,vm)
                              
c.....now call eigenvalue solver using rs subroutine of eispack.

      call rs(nx,nx,ar,er,npr,zr,fv1,fv2,ierr) 

c...... alterative option: using subroutine dsyev from lapack 
c......     call dsyev( 'v', 'u', nx, ar, nx, er, work, nx*3, ierr )
c......

      return 
      end

cccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccc
c.....      cmm stands for constructing meyer's matrix, ar, for 1d-hindered internal rotor using meyer's method, 
c.....      r. meyer j. chem. phys. 52 (1970) 2053.
c.....      both rotational energy profile and moment inertia are properly treated as functions of 
c.....      internal rotation angle, i.e. non-rigid rotor
c.....      input:      1) ncf, ncv = number of coefficients drived from fitting rotational energy profile and moment inertia
c.....            2) cv(ncv) = vector of coefficients from rotational energy profile, with ncv elements
c.....            3) cf(ncf) = vector of coefficients from moment inertia, with ncf elements
c.....            4) nsig = rotational symmetry number, it can be either one or two or three or whatever 
c.....      output: 1) ar(nx,nx) = meyer's symmetrical matrix, with a size of nx--number of grid points   


      subroutine cmm(ar,ncv,ncf,cv,cf,nsiv,nsif,vhr,bhr,phav,phab,vm)
      implicit double precision(a-h,o-z), integer(i-n)
      parameter (nx=501)
c
c.....declare arrays.
      dimension ar(nx,nx),xa(nx) 
      dimension d(nx,nx), f(nx), v(nx), dt(nx,nx),
     & t(nx,nx), ft(nx,nx)
      dimension  cf(ncf), cv(ncv)
      character vhr*5, bhr*5

c
c.....variable data input.

      pi=acos(-1.d0)
c
c.....set up grid
      rmin=0.0d0
      rmax=pi*2
      zl=(rmax-rmin)
      dx=zl/nx

c.....now compute meyer's symmetrical matrix:

      do i=1, nx
            xa(i)=rmin+dx*i
            call vx(v(i),xa(i)+phav,cv,ncv,nsiv,vhr,vm)
            call fx(f(i),xa(i)+phab,cf,ncf,nsif,bhr)
            do j=1, nx
                  imj=i-j
                  if(imj.eq.0) then
                        d(i,j)=0.0d0
                  else
                        tp1=pi*imj/dble(nx)
                        tp2=(-1.0d0)**imj
                        d(i,j)=tp2/2.0d0/dsin(tp1)
                  endif
                  dt(i,j)=-d(i,j)
                  ft(i,j)=f(i)*delta(i,j)
            enddo
      enddo

      do i=1, nx
            do j=1, nx
                  su=0.0d0
                  do k=1, nx
                        su=su + ft(i,k)*d(k,j)
                  enddo
                  t(i,j)=su
            enddo
      enddo

      do i=1, nx
            do j=1, nx
                  su=v(i)*delta(i,j)
                  do k=1, nx
                        su=su + dt(i,k)*t(k,j)
                  enddo
                  ar(i,j)=su
            enddo
      enddo 

      return 
      end
c..........
c..........
      subroutine vx(v,x,cv,n,nsig,vhr,vm)
      implicit double precision(a-h,o-z), integer(i-n)
      dimension cv(n)
      character vhr*5

      if((vhr.eq.'vhrd1').or.(vhr.eq.'vhrd1')) then
            v=0.0d0
            do i=1, n
                  v=v+cv(i)*(1.0d0-dcos(x*nsig*i))/2.0d0
            enddo
      elseif((vhr.eq.'vhrd2').or.(vhr.eq.'vhrd2')) then
            v=cv(1)
            do i=2, n
                  v=v+cv(i)*dcos(x*nsig*(i-1))
            enddo
        elseif((vhr.eq.'vhrd3').or.(vhr.eq.'vhrd3')) then
            v=cv(1)
            nc=(n+1)/2
            do i=2, n
                  if(i.le.nc) then
                        v=v+cv(i)*dcos(x*nsig*(i-1))
                  else
                        v=v+cv(i)*dsin(x*nsig*(i-nc))
                  endif
            enddo
      else
            write(*,*) "error at input for torsional function = vhrd"
            stop
      endif

      v=v-vm

      return
      end
c..........
c..........
      subroutine fx(f,x,cf,n,nsig,bhr)
      implicit double precision(a-h,o-z), integer(i-n)
      dimension cf(n)
      parameter (con=16.85763d0)
      character bhr*5

      if((bhr.eq.'bhrd1').or.(bhr.eq.'bhrd1')) then
            f=cf(1)
            do i=2, n
                  f=f+cf(i)*dcos(x*nsig*(i-1))
            enddo
        elseif((bhr.eq.'ihrd1').or.(bhr.eq.'ihrd1')) then
                f=cf(1)
                do i=2, n
                        f=f+cf(i)*dcos(x*nsig*(i-1))
                enddo
            f=con/f
      else
                write(*,*) "error at input for bhrd or ihrd"
            stop
      endif

      return
      end
c.........
c.........      
c      double precision function delta(i,j)
c      integer i, j
c
c      if(i.eq.j) then
c            delta=1.0d0
c      else
c            delta=0.0d0
c      endif
c      
c      return
c      end
c.........
c.........

ccccccccccccccccccccc
ccccccccccccccccccccc

      subroutine calei(ev,nn,nmax,emax,b)
 
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c calei: a code to compute/select eigenvalues of torsional motion based on emax 
c copyright (c) 2009 lam t. nguyen and john r. barker
c
c john r. barker
c jrbarker@umich.edu
c university of michigan
c ann arbor, mi 48109-2143
c (734) 763 6239
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c     input:      1) ev(nmax) : a vector of eigenvalues with nmax elements
c           2) nmax           : number of elements of vector ev
c           3) emax           : maximum energy in master equation
c           4) b        : average rotational constant, defined as:
c                 b = (integral of b(x)dx from 0 to 2*pi) / (2*pi)
c     output:     1) ev(nmax)
c           2) nmax     
c
c
      implicit double precision(a-h,o-z), integer(i-n)
      dimension ev(nn)

      i=nmax

c.......if the largest eigenvalue is larger than the maximum energy, 
c.......then searching for an eigenvalue, which is <= emax
11    e=ev(i)
      if(e.gt.emax) then
            i=i-1
            goto 11
      endif
      nmax=i

c.......looking for higher-lying eigenvalues with double degeneracy
111   ne1=nint(ev(i))
      ne2=nint(ev(i-1))
      if(ne1.gt.ne2) then
            i=i-1
            goto 111
      endif

c.......switching from a hindered rotor to a pitzer rotor
c.......pitzer rotor is defined as e = eo + b*(jp**2)
c.......two variables eo and jp will be computed as below
      tp=ev(i-1) - ev(i-2)
      jp=nint((tp/b - 1.0d0)/2.0d0)
      eo=ev(i-2)-b*jp*jp

1111  ev(i+1)=eo + b*(jp+2)*(jp+2)
      ev(i+2)=ev(i+1)
      if(ev(i+2).le.emax) then
            i=i+2
            jp=jp+1
            goto 1111
      endif
      nmax=i

      return
 
      end

cccccccccccccccccccc
cccccccccccccccccccc
c
c calvo is used to compute vmax and vmin 
c of the torsional potential energy function 
c
        subroutine calvo(vo,vm,cv,ncv,nsiv,vhr,pha)

        implicit double precision(a-h,o-z), integer(i-n)
        dimension cv(ncv)
        character vhr*5

        pi=acos(-1.0d0)

        rmin=0.0d0
        rmax=pi*2
        zl=(rmax-rmin)
        dx=zl/3600.0d0
        call vx(vo,0.0d0+pha,cv,ncv,nsiv,vhr,0.0d0)
        vm=vo
        do i=1, 3600
                xa = rmin + dx*i
                call vx(vtp,xa+pha,cv,ncv,nsiv,vhr,0.0d0)
                if(vtp.gt.vo) vo=vtp
                if(vtp.lt.vm) vm=vtp
        enddo
        return
        end

cccccccccccccccccccc
cccccccccccccccccccc









