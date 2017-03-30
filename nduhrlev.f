c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c uhrlev: a code to calculate eigenvalues of unsymmetrical, non-rigid 1d-hindered internal rotation
c copyright (c) 2009 lam t. nguyen and john r. barker
c
c     date: feb. 13, 2009
c
c john r. barker
c jrbarker@umich.edu
c university of michigan
c ann arbor, mi 48109-2143
c (734) 763 6239
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine uhrlev(t,at,dele,jmax,ncb,ncv,cb,cv,n1,n2,imax,
     &          zpe,vhr,bhr,phav,phab,control,nsig)

      implicit double precision(a-h,o-z), integer(i-n)
      parameter (nn=2000)
      integer control
      dimension ev(nn), cv(ncv), cb(ncb), t(jmax), at(jmax)
       dimension ir(nn)
        character vhr*5, bhr*5

        pi=acos(-1.0d0)
        emax=(jmax-1)*dele + 1000.              ! *****raise max energy to account for later subtraction of zpe******

        if((bhr.eq.'bhrd2').or.(bhr.eq.'bhrd2')) then   ! average rotational constant
                b=0.0d0
                do i=1, n2
                        b=b + (cb(i)/i)*((pi*2)**(i-1))
                enddo
        elseif((bhr.eq.'ihrd2').or.(bhr.eq.'ihrd2')) then       ! average moment of inertia
                b=0.0d0
                do i=1, n2
                        b=b + (cb(i)/i)*((pi*2)**(i-1))
                enddo
                b=16.85763d0 / b        ! convert from i to b
        elseif((bhr.eq.'ihrd1').or.(bhr.eq.'ihrd1')) then       ! average moment of inertia
                b=16.85763d0/cb(1)      ! convert from i to b
        elseif((bhr.eq.'bhrd1').or.(bhr.eq.'bhrd1')) then       ! average rotational constant
                b=cb(1)
        else
                write(*,*) "error at input for bhr or ihr"
        endif

        call ghrlev(ev,nn,imax,emax,b,ncv,ncb,cv,cb,n1,n2,
     &          vhr,bhr,phav,phab)

      zpe=ev(1)

      if(control.eq.0)return

c*********************************************************************
c     the following lines are for output of hindered rotor eigenvalues
c*********************************************************************
c        write(12,*) '    i     e(i)      e(i)-zpe     [uhrlev]'

        do i=1, imax+1
c          write(12, 99) i, ev(i), (ev(i)-zpe)
        enddo
99    format (3x, i5, 1x, f10.3 , 1x, f10.3 )
c*********************************************************************

      do i=1, imax
            tp=ev(i)-zpe
            ir(i)= nint( tp/dele )
            ir(i)= nint( tp/dele )
c......           write(*,*) i, tp
      enddo 

      do j=1, imax
            ll=ir(j)
            do l=1, jmax - ll
                  at(l+ll) = at(l+ll) + t(l)
            enddo
      enddo
      do j=1, jmax
            t(j)=at(j)/nsig
            at(j)=0.0d0
      enddo 

      return
      end


