      implicit  none
      character aversion*6,adate*8,prog*6,cut*46,tstmp*10
      integer   i,j,k,l,m,n,iv,ib,etop,datpts,jtop
      parameter ( aversion='2017', adate='Feb 2017', prog='ktools' )
      parameter ( etop = 50000, datpts=100, jtop=1001   )
      parameter (cut='**************INPUT DATA SUMMARY**************')
      
c
c i/o units
c
      integer iunit,lunit,efunit,canunit,veffout,densunit,nminsunit
      integer junit,checkunit
      parameter (iunit=40,lunit=41,efunit=43,canunit=42,veffout=44)
      parameter (densunit=45,nminsunit=46,junit=47,checkunit=48)
c
c noise parameters
c
      real*8     canonh,microh
      parameter  ( canonh=2.0d0, microh=1.5d0 )

c
c densum parameters
c
      real*8 degrain1,demax2
      integer dimax1,disize
c
c found in read_input
c

      character eunits*4,sunits*3,whatdo*5,inputfile*80,backup*6
      character ATYPE*6,fileroot*80,names*80
      dimension names(10)
      integer   numnames
      dimension atype(datpts,10)
      real*8    de,temps(datpts),maxtemp
      integer   etopuser,jtopuser,dj,nt,binmax,rcnt,ntts,pcnt
      character reprod*4,molname*80,title*180
      dimension reprod(datpts),molname(datpts)
      real*8    delh(datpts),delhf(datpts),delhr(datpts)
      character formulal*49,title1*150,title2*150,title3*150
      dimension formulal(datpts),title1(datpts)
      dimension title2(datpts),title3(datpts)
      real*8    sym(datpts),sopt(datpts)
      integer   nele(datpts),nelement(datpts)
      real*8    elev(datpts,datpts)
      integer   gele(datpts,datpts)
      integer   ndof(datpts),Natom(datpts,10)
      character keyword*4
      dimension keyword(datpts,2)
      integer   mode(datpts,datpts)
      character idofl*3
      dimension idofl(datpts,datpts)
      real*8    wel(datpts,datpts),anhl(datpts,datpts)
      integer   ngl(datpts,datpts)
      real*8    vimag,vvr,distl(datpts)
      character*5  vhrl(datpts,datpts),bhrl(datpts,datpts)
      integer    nsvl(datpts,datpts)
      real*8     phavl(datpts,datpts),cvl(datpts,datpts,datpts)
      integer    nsbl(datpts,datpts)
      real*8     phabl(datpts,datpts),cbl(datpts,datpts,datpts)
      integer    ncvl(datpts,datpts),ncbl(datpts,datpts)

c
c super_mol
c
      real*8     amu(datpts),tamu(2),ramu(2)


c
c part of hrd
c
      real*8     cbb(datpts*10),cvv(datpts*10)
      real*8     ev(datpts*20)
      real*8     evhl(datpts,datpts,datpts*20),evh(datpts*10)
      integer    valmax(datpts,datpts)

c
c found in nrgconvert
c
      real*8     converge,h,c,kb,kbj,avo,kcaltoj,pi,amu2kg
      real*8     kcal2nu,amua2nu,gmcm2nu,ghz2nu,mhz2nu
      real*8     rgascal,rgasjol,rgascm1

      parameter ( converge = 1.0d-16                              )
      parameter ( h = 6.62606957d-34,c = 29979245800.0d0          )
      parameter ( kb = 0.69503476d0, kbj = 1.3806488d-23          )
      parameter ( avo = 6.02214129d+23, kcaltoj = 4.184d+03       )
      parameter ( pi = 2.0d0*dacos(0.0d0)                         )
      parameter ( amu2kg= (1.0d0/avo)*(1.0d0/1.0d+03)             )
      parameter ( kcal2nu=(1.0d0/avo)*kcaltoj*(1.0d0/h)*(1.0d0/c) )
      parameter ( amua2nu=((h*avo*1.0d+23)/(8.0d0*c*(pi**2)))  )
      parameter ( gmcm2nu=((h*1.0d+7)/(8.0d0*c*(pi**2)))          )
      parameter ( ghz2nu=((1.0d+09)/c)                            )
      parameter ( mhz2nu=((1.0d+06)/c)                            )
      parameter ( rgasjol = 8.3144621d+00                         )
      parameter ( rgascal = rgasjol/kcaltoj                       )
      parameter ( rgascm1 = rgasjol/(h*c*avo)                     )
c      parameter ( rgasjol = 8.314462d+00                         )
c      parameter ( rgascal = 8.314462D+00/4.184d+00               )
c      parameter ( rgascm1 = 0.69503475684138738                  )
c
c found in bsort
c
      real*8     barr(datpts),aarr(datpts)

c
c found in partition_funcs
c
      integer    minvec(datpts),starvec(datpts)
      real*8     ratekinf(4,datpts),interr(datpts)
      integer    divide(datpts)


c
c found in find_jmax
c
      integer    maxj

c
c found in jloop/veff
c
      real*8    vef(datpts),epsil(datpts)
      real*8    v1l(jtop+1),vlist(jtop+1),vefmax(jtop+1)
      real*8    vmax
      integer   vmaxi,epsmax

c
c found in sterab
c

      integer   nbinl(jtop+1)
      real*8    rpmat(jtop+1,etop),minpmat(jtop+1,etop+1)
      real*8    ppmat(jtop+1,etop),minppmat(jtop+1,etop+1)
      real*8    mingmat(jtop+1,etop+1),umingmat(jtop+1,etop+1)

c
c found in calc_rates
c

      real*8     kejmat(jtop+1,etop+1),ukejmat(jtop+1,etop+1)
      real*8     cratef(datpts),ucratef(datpts)
      real*8     crater(datpts),ucrater(datpts)
      real*8     mrate(datpts),umrate(datpts)
      real*8     queues(datpts)
      logical    manymins

      common  /reals/  temps,maxtemp,delh,sym,sopt,elev,delhf,
     $                 delhr,distl,wel,anhl,phavl,cvl,phabl,cbl,
     $                 evhl,evh,barr,aarr,ratekinf,interr,epsil,
     $                 vmax,v1l,vlist,vefmax,rpmat,minpmat,ppmat,
     $                 minppmat,mingmat,umingmat,kejmat,ukejmat,
     $                 degrain1,demax2,amu,tamu,ramu,cratef,ucratef,
     $                 mrate,umrate,queues,crater,ucrater,vef
      common  /ints/   de,dj,etopuser,jtopuser,nt,rcnt,ntts,pcnt,
     $                 nele,gele,ndof,ngl,nsvl,nsbl,ncvl,ncbl,
     $                 valmax,binmax,minvec,starvec,maxj,
     $                 vmaxi,epsmax,nbinl,divide,numnames,
     $                 dimax1,disize,Nelement,Natom
      common  /logs/   manymins
      common  /chars/  title,formulal,tstmp,title1,title2,title3,
     $                 eunits,sunits,inputfile,fileroot,whatdo,
     $                 backup,reprod,keyword,idofl,molname,vhrl,bhrl,
     $                 names,ATYPE

