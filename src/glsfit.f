c      $Id$
      subroutine glsfit(npts,mode,chisqr,varfit,xmean,ymean,sum,nz,
     +     wmax,lw,ddmch,
     +     buf,npmsav,ksav,
     +     resfile2)

	implicit real*8 (a-h,o-z)
	include 'dim.h'
	include 'acom.h'
	include 'bcom.h'
	include 'vcom.h'
	include 'orbit.h'
	include 'tz.h'

c This is a version of the original tempo fit() function, rewritten
c to use a generalized least squares (GLS) approach, which can use
c a non-diagonal 'weighting' matrx.  This uses the LAPACK routines
c for various things.  As much as possible this is consistent with
c the original fit() in terms of input/outputs.
c PBD 2014/04

c       moved declaration of real*8 array(NPA,NPA) to acom.h, djn, 8 Sep 98
	real*8 xmean(NPA),fctn(NPAP1)
	real*8 a(NPAP1),atmp(NPAP1),sigmaa(NPAP1),gcor(NPA)
	logical lw
	real*8 ddmch(*)
        real*8 buf(*)
        integer npmsav(*), ksav(*)
	character mark*1,adn*14,date*9,damoyr*9
        character*80 resfile2

        integer isetflag

        integer resn(20)  ! small buffer to hold one toa's worth of information
	real*8 resr(9)
        equivalence (resn(2),resr(1)) ! real*8 will not be on 8-byte boundaries

c Just define some huge matrices here for now.  Once it's working, try
c replacing using the "malloc" routines?
        real*8 Adm(NPTSDEF,NPAP1) ! weighted fit design matrix
        real*8 VTsvd(NPAP1,NPAP1) ! SVD "V-transpose" result
        real*8 sv(NPAP1) ! The singular values
        real*8 r(NPTSDEF) ! The weighted prefit resids
        real*8 cts(NPTSDEF) ! copy of the times
        real*8 work(10*NPAP1*NPAP1)
        integer lwork
        integer iwork(10*NPAP1)
        integer inforv

c packed cov matrix, to be malloced
        integer ncovpts  
        integer*8 dcovoff, idx
        real*8 dcov(1)
        real*8 jit_phs

	integer fd
	integer nwrt
	integer flags, filemode  
 	integer open, close, write

        lwork = 10*NPAP1*NPAP1
	mprt=10**nprnt
	sigma=0.
	chisq=0.

c nparam is total number of fit params (including mean)
c Zero out various matrices
	nterms=nparam-1
	do 29 j=1,nterms
          xmean(j) = xmean(j)/sum
          do 28 k=1,nterms
            array(j,k)=0.
 28       continue
 29     continue

        do 32 j=1,NPAP1
          a(j)=0.
          atmp(j)=0.
          sigmaa(j)=0.
          sv(j)=0.
          do 30 k=1,NPAP1
            VTsvd(k,j)=0.
 30       continue
          do 31 k=1,NPTSDEF
            Adm(k,j)=0.
 31       continue
 32     continue 

        write (*,'(''  glsfit Ntoa='',i6)') npts

c Alloc space for COV matrix.  To access array element i
c use dcov(i+dcovoff) and pass to calls as dcov(1+dcovoff).
c Packed upper triangular storage means element (i,j) is accessed
c using index i+j*(j-1)/2, with i<=j.
        print *,'  ... allocate matrices'
        ncovpts = npts*(npts-1)/2 + npts
        call mallocxd(dcov,ncovpts,8,dcovoff)
        do i=1,ncovpts 
          dcov(dcovoff+i)=0.0
        enddo

c This is the part where we read in the data
c vmemr call reads info for TOA number i from memory
c Inputs:
c   y is the (pre-fit) residual, in pulse phase (turns)
c   fctn are the fit basis funcs (param derivs) evaluated for this TOA
c   weight is the TOA's weight (units phase^-2)
c   terr is the TOA's uncertainty (us, including efac/equad)
c Computed here:
c   Adm, design matrix
c   r, pre-fit resids
c   dcov, cov matrix (diagonal part only)
        print *,'  ... fill design matrix'
	do 67 i=1,npts
          call vmemr(i,fctn,ct,y,weight,dn,terr,frq,fmjd,rfrq,nterms,
     +     buf,npmsav,ksav)
          if(ldesign) write(37) ct, weight, (fctn(j)-xmean(j), j=1,nterms)
          cts(i) = ct
          r(i) = y
          Adm(i,1) = 1d0 ! Constant phase term
          do 66 j=2,nparam
            Adm(i,j) = (fctn(j-1)-xmean(j-1))
 66       continue
          dcov(dcovoff+i+i*(i-1)/2) = 1d0/weight
 67     continue

c Add in extra cov matrix terms
c Test "jitter"
        jit_phs = 0.01
        print *,'  ... compute covariance'
        do 51 i=1,npts
          do 50 j=1,i
            idx=dcovoff+j+i*(i-1)/2
            if (abs(cts(i)-cts(j)).lt.1d0/1440.) then
              dcov(idx) = dcov(idx) + jit_phs**2
            endif
 50       continue
 51     continue

c Notes for GLS:
c Cholesky factorize: dpotrf() or dpptrf() for packed
        print *,'  ... Cholesky decomp'
        call dpptrf('U',npts,dcov(1+dcovoff),inforv)
c invert triangular matrix: dtrtri() or dtptri() for packed
        print *,'  ... matrix inverse'
        call dtptri('U','N',npts,dcov(1+dcovoff),inforv)

	if(ldesign) rewind(37)
	if(ldesign) write(37) npts, nterms

c Multiply by inv cov matrix TODO check transpose
        print *,'  ... matrix multiply'
        call dtpmv('U','T','N',npts,dcov(1+dcovoff),r,1)
        do 68 i=1,nparam
          call dtpmv('U','T','N',npts,dcov(1+dcovoff),Adm(1,i),1)
 68     continue

c Call SVD routine.  On output, Adm will be replaced by "U".
        print *,'  ... SVD'
        call dgesdd('O',npts,nparam,Adm,NPTSDEF,sv,
     +    Adm,NPTSDEF,VTsvd,NPAP1,work,lwork,iwork,inforv)

c TODO could test for low SV's here
        write(*,'(''      log10(cond number) = '',f5.2)')
     +    log10(sv(1)) - log10(sv(nparam))

c Fill in "array" param cov matrix values
        do 72 j=1,nterms
          do 71 k=1,nterms
            array(j,k) = 0.0
            do 70 i=1,nparam
              array(j,k)=array(j,k)+VTsvd(i,j+1)*VTsvd(i,k+1)/sv(i)/sv(i)
 70         continue
 71       continue
 72     continue


c Need to keep "gcor" results for output?
c what is the reference for this formula??
	do 80 i=1,nterms
c          gcor(i)=sqrt(abs(1.d0 - zzz/array(i,i)))
          gcor(i)=0.0
 80     continue

	if(lw) then
	  write(31,1170)
1170	  format(//'Normalized covariance matrix in "number of 9''s"',
     +      ' format:'/)
	  call mxprt(array,gcor,nterms,mfit,nbin,eclcoord)
	  unit=1.d-6
	  if(p0firs.gt.0.1) unit=1.d-3
	  nunit=dlog10(unit)-0.0001
	  write(31,1050) nunit
1050	  format(//'   N     TOBS         Date',
     +     '      FREQ    WGT      NPULSE     R1(P)  R2(',i2,') PH/DT'/)
	endif

c Compute post-fit resids in "whitened" basis
c r_post = r_pre - U U^t r_pre
        call dgemv('T',npts,nparam,1d0,Adm,NPTSDEF,r,1,0d0,atmp,1)
        call dgemv('N',npts,nparam,-1d0,Adm,NPTSDEF,atmp,1,1d0,r,1)
        chisq = ddot(npts,r,1,r,1)
        print *,"GLS chisq=", chisq
        chisq = 0.0

c Computes the post-fit param values in a
	do i=1,nparam
          atmp(i) = atmp(i)/sv(i)
        enddo
        call dgemv('T',nparam,nparam,1d0,VTsvd,NPAP1,atmp,1,0d0,a,1)

c Note a(1) is the const phase term, aa0 is var name from orig fit.f
        aa0 = a(1)

        flags = isetflag()
	filemode  = 6*64 + 6*8 + 2  ! octal 662 = rw-rw-r--
	if (lw) fd = open(resfile2,flags,filemode)
	resn(1) = 72
	resn(20) = 72

c re-read data, compute post-fit residuals in dt2, output 
c them to resid2.tmp and tempo.lis
        fnpts = npts-nz
        wmean = sum/fnpts
	do 108 i=1,npts
          call vmemr(i,fctn,ct,y,weight,dn,terr,frq,fmjd,rfrq,nterms,
     +     buf,npmsav,ksav)
          dt2=y-aa0
          weight=weight/wmean
          do 107 j=1,nterms
            dt2=dt2-a(j+1)*(fctn(j)-xmean(j)) ! j+1 due to a(1)==aa0
 107      continue
          nct=ct
          mark=char(mod(nct,26)+65)
          date=damoyr(int(fmjd))
          dt2sec=dt2*p0firs
          phase=0.
          if(a1(1).ne.0.) then
            if (nbin.ne.9) then
               phase=dmod((ct-t0(1))*86400.d0/pb(1)+1000000.d0,1.d0)
            else
               phase=dmod((ct-t0asc)*86400.d0/pb(1)+1000000.d0,1.d0)
            endif
          endif

          if(nprnt.eq.0.or.mod(i,mprt).eq.1) then
            write(adn,1106) dn
 1106       format(f14.0)

            if(lw) write(31,1107) i,fmjd,mark,date,
     +           rfrq,weight,adn(1:13),y,dt2sec/unit,phase
 1107       format(i4,f12.5,1x,a1,1x,a9,f9.3,f6.3,a13,
     +           f8.3,f8.2,f7.3)
          endif

C Correct tz ref TOA
          if(i.eq.ntzref) then
            ftzrmjd=ftzrmjd-dt2sec/8.64d4
            if (ftzrmjd.lt.0.) then
              ftzrmjd = ftzrmjd + 1
              ntzrmjd = ntzrmjd - 1
            elseif (ftzrmjd.ge.1.) then
              ftzrmjd = ftzrmjd - 1
              ntzrmjd = ntzrmjd + 1
            endif
          endif

	  resr(1) = ct
	  resr(2) = dt2
	  resr(3) = dt2sec
	  resr(4) = phase
	  resr(5) = frq
	  resr(6) = weight
	  resr(7) = terr
	  resr(8) = y
	  resr(9) = ddmch(i)
          if (lw) nwrt = write(fd,resn,80)
          wmax=max(wmax,weight)
          chisq=chisq+weight*dt2**2
 108    continue
	if (lw) nwrt = close(fd)

	freen=fnpts-nterms-1
	chisqr=chisq*wmean/freen
	if(mode.eq.0) varnce=chisqr
	if(mode.ge.1) varnce=1./wmean
	varfit=chisq/fnpts
	if(mode.eq.0) chisqr=0.

	do 133 j=1,nterms
          sigmaa(j+1)=dsqrt(array(j,j))
 133    continue

c TODO fix sigma0 calculation? does it matter?
	sigma0=varnce/fnpts
	do 145 j=1,nterms
          do 144 k=1,nterms
c            sigma0=sigma0+xmean(j)*xmean(k)*array(j,k) 
            sigma0=sigma0
 144      continue
 145    continue

	sigmaa(1)=dsqrt(sigma0)

	do 146 j=1,NPAP1
          freq(j)=0.
          ferr(j)=0.
 146      continue

	do 150 j=1,nparam
          freq(mfit(j))=a(j)
          ferr(mfit(j))=sigmaa(j)
 150    continue

	return
	end
