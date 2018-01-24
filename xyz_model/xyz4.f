	program  xyz4
* 18/III/16 - Use gaussian random number of collisions rather than poisson in 'trvers'
**** 'trvers' has been changed..... !!!!!
* 21/XII/15
* Reproduce TC's multiple scattering code...

	common/hists/iha(402),ihb(402),ihc(402),ihd(402)

	write(6,*),'xyz4'
* Set up histograms
	call hist(-1)
* Set up elements
	call setel
* Set up materials
	call setmat
*
	dum = 999.
	call muscat(dum,-1)

	do 100 mu = 1,50000
	if(mod(mu,5000).eq.0)write(6,*),mu
	pmu = 240.0
	thx = 0.0
	thy = 0.0
*.............................SCATTERERS HERE ...........................
* 2 x 0.165 cm tracker planes
*	mat = 6
*	thick = 0.33
*	call trvers(pmu,mat,thick,thx,thy)
* Windows ~ 0.12 cm of Al
*	mat = 3
*	thick = 0.12
*	call trvers(pmu,mat,thick,thx,thy)	
* Absorber (or not) - Xe
*	thick = 35 ! thickness in cm
*	mat = 5 ! Xe
* Absorber - Li6-H
	thick = 6.5 ! thickness in cm
	mat = 12 ! Li6-H
	call trvers(pmu,mat,thick,thx,thy)
*........................................................

	call muscat(thx,0)
	thesq = thx**2 + thy**2
	themu = sqrt(thesq)
* Histogram theta
	i = 1 + themu/0.00025
	if(i.gt.400)i = 401
	if(i.lt.1)i = 402
	iha(i) = iha(i) + 1
* Histogram theta^2
	i = 1 + thesq/2.5e-5
	if(i.gt.400)i = 401
	if(i.lt.1)i = 402
	ihb(i) = ihb(i) + 1
* Theta-x
	i = 201+ thx/0.0005
	if(i.gt.400)i = 401
	if(i.lt.1)i = 402
	ihc(i) = ihc(i) + 1
* Theta-y
	i = 201+ thy/0.0005
	if(i.gt.400)i = 401
	if(i.lt.1)i = 402
	ihd(i) = ihd(i) + 1

100	continue
	call muscat(dum,1)
* write out histograms
	call hist(1)
	stop
	end
* -------------------------------------------------------------
	subroutine hist(iflag)

	common/hists/iha(402),ihb(402),ihc(402),ihd(402)

	if(iflag.eq.-1)then
	  do k = 1,402
	    iha(k) = 0
	    ihb(k) = 0
	    ihc(k) = 0
	    ihd(k) = 0
	  enddo
	  return
	elseif(iflag.eq.1)then
	  do k = 1,402
	     write(21,*) iha(k),ihb(k),ihc(k),ihd(k) 
	  enddo
	endif
	return
	end

* -------------------------------------------------------------

	subroutine trvers(pmu,mat,thick,thx,thy)
* Scatter this muon 
	common/element/izel(200),awtel(200)
	common/material/nmats,naspec(100),iza(100,20),mza(100,20),rho(100),amw(100),iprmat(100)

	real*4 pcn(20),pce(20),th1n(20),th2n(20),th1e(20),th2e(20)
	real*4 th1ns(20),th2ns(20),th1es(20),th2es(20)
	real*4 spcoll(40)
	real*8 thxd, thyd
	data ame,ammu,re,ana/0.510998,105.658,2.817940e-13,6.022142e23/

	data istart/0/

* calculate collision probabilities in this sample of material
	if(iprmat(mat).eq.0)write(6,*)'Material ',mat
	beta = pmu/sqrt(pmu**2+ammu**2)
	pi = 2.0*asin(1.0)
	const = 4.0*ana*pi*re*re*(ame/(beta*pmu))**2
*	write(6,*)const
	nsp = naspec(mat)
	pctot = 0.0
	do k = 1,nsp
	  iz = iza(mat,k)
	  m = mza(mat,k)
	  a = awtel(iz)
	  z = izel(iz)
	  th1n(k) = ame*z**(1.0/3.0)/(137*pmu)
	  th1e(k) = th1n(k) ! (don't know any better)
	  th2n(k) = 280*(ame/pmu)*a**(-1.0/3.0)
	  th2e(k) = 4.84e-3
	  if(iprmat(mat).eq.0)write(6,*)k,th1n(k),th1e(k),th2n(k),th2e(k)
	  th1ns(k) = th1n(k)**2
	  th1es(k) = th1e(k)**2
	  th2ns(k) = th2n(k)**2
	  th2es(k) = th2e(k)**2
* Nucleus
	  pcn(k) = const * z*z * m /amw(mat)
*	  pcn(k) = pcn(k)*(th2n(k)**2-th1n(k)**2)/(th1n(k)**2+th2n(k)**2)/(th1n(k)**2)
	  pcn(k) = pcn(k)*(th2ns(k)-th1ns(k))/(th1ns(k)+th2ns(k))&
                   /(th1ns(k))
	  pcn(k) = pcn(k)*rho(mat)*thick
* Electrons
	  pce(k) = const * z * m /amw(mat)
*	  pce(k) = pce(k)*(th2e(k)**2-th1e(k)**2)/(th1e(k)**2+th2e(k)**2)/(th1e(k)**2)
	  pce(k) = pce(k)*(th2es(k)-th1es(k))/(th1es(k)+th2es(k))/(th1es(k))
	  pce(k) = pce(k)*rho(mat)*thick
	  if(iprmat(mat).eq.0)write(6,*)'N ',mat,pcn(k)
	  if(iprmat(mat).eq.0)write(6,*)'e ',mat,pce(k)
	  pctot = pctot + pcn(k) 
	  spcoll(2*k-1) = pctot
	  pctot = pctot + pce(k)
	  spcoll(2*k) = pctot
	enddo
* Normalise probabilities
	do k=1,2*nsp
	  spcoll(k) = spcoll(k)/pctot
	  if(iprmat(mat).eq.0)write(6,*),k,spcoll(k)
	enddo
* Mean free path
	fpath = thick/pctot
	if(iprmat(mat).eq.0)then
	   write(6,*)'Material, Thickness', mat, thick
	   write(6,*)'pctot, fpath ',pctot, fpath
	   iprmat(mat) = 1
	endif
* Number of collisions fluctutated
	call gran(gfluc)
	ncoll = pctot + sqrt(pctot)*gfluc 
*----------------------------
	thxd = thx
	thyd = thy
	idum = 0
	nce = 0
	ncn = 0
* make ncoll collisions:
	do 100 mc = 1,ncoll
* Get type of collision
	  r = rand(idum)
	  do k = 1, 2*nsp
	    if(r.lt.spcoll(k))then
	      kt = k
	      go to 2
	    endif
	  enddo
2	  continue
* Species
	  ksp = (kt+1)/2
*	write(6,*),ksp,kt
* N or e?
	  if(mod(kt,2).eq.0)then
* electron
	    nce = nce + 1
	    r = rand(idum)
	    tsq = r*th1es(ksp)*th2es(ksp)/(th1es(ksp)+(1-r)*th2es(ksp))
	  else
* Nucleus
	    ncn = ncn + 1
	    r = rand(idum)
	    tsq = r*th1ns(ksp)*th2ns(ksp)/(th1ns(ksp)+(1-r)*th2ns(ksp))
          endif
	  th = sqrt(tsq)
	  phi = 2.0*pi*rand(idum)
	  thxd = thxd + th*cos(phi)
	  thyd = thyd + th*sin(phi)
*	write(6,*)tx,ty
100	continue
*	
	if(istart.lt.10)write(6,*)ncoll,nce,ncn
	istart = istart+1
	thx = thxd
	thy = thyd
	return
	end
* -------------------------------------------------------------
	subroutine setmat

	common/element/izel(200),awtel(200)
	common/material/nmats,naspec(100),iza(100,20),mza(100,20),rho(100),amw(100),
     1	iprmat(100)

* Up to 100 materials
* naspec(mat) = # atomic species in material mat
* iza(mat, i)  = z of atomic species i in molecule / mixture
* mza(mat, i)  = # of atomic species i in molecule
* rho(mat)    = density, gms/cm^3

* Set print flag
	do n = 1,100
	  iprmat(n) = 0
	enddo

* LH2 
	mat = 1
	naspec(mat) = 1
	iza(mat,1) = 1
	mza(mat,1) = 2
	rho(mat) = 0.0708
	rho(mat) = 0.0755 ! Muscat LH2 density
* He gas 
	mat = 2
	naspec(mat) = 1
	iza(mat,1) = 2
	mza(mat,1) = 1
	rho(mat) = 0.1786e-3 * (273.15/293.15)*(101.325/100.00) ! at NTP
* Aluminium
	mat = 3
	naspec(mat) = 1
	iza(mat,1) = 13
	mza(mat,1) = 1
	rho(mat) = 2.70
* LiH

	mat = 4
	naspec(mat) = 2
	iza(mat,1) = 3
	mza(mat,1) = 1
	iza(mat,2) = 1
	mza(mat,2) = 1
	rho(mat) = 0.82

* Xenon
	mat = 5
	naspec(mat) = 1
	iza(mat,1) = 54
	mza(mat,1) = 1
	rho(mat) = 5.858e-3 * (273.15/293.15)*(101.325/100.00) ! at NTP

* Polystyrene = (C6H5CHCH2)n == C8H8
	mat = 6
	naspec(mat) = 2
	iza(mat,1) = 6
	mza(mat,1) = 8
	iza(mat,2) = 1
	mza(mat,2) = 8
	rho(mat) = 1.060

* Lithium
	mat = 7
	naspec(mat) = 1
	iza(mat,1) = 3
	mza(mat,1) = 1
	rho(mat) = 0.53 ! Muscat Li density

* Berylium
	mat = 8
	naspec(mat) = 1
	iza(mat,1) = 4
	mza(mat,1) = 1
	rho(mat) = 1.85 ! Muscat density

* Polyethylene = (CH2)n
	mat = 9
	naspec(mat) = 2
	iza(mat,1) = 6
	mza(mat,1) = 1
	iza(mat,2) = 1
	mza(mat,2) = 2
	rho(mat) = 0.93 ! Muscat density

* Carbon
	mat = 10
	naspec(mat) = 1
	iza(mat,1) = 6
	mza(mat,1) = 1
	rho(mat) = 1.69 ! Muscat density

* Iron
	mat = 11
	naspec(mat) = 1
	iza(mat,1) = 26
	mza(mat,1) = 1
	rho(mat) = 7.86 ! Muscat density

* Li6H (i.e. LiH with Lithium-6)

	mat = 12
	naspec(mat) = 2
	iza(mat,1) = 103
	mza(mat,1) = 1
	iza(mat,2) = 1
	mza(mat,2) = 1
	rho(mat) = 0.692

* List everything
	nmats = 12
	write(6,*)nmats, ' Materials'
	do mat = 1, nmats
	  nsp = naspec(mat)
	  amw(mat) = 0.0
	  do k = 1,nsp
	    iz = iza(mat,k)
	    m = mza(mat,k)
	    a = awtel(iz)
	    z = izel(iz)
	    amw(mat) = amw(mat) + m*a
	    write(6,*),mat, m,z,a
	  enddo
	  write(6,*)'Molecular weight ', amw(mat)
	enddo
	write(6,*),'-------------'


	return
	end
* -------------------------------------------------------------
	subroutine setel
*
	common/element/izel(200),awtel(200)

	do j = 1,200
	   izel(j) = -999
	  awtel(j) = -999.
	enddo

* Hydrogen
	 izel(1) = 1
	awtel(1) = 1.00794
* d
	 izel(101) = 1
	awtel(101) = 2.014

* He
	 izel(2) = 2
	awtel(2) = 4.002602
* Li
	 izel(3) = 3
	awtel(3) = 6.941
* Li-6
	 izel(103) = 3
	awtel(103) = 6.0
* Be
	 izel(4) = 4
	awtel(4) = 9.0122
* C
	 izel(6) = 6
	awtel(6) = 12.011
* N
	 izel(7) = 7
	awtel(7) = 14.00674
* O
	 izel(8) = 8
	awtel(8) = 15.9994
* Al 
	 izel(13) = 13
	awtel(13) = 26.982
* Fe	
	 izel(26) = 26
	awtel(26) = 55.845
* Ar
	 izel(18) = 18
	awtel(18) = 39.948
* Xe
	 izel(54) = 54
	awtel(54) =  131.29
* Pb
	 izel(82) = 82
	awtel(82) = 207.2

	return
	end

*---------------------------------------------------------------------
	subroutine muscat(thx,iflag)

	save

	real*4 tbe(12)
	integer*4 nbin(11)

	data tbe/0.0,0.00269,0.00895,0.0162,0.0248,0.0347,0.0463,0.0597,
     1  0.0754,0.0938,0.1151,3.141/

	if(iflag.lt.0)then
	  do k=1,11
	    nbin(k) = 0
	  enddo
	return
	endif

	if(iflag.eq.0)then
  	  t = abs(thx)
	  if(t.gt.3)t=3.0
	  do k= 1,11
	     if(t.lt.tbe(k+1))then
	       nbin(k) = nbin(k) + 1
	       return
	     endif
	  enddo
	stop '????????????'
	endif


	if(iflag.eq.1)then
	
	  ntot = 0
	  write(6,*),iflag,ntot,(nbin(j),j=1,11)
	  do k = 1,11
	    ntot = ntot+nbin(k)
	  enddo
	  do k = 1,11
	    dt = tbe(k+1) - tbe(k)
	    p = nbin(k)/(dt*ntot)
	write(6,*),k,ntot,nbin(k),p
	  enddo
* Write out so can plot with paw g/h/e
	  do k = 1,10
	    dt = tbe(k+1) - tbe(k)
* Prob per radian and error
	    p = nbin(k)/(dt*ntot)/2.0
	    ep =sqrt(float(nbin(k)))/(dt*ntot)/2.0
* Bin centre & half bin width
            ab = (tbe(k+1)+tbe(k))/2.
	   eab = dt/2.
	write(7,*)ab,p,eab,ep
	  enddo
	endif

	return
	end

*-----------------------------------------------------------------------------
		subroutine gran(g)
* ~ Gaussian random number with mean 0, rms = 1.0
	do mu = 1,1000
	  g = 0.0	
	  do k = 1,6
	    idum = 0
	    g = g + rand(idum)
	  enddo
	  g = (g-3.)*sqrt(2.0) 
*	  write(12,*),g
	enddo
	return
	end
