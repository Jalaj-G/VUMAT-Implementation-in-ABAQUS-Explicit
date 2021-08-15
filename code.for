subroutine vumat(
C Read only (unmodifiable)variables -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     5  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
C Write only (modifiable) variables -
     7  stressNew, stateNew, enerInternNew, enerInelasNew )
C
      include 'vaba_param.inc'
C
      dimension props(nprops), density(nblock), coordMp(nblock,*),
     1  charLength(nblock), strainInc(nblock,ndir+nshr),
     2  relSpinInc(nblock,nshr), tempOld(nblock),
     3  stretchOld(nblock,ndir+nshr),
     4  defgradOld(nblock,ndir+nshr+nshr),
     5  fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6  stateOld(nblock,nstatev), enerInternOld(nblock),
     7  enerInelasOld(nblock), tempNew(nblock),
     8  stretchNew(nblock,ndir+nshr),
     8  defgradNew(nblock,ndir+nshr+nshr),
     9  fieldNew(nblock,nfieldv),
     1  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     2  enerInternNew(nblock), enerInelasNew(nblock),
C
      character*80 cmname
C
	parameter ( zero = 0.d0, one = 1.d0, two = 2.d0,
	* third = 1.d0 / 3.d0, half = 0.5d0, op5 = 1.5d0)
C
C 	For plane strain, axisymmetric, and 3D cases using
C 	the J2 Mises Plasticity with piecewise-linear isotropic hardening.
C
C 	The state variable is stored as:
C
C 	STATE(*,1) = equivalent plastic strain
C
C 	User needs to input
C 	props(1) Young’s modulus
C 	props(2) Poisson’s ratio
C 	props(3) syield0
C	props(4) hard

	e = props(1)
	xnu = props(2)
	hard = props(3)
	syield0 = props(4)

	twomu = e / ( one + xnu )
	alamda = xnu * twomu / ( one - two * xnu )
	thremu = op5 * twomu
C	nvalue = nprops/2-1
C
	if ( stepTime .eq. zero ) then
	do k = 1, nblock
C	
	trace = strainInc(k,1) + strainInc(k,2) + strainInc(k,3)
	stressNew(k,1) = stressOld(k,1)
	* + twomu * strainInc(k,1) + alamda * trace
	stressNew(k,2) = stressOld(k,2)
	* + twomu * strainInc(k,2) + alamda * trace
	stressNew(k,3) = stressOld(k,3)
	* + twomu * strainInc(k,3) + alamda * trace
	stressNew(k,4)=stressOld(k,4) + twomu * strainInc(k,4)
C
	if ( nshr .gt. 1 ) then
	stressNew(k,5)=stressOld(k,5) + twomu * strainInc(k,5)
	stressNew(k,6)=stressOld(k,6) + twomu * strainInc(k,6)
	end if
	end do
	else
C
	do k = 1, nblock
	peeqOld=stateOld(k,1)
	trace = strainInc(k,1) + strainInc(k,2) + strainInc(k,3)
	s11 = stressOld(k,1) + twomu * strainInc(k,1) + alamda * trace
	s22 = stressOld(k,2) + twomu * strainInc(k,2) + alamda * trace
	s33 = stressOld(k,3) + twomu * strainInc(k,3) + alamda * trace
	s12 = stressOld(k,4) + twomu * strainInc(k,4)
	if ( nshr .gt. 1 ) then
	s13 = stressOld(k,5) + twomu * strainInc(k,5)
	s23 = stressOld(k,6) + twomu * strainInc(k,6)
	end if
C
	smean = third * ( s11 + s22 + s33 )
	s11 = s11 - smean
	s22 = s22 - smean
	s33 = s33 - smean
	if ( nshr .eq. 1 ) then
	vmises = sqrt( op5*(s11*s11+s22*s22+s33*s33+two*s12*s12) )
	else
	vmises = sqrt( op5 * ( s11 * s11 + s22 * s22 + s33 * s33 +
	* two * s12 * s12 + two * s13 * s13 + two * s23 * s23 ) )
	end if
C
	sigdif = vmises - yieldOld
	facyld = zero
	if ( sigdif .gt. zero ) facyld = one
	deqps = facyld * sigdif / ( thremu + hard )
C
	syield =syield0 + hard*deqps	
C
C	 Update the stress
C
	yieldNew = yieldOld + hard * deqps
	factor = yieldNew / ( yieldNew + thremu * deqps )
	stressNew(k,1) = s11 * factor + smean
	stressNew(k,2) = s22 * factor + smean
	stressNew(k,3) = s33 * factor + smean
	stressNew(k,4) = s12 * factor
	if ( nshr .gt. 1 ) then
	stressNew(k,5) = s13 * factor
	stressNew(k,6) = s23 * factor
	end if
C
C	 Update the state variables
C
	stateNew(k,1) = stateOld(k,1) + deqps
C
C Update the specific internal energy -
C
 if ( nshr .eq. 1 ) then
	stressPower = half * (
	* ( stressOld(k,1) + stressNew(k,1) ) * strainInc(k,1) +
	* ( stressOld(k,2) + stressNew(k,2) ) * strainInc(k,2) +
	* ( stressOld(k,3) + stressNew(k,3) ) * strainInc(k,3) ) +
	* ( stressOld(k,4) + stressNew(k,4) ) * strainInc(k,4)
	else
	stressPower = half * (
	* ( stressOld(k,1) + stressNew(k,1) ) * strainInc(k,1) +
	* ( stressOld(k,2) + stressNew(k,2) ) * strainInc(k,2) +
	* ( stressOld(k,3) + stressNew(k,3) ) * strainInc(k,3) ) +
	* ( stressOld(k,4) + stressNew(k,4) ) * strainInc(k,4) +
	* ( stressOld(k,5) + stressNew(k,5) ) * strainInc(k,5) +
	* ( stressOld(k,6) + stressNew(k,6) ) * strainInc(k,6)
	end if
	enerInternNew(k) = enerInternOld(k) + stressPower / density(k)
	
C
C	 Update the dissipated inelastic specific energy -
C
	plasticWorkInc = half * ( yieldOld + yieldNew ) * deqps
	enerInelasNew(k) = enerInelasOld(k)
	* + plasticWorkInc / density(k)
	end do
	end if
C
	return
	end