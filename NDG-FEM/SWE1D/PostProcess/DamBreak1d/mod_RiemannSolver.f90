module mod_RiemannSolver
implicit none

contains
!-------------------------------------------------------------------------!
subroutine Output
implicit none
integer :: i, mCells
real :: chaLen, gate, timeout, xCoord
integer, parameter :: mx = 3000
real :: d(mx), u(mx)

common /soluti/ d, u
common /domain/ chaLen, gate, mCells, timeout

open(unit=1, file = 'result.out', status='unknown')
! write(1, *)'result of Riemann Solver: '
do i = 1,mCells
	xCoord = real(i)*chaLen/real(mCells)
	write(1,*)xCoord, d(i), u(i), u(i)*d(i)
enddo
close(1)
end subroutine
!-------------------------------------------------------------------------!
subroutine WetBed
implicit none

integer :: i, it, mCells, niter
real :: cha, chaLen, cl, cr, cs, d0, dl, dr, ds, dsam, fl, fr, &
	& fld, frd, gate, gravity, s, timeout, TOL, ul, ur, us, usam, xCoord
integer, parameter :: mx = 3000
real :: u(mx), d(mx)
common /soluti/ d, u
common /states/ cl, dl, ul, cr, dr, ur
common /starso/ cs, ds, us
common /accele/ gravity
common /tolera/ niter, TOL
common /domain/ chalen, gate, mCells, timeout

write(*,*)
write(*,*)'Exact Solution in Start Region'
write(*,*)'=============================='
write(*,*)

call Starte

d0 = ds
write(*,*)'   IT   ', '   DS   ', '   CHA'
write(*,*)
do it = 1, niter
	call GeoFun(fl, fld, ds, dl, cl)
	call GeoFun(fr, frd, ds, dr, cr)
	ds = ds - (fl + fr + ur-ul)/(fld + frd)
	cha = abs(ds - d0)/ (0.5*(ds+d0) )
	write(*,*)it, ds, cha
	if(cha .le. TOL) exit
	if(ds .lt. 0.0)ds = TOL
	d0 = ds
enddo

if (it > niter) then
	write(*,*)'Number of niter iterations exceeded, stop'
	stop
endif

30 format(i6, 2x, 2(f12.7, 2x))

us = 0.5*(ur+ul) +0.5*(fr- fl)

write(*,*)
write(*,*)'Depth in Star Region  	h* =', ds
write(*,*)'Velocity in Star Region  u* =', us
write(*,*)

cs = sqrt(gravity*ds)

do i =1,mCells
	xCoord = real(i)*chaLen/real(mCells) - gate
	s = xCoord/timeout

	call SamWet(dsam, usam, s)

	d(i) = dsam
	u(i) = usam
enddo
end subroutine WetBed
!-------------------------------------------------------------------------!
subroutine GeoFun(f, fd, d, dk, ck)
! Purpose: to evaluate functions FL, FR and their derivatives in
!			iterative Riemann solver for wet-bed case
implicit none

real :: c, ck, d, dk, f, fd, ges, gravity

common /accele/ gravity
if (d <= dk) then
 	! rarefaction wave
 	c = sqrt(gravity*d)
 	f = 2.0*(c - ck)
 	fd = gravity/c
else
 	! shock wave
 	ges = sqrt(0.5*gravity*(d+dk)/(d*dk))
 	f = (d-dk)*ges
 	fd = ges - 0.25*gravity*(d - dk)/(ges*d*d)
endif
end subroutine GeoFun
!-------------------------------------------------------------------------!
subroutine Starte
! Purpose: to provide starting value for Newton-Raphson iteration.
! 			the Two-Rarefaction Riemann Solver (TRRS) and Two-Shock 
! 			Riemann Solver (TSRS) are used adaptively
implicit none
real :: cl, cr, cs, dl, dmin, dr, ds, gel, ger, gravity, &
	& ul, ur, us
common /states/ cl, dl, ul, cr, dr, ur
common /starso/ cs, ds, us
common /accele/ gravity

dmin = min(dl, dr)

! use Two-Rarefaction (TRRS) solution as starting value
ds = (1.0/gravity)*(0.5*(cl+cr) - 0.25*(ur-ul))**2

if ( ds .le. dmin) then
	! use Tow-Rarefaction (TRRS) solution as starting value
	! with DS as computed from TRRS as estimate
	write(*,*)'TS approximation, h* = ', ds
else
	! use Two-Shock (TSRS) solution as starting value with 
	! DS as computed from TRRS as estimate
	write(*,*)'TS approximation, h* = ', ds

	gel = sqrt(0.5*gravity*(ds+dl)/(ds*dl))
	ger = sqrt(0.5*gravity*(ds+dl)/(ds*dl))
	ds = (gel*dl + ger*dr -(ur - ul))/(gel + ger)
endif
write(*,*)
end subroutine Starte
!-------------------------------------------------------------------------!
subroutine SamWet(d, u, s)
! Purpose: to sample solution through wave structure at TIMEOUT for 
!  			wet-bed case
implicit none
real :: c, cl, cr, cs, d, dl, dr, ds, gravity, ql, qr, s, &
		shl, shr, sl, sr, stl, str, u, ul, ur, us
common /states/ cl, dl, ul, cr, dr, ur
common /starso/ cs, ds, us
common /accele/ gravity

if (s .le. us) then
	! sample left wave
	if (ds .ge. dl) then
		! left shock
		ql = sqrt((ds+dl)*ds/(2.0*dl*dl))
		sl = ul - cl*ql

		if (s .le. sl) then
			d = dl
			u = ul
		else
			d = ds
			u = us 
		endif
	else
		! left rarefaction
		shl = ul - cl
		if (s .le. shl) then
			d = dl
			u = ul
		else 
			stl = us -cs
			if (s .le. stl) then
				u = (ul + 2.0*cl +2.0*s)/3.0
				c = (ul + 2.0*cl -s)/3
				d = c*c/gravity
			else
				d = ds
				u = us
			endif
		endif
	endif
else
	! sample right wave
	if (ds .ge. dr) then
		! right shock
		qr = sqrt((ds+dr)*ds/(2.0*dr*dr))
		sr = ur + cr*qr

		if (s .ge. sr) then
			! sample point line to the right of the shcok
			d = dr
			u = ur
		else
		 	! simple point line to the left of the shcok
		 	d = ds
		 	u = us
		endif
	else
		! right rarefaction
		shr = ur + cr
		if (s .ge. shr) then
			d = dr
			u = ur
		else
			str = us + cs
			if (s .ge. str)then
				u = (ur - 2.0*cr + 2.0*s)/3
				c = (-ur + 2.0*cr + s)/3
				d = c*c/gravity
			else
				d = ds
				u = us
			endif
		endif
	endif
endif
end subroutine SamWet
!-------------------------------------------------------------------------!
subroutine DryBed
! Purpose: to compute the exact solution in the case in which 
! 			a portion of dry bed is present
implicit none

integer :: i, mCells
real :: chaLen, cl, cr, dl, dr, dsam, gate, s, timeout, &
 		ul, ur, usam, xCoord
integer, parameter :: mx = 3000
real :: d(mx), u(mx)

common /soluti/ d, u
common /states/ cl, dl, ul, cr, dr, ur
common /domain/ chaLen, gate, mCells, timeout

do i = 1, mCells
	xCoord = real(i)*chaLen/real(mCells) - gate
	s = xCoord/timeout
	if (dl .le. 0.0) then
		call SamLef(dsam, usam, s)
	else
	 	if (dr .le. 0.0) then
	 		call SamRig(dsam, usam, s)
	 	else
	 		call SamMid(dsam, usam, s)
	 	endif
	endif
	d(i) = dsam
	u(i) = usam
enddo

end subroutine DryBed
!-------------------------------------------------------------------------!
subroutine SamLef(d, u, s)
! Purpose: to sample the solution through the wave structure at TimeOut,
! 			for the case in which the left state is dry. Solution consists
! 			of single right rarefaction
implicit none

real :: c, cl, cr, d, dl, dr, gravity, s, shr, str, u, ul, ur

common /states/ cl, dl, ul, cr, dr, ur
common /accele/ gravity

shr = ur + cr
if (s .ge. shr) then
	! sample point lies to the right of the rarefaction
	d = dr
	u = ur

else
	str = ur - 2.0*cr
	if(s .ge. str) then
		u = (ur - 2.0*cr + 2.0*s)/3.0
		c = (-ur + 2.0*cr +s)/3.0
		d = c*c/gravity
	else
		d = dl
		u = ul
	endif
endif
end subroutine SamLef
!-------------------------------------------------------------------------!
subroutine SamMid(d, u, s)
! Purpose: to sample the solution through the wave structure at time 
! 			TimeOut, for the case in which the middle state is dry.
! 			Solution consists of a left and a right rarefaction with
! 			a dry portion in the middle
implicit none
real :: c, cl, cr, d, dl, dr, gravity, s, shl, shr, ssl, ssr, &
		u, ul, ur
common /states/ cl, dl, ul, cr, dr, ur
common /accele/ gravity

shl = ul - cl
ssl = ul + 2.0*cl
ssr = ur - 2.0*cr
shr = ur + cr

if (s .le. shl) then
	d = dl
	u = ul
endif

if ((s .gt. shl).and.(s .le. ssl)) then
	u = (ul + 2.0*cl + 2.0*s)/3.0
	c = (ul + 2.0*cl - s)/3.0
	d = c*c/gravity
endif

if ((s .ge. shl).and. (s .le. ssr)) then
	d = 0.0
	u = 0.0
endif

if ((s .gt. ssr).and. (s .le. shr)) then
	u = (ur - 2.0*cr + 2.0*s)/3.0
	c = (-ur + 2.0*cr + s)/3.0
	d = c*c/gravity
endif

if(s .ge. shr) then
	d = dr
	u = ur
endif

end subroutine SamMid
!-------------------------------------------------------------------------!
subroutine SamRig(d, u, s)
implicit none
real :: c, cl, cr, d, dl, dr, gravity, s, shl, stl, u, ul, ur
common /states/ cl, dl, ul, cr, dr, ur
common /accele/ gravity

shl = ul - cl
if (s .le. shl) then
	d = dl
	u = ul
else
	stl = ul + 2.0*cl
	if (s .le. stl) then
		u = (ul + 2.0*cl + 2.0*s)/3.0
		c = (ul + 2.0*cl - s)/3.0
		d = c*c/gravity
	else
		d = dr
		u = ur
	endif
endif
end subroutine SamRig

end module mod_RiemannSolver

