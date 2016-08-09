program main
use mod_RiemannSolver
implicit none
real :: chaLen, cl, cr, dcrit, dl, dr, gate, gravity, timeout, TOL, ul, ur
integer :: mCells, nIter
character(len = 100):: filename

common /states/ cl, dl, ul, cr, dr, ur
common /accele/ gravity
common /tolera/ nIter, TOL
common /domain/ chaLen, gate, mCells, timeout

call getarg(1, filename)

open(unit=1, file=trim(filename), status='unknown')

read(1,*)chaLen 	! length of channel
write(*,*)"chnnel length = ", chaLen
read(1,*)gate 		! position of gate
read(1,*)gravity 	! 
read(1,*)mCells 	! number of cells
read(1,*)TOL
read(1,*)nIter 		! iterations in exact solver
read(1,*)timeout 	!
read(1,*)dl
read(1,*)ul
read(1,*)dr
read(1,*)ur

close(1)

cl = sqrt(gravity*dl)
cr = sqrt(gravity*dr)

dcrit = (ur - ul) - 2.0*(cl*cr)

if ((dl < 0.0) .or. (dr < 0.0) .or. (dcrit>0.0) ) then
	call DryBed
else
	call WetBed
endif
call Output

end program