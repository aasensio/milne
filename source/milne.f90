module milneMod
use vars
use atomic_functions
use math_functions, only : init_maths
implicit none

contains

	subroutine setLine(lineData)
	real(kind=8) :: lineData(8)
	integer :: k
	
!f2py integer, intent(in) :: nLines
!f2py real(8), intent(in), dimension(nLines,8) :: lineData
			
		line%wave0 = lineData(1)
		line%Jup = lineData(2)
		line%Jlow = lineData(3)
		line%gup = lineData(4)
		line%glow = lineData(5)
		line%lambdaInit = lineData(6)
		line%lambdaStep = lineData(7)
		line%nLambda = lineData(8)
		
		Stokes_Syn%nlambda = line%nLambda
						
		if (.not.associated(Stokes_Syn%lambda)) allocate(Stokes_Syn%lambda(line%nLambda))
		if (.not.associated(Stokes_Syn%stokes)) allocate(Stokes_Syn%stokes(4,line%nLambda))
		
		Stokes_Syn%stokes = 0.d0
		
		do k = 1, line%nLambda
			Stokes_Syn%lambda(k) = line%lambdaInit + line%lambdaStep * (k-1)			
		enddo
		
		call init_maths
		
	end subroutine setLine
	
	subroutine milneSynth(modelIn,muIn,waveOut, stokesOut,n)
	integer :: n
	real(kind=8) :: modelIn(8), muIn, waveOut(n), stokesOut(4,n)
	
!f2py integer, optional, intent(in) :: n
!f2py real(8), intent(in), dimension(8) :: modelIn
!f2py real(8), intent(in) :: muIn
!f2py real(8), intent(out), dimension(4,n) :: stokesOut
!f2py real(8), intent(out), dimension(n) :: waveOut
	
 		model%Bfield = modelIn(1)
 		model%theta = modelIn(2)
 		model%chi = modelIn(3)
 		model%vmac = modelIn(4)
 		model%damping = modelIn(5)
 		model%beta = modelIn(6) 		
 		model%doppler = modelIn(7)
 		model%kl = modelIn(8)
 		model%mu = muIn
 				
 		call synthesize(model,line,Stokes_Syn)
 		 		
 		waveOut = Stokes_Syn%lambda
 		StokesOut = Stokes_Syn%stokes
 		
	end subroutine milneSynth
	
end module milneMod