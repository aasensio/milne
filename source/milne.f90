module milneMod
implicit none
use vars

	subroutine setLine(lineData)
	real(kind=8) :: lineData(nLines,8)
	integer :: k
	
	!f2py integer, intent(in) :: nLines
	!f2py real(8), intent(in), dimension(nLines,8) :: lineData
			
		line%linewave0 = lineData(i,1)
		line%Jup = lineData(i,2)
		line%Jlow = lineData(i,3)
		line%gup = lineData(i,4)
		line%glow = lineData(i,5)
		line%lambdaInit = lineData(i,6)
		line%lambdaStep = lineData(i,7)
		line%nLambda = lineData(i,8)
		
		Stokes_Syn%nlambda = line%nLambda
						
		if (.not.associated(Stokes_Syn%lambda)) allocate(Stokes_Syn%lambda(line%nLambda))
		if (.not.associated(Stokes_Syn%stokes)) allocate(Stokes_Syn%stokes(4,line%nLambda))
		
		Stokes_Syn%stokes = 0.d0
		
		do k = 1, line%nLambda
			Stokes_Syn%lambda(k) = line%lambdaInit + line%lambdaStep * (k-1)			
		enddo
		
	end subroutine setLine


	subroutine milneInit(nLines)
	!f2py integer, intent(in) :: nLines
		
	end subroutine milneInit
	
end module milneMod