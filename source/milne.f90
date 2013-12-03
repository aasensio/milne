module milneMod
implicit none
use vars

	subroutine setLine(nLines, lineData)
	integer :: nLines
	real(kind=8) :: lineData(nLines,8)
	
	!f2py integer, intent(in) :: nLines
	!f2py real(8), intent(in), dimension(nLines,8) :: lineData
	
		allocate(lines(nLines))
		
		do i = 1, nLines
			lines(i)%linewave0 = lineData(i,1)
			lines(i)%Jup = lineData(i,2)
			lines(i)%Jlow = lineData(i,3)
			lines(i)%gup = lineData(i,4)
			lines(i)%glow = lineData(i,5)
			lines(i)%lambdaInit = lineData(i,6)
			lines(i)%lambdaStep = lineData(i,7)
			lines(i)%nLambda = lineData(i,8)
		enddo
						
		
	end subroutine setLine


	subroutine milneInit(nLines)
	!f2py integer, intent(in) :: nLines
		
	end subroutine milneInit
	
end module milneMod