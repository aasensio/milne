module milneMod
use vars
use atomic_functions
use math_functions, only : init_maths
use iso_c_binding, only: c_int, c_double
implicit none

contains

	subroutine c_setline(nLambda, lineData, waveOut) bind(c)
	integer(c_int), intent(in) :: nLambda
	real(c_double), intent(in) :: lineData(7)
	real(c_double), intent(out), dimension(nLambda) :: waveOut
	integer :: k
				
		line%wave0 = lineData(1)
		line%Jup = lineData(2)
		line%Jlow = lineData(3)
		line%gup = lineData(4)
		line%glow = lineData(5)
		line%lambdaInit = lineData(6)
		line%lambdaStep = lineData(7)
		line%nLambda = nLambda
		
		Stokes_Syn%nlambda = line%nLambda
						
		if (.not.associated(Stokes_Syn%lambda)) allocate(Stokes_Syn%lambda(line%nLambda))
		if (.not.associated(Stokes_Syn%stokes)) allocate(Stokes_Syn%stokes(4,line%nLambda))
		
		Stokes_Syn%stokes = 0.d0
		
		do k = 1, line%nLambda
			Stokes_Syn%lambda(k) = line%lambdaInit + line%lambdaStep * (k-1)						
		enddo
		
		waveOut = Stokes_Syn%lambda
		
		call init_maths
		
		if (allocated(zeeman_voigt)) deallocate(zeeman_voigt)
		if (allocated(zeeman_faraday)) deallocate(zeeman_faraday)
		
		if (allocated(ki)) deallocate(ki)
		if (allocated(kq)) deallocate(kq)
		if (allocated(ku)) deallocate(ku)
		if (allocated(kv)) deallocate(kv)
		if (allocated(fq)) deallocate(fq)
		if (allocated(fu)) deallocate(fu)
		if (allocated(fv)) deallocate(fv)
		
		if (allocated(delta)) deallocate(delta)
		if (allocated(stokes)) deallocate(stokes)
		
		if (allocated(profile)) deallocate(profile)
		if (allocated(v)) deallocate(v)
		
		allocate(profile(2,line%nLambda))
		allocate(v(line%nLambda))
				
		allocate(zeeman_voigt(3,line%nLambda))
		allocate(zeeman_faraday(3,line%nLambda))
		
		allocate(ki(line%nLambda))
		allocate(kq(line%nLambda))
		allocate(ku(line%nLambda))
		allocate(kv(line%nLambda))
		allocate(fq(line%nLambda))
		allocate(fu(line%nLambda))
		allocate(fv(line%nLambda))
		
		allocate(delta(line%nLambda))
		allocate(stokes(4,line%nLambda))

		
	end subroutine c_setline
	
	subroutine c_milnesynth(nLambda, modelIn, muIn, stokesOut) bind(c)
	integer(c_int), intent(in) :: nLambda
	real(c_double), intent(in) :: modelIn(9), muIn
	real(c_double), intent(out), dimension(4,nLambda) :: stokesOut
	
 		model%Bfield = modelIn(1)
 		model%theta = modelIn(2)
 		model%chi = modelIn(3)
 		model%vmac = modelIn(4)
 		model%damping = modelIn(5)
 		model%B0 = modelIn(6)
 		model%B1 = modelIn(7)
 		model%doppler = modelIn(8)
 		model%kl = modelIn(9)
 		
 		model%mu = muIn
 				
 		call synthesize(model,line,Stokes_Syn)
 		 		
 		StokesOut = Stokes_Syn%stokes
 		
	end subroutine c_milnesynth

	subroutine c_milnesynthmany(nLambda, nModels, modelIn, muIn, stokesOut) bind(c)
	integer(c_int), intent(in) :: nLambda, nModels
	real(c_double) :: muIn
	real(c_double), intent(in), dimension(nModels,9) :: modelIn
	real(c_double), intent(out), dimension(nModels,4,nLambda) :: stokesOut
	integer(c_int) :: i
	
		do i = 1, nModels
			model%Bfield = modelIn(i,1)
			model%theta = modelIn(i,2)
			model%chi = modelIn(i,3)
			model%vmac = modelIn(i,4)
			model%damping = modelIn(i,5)
			model%B0 = modelIn(i,6)
			model%B1 = modelIn(i,7)
			model%doppler = modelIn(i,8)
			model%kl = modelIn(i,9)
			model%mu = muIn
					
			call synthesize(model,line,Stokes_Syn)
					
			StokesOut(i,:,:) = Stokes_Syn%stokes
		enddo
 		
	end subroutine c_milnesynthmany

end module milneMod