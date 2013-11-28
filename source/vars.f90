!********************************************************
! Global variables
!********************************************************
module vars
implicit none

	real(kind=8), parameter :: PK = 1.3806503d-16, UMA = 1.66053873d-24, PC = 2.99792458d10
	real(kind=8), parameter :: PH = 6.62606876d-27, PHK = PH / PK, PHC = 2.d0 * PH / PC**2
	real(kind=8), parameter :: PME = 9.10938188d-28, PE = 4.8032d-10, PI = 3.14159265359d0
	real(kind=8), parameter :: PHC2 = 2.d0 * PH * PC**2, OPA = PI * PE**2 / (PME * PC)
	real(kind=8), parameter :: SQRTPI = 1.77245385091d0, EARTH_SUN_DISTANCE = 1.495979d13
	real(kind=8), parameter :: LARMOR = PE/(4.d0*PI*PME*PC), PMP = 1.67262158d-24
	real(kind=8), parameter :: NUCLEAR_MAGNETON = PE*PH/(2.d0*PI*PMP)
	real(kind=8), parameter :: BOHR_MAGNETON = PE*PH/(2.d0*PI*PME)
	real(kind=8), parameter :: PH2C3 = 2.d0 * PH**2 * PC**3, PERTURBATION = 0.005d0	
	
	real(kind=8), allocatable :: zeeman_voigt(:,:), zeeman_faraday(:,:)
	
	integer :: n_lineas
	character(len=80) :: fich_malla, fich_lineas, file_output_profile, file_output_model
	character(len=80) :: file_stray_light
	
	type stokes_type
		integer :: nlambda
		real(kind=8), pointer :: stokes(:,:), lambda(:), sigma(:,:)
	end type stokes_type
	
	type modelo_type
		logical :: stray_light_component
		integer :: stray_light_nlambda
		real(kind=8) :: Bfield, theta, chi, ff, vmac, damping, beta, macrot, mu
		real(kind=8) :: doppler, filling_factor
		real(kind=8), dimension(2) :: Bfield_range, theta_range, chi_range, ff_range
		real(kind=8), dimension(2) :: vmac_range, damping_range, beta_range, macrot_range
		real(kind=8), dimension(2) :: doppler_range, filling_factor_range
		real(kind=8), pointer :: isotopic_abun(:), kl_range(:,:)
		real(kind=8), pointer :: straylight(:,:)
		real(kind=8) :: kl(10)
	end type modelo_type
	
	type inversion_type
		integer :: n_cycles, n_params, nparams_invert
		integer, pointer, dimension(:) :: which_to_invert, filling_factor_position
		real(kind=8), pointer, dimension(:,:) :: range
		real(kind=8), pointer, dimension(:,:,:) :: dydx
		real(kind=8) :: lambda, chisq
		real(kind=8) :: stokes_weights(4)
	end type inversion_type
	
	type line_type
		character(len=2) :: theory
		real(kind=8) :: lambda_init, lambda_end, lambda_step, lambda
		real(kind=8) :: wave0, Jup, Jlow, gup, glow
		real(kind=8) :: Lu, Su, Au, Bu, Ll, Sl, Al, Bl, I
	end type line_type
	
	type(stokes_type) :: Observation, Emergent, Stokes_unperturbed, Stokes_perturbed
	type(stokes_type), pointer :: Stokes_Syn(:)
	
	type(modelo_type), pointer :: model(:)
	
	type(inversion_type) :: inversion
	
	type(line_type), pointer :: linea(:)
	
	integer :: what_to_do, number_of_components, verbose
	
	character(len=10), allocatable :: parameters_name(:)
end module vars
