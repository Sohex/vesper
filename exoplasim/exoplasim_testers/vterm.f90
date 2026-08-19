! Calculate terminal velocity of aerosol particles
! Sources:
! 
! Parmentier, V. et al. (2013) 3D mixing in hot Jupiter atmospheres. I: Application to the day/night cold trap in HD 209458b
! Steinrueck, M. et al. (2021) 3D simulations of photochemical hazes in the atmosphere of hot Jupiter HD 189733b
!

PROGRAM generate_data
	IMPLICIT NONE

	INTEGER,PARAMETER :: im=3
	INTEGER,PARAMETER :: jm=3
	INTEGER,PARAMETER :: nl=6
	INTEGER,PARAMETER :: ts=10
	INTEGER :: i,j,k,t


	REAL :: ta(im,jm,nl)
	REAL :: sigma(im,jm,nl)
	REAL :: ps(im,jm,nl)
	REAL :: beta(im,jm,nl)
	REAL :: rho(im,jm,nl)
	REAL :: mu(im,jm,nl)
	REAL :: vels(im,jm,nl)
	
	REAL :: apart
	REAL :: rhop

	DO k=1,nl
		DO i=1,im
			DO j=1,jm
				ta(i,j,k) = 300.-i-j-k
				ps(i,j,k) = 101325
				sigma(i,j,k) = 101325-i-j-10*k
				rho(i,j,k) = 1.255-0.1*k
			END DO
		END DO
	END DO
	
	apart = 50E-06
	rhop = 1000
	
	call viscos(ta,im,jm,nl,mu)
	call chamfac(ta,im,jm,nl,sigma,ps,apart,beta)
	call vterm(im,jm,nl,beta,apart,rho,rhop,mu,ta,sigma,ps,vels)
	
	WRITE(*,*) mu,vels
	

END PROGRAM

!****6***0*********0*********0*********0*********0*********0**********72
	SUBROUTINE viscos(ta,im,jm,nl,	 &
				mu)
	
! 	Calculate the viscosity for molecular nitrogen, N2
!****6***0*********0*********0*********0*********0*********0**********72
	
	IMPLICIT NONE
	
! Dimensions of arrays

	INTEGER,INTENT(IN) :: im,jm,nl ! Length of x, y, and z dimensions
	
! Define arrays

	REAL,INTENT(IN) :: ta(im,jm,nl) ! Temperature array
	REAL,INTENT(OUT) :: mu(im,jm,nl) ! Viscosity array
	
! Define constants

	REAL :: kb ! Boltzmann constant
	REAL :: m ! Molecular mass of N2
	REAL :: d ! Molecular diameter of N2
	REAL :: eps ! Lennard-Jones potential well of N2
	REAL :: PI ! Value of pi
	
	kb = 1.380649E-23 ! m2 kg s-2 K-1
	m = 4.652E-26 ! kg
	d = 3.64E-10 ! m
	eps = 95.5 ! K
	PI = 3.14159
	
!---------------------------------------------------------------
! Calculation
	
	mu = (5./16.)*(1./1.22)*(1./PI)*SQRT(PI*ta*m)*(SQRT(kb)/(d**2))*((ta/eps)**(0.16))
			
	RETURN
	END


!****6***0*********0*********0*********0*********0*********0**********72
	SUBROUTINE chamfac(ta,im,jm,nl,sigma,ps,apart,	 &
						beta)
	
! 	Calculate the Cunningham factor used in terminal velocity calculation
!****6***0*********0*********0*********0*********0*********0**********72

	IMPLICIT NONE

! Dimensions of arrays

	INTEGER,INTENT(IN) :: im,jm,nl ! Length of x, y, and z dimensions
	
! Define arrays

	REAL,INTENT(IN) :: ta(im,jm,nl) ! Temperature array
	REAL,INTENT(IN) :: sigma(im,jm,nl) ! Sigma array
	REAL,INTENT(IN) :: ps(im,jm,nl) ! Surface air pressure
	REAL,INTENT(OUT) :: beta(im,jm,nl) ! Cunningham factor array
	REAL :: lambda(im,jm,nl) ! Mean free path array
	REAL :: air_pres(im,jm,nl)
	REAL :: kn(im,jm,nl)
	
! Define constants

	REAL :: kb ! Boltzmann constant
	REAL :: rad ! Molecular radius of N2
	REAL :: PI ! Value of pi
	REAL :: apart ! Particle radius
	
	kb = 1.380649E-23 ! J/K
	rad = 0.5*3.64E-10 ! m
	PI = 3.14159


	!---------------------------------------------------------------
! Calculation of mean free path of N2

	air_pres = ps*sigma
	lambda = (kb/rad)*(ta/air_pres)*(1/(4.*PI*rad))
	
! Calculation of Cunningham factor

	kn = lambda/apart ! Knudsen number
	beta = 1 + kn*(1.256 + 0.4*EXP(-1.1/kn))
	
	
	RETURN
	END

	
!****6***0*********0*********0*********0*********0*********0**********72
	SUBROUTINE vterm(im,jm,nl,beta,apart,rho,rhop,mu,	 &
					ta,sigma,ps,vels)
	
! 	Calculate the terminal velocity of aerosol particles in the vertical direction
!****6***0*********0*********0*********0*********0*********0**********72

	IMPLICIT NONE

! Dimensions of arrays
	
	INTEGER,INTENT(IN) :: im,jm,nl ! Length of x, y, and z dimensions
	
! Define arrays

	REAL,INTENT(INOUT) :: beta(im,jm,nl) ! Cunningham factor
	REAL,INTENT(IN) :: rho(im,jm,nl) ! Air density
	REAL,INTENT(IN) :: mu(im,jm,nl) ! Viscosity
	REAL,INTENT(IN) :: ta(im,jm,nl) ! Temperature array
	REAL,INTENT(IN) :: sigma(im,jm,nl) ! Sigma array
	REAL,INTENT(IN) :: ps(im,jm,nl) ! Surface air pressure
	REAL,INTENT(OUT) :: vels(im,jm,nl) ! Terminal velocity

! Define constants

	REAL :: apart ! Aerosol particle radius
	REAL :: rhop ! Density of aerosol particle
	REAL :: g ! Gravity
	
	g = 10.9 ! m/s2, gravitational constant for ProxB
	
	call chamfac(ta,im,jm,nl,sigma,ps,apart,beta)
	CONTINUE
	
	vels = 2*beta*apart**2*g*(rhop - rho)/(9*mu)
	
	RETURN
	END

	
!****6***0*********0*********0*********0*********0*********0**********72
	SUBROUTINE gsettle(t,rho,im,jm,nl,q,vterm)
	
! 	Calculate mass mixing ratio of haze particles at each timestep
!****6***0*********0*********0*********0*********0*********0**********72

	IMPLICIT NONE

! Dimensions of arrays
	
	INTEGER,INTENT(IN) :: im,jm,nl ! Length of x, y, and z dimensions
	
! Define arrays
	
	REAL,INTENT(IN) :: rho(im,jm,nl)
	REAL,INTENT(INOUT) :: vterm(im,jm,nl)
	REAL,INTENT(INOUT) :: q(im,jm,nl)
	
	DO k=1,nl
		q(:,:,k) = rho(:,:,k+1)*q(:,:,k+1)*vterm(:,:,k+1) -
	
	
	
	RETURN
	END

	