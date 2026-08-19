! Cunningham factor
! Sources:
! 
! Parmentier, V. et al. (2013) 3D mixing in hot Jupiter atmospheres. I: Application to the day/night cold trap in HD 209458b
! Steinrueck, M. et al. (2021) 3D simulations of photochemical hazes in the atmosphere of hot Jupiter HD 189733b
!

PROGRAM generate_data
	IMPLICIT NONE

	INTEGER,PARAMETER :: im=3
	INTEGER,PARAMETER :: jm=3
	INTEGER,PARAMETER :: nl=3
	INTEGER :: i,j,k


	REAL :: ta(im,jm,nl)
	REAL :: sigma(im,jm,nl)
	REAL :: ps(im,jm,nl)
	REAL :: beta(im,jm,nl)
	
	REAL :: apart

	DO k=1,nl
		DO i=1,im
			DO j=1,jm
				ta(i,j,k) = 280.+i+j+k
				ps(i,j,k) = 101325
				sigma(i,j,k) = 101325-i-j-k
			END DO
		END DO
	END DO
	
	apart = 50E-06
	
	WRITE(*,*) ta

	call chamfac(ta,im,jm,nl,sigma,ps,apart,beta)
	
	WRITE(*,*) beta
	

END PROGRAM

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
	REAL,INTENT(OUT) :: beta(im,jm,nl) ! Viscosity array
	REAL :: lambda(im,jm,nl) ! Mean free path array
	REAL :: air_pres(im,jm,nl)
	REAL :: kn(im,jm,nl)
	
! Define constants

	REAL :: kb ! Boltzmann constant
	REAL :: rad ! Molecular radius of N2
	REAL :: PI ! Value of pi
	REAL :: apart
	
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
	