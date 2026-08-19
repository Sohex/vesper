! Parametrization of viscosity
! Sources:
! Rosner, D. (2000) Transport Processes in Chemically Reacting Flow Systems
! Parmentier, V. et al. (2013) 3D mixing in hot Jupiter atmospheres. I: Application to the day/night cold trap in HD 209458b
! Steinrueck, M. et al. (2021) 3D simulations of photochemical hazes in the atmosphere of hot Jupiter HD 189733b
!

PROGRAM calc_vis
	IMPLICIT NONE

	INTEGER,PARAMETER :: im=3
	INTEGER,PARAMETER :: jm=3
	INTEGER,PARAMETER :: nl=3
	INTEGER :: i,j,k


	REAL :: ta(im,jm,nl)
	REAL :: mu(im,jm,nl)

	DO k=1,nl
		DO i=1,im
			DO j=1,jm
				ta(i,j,k) = 280.+i+j+k
			END DO
		END DO
	END DO
	
	WRITE(*,*) ta

	call viscos(ta,im,jm,nl,mu)
	
	WRITE(*,*) ta, mu
	

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
