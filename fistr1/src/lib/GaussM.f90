!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module provides data for gauss quadrature
module gauss_integration
  use hecmw
  implicit none
  real(kind=kreal) :: XG(3, 3) !< abscissa of gauss points
  real(kind=kreal) :: WGT(3, 3) !< weight of gauss points
  !****************************
  !* Gauss Integration Table **
  !****************************
  !** 1st ***
  data XG (1,1)/0.0/
  data WGT(1,1)/2.0/
  !** 2nd ***
  data XG(2,1),XG(2,2)/-0.577350269189626,0.577350269189626/
  data WGT(2,1),WGT(2,2)/1.0,1.0/
  !** 3rd ***
  data XG(3,1),XG(3,2),XG(3,3)/  &
    -0.7745966692,               &
    0.0,                         &
    0.7745966692/
  data WGT(3,1),WGT(3,2),WGT(3,3)/ &
    0.5555555555,                  &
    0.8888888888,                  &
    0.5555555555/
  ! end of this module

  CONTAINS
!********************************************************************************
!* Generates gaussian points and weights to be used in Gauss-Legendre quadrature.
!**********************************************
	SUBROUTINE  gauleg(n, x, w)
      INTEGER, INTENT(IN) :: n            ! # of Gauss Points
      REAL(kind=kreal), INTENT(OUT) :: x(n), w(n)

      REAL(kind=kreal), PARAMETER :: M_PI=3.141592654d0
      REAL(kind=kreal), PARAMETER :: EPS=3.0d-14

      INTEGER  i, j, m
      REAL(kind=kreal)  p1, p2, p3, pp, z, z1

	  m = (n + 1) / 2
	  do i = 1, m

     	z = cos( M_PI * (i-0.25d0) / (n+0.5d0) )
        do
     	    p1 = 1.0d0
        	p2 = 0.0d0
        	do j = 1, n
           	  p3 = p2
           	  p2 = p1
           	  p1 = ((2.0d0*j-1.0d0) * z * p2 - (j-1.0d0)*p3) / j
        	enddo
        	pp = n*(z*p1-p2)/(z*z-1.0d0)
        	z1 = z
        	z = z1 - p1/pp   

        	if(abs(z-z1) <= eps) exit
	    enddo

      	x(i) =  -z  
      	x(n+1-i) = z
      	w(i) = 2.0d0/((1.0d0-z*z)*pp*pp)
      	w(n+1-i) = w(i)

      end do

   End subroutine gauleg
end module gauss_integration
