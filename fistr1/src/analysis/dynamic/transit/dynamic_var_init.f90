!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module provides functions to initialize variables
!> when initial velocity or acceleration boundary conditions are given.
!> attention : just for rigid motion in the initial state.

module m_dynamic_init_variables

  use m_fstr
  use m_dynamic_mat_ass_load

contains

  subroutine dynamic_init_varibles( hecMESH, hecMAT, fstrSOLID, fstrEIG, fstrDYNAMIC, fstrPARAM )

    implicit none

    type(hecmwST_local_mesh) :: hecMESH
    type(hecmwST_matrix)     :: hecMAT
    type(fstr_eigen)         :: fstrEIG
    type(fstr_solid)         :: fstrSOLID
    type(fstr_dynamic)       :: fstrDYNAMIC
    type(fstr_param)         :: fstrPARAM

    integer(kind=kint) :: j

    call dynamic_mat_ass_load (hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC, fstrPARAM)

    if( fstrSOLID%VELOCITY_type == kbcInitial ) then
      do j = 1, hecMESH%n_node*hecMESH%n_dof
        fstrDYNAMIC%ACC(j,1)=(hecMAT%B(j)-fstrDYNAMIC%ray_m*fstrEIG%mass(j)*fstrDYNAMIC%VEL(j,1))/&
          fstrEIG%mass(j)
      enddo
    elseif( fstrSOLID%ACCELERATION_type == kbcInitial ) then
      do j = 1, hecMESH%n_node*hecMESH%n_dof
        fstrDYNAMIC%VEL(j,1)=(hecMAT%B(j)-fstrEIG%mass(j)*fstrDYNAMIC%ACC(j,1))/&
          (fstrDYNAMIC%ray_m*fstrEIG%mass(j))
      enddo
    endif

  end subroutine dynamic_init_varibles

end module m_dynamic_init_variables
