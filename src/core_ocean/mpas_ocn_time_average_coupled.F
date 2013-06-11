! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_average_coupled
!
!> \brief MPAS ocean coupled time averager
!> \author Doug Jacobsen
!> \date   06/08/2013
!> \details
!>  This module contains the routines for time averaging
!>  coupling fields for the ocean core.
!
!-----------------------------------------------------------------------

module ocn_time_average_coupled

    use mpas_grid_types
    use ocn_constants

    implicit none
    save
    public

    contains 

!***********************************************************************
!
!  routine ocn_time_average_coupled_init
!
!> \brief   Coupled time averager initialization
!> \author  Doug Jacobsen
!> \date    06/08/2013
!> \details 
!>  This routine initializes the coupled time averaging fields
!
!-----------------------------------------------------------------------
    subroutine ocn_time_average_coupled_init(forcing)!{{{
        type (forcing_type), intent(inout) :: forcing

#ifdef MPAS_CESM
        real (kind=RKIND), dimension(:), pointer :: temperatureState, salinityState, zonalVelocityState, meridionalVelocityState, zonalSSHGradientState, meridionalSSHGradientState

        temperatureState => forcing % temperatureState % array
        salinityState => forcing % salinityState % array
        zonalVelocityState => forcing % zonalVelocityState % array
        meridionalVelocityState => forcing % meridionalVelocityState % array

        temperatureState(:) = 0.0_RKIND
        salinityState(:) = 0.0_RKIND
        zonalVelocityState(:) = 0.0_RKIND
        meridionalVelocityState(:) = 0.0_RKIND
        forcing % nAccumulatedCoupled % scalar = 0
#endif

    end subroutine ocn_time_average_coupled_init!}}}

!***********************************************************************
!
!  routine ocn_time_average_coupled_accumulate
!
!> \brief   Coupled time averager accumulation
!> \author  Doug Jacobsen
!> \date    06/08/2013
!> \details 
!>  This routine accumulated the coupled time averaging fields
!
!-----------------------------------------------------------------------
    subroutine ocn_time_average_coupled_accumulate(state, forcing)!{{{
        type (state_type), intent(in) :: state
        type (forcing_type), intent(in) :: forcing

#ifdef MPAS_CESM
        real (kind=RKIND), dimension(:,:), pointer :: zonalVelocity, meridionalVelocity
        real (kind=RKIND), dimension(:,:,:), pointer :: tracers
        integer :: index_temperature, index_salinity

        real (kind=RKIND), dimension(:), pointer :: temperatureState, salinityState, zonalVelocityState, meridionalVelocityState, zonalSSHGradientState, meridionalSSHGradientState


        index_temperature = state % index_temperature
        index_salinity = state % index_salinity
        zonalVelocity => state % normalVelocityZonal % array
        meridionalVelocity => state % normalVelocityMeridional % array
        tracers => state % tracers % array

        temperatureState => forcing % temperatureState % array
        salinityState => forcing % salinityState % array
        zonalVelocityState => forcing % zonalVelocityState % array
        meridionalVelocityState => forcing % meridionalVelocityState % array

        write(0,*) 'Tk = ', T0_Kelvin
        temperatureState(:) = temperatureState(:) + tracers(index_temperature, 1, :) + T0_Kelvin
        salinityState(:) = salinityState(:) + tracers(index_salinity, 1, :)
        zonalVelocityState(:) = zonalVelocityState(:) + zonalVelocity(1, :)
        meridionalVelocityState(:) = meridionalVelocityState(:) + meridionalVelocity(1, :)
        forcing % nAccumulatedCoupled % scalar = forcing % nAccumulatedCoupled % scalar + 1
#endif

    end subroutine ocn_time_average_coupled_accumulate!}}}

!***********************************************************************
!
!  routine ocn_time_average_coupled_normalize
!
!> \brief   Coupled time averager normalization
!> \author  Doug Jacobsen
!> \date    06/08/2013
!> \details 
!>  This routine normalizes the coupled time averaging fields
!
!-----------------------------------------------------------------------
    subroutine ocn_time_average_coupled_normalize(forcing)!{{{

        type (forcing_type), intent(inout) :: forcing

#ifdef MPAS_CESM
        real (kind=RKIND), dimension(:), pointer :: temperatureState, salinityState, zonalVelocityState, meridionalVelocityState, zonalSSHGradientState, meridionalSSHGradientState

        temperatureState => forcing % temperatureState % array
        salinityState => forcing % salinityState % array
        zonalVelocityState => forcing % zonalVelocityState % array
        meridionalVelocityState => forcing % meridionalVelocityState % array

        if(forcing % nAccumulatedCoupled % scalar > 0) then
           temperatureState(:) = temperatureState(:) / forcing % nAccumulatedCoupled % scalar
           salinityState(:) = salinityState(:) / forcing % nAccumulatedCoupled % scalar
           zonalVelocityState(:) = zonalVelocityState(:) / forcing % nAccumulatedCoupled % scalar
           meridionalVelocityState(:) = meridionalVelocityState(:) / forcing % nAccumulatedCoupled % scalar
           forcing % nAccumulatedCoupled % scalar = 0
        end if
#endif

    end subroutine ocn_time_average_coupled_normalize!}}}

end module ocn_time_average_coupled
