!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
! Three band model (3D)
module global_variables
! mathematical parameters
  real(8),parameter :: pi = 4d0*atan(1d0)
  complex(8),parameter :: zI = (0d0,1d0)
  real(8),parameter :: a_B=0.529177d0,Ry=13.6058d0

! Material parameter
  integer :: nkx, nky, nkz, nk, nk_s, nk_e
  integer :: nk_ave, nk_remainder
  integer,allocatable :: ikx_table(:),iky_table(:),ikz_table(:)
  real(8) :: dkx, dky, dkz
  real(8),allocatable :: kx(:), ky(:), kz(:)
  real(8),parameter :: fact_intra = 1d0, fact_inter = 1d0
! d: (semi) core , v: valence, c: conduction
  real(8),parameter :: eps_d = -(0.1440380881432848d+001 - 0.2586486109145065d+00)
  real(8),parameter :: eps_c1 = 0d0
  real(8),parameter :: eps_c2 = eps_c1 + 0.2d0/(2d0*Ry) 
  real(8),parameter :: mass_c1 = 2.5d0, mass_c2 = mass_c1
  real(8),parameter :: piz_dc1 = 1d0,piz_dc2 = 0d0
  real(8),parameter :: piz_dcc = 0d0
  real(8),parameter :: lattice_constant

! Time-propagation
  integer :: Nt
  real(8) :: dt
  real(8),allocatable :: eps(:,:)
  complex(8),allocatable :: zCt(:,:,:)
  real(8),allocatable :: Act(:),Act_dt2(:),jtz(:),jtz_intra(:),jtz_inter(:)
  character(20) :: envelope_1,envelope_2

! Laser parameters
  real(8) :: E0_1,omega_1,tpulse_1,omega_ev_1,tpulse_fs_1,Iwcm2_1
  real(8) :: E0_2,omega_2,tpulse_2,omega_ev_2,tpulse_fs_2,Iwcm2_2
  real(8) :: Tdelay_fs,Tdelay

! MPI
  include 'mpif.h'
  integer :: Myrank,Nprocs,ierr

! I/O
  integer,parameter :: nf_current = 21
  character(99),parameter :: cf_current = "current.out"

  contains
    subroutine calc_single_particle_energy(eps_t, kx_t, ky_t, kz_t)
      implicit none
      real(8),intent(out) :: eps_t(3)
      real(8),intent(in):: kx_t, ky_t, kz_t

      eps_t = 0d0

      eps_t(1) = eps_d
      eps_t(2) = eps_c1 + 0.5d0/mass_c1*(1d0 &
          -cos(0.5d0*lattice_const*kx_t)**2 &
          *cos(0.5d0*lattice_const*ky_t)**2 &
          *cos(0.5d0*lattice_const*kz_t)**2)

      eps_t(3) = eps_c2 + 0.5d0/mass_c2*(1d0 &
          -cos(0.5d0*lattice_const*kx_t)**2 &
          *cos(0.5d0*lattice_const*ky_t)**2 &
          *cos(0.5d0*lattice_const*kz_t)**2)

    end subroutine calc_single_particle_energy
end module global_variables