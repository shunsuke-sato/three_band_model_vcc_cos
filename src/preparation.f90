!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine preparation
  use global_variables
  implicit none
  integer :: ikx,iky,ikz,ik


  nk = nkx*nky*nkz
  nk_ave=nk/Nprocs; nk_remainder=mod(nk,Nprocs)
  if (nk_remainder == 0) then
    nk_s=nk_ave*Myrank+1
    nk_e=nk_ave*(Myrank+1)
  else
    if (Myrank +1 <= nk_remainder) then
      nk_s=(nk_ave+1)*Myrank+1
      nk_e=nk_s + (nk_ave+1)-1 !(NKrz_ave+1)*(Myrank+1)
    else
      nk_s=(nk_ave+1)*nk_remainder + nk_ave*(Myrank-nk_remainder)+1
      nk_e=nk_s + nk_ave-1
    end if
  end if

  if(myrank==0)then
    write(*,"(A,2x,I7)")"nk=",nk
    write(*,"(A,2x,I7)")"nk_ave",nk_ave
    write(*,"(A,2x,I7)")"nk_remainder",nk_remainder
  end if
!  write(*,"(9I9)")myrank,NKrz_s,NKrz_e

  allocate(zCt(3,1,nk_s:nk_e),eps(3,nk_s:nk_e))
  allocate(kx0(nk),kx(nk))
  allocate(ky0(nk),ky(nk))
  allocate(kz0(nk),kz(nk))

  allocate(ikx_table(nk))
  allocate(iky_table(nk))
  allocate(ikz_table(nk))
  zCt = 0d0; zCt(1,1,:) = 1d0
  eps = 0d0

! table
  ik = 0
  do ikx = 0,nkx-1
    do iky = 0,nky-1
      do ikz = 0,nkz-1

        ik = ik + 1
        ikx_table(ik) = ikx
        iky_table(ik) = iky
        ikz_table(ik) = ikz

      end do
    end do
  end do
  

  dkx  = ((2d0*pi)/lattice_constant)/nkx
  dky  = ((2d0*pi)/lattice_constant)/nky
  dkz  = ((2d0*pi)/lattice_constant)/nkz

  
  ik = 0
  do ikx = 0,nkx-1
    do iky = 0,nky-1
      do ikz = 0,nkz-1

        ik = ik + 1
        kx0(ik) = ikx_table(ik)*dkx
        ky0(ik) = iky_table(ik)*dky
        kz0(ik) = ikz_table(ik)*dkz

      end do
    end do
  end do

  kx = kx0
  ky = ky0
  kz = kz0
  
  do ik = nk_s, nk_e
    call calc_single_particle_energy(eps(:,ik), kx(ik), ky(ik), kz(ik))
  end do

end subroutine preparation
