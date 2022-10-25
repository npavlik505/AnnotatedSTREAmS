subroutine wrap_startmpi()
    call startmpi()
end subroutine

subroutine wrap_setup()
    call setup()
end subroutine

subroutine wrap_init_solver()
    call init_solver()
end subroutine

subroutine wrap_step_solver()
    call step_solver()
end subroutine

subroutine wrap_finalize_solver()
    call finalize_solver()
end subroutine

subroutine wrap_finalize()
    use mod_streams

    call finalize()
    call mpi_finalize(iermpi)
end subroutine

subroutine wrap_copy_gpu_to_cpu()
    call copy_gpu_to_cpu()
end subroutine

subroutine wrap_copy_cpu_to_gpu()
    call copy_cpu_to_gpu()
end subroutine

subroutine wrap_tauw_calculate() 
 use mod_streams
 implicit none
!
 integer :: j,m
 real(mykind), dimension(nvmean,ny/2) :: w1dh
 real(mykind), dimension(ny/2) :: ufav,vfav,wfav
 real(mykind), dimension(ny/2) :: uvd,ut,yt,dft,ft,gt
 real(mykind), dimension(3)    :: cc
 real(mykind) :: rhow,rmuw,rnuw,deltav
 real(mykind) :: d1,d2
 real(mykind) :: dudyw,tauw,utau,retau
 real(mykind) :: rr,rn
 real(mykind) :: yy,uu2,vv2,ww2,uv,rho2,tt2
 logical, dimension(nvmean) :: symm
!
 cc(1) =  1.83333333333333_mykind
 cc(2) = -1.16666666666667_mykind
 cc(3) =  0.33333333333333_mykind
!
 if (masterproc) then
!
  symm     = .true.
  symm(3 ) = .false.
  symm(14) = .false.
  symm(19) = .false.
!
  do j=1,ny/2
   do m=1,nvmean
   if (symm(m)) then
    w1dh(m,j) = 0.5_mykind*(w_av_1d(m,j)+w_av_1d(m,ny-j+1))
   else
    w1dh(m,j) = 0.5_mykind*(w_av_1d(m,j)-w_av_1d(m,ny-j+1))
   endif
  enddo
 enddo
!
! Wall density and viscosity
  rhow = cc(1)*w1dh( 1,1)+cc(2)*w1dh( 1,2)+cc(3)*w1dh( 1,3)
  rmuw = cc(1)*w1dh(20,1)+cc(2)*w1dh(20,2)+cc(3)*w1dh(20,3)
  rnuw = rmuw/rhow
!
  do j=1,ny/2
   ufav(j) = w1dh(13,j)/w1dh(1,j)
   vfav(j) = w1dh(14,j)/w1dh(1,j)
   wfav(j) = w1dh(15,j)/w1dh(1,j)
  enddo
  d1 = y(1)+1._mykind
  d2 = y(2)+1._mykind
  dudyw = (ufav(1)*d2**2-ufav(2)*d1**2)/(d1*d2*(d2-d1))
  tauw  = rmuw*dudyw
 endif
end subroutine
