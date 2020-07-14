program main
    implicit none
    real(8), parameter :: Ra = 5d9, Pr = 4.4d0
    integer :: i,j,k,e
    integer,parameter :: lx = 8, ly = 8, lz = 8, le = 59
    real(8), dimension(1:lx,1:ly,1:lz,1:le) :: x_in, y_in, z_in
    real(8), dimension(1:lx,1:ly,1:lz,1:le) :: u_in, v_in, w_in, th_in
    real(8), dimension(1:lx,1:ly,1:lz,1:le) :: uu_in, vv_in, ww_in, tt_in
    real(8), dimension(1:lx,1:ly,1:lz,1:le) :: ut_in, vt_in, wt_in, uw_in
    real(8), dimension(1:lx,1:ly,1:lz,1:le) :: utt_in, vtt_in, wtt_in
    real(8), dimension(1:lx,1:ly,1:lz,1:le) :: tx2_in, ty2_in, tz2_in
    real(8) :: x,y,z
    real(8) :: t_mean
    real(8) :: t_var
    real(8) :: kappa_t,kappa_f
    real(8) :: mc,tc,prod,diss

    x_in = 0.d0; y_in = 0.d0; z_in = 0.d0
    u_in = 0.d0; v_in = 0.d0; w_in = 0.d0; th_in = 0.d0
    uu_in = 0.d0; vv_in = 0.d0; ww_in = 0.d0; tt_in = 0.d0
    ut_in = 0.d0; vt_in = 0.d0; wt_in = 0.d0; uw_in = 0.d0
    utt_in = 0.d0; vtt_in = 0.d0; wtt_in = 0.d0
    open(unit=10, file="c00_x_y_z.bin", status='unknown', form='unformatted')
    read(10) x_in,y_in,z_in
    close(unit = 10)

    open(unit=11, file="c01_u_v_w_t.bin", status='unknown', form='unformatted')
    read(11) u_in,v_in,w_in,th_in
    close(unit = 11)

    open(12, file = "c02_uu_vv_ww_tt.bin",status = "unknown",form = "unformatted")
    read(12) uu_in,vv_in,ww_in,tt_in
    close(12)

    open(13, file = "c03_ut_vt_wt_uw.bin",status = "unknown",form = "unformatted")
    read(13) ut_in,vt_in,wt_in,uw_in
    close(13)

    open(14, file = "c04_utt_vtt_wtt.bin",status = "unknown",form = "unformatted")
    read(14) utt_in,vtt_in,wtt_in
    close(14)
        
    open(15, file = "c14_tx2_ty2_tz2.bin",status = "unknown",form = "unformatted")
    read(15) tx2_in,ty2_in,tz2_in
    close(15)
    
!====calculate temperature preofile====
    z = 0.d0
    t_mean = 0.d0
    open(20, file = "mean_temperature.txt", status = 'unknown')
    write(20,*)' z ',' t_mean'
    do e = 1, le, 1
    do j = 1, lz-1, 1
    do k = ly-1, 1, -1
    do i = lx-1, 1, -1
    t_mean = t_mean + th_in(i,j,k,e)
    z  = z  + z_in(i,j,k,e)
    enddo
    enddo
    z = z/((lx-1)*(ly-1))
    t_mean = t_mean/((lx-1)*(ly-1))
    write(20,'(2es14.6)') z, t_mean
    z = 0.d0
    t_mean = 0.d0
    enddo
    enddo
    close(20)
!====calculate variance temperature profile====
    z = 0.d0
    t_var = 0.d0
    open(20, file = "variance_temperature.txt", status = 'unknown')
    write(20,*)' z ',' t_var'
    do e = 1, le, 1
    do j = 1, lz-1, 1
    do k = ly-1, 1, -1
    do i = lx-1, 1, -1
    t_var = t_var + dabs(tt_in(i,j,k,e)-th_in(i,j,k,e)*th_in(i,j,k,e))
    z  = z  + z_in(i,j,k,e)
    enddo
    enddo
    z = z/((lx-1)*(ly-1))
    t_var = t_var/((lx-1)*(ly-1))
    write(20,'(2es14.6)') z, t_var
    z = 0.d0
    t_var = 0.d0
    enddo
    enddo
    close(20)
!====calculate kappa_t====
    z = 0.0d0
    kappa_t = 0.d0
    open(20, file = "kappa_t.txt", status = 'unknown')
    write(20,*)' z ',' kappa_t'
    do e = 1, le, 1
    do j = 1, lz-1, 1
    do k = ly, 1, -1
    do i = lx, 1, -1
    kappa_t = kappa_t - (wt_in(i,j,k,e)-w_in(i,j,k,e)*th_in(i,j,k,e))/((th_in(i,j+1,k,e)-th_in(i,j,k,e))/(z_in(i,j+1,k,e)-z_in(i,j,k,e)))
    z  = z  + z_in(i,j,k,e)
    enddo
    enddo
    z = z/(lx*ly)
    kappa_t = kappa_t/(lx*ly)*dsqrt(Ra*Pr)
    write(20,'(2es14.6)') z, kappa_t
    z = 0.d0
    kappa_t = 0.d0
    enddo
    enddo
    close(20)
!====calculate kappa_f====
    z = 0.0d0
    kappa_f = 0.d0
    open(20, file = "kappa_f.txt", status = 'unknown')
    write(20,*)' z ',' kappa_f'
    do e = 1, le, 1
    do j = 1, lz-1, 1
    do k = ly, 1, -1
    do i = lx, 1, -1
    kappa_f = kappa_f - (wtt_in(i,j,k,e)-2*wt_in(i,j,k,e)*th_in(i,j,k,e)-w_in(i,j,k,e)*tt_in(i,j,k,e)+2*w_in(i,j,k,e)*th_in(i,j,k,e)*th_in(i,j,k,e)) &
              /((tt_in(i,j+1,k,e)-th_in(i,j+1,k,e)*th_in(i,j+1,k,e))-(tt_in(i,j,k,e)-th_in(i,j,k,e)*th_in(i,j,k,e))/(z_in(i,j+1,k,e)-z_in(i,j,k,e)))
    z  = z  + z_in(i,j,k,e)
    enddo
    enddo
    z = z/(lx*ly)
    kappa_f = kappa_f/(lx*ly)*dsqrt(Ra*Pr)
    write(20,'(2es14.6)') z, kappa_f
    z = 0.d0
    kappa_f = 0.d0
    enddo
    enddo
    close(20)
!====calculate mean convection term====
    z = 0.0d0
    mc = 0.d0
    open(20, file = "mean_convection.txt", status = 'unknown')
    write(20,*)' z ',' mc'
    do e = 1, le, 1
    do j = 1, lz-1, 1
    do k = ly, 1, -1
    do i = lx-1, 1, -1
    mc = mc + u_in(i,j,k,e)*(((tt_in(i+1,j,k,e)-th_in(i+1,j,k,e)**2)-(tt_in(i,j,k,e)-th_in(i,j,k,e)**2))/(x_in(i+1,j,k,e)-x_in(i,j,k,e))) &
            + v_in(i,j,k,e)*(((tt_in(i,j+1,k,e-th_in(i,j+1,k,e)**2))-(tt_in(i,j,k,e)-th_in(i,j,k,e)**2))/(z_in(i,j+1,k,e)-z_in(i,j,k,e)))
    z  = z  + z_in(i,j,k,e)
    enddo
    enddo
    z = z/((lx-1)*ly)
    mc = mc/((lx-1)*ly)
    write(20,'(2es14.6)') z, mc
    z = 0.d0
    mc = 0.d0
    enddo
    enddo
    close(20)
!====calculate turbulent convection term====
    z = 0.d0
    tc = 0.d0
    open(20, file = "turbulent_convection.txt", status = 'unknown')
    write(20,*)' z ',' tc'
    do e = 1, le, 1
    do j = 1, lz-1, 1
    do k = ly, 1, -1
    do i = lx, 1, -1
    tc = tc + ((wtt_in(i,j+1,k,e)-2*wt_in(i,j+1,k,e)*th_in(i,j+1,k,e)-w_in(i,j+1,k,e)*th_in(i,j+1,k,e)**2+2*w_in(i,j+1,k,e)*tt_in(i,j+1,k,e))&
            -(wtt_in(i,j,k,e)-2*wt_in(i,j,k,e)*th_in(i,j,k,e)-w_in(i,j,k,e)*th_in(i,j,k,e)**2+2*w_in(i,j,k,e)*tt_in(i,j,k,e)))/(z_in(i,j+1,k,e)-z_in(i,j,k,e))
    z  = z  + z_in(i,j,k,e)
    enddo
    enddo
    z = z/(lx*ly)
    tc = tc/(lx*ly)
    write(20,'(2es14.6)') z, tc
    z = 0.d0
    tc = 0.d0
    enddo
    enddo
    close(20)
!====calculate production term====
    z = 0.0d0
    prod = 0.d0
    open(20, file = "production.txt", status = 'unknown')
    write(20,*)' z ',' mc'
    do e = 1, le, 1
    do j = 1, lz-1, 1
    do k = ly, 1, -1
    do i = lx, 1, -1
    prod = prod + abs(2*(wt_in(i,j,k,e)-w_in(i,j,k,e)*th_in(i,j,k,e))*((th_in(i,j+1,k,e)-th_in(i,j,k,e)))/(z_in(i,j+1,k,e)-z_in(i,j,k,e)))
    z  = z  + z_in(i,j,k,e)
    enddo
    enddo
    z = z/(lx*ly)
    prod = prod/(lx*ly)
    write(20,'(2es14.6)') z, prod
    z = 0.d0
    prod = 0.d0
    enddo
    enddo
    close(20)
!====calculate diffusion term====
!diff

!====calculate thermal dissipation term==== 
    z = 0.0d0
    diss = 0.d0
    open(20, file = "thermal_dissipation.txt", status = 'unknown')
    write(20,*)' z ',' diss'
    do e = 1, le, 1
    do j = 1, lz-1, 1
    do k = ly, 1, -1
    do i = lx-1, 1, -1
    diss = diss + dabs(tx2_in(i,j,k,e)-((th_in(i+1,j,k,e)-th_in(i,j,k,e))/(x_in(i+1,j,k,e)-x_in(i,j,k,e)))**2 &
                + tz2_in(i,j,k,e)-((th_in(i,j+1,k,e)-th_in(i,j,k,e))/(z_in(i,j+1,k,e)-z_in(i,j,k,e)))**2)
    z  = z  + z_in(i,j,k,e)
    enddo
    enddo
    z = z/((lx-1)*ly)
    diss = diss/((lx-1)*ly)*(2/dsqrt(Ra*Pr))
    write(20,'(2es14.6)') z, diss
    z = 0.d0
    diss = 0.d0
    enddo
    enddo
    close(20)

    end program main
