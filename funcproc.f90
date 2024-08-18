module funcproc

contains
  
  subroutine param(x_max, x_min, y_max, y_min, num, r_fit, al, rough)
    
    real, intent(in) :: x_max, x_min, y_max, y_min, rough
    real :: r_fit, al
    integer, intent(in) :: num
    
    r_fit = 3*sqrt((x_max - x_min)**2 + (y_max - y_min)**2)/sqrt(real(num)-1)
    
    al = (num-1)*100/((x_max - x_min)*(y_max - y_min))
    al = exp(0.4247*rough*log(al)-1.3501/rough)
    
  end subroutine param
  !==============================================================================

  subroutine param_esf(r_med, num, r_fit, al, rough)
    
    real, parameter :: pi = 3.141592654
    real, intent(in) :: r_med, rough
    real :: r_fit, al
    integer :: num
    
    r_fit = 6*r_med*pi/sqrt(real(num-1))
    
    al = (num-1)*100/(4*pi*r_med*r_med)
    al = exp(0.4984*rough*log(al) - 1.06016110229/rough)
    
  end subroutine param_esf
  !==============================================================================
  
  subroutine abre(arq, file_number, ext, back)
    
    logical :: ex, back
    integer :: ierr, file_number
    character(len=10) :: arq
    character(len=3) :: ext

    inquire(file=trim(arq)//'.'//ext, exist=ex)
    
    if (ex) then
       
       call execute_command_line("mv "//trim(arq)//"."//ext//" "//trim(arq)//"_p."//ext)
       back = .true.

    end if
    
    open(file_number, file= trim(arq)//'.'//ext, status='new', iostat=ierr)
    
    if(ierr /= 0) then
       
       write(*, *)
       write(*, *) 'Unable to open file '//trim(arq)//'.'//ext
       write(*, *)
       stop
       
    endif
    
  end subroutine abre
  !==========================================================================

  subroutine def_frame(frame, fr_in, fr_end, skip, n_skip, end, begin)

    logical :: begin, end, skip
    
    integer :: fr_in, fr_end, n_skip, frame
    
    frame = 0
    
    if (.not.begin) then
       
       fr_in = 1
       
    end if
    
    if (.not.end) then
       
       fr_end = fr_in + 2 !==para que fr_end seja maior que frame                                                                                  
    end if
    
    if (.not.skip) then
       
       n_skip = 0
       
    end if
    
  end subroutine def_frame
  !=================================================================================

  subroutine def_bin(outer, bin_coarse, n_index, bin, n_grid, bin_out)

    logical :: bin, outer    
    integer :: n_grid, bin_out, bin_coarse, n_index
    
    bin_coarse = nint(sqrt(real(n_index-1)) - 1)
    
    if (.not.bin)then
       
       n_grid = nint(sqrt(real(n_index-1)) - 1)
       write(*, *)
       write(*, *) 'STD_BIN = ', n_grid
       
    end if
    
  end subroutine def_bin
  !====================================================================

  subroutine def_bin_sph(bin_coarse, n_index, bin, n_grid)

    logical :: bin
    integer :: n_grid, bin_coarse, n_index

    
    bin_coarse = nint(sqrt(real(2*(n_index-1))))
    
    if (.not.bin)then
       
       n_grid = nint(sqrt(real(2*(n_index-1))))
       write(*, *)
       write(*, *) 'STD_BIN = ', n_grid
       
    end if

  end subroutine def_bin_sph
  !====================================================================
  
  subroutine abre_ndx(ind, in_num, n_index)

    character(len=30), intent(in) :: ind

    integer :: ierr, n_index, aux
    integer, dimension(500000) :: in_num

    open(2, file=ind, status='old', iostat=ierr)
    
    if(ierr /= 0) then
       
       write(*, *)
       write(*, *) 'Unable to open file ', ind
       write(*, *)
       stop
       
    endif
     
    ierr = 0
    n_index = 1
    
    do while (ierr==0)
       
       read(2, *, iostat=ierr) aux
       
       if (ierr==0) then
          
          in_num(n_index) = aux
          n_index = n_index + 1
          
       end if
       
    end do
    
    if(ierr>0) then
       
       write(*, *)
       write(*, *) 'Problem by reading ', ind
       write(*, *)
       stop
       
    end if
    
    do i=n_index, 500000
       
       in_num(i) = -1
       
    end do
    
    write(*, *)
    write(*, *) ind, 'input file is OK'

    close(2)

  end subroutine abre_ndx
  !====================================================================

  subroutine abre_trj(n_file, name)
    
    character(len=30) :: name 
    integer ::n_file

    open(n_file, file=name, status='old', iostat=ierr)
    
    if(ierr /= 0) then
       
       write(*, *)
       write(*, *) 'Unable to open file',name
       write(*, *)
       stop
       
    endif
  
  end subroutine abre_trj
  !====================================================================

  real function calc_rmsd(store, store2, grid, grid2, x_min, y_min, dx, dy, noi1, noi2)

    use types

    integer :: i, noi1, noi2, a, b
    
    real :: dist_z, noir
    real :: x_min, y_min, dx, dy
    
    type(vet1), dimension(50000) :: store, store2
    type(vet2), dimension(1001,1001) :: grid, grid2

    noir = 0
    
    do i=1, noi1
       
       a = nint((store(i)%x-x_min)/dx) + 1
       b = nint((store(i)%y-y_min)/dy) + 1
       dist_z = store(i)%z-grid(a,b)%z
       noir = noir + dist_z*dist_z/(noi1+noi2)
       
    end do
    
    do i=1, noi2
       
       a = nint((store2(i)%x-x_min)/dx) + 1
       b = nint((store2(i)%y-y_min)/dy) + 1
       dist_z = store2(i)%z-grid2(a,b)%z
       noir = noir + dist_z*dist_z/(noi1+noi2)
       
    end do
 
    calc_rmsd = sqrt(noir)

  end function calc_rmsd
  !====================================================================

  function calc_rmsd_sph(noi1, noi2, store, store2, grid, grid2, dph, dth)

    use types
    
    integer :: noi1, noi2, i, j, a, b
    real :: noir, dist_r, dth, dph
    type(vet3) :: store(50000), store2(50000)
    type(vet3) :: grid(1001, 1001), grid2(1001, 1001)
    
    noir = 0
    
    do i=1, noi1
       
       a = nint((store(i)%phi)/dph) + 1
       b = nint((store(i)%theta)/dth) + 1
       dist_r = store(i)%rho-grid(a,b)%rho
       noir = noir + dist_r*dist_r/(noi1+noi2)
       
    end do
    
    do i=1, noi2
       
       a = nint((store2(i)%phi)/dph) + 1
       b = nint((store2(i)%theta)/dth) + 1
       dist_r = store2(i)%rho-grid2(a,b)%rho
       noir = noir + dist_r*dist_r/(noi1+noi2)
       
    end do
    
    calc_rmsd_sph = sqrt(noir)
    
  end function calc_rmsd_sph
  !====================================================================

  real function calc_rmsd_inert(noi1, store, dz, dth, grid)

    use types

    integer :: i, noi1, a, b
    real :: dist_z, noir, dth, dz
    type(vet3), dimension(50000) :: store
    type(vet3), dimension(1001,1001) :: grid
    
    noir = 0
    
    do i=1, noi1
       
       !esse calculo é único para este programa                                                                      
       a = nint((cos(store(i)%phi)+1)/dz + 1/2)
       b = nint((store(i)%theta)/dth) + 1
       
       dist_z = store(i)%rho-grid(a,b)%rho
       
       noir = noir + dist_z*dist_z/(noi1)
       
    end do
    
    calc_rmsd_inert = sqrt(noir)
    
  end function calc_rmsd_inert
  !====================================================================
  
  subroutine calc_dens_sph(a_dens, hist, s_vol, dens, grid, grid2, div, del, dth, dph)

    use types

    integer :: a, b, i, bini, a_dens, div
    real :: k_inf, k_sup, del, rad_aver, s_vol, hist(1000)
    type(vet3) :: grid(1001, 1001), grid2(1001, 1001)
    type(vet3), dimension(1000000) :: dens
    
    do i=1, a_dens - 1
       
       a = nint((dens(i)%phi)/dph) + 1
       b = nint((dens(i)%theta)/dth) + 1
       
       rad_aver = (grid(a,b)%rho + grid2(a,b)%rho)/2
       del = 2.0/div
       
       bini = int((dens(i)%rho/rad_aver)/del)
       k_inf = bini*del
       k_sup = bini*del + del

       hist(bini+500) = hist(bini+500) + 1*1000/(s_vol*(k_sup**3 - k_inf**3))
       
    end do
    
  end subroutine calc_dens_sph
  !====================================================================
  
  real function calc_area(n_grid, grid3)

    use types
    
    integer ::i, j, n_grid
    real :: dist_x, dist_y, dist_z, la, lb, lc, s_soma, s_area
    type(vet2) :: grid3(1001, 1001)

    ! Calculando a área pela fórmula de Heron
    s_area = 0
    
    do i=2, n_grid+1
       
       do j=2, n_grid+1
          
          dist_x = abs(grid3(i-1,j-1)%x - grid3(i-1,j)%x)
          dist_y = abs(grid3(i-1,j-1)%y - grid3(i-1,j)%y)
          dist_z = abs(grid3(i-1,j-1)%z - grid3(i-1,j)%z)
          
          la = sqrt(dist_x**2 + dist_y**2 + dist_z**2)
          
          dist_x = abs(grid3(i-1,j-1)%x - grid3(i,j-1)%x)
          dist_y = abs(grid3(i-1,j-1)%y - grid3(i,j-1)%y)
          dist_z = abs(grid3(i-1,j-1)%z - grid3(i,j-1)%z)
          
          lb = sqrt(dist_x**2 + dist_y**2 + dist_z**2)
          
          dist_x = abs(grid3(i-1,j)%x - grid3(i,j-1)%x)
          dist_y = abs(grid3(i-1,j)%y - grid3(i,j-1)%y)
          dist_z = abs(grid3(i-1,j)%z - grid3(i,j-1)%z)
          
          lc = sqrt(dist_x**2 + dist_y**2 + dist_z**2)
          
          s_soma = (la + lb + lc)/2
          
          s_area = s_area + sqrt(s_soma*(s_soma-la)*(s_soma-lb)*(s_soma-lc))
          
          ! segundo triângulo
          
          dist_x = abs(grid3(i-1,j)%x - grid3(i,j)%x)
          dist_y = abs(grid3(i-1,j)%y - grid3(i,j)%y)
          dist_z = abs(grid3(i-1,j)%z - grid3(i,j)%z)
          
          la = sqrt(dist_x**2 + dist_y**2 + dist_z**2)
          
          dist_x = abs(grid3(i,j-1)%x - grid3(i,j)%x)
          dist_y = abs(grid3(i,j-1)%y - grid3(i,j)%y)
          dist_z = abs(grid3(i,j-1)%z - grid3(i,j)%z)
          
          lb = sqrt(dist_x**2 + dist_y**2 + dist_z**2)
          
          s_soma = (la + lb + lc)/2
          
          s_area = s_area + sqrt(s_soma*(s_soma-la)*(s_soma-lb)*(s_soma-lc))
          
       end do
       
    end do

    calc_area = s_area
    
  end function calc_area
  !====================================================================

  subroutine calc_area_sph(lim_i, lim_j, s_area, s_vol, dph, dth, grid, grid3)

    use types
    
    integer :: lim_i, lim_j
    real :: s_area, s_vol, dist_x, dist_y, dist_z
    real :: la, lb, lc, s_soma, c_angle, aux2, dph, dth
    type(vet2) :: grid3(1001, 1001)
    type(vet3) :: grid(1001, 1001)
    
    s_area = 0
    s_vol  = 0
    
    do i=2, lim_i
       
       do j=2, lim_j
          
          dist_x = abs(grid3(i-1,j-1)%x - grid3(i-1,j)%x)
          dist_y = abs(grid3(i-1,j-1)%y - grid3(i-1,j)%y)
          dist_z = abs(grid3(i-1,j-1)%z - grid3(i-1,j)%z)
          
          la = sqrt(dist_x**2 + dist_y**2 + dist_z**2)
          
          dist_x = abs(grid3(i-1,j-1)%x - grid3(i,j-1)%x)
          dist_y = abs(grid3(i-1,j-1)%y - grid3(i,j-1)%y)
          dist_z = abs(grid3(i-1,j-1)%z - grid3(i,j-1)%z)
          
          lb = sqrt(dist_x**2 + dist_y**2 + dist_z**2)
          
          dist_x = abs(grid3(i-1,j)%x - grid3(i,j-1)%x)
          dist_y = abs(grid3(i-1,j)%y - grid3(i,j-1)%y)
          dist_z = abs(grid3(i-1,j)%z - grid3(i,j-1)%z)
          
          lc = sqrt(dist_x**2 + dist_y**2 + dist_z**2)
          
          s_soma = (la + lb + lc)/2
          
          aux2 = sqrt(s_soma*abs(s_soma-la)*abs(s_soma-lb)*abs(s_soma-lc))
          c_angle = ang(grid3(i-1,j-1), grid3(i-1,j), grid3(i,j-1), (i-1-0.5)*dph, (j-1-0.5)*dth)
          s_area = s_area + aux2
          
          s_vol = s_vol + c_angle*aux2*grid(i-1,j-1)%rho/3
          
          ! segundo triângulo
          
          dist_x = abs(grid3(i-1,j)%x - grid3(i,j)%x)
          dist_y = abs(grid3(i-1,j)%y - grid3(i,j)%y)
          dist_z = abs(grid3(i-1,j)%z - grid3(i,j)%z)
          
          la = sqrt(dist_x**2 + dist_y**2 + dist_z**2)
          
          dist_x = abs(grid3(i,j-1)%x - grid3(i,j)%x)
          dist_y = abs(grid3(i,j-1)%y - grid3(i,j)%y)
          dist_z = abs(grid3(i,j-1)%z - grid3(i,j)%z)
          
          lb = sqrt(dist_x**2 + dist_y**2 + dist_z**2)
          
          s_soma = (la + lb + lc)/2
          
          aux2 = sqrt(s_soma*abs(s_soma-la)*abs(s_soma-lb)*abs(s_soma-lc))
          c_angle = ang(grid3(i-1,j), grid3(i,j-1), grid3(i,j), (i-1-0.5)*dph, (j-1-0.5)*dth)
          s_area = s_area + aux2
          
          s_vol = s_vol + c_angle*aux2*grid(i,j)%rho/3
          
       end do
       
    end do

  end subroutine calc_area_sph
  !====================================================================
  
  subroutine calc_order(n_grid, grid, r_xpm, aver, aver2, hist)

    use types
    real, parameter :: pi = 3.141592654
    integer i, j, n_grid, bini
    real :: aver, aver2, r_xpm(n_grid, n_grid), la, aux2, hist(100)
    type(vet2) :: v1, v2, v3
    type(vet2), dimension(1001,1001) :: grid
    
    aver = 0
    aver2 = 0
    
    do i=2, n_grid+1
       
       do j=2, n_grid+1
          
          v1%x = grid(i,j)%x-grid(i-1,j-1)%x
          v1%y = grid(i,j)%y-grid(i-1,j-1)%y
          v1%z = grid(i,j)%z-grid(i-1,j-1)%z
          
          v2%x = grid(i-1,j)%x-grid(i,j-1)%x
          v2%y = grid(i-1,j)%y-grid(i,j-1)%y
          v2%z = grid(i-1,j)%z-grid(i,j-1)%z
          
          v3%x = - v2%y*v1%z + v2%z*v1%y
          v3%y = - v2%z*v1%x + v2%x*v1%z
          v3%z = - v2%x*v1%y + v2%y*v1%x
          
          la = acos(v3%z/sqrt(v3%x**2 + v3%y**2 + v3%z**2))*180/pi
          
          aux2 = 0.5*(3*cos(la*pi/180)**2 - 1)
          
          r_xpm(i-1, j-1) = r_xpm(i-1, j-1) + aux2
          
          aver = aver + aux2/(n_grid*n_grid)
          aver2 = aver2 + aux2*aux2/(n_grid*n_grid)
          
          if ((la<0).or.(la>90)) then
             
             write(*, *) "problem ... la = ", la
             stop
             
          end if
          
          bini = nint(la) + 1
          
          hist(bini) = hist(bini) + 1
          
       end do
       
    end do
    
  end subroutine calc_order
  !====================================================================

  subroutine calc_order_sph(lim_i, lim_j, dph, dth, grid3, r_xpm1, aver, aver2, hist)

    use types
    real, parameter :: pi = 3.141592654
    integer i, j, lim_i, lim_j, bini
    real :: aver, aver2, r_xpm1(lim_i-1, lim_j-1), la, aux2, hist(1000)
    real :: dth, dph
    type(vet2) :: v1, v2, v3
    type(vet2), dimension(1001,1001) :: grid3
    
    aver = 0
    aver2 = 0
    
    do i=2, lim_i
       
       do j=2, lim_j
          
          v1%x = grid3(i,j)%x-grid3(i-1,j-1)%x
          v1%y = grid3(i,j)%y-grid3(i-1,j-1)%y
          v1%z = grid3(i,j)%z-grid3(i-1,j-1)%z
          
          v2%x = grid3(i-1,j)%x-grid3(i,j-1)%x
          v2%y = grid3(i-1,j)%y-grid3(i,j-1)%y
          v2%z = grid3(i-1,j)%z-grid3(i,j-1)%z
          
          v3%x = - v2%y*v1%z + v2%z*v1%y
          v3%y = - v2%z*v1%x + v2%x*v1%z
          v3%z = - v2%x*v1%y + v2%y*v1%x
          
          !reutilizando a variável v1 para facilitar o cálculo ======
          ! v1 é o gradiente definido entre os quatro pontos avaliados
          
          v1%x = sin((i-1-0.5)*dph)*cos((j-1-0.5)*dth)
          v1%y = sin((i-1-0.5)*dph)*sin((j-1-0.5)*dth)
          v1%z = cos((i-1-0.5)*dph)
          
          !==========================================================
          
          la = v1%x*v3%x + v1%y*v3%y + v1%z*v3%z
          la = la/sqrt(v3%x**2 + v3%y**2 + v3%z**2)
          la = la/sqrt(v1%x**2 + v1%y**2 + v1%z**2)
          
          !erros de arredondamento levam a valores 
          !levemente maiores que 1.00 (1.00000012)
          la = min(la, 1.00)
          
          aux2 = 0.5*(3*(la)**2 - 1)
          
          r_xpm1(i-1, j-1) = r_xpm1(i-1, j-1) + aux2
          
          aver = aver + aux2/((lim_i-1)*(lim_j-1))
          aver2 = aver2 + aux2*aux2/((lim_i-1)*(lim_j-1))
          
          la = acos(la)*180/pi !colocando em graus

          if ((la<0).or.(la>90)) then
             
             write(*, *) "problem ... la = ", la
             stop
             
          end if
          
          bini = nint(la) + 1
          
          hist(bini) = hist(bini) + 1
          
       end do
       
    end do
   
  end subroutine calc_order_sph
  !====================================================================
  
  subroutine calc_thick(n_grid, aver, aver2, s_v, dx, dy, grid, grid2, grid3, r_xpm)

    use types
    
    integer :: i, j, n_grid 
    real :: aver, aver2, r_xpm(n_grid, n_grid)
    real :: la, lb, aux2, s_v, dx, dy
    type(vet2) :: v1, v2, v3
    type(vet2), dimension(1001,1001) :: grid, grid2, grid3
    
    s_v = 0
    aver = 0
    aver2 = 0
    
    do i=2, n_grid+1
       
       do j=2, n_grid+1
          
          v1%x = grid3(i,j)%x-grid3(i-1,j-1)%x
          v1%y = grid3(i,j)%y-grid3(i-1,j-1)%y
          v1%z = grid3(i,j)%z-grid3(i-1,j-1)%z
          
          v2%x = grid3(i-1,j)%x-grid3(i,j-1)%x
          v2%y = grid3(i-1,j)%y-grid3(i,j-1)%y
          v2%z = grid3(i-1,j)%z-grid3(i,j-1)%z
          
          v3%x = - v2%y*v1%z + v2%z*v1%y
          v3%y = - v2%z*v1%x + v2%x*v1%z
          v3%z = - v2%x*v1%y + v2%y*v1%x
          
          la = (grid(i-1,j-1)%z +grid(i-1,j)%z + grid(i,j-1)%z + grid(i,j)%z)/4
          lb = (grid2(i-1,j-1)%z +grid2(i-1,j)%z + grid2(i,j-1)%z + grid2(i,j)%z)/4
          r_xpm(i-1,j-1) = r_xpm(i-1, j-1) + abs(v3%z*(la-lb))/sqrt(v3%x**2 + v3%y**2 + v3%z**2)
          
          s_v = s_v + abs((la-lb)*dx*dy)
          
          aux2 = abs(v3%z*(la-lb))/sqrt(v3%x**2 + v3%y**2 + v3%z**2)/10
          aver = aver + aux2/((n_grid+1)*(n_grid+1))
          aver2 = aver2 + aux2*aux2/((n_grid+1)*(n_grid+1))
          
       end do
       
    end do

  end subroutine calc_thick
  !====================================================================

  subroutine calc_thick_sph(lim_i, lim_j, aver, aver2, r_xpm1, grid, grid2)

    use types

    integer :: i, j, lim_i, lim_j
    real :: aver, aver2, aux2, r_xpm1(lim_i-1, lim_j-1)
    type(vet3) :: grid(1001, 1001), grid2(1001, 1001)
    
    aver = 0
    aver2 = 0
    
    do i=2, lim_i
       
       do j=2, lim_j
          
          aux2 = grid(i,j)%rho + grid(i-1,j)%rho + grid(i,j-1)%rho + grid(i-1,j-1)%rho
          aux2 = aux2 - (grid2(i,j)%rho + grid2(i-1,j)%rho + grid2(i,j-1)%rho + grid2(i-1,j-1)%rho)
          aux2 = abs(aux2/4)/10
          aver = aver + aux2/(lim_i*lim_j)
          aver2 = aver2 + aux2*aux2/(lim_i*lim_j)
          
          r_xpm1(i-1,j-1) = r_xpm1(i-1,j-1) + aux2*10
          
       end do
       
    end do
    
  end subroutine calc_thick_sph
  !====================================================================
  
  subroutine calc_topog(n_grid, aver, aver2, r_xpm, grid, grid2)

    use types

    integer :: n_grid, i, j
    real :: aver, aver2, r_xpm(n_grid, n_grid), la, lb, lc
    type(vet2), dimension(1001,1001) :: grid, grid2
    
    lc = 0
    aver = 0
    aver2 = 0
    
    do i=2, n_grid+1
       
       do j=2, n_grid+1
          
          la = (grid(i-1,j-1)%z +grid(i-1,j)%z + grid(i,j-1)%z + grid(i,j)%z)/4
          lb = (grid2(i-1,j-1)%z + grid2(i-1,j)%z + grid2(i,j-1)%z + grid2(i,j)%z)/4
          lc = (la+lb)
          
          r_xpm(i-1,j-1) = r_xpm(i-1, j-1) +lc/10
          
          aver = aver + lc/10/(n_grid*n_grid)
          aver2 = aver2 + lc*lc/100/(n_grid*n_grid)
          
       end do

    end do
    
  end subroutine calc_topog
  !====================================================================
  
  subroutine calc_inertia(n_grid, grid, grid3, MI)

    use types

    integer :: i, j, n_grid
    real :: averx, avery, averz
    double precision :: MI(3,3)
    type(vet3), dimension(1001,1001) :: grid
    type(vet2), dimension(1001,1001) :: grid3
    
    do i=1, 3

       do j=1, 3

          MI(i,j) = 0
    
       end do

    end do
    
    averx = 0
    avery = 0
    averz = 0
    
    do i=1, (n_grid) 
       
       do j=1, n_grid 
          
          grid3(i,j)%x = grid(i,j)%rho*sin(grid(i,j)%phi)*cos(grid(i,j)%theta)/10
          grid3(i,j)%y = grid(i,j)%rho*sin(grid(i,j)%phi)*sin(grid(i,j)%theta)/10
          grid3(i,j)%z = grid(i,j)%rho*cos(grid(i,j)%phi)/10
          
          averx = averx + grid3(i,j)%x
          avery = avery + grid3(i,j)%y
          averz = averz + grid3(i,j)%z
          
       end do
       
    end do
    
    averx = averx/(n_grid*n_grid)
    avery = avery/(n_grid*n_grid)
    averz = averz/(n_grid*n_grid)
    
    do i=1, (n_grid) 
       
       do j=1, n_grid 
          
          MI(1,1) = MI(1,1) + (grid3(i,j)%y-avery)**2 + (grid3(i,j)%z-averz)**2
          MI(2,2) = MI(2,2) + (grid3(i,j)%x-averx)**2 + (grid3(i,j)%z-averz)**2
          MI(3,3) = MI(3,3) + (grid3(i,j)%x-averx)**2 + (grid3(i,j)%y-avery)**2
          MI(1,2) = MI(1,2) - (grid3(i,j)%x-averx)*(grid3(i,j)%y-avery)
          MI(1,3) = MI(1,3) - (grid3(i,j)%x-averx)*(grid3(i,j)%z-averz)
          MI(2,3) = MI(2,3) - (grid3(i,j)%y-avery)*(grid3(i,j)%z-averz)
          
       end do
       
    end do
    
    MI(1,1) = MI(1,1)/(n_grid*n_grid)
    MI(2,2) = MI(2,2)/(n_grid*n_grid)
    MI(3,3) = MI(3,3)/(n_grid*n_grid)
    MI(1,2) = MI(1,2)/(n_grid*n_grid)
    MI(1,3) = MI(1,3)/(n_grid*n_grid)
    MI(2,3) = MI(2,3)/(n_grid*n_grid)
    
    MI(2,1) = MI(1,2)
    MI(3,1) = MI(1,3)
    MI(3,2) = MI(2,3)
    
  end subroutine calc_inertia
  !====================================================================
  
  subroutine print_pdb(grid, lim_i, lim_j, file, boxx, boxy, boxz, frame, back)
    
    use types
    
12  format(a6, i5.1, a5, a5, a1, i4.1, a4, 3f8.3)
3   format(a8, f7.2, a2, f7.2, a2, f7.2, a35)
    
    logical :: back
    integer :: ierr, i, j, lim_i, lim_j, file, frame
    real :: boxx, boxy, boxz
    type(vet2), dimension(1001,1001) :: grid

    write(4, '(a26)') 'REMARK Generated by s_grid'
    write(4, '(a34)') 'TITLE     Rectangular Regular GRID'
    write(4, '(a22)') 'REMARK trajectory file'
    write(4, 3) 'CRYST1  ', boxx, '  ', boxy, '  ', boxz, '  90.00  90.00  90.00 P  1       1 '
    write(4, '(a6, i4)') 'MODEL ', frame
    
    do i=1, lim_i
       
       do j=1, lim_j
          
          write(file, 12, iostat=ierr) 'ATOM  ', (i-1)*(lim_j) + j, '  DOT', &
               '  DOT',' ', 1,'    ', grid(i,j)%x, &
               grid(i,j)%y, grid(i,j)%z
          
       end do
       
    end do
    
    write(4, '(a3)') 'TER'
    write(4, '(a6)') 'ENDMDL'
     
  end subroutine print_pdb
  !====================================================================

  subroutine print_grid(gridc, lim_i, lim_j, name, back)

    use types
    
12  format(a6, i5.1, a5, a5, a1, i4.1, a4, 3f8.3)
    
    logical :: ex, back
    integer :: ierr, i, j, lim_i, lim_j
    character(len=5) :: name
    character(len=30) :: junta
    type(vet2), dimension(1001,1001) :: gridc
    
    junta = name//'.pdb'
    inquire(file=junta, exist=ex)
    
    if (ex) then
       
       junta = 'mv  '//name//'.pdb '//name//'_p.pdb'
       call execute_command_line(junta)
       back = .true.

    end if
    
    open(3, file=name//'.pdb', status='new', iostat=ierr)
    
    if(ierr /= 0) then
       
       write(*, *)
       write(*, *) 'Unable to open file '//name//'.pdb'
       write(*, *)
       stop
       
    end if
    
    write(3, *) 'TITLE     Rectangular Regular GRID'
    write(3, *) 'REMARK    THIS IS A SIMULATION BOX'
    write(3, *) 'MODEL        1'
    
    do i=1, lim_i
       
       do j=1, lim_j
          
          write(3, 12, iostat=ierr) 'ATOM  ', (i-1)*(lim_j) + j, '  DOT', &
               '  DOT',' ', 1,'    ', gridc(i,j)%x, &
               gridc(i,j)%y, gridc(i,j)%z
          
       end do
       
    end do
    
    write(3, *) 'TER'
    write(3, *) 'ENDMDL'
    
    close(3)
    
  end subroutine print_grid
  !====================================================================

  subroutine print_grid_xpm(gridc, r_xpm, lim_i, lim_j, name, back)

    use types
    
13  format(a6, i5.1, a5, a5, a1, i4.1, a4, 3f8.3, f5.1)
    
    logical :: ex, back
    integer :: ierr, i, j, lim_i, lim_j
    real :: r_xpm(lim_i, lim_j)
    character(len=5) :: name
    character(len=30) :: junta
    type(vet2), dimension(1001,1001) :: gridc
    
    junta = name//'.pdb'
    inquire(file=junta, exist=ex)
    
    if (ex) then
       
       junta = 'mv  '//name//'.pdb '//name//'_p.pdb'
       call execute_command_line(junta)
       back = .true.

    end if
    
    open(3, file=name//'.pdb', status='new', iostat=ierr)
    
    if(ierr /= 0) then
       
       write(*, *)
       write(*, *) 'Unable to open file '//name//'.pdb'
       write(*, *)
       stop
       
    end if
    
    write(3, *) 'TITLE     Rectangular Regular GRID'
    write(3, *) 'REMARK    THIS IS A SIMULATION BOX'
    write(3, *) 'MODEL        1'
    
    do i=1, lim_i
       
       do j=1, lim_j

          write(3, 13, iostat=ierr) 'ATOM  ', (i-1)*(lim_j) + j, '  DOT', &
               '  DOT',' ', 1,'    ', gridc(i,j)%x, &
               gridc(i,j)%y, gridc(i,j)%z, r_xpm(i,j)
          
       end do
       
    end do
    
    write(3, *) 'TER'
    write(3, *) 'ENDMDL'
    
    close(3)
    
  end subroutine print_grid_xpm
  !====================================================================
  
  subroutine print_grid_xpm_xpm(gridc, r_xpm1, r_xpm2, lim_i, lim_j, name, back)

    use types
    
13  format(a6, i5.1, a5, a5, a1, i4.1, a4, 3f8.3, 2f5.1)
    
    logical :: ex, back
    integer :: ierr, i, j, lim_i, lim_j
    real :: r_xpm1(lim_i, lim_j), r_xpm2(lim_i, lim_j)
    character(len=5) :: name
    character(len=30) :: junta
    type(vet2), dimension(1001,1001) :: gridc
    
    junta = name//'.pdb'
    inquire(file=junta, exist=ex)
    
    if (ex) then
       
       junta = 'mv  '//name//'.pdb '//name//'_p.pdb'
       call execute_command_line(junta)
       back = .true.

    end if
    
    open(3, file=name//'.pdb', status='new', iostat=ierr)
    
    if(ierr /= 0) then
       
       write(*, *)
       write(*, *) 'Unable to open file '//name//'.pdb'
       write(*, *)
       stop
       
    end if
    
    write(3, *) 'TITLE     Rectangular Regular GRID'
    write(3, *) 'REMARK    THIS IS A SIMULATION BOX'
    write(3, *) 'MODEL        1'
    
    do i=2, lim_i, 2 !limites do s_gauss, onde essa rotina é utilizada
       
       do j=2, lim_j, 2
          
          write(3, 13, iostat=ierr) 'ATOM  ', (i-1)*(lim_j) + j, '  DOT', &
               '  DOT',' ', 1,'    ', gridc(i,j)%x, &
               gridc(i,j)%y, gridc(i,j)%z, r_xpm1(i,j), r_xpm2(i,j)
          
       end do
       
    end do
    
    write(3, *) 'TER'
    write(3, *) 'ENDMDL'
    
    close(3)
    
  end subroutine print_grid_xpm_xpm
  !====================================================================

  subroutine calcula_bin(n_grid, X, Y)
    
    integer :: n_grid
    real :: X, Y

    !o valor de 0.08 veio de um estudo com diferentes funções
    !onde verificou-se que valores de dx ou dy menores que 0.08
    !poderiam levar a erros nos calculos das curvaturas                                                                                                                                                       
    if (X>Y) then
       
       n_grid = int(Y/0.08)
       
    else
       
       n_grid = int(X/0.08)
       
    end if
    
    if (mod(n_grid,2) /= 0.0)  n_grid = n_grid + 1
    
    if (n_grid>500) n_grid = 500
    
  end subroutine calcula_bin
  !====================================================================
  
  subroutine ending(back, finish, start, clock_rate)

    logical :: back
    integer :: start, finish, clock_rate
    real ::  hour, minu, sec

1   format(a19, i6.3, a2, i6.3, a4, i6.3, a4)
   
    write(*, *)
    write(*, *) "Finished"
    write(*, *)
    if (back) write(*, '(a31)') ' Previous files were backed up!'
    
    hour = real(finish-start)/(3600*real(clock_rate))
    minu = (hour - int(hour))*60
    sec = (minu-int(minu))*60
    
    if (back) write(*, *)
    write(*, 1) " Processing time : ", int(hour), " h", int(minu), " min", int(sec), " sec"
    write(*, *)
    
  end subroutine ending
  !====================================================================
  
  subroutine print_xpm(n_grid, del, dx, dy, gx, gy, xpm, r_xpm1, minv, cor)
    
7   format(a20, f6.3, a5)
3   format(a10)
4   format(f6.2, a1)
5   format(a1)

    character :: xpm(n_grid, n_grid), cor
    integer :: i, j, n_grid 
    real :: del, r_xpm1(n_grid, n_grid), minv
    real :: dx, dy, gy, gx

    if (cor .eq. 'r') then
    
       !vamos considerar que o arquivo será sempre o 4
       write(4, 7) '"A  c #993333" /* "',minv,'" */,'
       write(4, 7) '"B  c #cc0000" /* "',(minv+del),'" */,'
       write(4, 7) '"C  c #ff5050" /* "',(minv+2*del),'" */,'
       write(4, 7) '"D  c #ff6666" /* "',(minv+3*del),'" */,'
       write(4, 7) '"E  c #ff9999" /* "',(minv+4*del),'" */,'
       write(4, 7) '"F  c #ffcccc" /* "',(minv+5*del),'" */,'
       write(4, 7) '"G  c #ffffff" /* "',(minv+6*del),'" */,'
       
    end if

    if(cor .eq. 'b') then

       write(4, 7) '"A  c #ffffff " /* "',minv,'" */,'
       write(4, 7) '"B  c #87cefa " /* "',(minv+del),'" */,'
       write(4, 7) '"C  c #00bfff " /* "',(minv+2*del),'" */,'
       write(4, 7) '"D  c #1e90ff " /* "',(minv+3*del),'" */,'
       write(4, 7) '"E  c #4169e1 " /* "',(minv+4*del),'" */,'
       write(4, 7) '"F  c #0000ff " /* "',(minv+5*del),'" */,'
       write(4, 7) '"G  c #00008b " /* "',(minv+6*del),'" */,'

    end if

    do i=1, n_grid
       
       do j=1, n_grid

          if (r_xpm1(i,j)<minv+del)then
             xpm(i,j)= 'A'
          else
             if (r_xpm1(i,j)<minv+2*del)then
                xpm(i,j) = 'B'
             else
                if (r_xpm1(i,j)<minv+3*del)then
                   xpm(i,j) = 'C'
                else
                   if (r_xpm1(i,j)<minv+4*del) then
                      xpm(i,j) = 'D'
                   else
                      if (r_xpm1(i,j)<minv+5*del) then
                         xpm(i,j) = 'E'
                      else
                         if (r_xpm1(i,j)<minv+6*del) then
                            xpm(i,j) = 'F'
                         else
                            xpm(i,j) = 'G'
                            
                         end if
                         
                      end if
                      
                   end if
                   
                end if
                
             end if
             
          end if

       end do
       
    end do

    k = 1
    
    do i=1, int(n_grid/20)+1
       
       write(4, 3, advance='no') '/* y-axis:'
       j = 1
       
       do while ((j<=20).and.(k<=n_grid))
          
          write(4, 4, advance='no') ((k-1/2)*dy + gy)/10, ' '
          k = k+1
          j = j+1
          
       end do
       write(4,*) '*/'
       
    end do
    
    k = 1
    do i=1, int(n_grid/20)+1
       
       write(4, 3, advance='no') '/* x-axis:'
       j = 1
       
       do while ((j<=20).and.(k<=n_grid))
          
          write(4, 4, advance='no') ((k-1/2)*dx + gx)/10, ' '
          k = k+1
          j = j+1
          
       end do
       write(4,*) '*/'
       
    end do
    
    do j=n_grid, 1, -1
       
       write(4, 5, advance='no') '"'
       
       do i=1, n_grid
          
          write(4, 5, advance='no') xpm(i,j)
          
       end do
       
       write(4, 5, advance='no')'"'
       write(4, 5) ','
       
    end do
    
  end subroutine print_xpm
  !====================================================================

  subroutine print_xpm_gauss(n_grid, del, dx, dy, gx, gy, xpm, r_xpm1, minv, cor)
              
7   format(a20, f6.3, a5)
3   format(a10)
4   format(f6.2, a1)
5   format(a1)
  
    character :: xpm(n_grid, n_grid), cor
    integer :: i, j, n_grid 
    real :: del, r_xpm1(n_grid, n_grid), minv
    real :: dx, dy, gy, gx

    if (cor .eq. 'r') then

       !vamos considerar que o arquivo será sempre o 4                                                                                                                                                      
       write(4, 7) '"A  c #993333" /* "',minv,'" */,'
       write(4, 7) '"B  c #cc0000" /* "',(minv+del),'" */,'
       write(4, 7) '"C  c #ff5050" /* "',(minv+2*del),'" */,'
       write(4, 7) '"D  c #ff6666" /* "',(minv+3*del),'" */,'
       write(4, 7) '"E  c #ff9999" /* "',(minv+4*del),'" */,'
       write(4, 7) '"F  c #ffcccc" /* "',(minv+5*del),'" */,'
       write(4, 7) '"G  c #ffffff" /* "',(minv+6*del),'" */,'

    end if

    if(cor .eq. 'b') then

       write(4, 7) '"A  c #ffffff " /* "',minv,'" */,'
       write(4, 7) '"B  c #87cefa " /* "',(minv+del),'" */,'
       write(4, 7) '"C  c #00bfff " /* "',(minv+2*del),'" */,'
       write(4, 7) '"D  c #1e90ff " /* "',(minv+3*del),'" */,'
       write(4, 7) '"E  c #4169e1 " /* "',(minv+4*del),'" */,'
       write(4, 7) '"F  c #0000ff " /* "',(minv+5*del),'" */,'
       write(4, 7) '"G  c #00008b " /* "',(minv+6*del),'" */,'

    end if
  
    do i=2, n_grid, 2
       
       do j=2, n_grid, 2
          
          if (r_xpm1(i,j)<minv+del)then
             xpm(i,j)= 'A'
          else
             if (r_xpm1(i,j)<minv+2*del)then
                xpm(i,j) = 'B'
             else
                if (r_xpm1(i,j)<minv+3*del)then
                   xpm(i,j) = 'C'
                else
                   if (r_xpm1(i,j)<minv+4*del) then
                      xpm(i,j) = 'D'
                   else
                      if (r_xpm1(i,j)<minv+5*del) then
                         xpm(i,j) = 'E'
                      else
                         if (r_xpm1(i,j)<minv+6*del) then
                            xpm(i,j) = 'F'
                         else
                            xpm(i,j) = 'G'
                            
                         end if
                         
                      end if
                      
                   end if
                   
                end if
                
             end if
             
          end if
          
       end do
       
    end do
    
    k = 1
    
    do i=1, int(n_grid/40)+1
       
       write(4, 3, advance='no') '/* y-axis:'
       j = 1
       
       do while ((j<=20).and.(k<=n_grid/2))
          
          write(4, 4, advance='no') ((k*2-1)*dy + gy)/10, ' '
          k = k+1
          j = j+1
          
       end do
       write(4,*) '*/'
       
    end do
    
    k = 1
    do i=1, int(n_grid/40)+1
       
       write(4, 3, advance='no') '/* x-axis:'
       j = 1

       do while ((j<=20).and.(k<=n_grid/2))
          
          write(4, 4, advance='no') ((k*2-1)*dx + gx)/10, ' '
          k = k+1
          j = j+1
          
       end do
       write(4,*) '*/'
       
    end do
    
    do j=n_grid, 2, -2
       
       write(4, 5, advance='no') '"'
       
       do i=2, n_grid, 2
          
          write(4, 5, advance='no') xpm(i,j)
          
       end do
       
       write(4, 5, advance='no')'"'
       write(4, 5) ','
       
    end do

  end subroutine print_xpm_gauss
  !====================================================================

  
  subroutine  print_xpm_sph(lim_i, lim_j, del, dth, dph, xpm, r_xpm1, minv, cor)
    
7   format(a20, f4.1, a5)
3   format(a10)
4   format(f6.2, a1)
5   format(a1)

    character :: xpm(lim_i, lim_j), cor
    integer :: i, j, k
    real :: del, r_xpm1(lim_i, lim_j), minv
    real :: dph, dth
  
    if (cor .eq. 'r') then

       !vamos considerar que o arquivo será sempre o 4                                                                              
       write(4, 7) '"A  c #993333" /* "',minv,'" */,'
       write(4, 7) '"B  c #cc0000" /* "',(minv+del),'" */,'
       write(4, 7) '"C  c #ff5050" /* "',(minv+2*del),'" */,'
       write(4, 7) '"D  c #ff6666" /* "',(minv+3*del),'" */,'
       write(4, 7) '"E  c #ff9999" /* "',(minv+4*del),'" */,'
       write(4, 7) '"F  c #ffcccc" /* "',(minv+5*del),'" */,'
       write(4, 7) '"G  c #ffffff" /* "',(minv+6*del),'" */,'

    end if

    if(cor .eq. 'b') then

       write(4, 7) '"A  c #ffffff " /* "',minv,'" */,'
       write(4, 7) '"B  c #87cefa " /* "',(minv+del),'" */,'
       write(4, 7) '"C  c #00bfff " /* "',(minv+2*del),'" */,'
       write(4, 7) '"D  c #1e90ff " /* "',(minv+3*del),'" */,'
       write(4, 7) '"E  c #4169e1 " /* "',(minv+4*del),'" */,'
       write(4, 7) '"F  c #0000ff " /* "',(minv+5*del),'" */,'
       write(4, 7) '"G  c #00008b " /* "',(minv+6*del),'" */,'

    end if
  
    do i=1, lim_i
       
       do j=1, lim_j
          
          if (r_xpm1(i,j)<minv+del)then
             xpm(i,j)= 'A'
          else
             if (r_xpm1(i,j)<minv+2*del)then
                xpm(i,j) = 'B'
             else
                if (r_xpm1(i,j)<minv+3*del)then
                   xpm(i,j) = 'C'
                else
                   if (r_xpm1(i,j)<minv+4*del) then
                      xpm(i,j) = 'D'
                   else
                      if (r_xpm1(i,j)<minv+5*del) then
                         xpm(i,j) = 'E'
                      else
                         if (r_xpm1(i,j)<minv+6*del) then
                            xpm(i,j) = 'F'
                         else
                            xpm(i,j) = 'G'
                            
                         end if
                         
                      end if
                      
                   end if
                   
                end if
                
             end if
             
          end if
          
       end do
       
    end do
    
    k = lim_i
    
    do i=1, int((lim_i)/20)+1
       
       write(4, 3, advance='no') '/* y-axis:'
       j = 1
       
       do while ((j<=20).and.(k>=1))
          
          write(4, 4, advance='no') ((k-1/2)*dph), ' '
          k = k-1
          j = j+1
          
       end do
       write(4,*) '*/'
       
    end do
    
    k = 1
    do i=1, int((lim_j)/20)+1
       
       write(4, 3, advance='no') '/* x-axis:'
       j = 1
       
       do while ((j<=20).and.(k<=lim_j))
          
          write(4, 4, advance='no') ((k-1/2)*dth), ' '
          k = k+1
          j = j+1
          
       end do
       write(4,*) '*/'
       
    end do
    
    do i=1, lim_i
       
       write(4, 5, advance='no') '"'
       
       do j=1, lim_j
          
          write(4, 5, advance='no') xpm(i,j)
          
       end do
       
       write(4, 5, advance='no')'"'
       write(4, 5) ','
       
    end do
    
  end subroutine print_xpm_sph
  !====================================================================
  
  subroutine imprime(vetor, indice)
    
    character(len=30), intent(in) :: vetor(500)
    integer, intent(in) :: indice

5   format(i4, a10)

    do i=1, int((indice+4)/5)

       do j=1, 5

          if ((i-1)*5+j<=indice)then

             write(*, 5, advance='no') (i-1)*5+j, vetor((i-1)*5+j)

          end if

       end do

       write(*,*)

    end do
    
  end subroutine imprime
  !====================================================================

  subroutine cart2sphe(num, r_med, cent_x, cent_y, cent_z, store, spher)

    use types
    
    real, parameter :: pi = 3.141592654
    integer :: num
    real :: cent_x, cent_y, cent_z, r_med
    type(vet3) :: store
    type(vet1):: spher
    
    store%rho = sqrt((spher%x-cent_x)**2 + (spher%y-cent_y)**2 + (spher%z-cent_z)**2)
    store%phi = acos((spher%z-cent_z)/store%rho)
    
    if ((spher%y-cent_y)==0)then
       
       store%theta = 0.000
       
    else
       
       store%theta = acos((spher%x-cent_x)/sqrt((spher%x-cent_x)**2 + (spher%y-cent_y)**2))
       
       if ((spher%y-cent_y)<=0)then
          
          store%theta = 2*pi - store%theta
          
       end if
       
    end if
    
    r_med = (r_med*(num-1) + store%rho)/num
    
    num = num + 1
    
  end subroutine cart2sphe
  !====================================================================

  subroutine sphe2cart(lim_i, lim_j, cent_x, cent_y, cent_z, grid, grid3)

    use types
    
    integer :: lim_i, lim_j, i, j
    real :: cent_x, cent_y, cent_z
    type(vet3) :: grid(1001, 1001)
    type(vet2) :: grid3(1001, 1001)
    
    do i=1, lim_i
       
       do j=1, lim_j
          
          grid3(i,j)%x = grid(i,j)%rho*sin(grid(i,j)%phi)*cos(grid(i,j)%theta) + cent_x
          grid3(i,j)%y = grid(i,j)%rho*sin(grid(i,j)%phi)*sin(grid(i,j)%theta) + cent_y
          grid3(i,j)%z = grid(i,j)%rho*cos(grid(i,j)%phi) + cent_z
          
       end do
       
    end do
    
  end subroutine sphe2cart
  !====================================================================
  
  subroutine calc_gyrat(n_grid, cent_x, cent_y, cent_z, girat, girat2, grid, grid2)

    use types

    integer :: i, j, n_grid
    real :: cent_x, cent_y, cent_z, girat, girat2
    type(vet3), dimension(1001,1001) :: grid, grid2
    type(vet2), dimension(1001,1001) :: grid3
    
    girat = 0
    girat2 = 0
    
    do i=1, (n_grid - int(n_grid/2))+1
       
       do j=1, n_grid+1
          
          girat = girat + grid(i,j)%rho*grid(i,j)%rho
          girat2 = girat2 + grid2(i,j)%rho*grid2(i,j)%rho
          
       end do
       
    end do
    
    girat = sqrt(girat/((n_grid+1)*(n_grid - int(n_grid/2)+1)))
    girat2 = sqrt(girat2/((n_grid+1)*(n_grid - int(n_grid/2)+1)))

  end subroutine calc_gyrat
  !====================================================================
  
  real function ang(p1, p2, p3, phi, theta)
    
    use types
    
    type(vet2), intent(in) :: p1, p2, p3
    type(vet2) :: v1, v2, v3
    
    real :: la, phi, theta
    
    v1%x = p2%x-p1%x
    v1%y = p2%y-p1%y
    v1%z = p2%z-p1%z
    
    v2%x = p3%x-p2%x
    v2%y = p3%y-p2%y
    v2%z = p3%z-p2%z
    
    v3%x = - v2%y*v1%z + v2%z*v1%y
    v3%y = - v2%z*v1%x + v2%x*v1%z
    v3%z = - v2%x*v1%y + v2%y*v1%x
    
    !reutilizando a variável v1 para facilitar o cálculo ====== 
    
    v1%x = sin(phi)*cos(theta)
    v1%y = sin(phi)*sin(theta)
    v1%z = cos(phi)
    
    !========================================================== 
    
    if (v3%x**2 + v3%y**2 + v3%z**2<0.000001) then
       
       !serve para os casos das extremidades de phi, onde a área
       !do triangulo é nula. Neste caso não importa o ângulo.
       !Por isso utilizamos o valor 1, pois a área já é nula
       
       ang = 1
       
    else
       
       la = v1%x*v3%x + v1%y*v3%y + v1%z*v3%z
       la = la/sqrt(v3%x**2 + v3%y**2 + v3%z**2)
       la = la/sqrt(v1%x**2 + v1%y**2 + v1%z**2)
       
       ang = abs(la)
       
    end if
    
  end function ang
  !====================================================================

  subroutine calc_stat_aver(aver, aver2, desv, skew, kurt, st_mom, n_index, func)

    double precision :: aver, aver2, desv, skew, kurt, st_mom
    integer :: i, n_index
    real, dimension(10000000) ::func
    
    aver = 0
    aver2 = 0
    desv = 0
    
    do i=1, n_index
       
       aver = aver + func(i)/n_index
       aver2 = aver2 + func(i)*func(i)/n_index
       
    end do
    
    desv = sqrt((aver2 - aver*aver))
    
    write(*, *)
    write(*, '(a45)')"****************** SUMMARY ******************"
    write(*, '(a45)')"*                                           *"
    write(*, '(a13, f10.3, a22)')"* Average  = ", aver, "                   *"
    write(*, '(a13, f10.3, a22)')"* S. dev.  = ", desv, "                   *"
    
    skew = 0
    kurt = 0
    do i=1, n_index
       
       st_mom = (func(i) - aver)/desv
       kurt = kurt + st_mom*st_mom*st_mom*st_mom/n_index
       skew = skew + st_mom*st_mom*st_mom/n_index
       
    end do
    
  end subroutine calc_stat_aver
  !====================================================================

  subroutine do_histogram(n_index, hist, minf, desv, aver, aux, del, func, n_file)

    real, parameter :: pi = 3.141592654
    integer :: i, n_index, bini, n_file
    double precision :: hist(1000), minf, desv, aver, aux, del
    real, dimension(10000000) ::func
    
    !criando histograma ===========================
    do i=1, n_index
       
       bini = nint((func(i)-minf)/del)+100
       hist(bini) = hist(bini) + 1/(n_index*del)
       
    end do
    
    !Original =========
    do i=1, 1000
       
       write(n_file, *) (i-100)*del+minf, hist(i)
       
    end do
    
    write(n_file, *)'&'
    
    !Modelo Gaussiano ========
    do i=1, 1000
       
       aux = (i-100)*del+minf
       aux = exp(-(aux-aver)*(aux-aver)/(2*desv*desv))
       aux = aux/(desv*sqrt(2*pi))
       write(n_file, *) (i-100)*del+minf, aux
       
    end do
        
  end subroutine do_histogram
  !====================================================================

  subroutine calc_stat_all(n_index, hist, del, minf, skew, kurt, aver, desv)
    
    integer :: n_index, i, class_modal
    integer :: class_med, class_q1, class_q3, class_d1, class_d9
    
    double precision :: aux, aver, desv, kurt, delta1, delta2
    double precision :: del, hist(1000), sum, acum, moda, mediana
    double precision :: acum25, acum75, quart1, quart3, acum10, acum90
    double precision :: decil1, decil9, skew, minf
    
    !Calculando mediana, moda, quartil e decil=====
    sum = 0
    class_modal = 1
    class_med = 1
    class_q1 = 1
    class_q3 = 1
    class_d1 = 1
    class_d9 = 1
    
    do i=1, 1000
       
       if (hist(i)>hist(class_modal)) class_modal = i
       
       sum = sum + hist(i)*del*n_index
       
       if (sum<n_index/2)then
          
          class_med = i+1
          acum = sum
          
       end if
       
       if (sum<n_index/4)then
          
          class_q1 = i+1
          acum25 = sum
          
       end if
       
       if (sum<3*n_index/4)then
          
          class_q3 = i+1
          acum75 = sum
          
       end if
       
       if (sum<n_index/10)then
          
          class_d1 = i+1
          acum10 = sum
          
       end if
       
       if (sum<9*n_index/10)then
          
          class_d9 = i+1
          acum90 = sum
          
       end if
       
    end do
    
    !moda de Czuber para dados tabulados====
    delta1 = hist(class_modal) - hist(class_modal-1)
    delta2 = hist(class_modal) - hist(class_modal+1)
    moda = (((class_modal-100)*del+minf)-del/2)
    moda = moda + del*(delta1/(delta1+delta2))
    
    !mediana para dados tabulados===========
    mediana = (((class_med-100)*del+minf)-del/2)
    mediana = mediana + (n_index/2 - acum)/(hist(class_med)*n_index)
    
    !primeiro Quartil para dados tabulados =
    quart1 = (((class_q1-100)*del+minf)-del/2)
    quart1 = quart1 + (n_index/4 - acum25)/(hist(class_q1)*n_index)
    
    !terceiro Quartil para dados tabulados =
    quart3 = (((class_q3-100)*del+minf)-del/2)
    quart3 = quart3 + (3*n_index/4 - acum75)/(hist(class_q3)*n_index)
    
    !primeiro Decil para dados tabulados====
    decil1 = (((class_d1-100)*del+minf)-del/2)
    decil1 = decil1 + (n_index/10 - acum10)/(hist(class_d1)*n_index)
    
    !nono Decil para dados tabulados========
    decil9 = (((class_d9-100)*del+minf)-del/2)
    decil9 = decil9 + (9*n_index/10 - acum90)/(hist(class_d9)*n_index)
    
    write(2, *) '&'
    write(2, *) quart1, hist(class_q1)
    write(2, *) '&'
    write(2, *) quart3, hist(class_q3)
    write(2, *) '&'
    write(2, *) decil1, hist(class_d1)
    write(2, *) '&'
    write(2, *) decil9, hist(class_d9)
    
    
    write(*, '(a13, f10.3, a22)')"* Mode     = ", moda, "                   *"
    write(*, '(a13, f10.3, a22)')"* Median   = ", mediana, "                   *"
    write(*, '(a45)') "*                                           *"
    write(*, '(a20, f10.3, a15)')"* Lower Quartile  = ", quart1, "          *" 
    write(*, '(a20, f10.3, a15)')"* Upper Quartile  = ", quart3, "          *"  
    write(*, '(a20, f10.3, a15)')"* Lower Decile    = ", decil1, "          *"
    write(*, '(a20, f10.3, a15)')"* Upper Decile    = ", decil9, "          *"
    write(*, '(a45)') "*                                           *"
    write(*, '(a30, f10.3, a5)')"* Coefficient of Variation  = ", desv/aver, "    *"
    write(*, '(a45)')"*===================================        *"
    
    if (desv*desv/aver<1) write(*, '(a45)') "* Under-Dispersed Distribution              *"
    if (desv*desv/aver>1) write(*, '(a45)') "* Over-Dispersed Distribution               *"
    
    write(*, '(a45)') "*===================================        *"
    write(*, '(a45)') "*                                           *"
    write(*, '(a14, f10.3, a21)') "* Skewness  = ", skew, "                    *" 
    write(*, '(a45)') "*===================================        *"
    
    if (skew<0) write(*, '(a45)') "* Negative Skew                             *"
    if (skew==0) write(*, '(a45)') "* Symmetrical Distribution                  *"
    if (skew>0) write(*, '(a45)') "* Positive Skew                             *"
    
    write(*, '(a45)') "*===================================        *"
    write(*, '(a45)') "*                                           *"
    write(*, '(a21, f10.3, a14)')"* Excess Kurtosis  = ", kurt-3, "             *"
    write(*, '(a45)') "*===================================        *"
    
    if (kurt-3<0) write(*, '(a45)') "* Platykurtic Distribution                  *"
    if (kurt-3==0) write(*, '(a45)') "* Mesokurtic Distribution                   *"
    if (kurt-3>0) write(*, '(a45)') "* Leptokurtic Distribution                  *"
    
    write(*, '(a45)') "*===================================        *"
    write(*, '(a45)') "*                                           *"
    write(*, '(a45)') "****************** SUMMARY ******************"
    
  end subroutine calc_stat_all
  !====================================================================

  subroutine calc_acf(n_index, aver, desv, func)

    integer :: i, j, n_index
    double precision :: acf, aver, desv
    real, dimension(10000000) ::func
    
    write(*, *)
    write(*, *) "Calculating Autocorrelation Function ......"
    
    do j=1,  int(n_index/2)
       
       acf = 0
       do i=1, int(n_index/2)
          
          acf = acf + (func(i)-aver)*(func(i+j)-aver)/(int(n_index/2)*desv*desv)
          
       end do
       
       write(3, *) j, acf
       
    end do
    
  end subroutine calc_acf
  !====================================================================

  subroutine calc_gauss(n_grid, r_xpmg, r_xpmh, grid, hist_g, hist_h, maxg, ming, maxh, minh, dx, dy)

    use types
    
    integer :: bini_g, bini_h, n_grid, i, j
    
    real :: maxg, ming, maxh, minh, dx, dy  
    real, dimension(n_grid, n_grid) :: r_xpmg, r_xpmh
    real :: aveg, ave2g, aveh, ave2h, desvg, desvh
    real, dimension(1000) :: hist_g, hist_h
    
    double precision :: fx, fy, fxx, fyy, fxy, kg, hm
    type(vet2), dimension(1001,1001) :: grid
    
    !==================================================================
    ! Cálculo das curvaturas gaussiana e média
    ! Pela metodologia implementada n_grid deve ser PAR!!!!
    ! Pois estamos pulando de 2 em 2
    
    aveg = 0
    ave2g = 0
    aveh = 0
    ave2h = 0
    
    do i=2, n_grid, 2
       
       do j=2, n_grid, 2
          
          fx = (grid(i+1,j)%z-grid(i-1,j)%z)/(2*dx)
          fy = (grid(i,j+1)%z-grid(i,j-1)%z)/(2*dy)
          
          fxx = (grid(i+1,j)%z-2*grid(i,j)%z+grid(i-1,j)%z)/(dx*dx)
          fyy = (grid(i,j+1)%z-2*grid(i,j)%z+grid(i,j-1)%z)/(dy*dy)
          fxy = (grid(i+1,j+1)%z-grid(i-1,j+1)%z-grid(i+1,j-1)%z+grid(i-1,j-1)%z)/(4*dy*dx)
          
          kg = (fxx*fyy-fxy*fxy)/((1+fx*fx+fy*fy)**2)
          hm = ((1+fy*fy)*fxx-fx*fy*fxy+(1+fx*fx)*fyy)/(2*(1+fx*fx+fy*fy)*sqrt(1+fx*fx+fy*fy))
          
          r_xpmg(i, j) = r_xpmg(i, j) + kg
          r_xpmh(i, j) = r_xpmh(i, j) + hm
          
          bini_g = nint(kg/0.0005) + 500 ! parametros ad hoc
          hist_g(bini_g) = hist_g(bini_g) + 1
          bini_h = nint(hm/0.0005) + 500
          hist_h(bini_h) = hist_h(bini_h) + 1
          aveg = aveg + kg/(n_grid*n_grid/4)
          ave2g = ave2g + kg*kg/(n_grid*n_grid/4)
          aveh = aveh + hm/(n_grid*n_grid/4)
          ave2h = ave2h + hm*hm/(n_grid*n_grid/4)
          
       end do
       
    end do
    
    desvg = sqrt(ave2g - aveg*aveg)
    desvh = sqrt(ave2h - aveh*aveh)
    
    maxg = aveg + 2*desvg
    ming = aveg - 2*desvg
    maxh = aveh + 2*desvh
    minh = aveh - 2*desvh
    
  end subroutine calc_gauss
  !====================================================================
  
  subroutine filter_suave(n_index, paramet, func, dt, inicial)

    real, dimension(3000000) ::func
    integer :: n_index, i, j
    real :: paramet, aux, sum, dt, inicial
    
     do i=1, n_index

        aux = 0
        sum = 0

        do j=1, n_index

           aux = aux + func(j)*exp(-paramet*paramet*(j-i)*(j-i))
           sum = sum + exp(-paramet*paramet*(j-i)*(j-i))

        end do

        aux = aux/sum
        write(4, *) (i-1)*dt+inicial, aux

     end do

   end subroutine filter_suave
   !====================================================================

   subroutine filter_ft(n_index, dt, dw, inicial, func, re, im)

     real, parameter :: pi = 3.141592654
     integer :: n_index, i, j, ierr
     real :: dt, dw, t, w, inicial
     real, dimension(3000000) ::func, re, im

     
     dw = (2*pi/dt)/(n_index) !utilizando 3 vezes mais pontos que na dimensão temporal
     
     do i=1, n_index

        w = -pi/dt + (i-1)*dw ! com isso w vai de 0 até 2pi/dt
        
        re(i) = 0
        im(i) = 0

        do j=1, n_index

           t = inicial + (j-1)*dt

           re(i) = re(i) + func(j)*cos(w*t)
           im(i) = im(i) - func(j)*sin(w*t)

        end do

     end do

     do i=1, n_index
        
	w = -pi/dt + (i-1)*dw
 
        write(2, *, iostat=ierr) w, sqrt(re(i)*re(i) + im(i)*im(i))
        
        if(ierr>0) then
           
           write(*, *)
           write(*, *) ' Problem by writing output '
           write(*, *)
           stop
           
        end if
        
     end do
          
   end subroutine filter_ft
   !====================================================================

   real function sinal(re)
     
     real :: re
     
     sinal = 1
     
     if (re<0) then
        
        sinal = -1
        
     end if
     
   end function sinal
   !====================================================================
   
   subroutine filter_ift(n_index, dt, inicial, func, re, im, ref, imf)
     
     real, parameter :: pi = 3.141592654
     integer :: n_index, i, j, ierr
     real :: dt, dw, t, w, inicial
     real, dimension(3000000) ::func, re, im, ref, imf

     dw = (2*pi/dt)/(n_index)
     
     do j=1, n_index
        
        t = inicial + (j-1)*dt

        ref(j) = 0
        imf(j) = 0

        do i=1, n_index

           w =  -pi/dt + (i-1)*dw ! com isso w vai de -2pi/dt até 2pi/dt
           
           ref(j) = ref(j) + re(i)*cos(w*t) - im(i)*sin(w*t)
           imf(j) = imf(j) + re(i)*sin(w*t) + im(i)*cos(w*t)

        end do

     end do

     do j=1, n_index
        
	t = inicial + (j-1)*dt
        write(3, *, iostat=ierr) t, sinal(ref(j))*sqrt(ref(j)*ref(j) + imf(j)*imf(j))/(n_index)
 
        if(ierr>0) then

           write(*, *)
           write(*, *) 'Problem by writing output '
           write(*, *)
           stop

        end if

     end do
     
   end subroutine filter_ift
   !==================================================================== 
   
   subroutine ordem(W, ORD)
     
     integer :: i
     double precision :: W(3), ORD(3)
     
     ORD(1) = maxval(W)
     ORD(3) = minval(W)
     
     do i=1, 3
        
        if (W(i) < ORD(1)) then
           
           if (W(i) > ORD(3)) then
              
              ORD(2) = W(i)

           end if
           
        end if
        
     end do
     
   end subroutine ordem
   !====================================================================
   

   



 end module funcproc
