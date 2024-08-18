  !===============================================================!
  !Arquivo contendo todas as variaveis utilizadas nos programas   !
  !===============================================================!

module variables

  use types

#ifdef INDEX

  character(len=7), parameter :: version = "2.24.07"
  logical :: ex, res, sph, igual, join, gro, bound, back

  integer :: i, j, ierr, res_i, res_j, id(500)
  integer :: n_index, i_atom, n_sph, aux2, n_index2
  integer, dimension(15)  :: ax

  real :: cent_x, cent_y, cent_z, rho, dist_z, aux3

  character(len=30) :: coord, cond, index, sai, atom, resid, gro_file
  character(len=30), dimension(30) :: get
  character(len=30) :: atomid(500,500), residue(500), aux

  type(vet1) :: buff 
  type(vet1), dimension(1000000) :: ind_store, ind_store2, sphstore

#elif CART

  real, parameter :: pi = 3.141592654
  character(len=7), parameter :: version = "2.24.07"

  logical :: ex, bin, outer, p_grid, rmsd, l_coarse, begin, end, skip
  logical :: eval_skip, lipid, slices, inside, map, back, range, next_frame

  integer :: i, j, k, l, ierr, n_grid, frame, aux, bin_out, noi1, noi2
  integer :: n_index, n_lipid, num, num2, i_atom, bin_coarse, a, b
  integer :: fr_in, fr_end, n_skip, tot_frame, inter, a_dens, tot_dens
  integer :: start, finish, clock_rate, clock_max, div, bini
  integer, dimension(500000) :: in_num, in_num2, in_dens

  real :: dist_x, dist_y, dist_z, s_soma, hour, minu, sec, noir, minv, del
  real :: x_max, x_min, y_max, y_min, dx, dy, dist, s_grid, al, peso, area 
  real :: s_area, la, lb, lc, r_fit, aux2, gridx, gridy, tempo, rough, gz, n_inside
  real :: z_max, z_min,  maxv, hist(1000)
  real, dimension(: , :), allocatable :: r_xpm1, r_xpm2

  character(len=30) :: coord,  ind, atom, ind2, ind3, name_in, traj_type
  character(len=30), dimension(30) :: get
  character(len=1), dimension(:,:), allocatable :: xpm

  type(vet1), dimension(50000) :: store, store2, coarse, coarse2
  type(vet1) :: buff
  type(vet1), dimension(500000) :: out
  type(vet2), dimension(500000) :: dens
  type(vet2), dimension(1001,1001) :: grid, grid2, grid3

#elif SPHE

  real, parameter :: pi = 3.141592654
  character(len=7), parameter :: version = "2.24.07"
  
  logical :: ex, bin, p_grid, rmsd, l_coarse, begin, end, skip
  logical :: eval_skip, lipid, help, outer, back, range, slices, next_frame
  
  integer :: i, j, ierr, n_grid, frame, aux, noi1, noi2
  integer :: n_index, n_lipid, num, num2, i_atom, bin_coarse, a, b
  integer :: start, finish, clock_rate, clock_max, a_dens, div
  integer :: fr_in, fr_end, n_skip
  integer :: in_dens(1000000)
  integer, dimension(500000) :: in_num, in_num2

  real :: dist_x, dist_y, dist_z, s_soma, hour, minu, sec, noir
  real :: dist, s_grid, r_med1, r_med2, girat, girat2, aux2, deq, deq2
  real :: s_area, s_area2, la, lb, lc, r_fit, al, rough, c_angle
  real :: cent_x, cent_y, cent_z, dph, dth, s_vol, s_vol2
  real :: maxv, minv, hist(1000)
  real, dimension(: , :), allocatable :: r_xpm1, r_xpm2
  
  character(len=30) :: coord,  ind, atom, ind2, ind3, traj_type
  character(len=30), dimension(20) :: get
  character(len=1), dimension(:,:), allocatable :: xpm1, xpm2
  
  type(vet3), dimension(50000) :: store, store2, coarse, coarse2
  type(vet1) :: buff
  type(vet1), dimension(500000) :: spher
  type(vet1), dimension(1000000) :: dens_spher
  type(vet3), dimension(1000000) :: dens
  type(vet3), dimension(1001,1001) :: grid, grid2
  type(vet2), dimension(1001,1001) :: grid3
  type(vet2) :: center  

#elif STAT

  real, parameter :: pi = 3.141592654
  character(len=7), parameter :: version = "2.24.07"
  
  logical :: ex, help, back
  
  character(len=30) :: get(20), signal
  
  integer :: n_index, i, j, ierr, bini, class_modal
  integer :: class_med, class_q1, class_q3, class_d1, class_d9
  integer :: start, finish, clock_rate, clock_max
  
  double precision :: aux, aux2, aver, aver2, desv, kurt, delta1, delta2
  double precision :: del, hist(1000), sum, acum, moda, mediana, st_mom
  double precision :: acum25, acum75, quart1, quart3, acum10, acum90
  double precision :: decil1, decil9, skew, acf, maxf, minf
  
  real, dimension(10000000) ::func
  
#endif
  
end module variables

