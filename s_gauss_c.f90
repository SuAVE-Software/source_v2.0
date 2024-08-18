program suave_gauss
  
  use types
  use variables
  use funcproc

  use gmxfort_trajectory
  use gmxfort_utils
 
  !============================
  type(Trajectory) :: trj
  real :: copy(3)
  !============================

  real :: maxg, ming, maxh, minh, Lx, Ly
  real, dimension( : , : ), allocatable :: r_xpmg, r_xpmh
  real, dimension(1000) :: hist_g, hist_h
  
  ! ==================================
  
  call startup(outer, bin, p_grid, coord, ind, ind2, rmsd, map, ind3, &
     l_coarse, begin, end, skip, lipid, rough, slices, inside, range, &
     n_grid, bin_out, fr_in, fr_end, n_skip, n_lipid, get, div, 's_gauss   ', version)

  !Leitura do PDB
12 format(a6, i5.1, a5, a5, a1, i4.1, a4, 3f8.3)

  !escrita do PDB com XPM
!13 format(a6, i5.1, a5, a5, a1, i4.1, a4, 3f8.3, f5.1)

  !==Lendo os arquivos de index=========================
  if (.not.outer)then
     
     call abre_ndx(ind, in_num, n_index)
     
  end if

  !=================definindo frames para inicio e fim========

  call def_frame(frame, fr_in, fr_end, skip, n_skip, end, begin)
  
  !=================calculando o espaçamento==================

  call def_bin(outer, bin_coarse, n_index, bin, n_grid, bin_out)
  
  !=================== zerando memoria do XPM====================
  !============== e outras variáveis ============================
  
  do i=1, 1000
  
     hist_g(i) = 0
     hist_h(i) = 0

  end do

  maxg = -10000000
  ming = 10000000
  maxh = -100000000
  minh = 10000000

  !==================== abrindo trajetórias =====================
 
   traj_type = coord(len(trim(coord))-2:len(trim(coord)))

   if (traj_type=='xtc') call trj%open(coord)
   if (traj_type=='pdb') call abre_trj(1, coord)
 
   !==================== abrindo outputs =========================

  if (rmsd)then

     call abre('rmsd      ', 3, 'xvg', back)

     write(3, '(a7, a7)') "#SuAVE ", version
     write(3, '(a14)') '#Command Line:'
     write(3, '(a10)', advance='no') '#s_gauss  '
     write(3, *) (trim(get(i)),"  ", i=1, 30)
     write(3, *) '@    title "RMSD X Frame"'
     write(3, *) '@    xaxis  label "#Frame"'
     write(3, *) '@    yaxis  label "RMSD [nm]"'

  end if

  !Alocando XPM ======================================================
  !===================================================================
  allocate (xpm(n_grid, n_grid), stat = ierr)

  allocate (r_xpmg(n_grid, n_grid), stat = ierr)
  allocate (r_xpmh(n_grid, n_grid), stat = ierr)

  if (ierr /= 0) stop 'Not enough memory to initialize xpm matrix'
  !===================================================================
  
  do i=1, n_grid

     do j=1, n_grid

        r_xpmg(i,j) = 0
        r_xpmh(i,j) = 0

     end do

  end do
  
  call system_clock(start, clock_rate, clock_max)

! Leitura do arquivo .gro
  
  i_atom = 0
  n_index = 1
  frame = 0
  ierr = 0
  num = 1
  x_min = 1000
  y_min = 1000
  x_max = 0
  y_max = 0
  tot_frame = 0

  do while (ierr>=0)

     if (traj_type=='pdb') then

        next_frame = .false.

        do while (.not.next_frame)

           read(1, 12, iostat=ierr) atom, buff%n_atom, buff%atom, &
                buff%resid, buff%ident, buff%n_resid, buff%code, buff%x, &
                buff%y, buff%z

           if (ierr<0) exit

           if ((atom.eq.'ATOM  '))then

              if(ierr > 0) then

                 write(*, *)
                 write(*, *) 'Problem reading atomic position!'
                 write(*, *)
                 stop

              endif

              i_atom = i_atom + 1

              buff%n_atom = i_atom

              if (outer)then

                 out(i_atom) = buff
                 n_index = n_index + 1

                 x_max = max(x_max, buff%x)
                 x_min = min(x_min, buff%x)
                 y_max = max(y_max, buff%y)
                 y_min = min(y_min, buff%y)
           
              else

                 if (i_atom == in_num(n_index)) then

                    store(n_index) = buff
                    n_index = n_index + 1

                    x_max = max(x_max, buff%x)
                    x_min = min(x_min, buff%x)
                    y_max = max(y_max, buff%y)
                    y_min = min(y_min, buff%y)
              
                 end if

              end if
        
           end if

           next_frame = ((atom.ne.'ATOM  ').and.(n_index > 1))

        end do

        frame = frame + 1

     end if
     ! ===================================

     if (traj_type=='xtc') then

        frame = frame + trj%read_next(1)

        do while (i_atom < trj%natoms())

           i_atom  = i_atom + 1

           if (outer) then
              copy = trj%x(1, i_atom)
              store(num)%x = copy(1)*10
              store(num)%y = copy(2)*10
              store(num)%z = copy(3)*10
              n_index = n_index + 1
              num = num + 1

              x_max = max(x_max, store(num-1)%x)
              x_min = min(x_min, store(num-1)%x)
              y_max = max(y_max, store(num-1)%y)
              y_min = min(y_min, store(num-1)%y)

           else

              if (i_atom == in_num(n_index)) then

                 copy = trj%x(1, i_atom)
                 store(num)%x = copy(1)*10
                 store(num)%y = copy(2)*10
                 store(num)%z = copy(3)*10
                 n_index = n_index + 1
                 num = num + 1

                 x_max = max(x_max, store(num-1)%x)
                 x_min = min(x_min, store(num-1)%x)
                 y_max = max(y_max, store(num-1)%y)
                 y_min = min(y_min, store(num-1)%y)

              end if

           end if

        end do

     end if

     ! Fitting process =============================================================                                                                                                                        

     eval_skip = ((mod((frame-fr_in),(n_skip+1))==0))

     if ((frame>fr_in-1).and.(frame<fr_end+1).and.(eval_skip).and.(ierr>=0)) then

        tot_frame = tot_frame + 1 ! para programas que calculam propriedades medias
              
        noi1 = n_index - 1
              
        if (l_coarse) then
                 
           ! Estruturação do grid
                 
           dx = (x_max - x_min)/bin_coarse
           dy = (y_max - y_min)/bin_coarse
                 
           ! adaptação dos valores de r_fit e alfa

           call param(x_max, x_min, y_max, y_min, n_index, r_fit, al, rough)

           !$OMP parallel do private(s_grid, dist, peso, k, aux, j)
           do i=1, bin_coarse + 1
                    
              do j=1, bin_coarse + 1
                       
                 aux = (i-1)*(bin_coarse + 1) + j
                 coarse(aux)%x = (i-1)*dx + x_min
                 coarse(aux)%y = (j-1)*dy + y_min
                       
                 s_grid = 0
                 coarse(aux)%z = 0
                       
                 do k=1, n_index - 1

                    dist = (coarse(aux)%x - store(k)%x)**2
                    dist = dist + (coarse(aux)%y - store(k)%y)**2

                    if (dist<r_fit*r_fit)then

                       peso = exp(-(dist*al*al)/pi)
                       s_grid = s_grid + peso
                       coarse(aux)%z = coarse(aux)%z + peso*store(k)%z

                    end if

                 end do
                       
                 coarse(aux)%z = coarse(aux)%z/s_grid

              end do
                    
           end do
           !$OMP end parallel do
                 
           n_index = (bin_coarse + 1)**2 + 1
                 
        else
                 
           coarse = store
                 
        end if
              
        ! estruturação do prmeiro grid de alta resolução
        dx = (x_max - x_min)/n_grid
        dy = (y_max - y_min)/n_grid

        call param(x_max, x_min, y_max, y_min, n_index, r_fit, al, rough)

        !$OMP parallel do private(s_grid, dist, peso, k, j)
        do i=1, n_grid+1
                 
           do j=1, n_grid+1
                    
              grid(i,j)%x = (i-1)*dx + x_min
              grid(i,j)%y = (j-1)*dy + y_min
                    
              s_grid = 0
              grid(i,j)%z = 0
                    
              do k=1, n_index - 1

                 dist = (grid(i,j)%x - coarse(k)%x)**2
                 dist = dist + (grid(i,j)%y - coarse(k)%y)**2

                 if (dist<r_fit*r_fit) then

                    peso = exp(-(dist*al*al)/pi)
                    s_grid = s_grid + peso
                    grid(i,j)%z = grid(i,j)%z + peso*coarse(k)%z

                 end if

              end do
                    
              grid(i,j)%z = grid(i,j)%z/s_grid
                    
           end do
                 
        end do
        !$OMP end parallel do
        !=================== Fim da estruturação ==========================

        ! Calculando as curvaturas pelo método das diferenças finitas =====

        call calc_gauss(n_grid, r_xpmg, r_xpmh, grid, hist_g, hist_h, maxg, ming, maxh, minh, dx, dy)
                 
        !==================================================================
              
        if (rmsd)then
                 
           noir = 0
                 
           do i=1, noi1
                    
              a = nint((store(i)%x-x_min)/dx) + 1
              b = nint((store(i)%y-y_min)/dy) + 1
              dist_z = store(i)%z-grid(a,b)%z
              noir = noir + dist_z*dist_z/noi1
                    
           end do
                 
           noir = sqrt(noir)
                 
           write(3, *) noir/10
                 
        end if

        !! ******* Fim do cálculo

        Lx = x_max - x_min
        Ly = y_max - y_min
        gx = x_min
        gy = y_min
        n_index = 1
        num = 1
        i_atom = 0
        x_min = 1000
        y_min = 1000
        x_max = 0
        y_max = 0

     end if !======((frame<fr_in-1).and.(frame>fr_end+1))

     !====garante que fr_end sempre seja maior que frame ===
     !====caso essa variável não tenha sido fixada==========
     if (.not.end) then
              
        fr_end = frame + 1
              
     end if
     !======================================================

     if ((frame>=trj%NFRAMES).and.(traj_type=='xtc')) ierr = -1 ! apenas quando lendo o tipo XTC

  end do

  call system_clock (finish, clock_rate, clock_max)

  ! Fim da leitura do arquivo .pdb                                                                                                                                                                          

  if (traj_type=='xtc') call trj%close()
  if (traj_type=='pdb') close(1)
  
  if (rmsd) then
     
     close(3)

  end if
  
  ! escrita do arquivo XPMG========================================

  call abre('gaussian  ', 4, 'xpm', back)
  
  write(4, '(a7, a7)') "#SuAVE ", version
  write(4, '(a14)') '#Command Line:'
  write(4, '(a9)', advance='no') '#s_gauss  '
  write(4, *) (trim(get(i)),"  ", i=1, 30)
  write(4, *) '/* XPM */'
  write(4, *)'/* This matrix is generated by s_gauss.f90 */'
  write(4, *)'/* title:   "Gaussian Curvature" */'
  write(4, *)'/* x-label: "x axis [nm]" */'
  write(4, *)'/* y-label: "y axis [nm]" */'
  write(4, *)'/* type:    "Continuous" */'
  write(4, *)'static char * gv_xpm[] = {'
  write(4, *) '"',n_grid/2,n_grid/2,' 7 1",'

  do i=2, n_grid, 2
     
     do j=2, n_grid, 2

        r_xpmg(i,j) = r_xpmg(i,j)/real(tot_frame)

     end do
     
  end do
  
  del = (maxg - ming)/6

  call print_xpm_gauss(n_grid, del, dx, dy, gx, gy, xpm, r_xpmg, ming, 'r')

  close(4)
  
  ! Final da escrita do XPMG ========================================

  ! escrita do arquivo XPMH ========================================

  call abre('mean      ', 4, 'xpm', back)
  
  write(4, '(a7, a7)') "#SuAVE ", version
  write(4, '(a14)') '#Command Line:'
  write(4, '(a9)', advance='no') '#s_gauss  '
  write(4, *) (trim(get(i)),"  ", i=1, 30)
  write(4, *) '/* XPM */'
  write(4, *)'/* This matrix is generated by s_gauss.f90 */'
  write(4, *)'/* title:   "Mean Curvature" */'
  write(4, *)'/* x-label: "x axis [nm]" */'
  write(4, *)'/* y-label: "y axis [nm]" */'
  write(4, *)'/* type:    "Continuous" */'
  write(4, *)'static char * gv_xpm[] = {'
  write(4, *) '"',n_grid/2,n_grid/2,' 7 1",'

  do i=2, n_grid, 2
     
     do j=2, n_grid, 2
        
        r_xpmh(i,j) = r_xpmh(i,j)/real(tot_frame)
                
     end do
     
  end do

  del = (maxh - minh)/6
  
  call print_xpm_gauss(n_grid, del, dx, dy, gx, gy, xpm, r_xpmh, minh, 'r')
  
  close(4)
  
  ! Final da escrita do XPMH ================================
  
  ! Escrita das distribuições de probabilidade ==============

  call abre('hist_gauss', 5, 'xvg', back)
  
  write(5, '(a7, a7)') "#SuAVE ", version
  write(5, '(a14)') '#Command Line:'
  write(5, '(a9)', advance='no') '#s_gauss  '
  write(5, *) (trim(get(i)),"  ", i=1, 30)
  write(5, *) '@    title "Gaussian Curvature Distribution"'
  write(5, *) '@    xaxis  label "Gaussian Curvature"'
  write(5, *) '@    yaxis  label "PDF"'

  do i=200, 800

     write(5, *) ((i-500)*0.0005), hist_g(i)/(tot_frame*0.0005*(n_grid*n_grid/4))

  end do

  close(5)
  
  ! Mean Curvature ==========================================

  call abre('hist_mean ', 5, 'xvg', back)
  
  write(5, '(a7, a7)') "#SuAVE ", version
  write(5, '(a14)') '#Command Line:'
  write(5, '(a9)', advance='no') '#s_gauss  '
  write(5, *) (trim(get(i)),"  ", i=1, 30)
  write(5, *) '@    title "Mean Curvature Distribution"'
  write(5, *) '@    xaxis  label "Mean Curvature"'
  write(5, *) '@    yaxis  label "PDF"'

  do i=1, 1000

     write(5, *) ((i-500)*0.0005), hist_h(i)/(tot_frame*0.0005*(n_grid*n_grid/4))

  end do

  close(5)
  
  ! Fim da escrita ==========================================
  
  !== escrita do GRID========================================
  if (p_grid) then
     
     do i=2, n_grid, 2
        
        do j=2, n_grid, 2
    
           r_xpmg(i,j) = (r_xpmg(i,j)-ming)/(maxg-ming)
           r_xpmh(i,j) = (r_xpmh(i,j)-minh)/(maxh-minh)
    
        end do

     end do
       
     call print_grid_xpm_xpm(grid, r_xpmg, r_xpmh, n_grid, n_grid, 'grid1', back)
     
  end if
  !== escrita do GRID========================================

  deallocate(r_xpmg)
  deallocate(r_xpmh)
  deallocate(xpm)

  call calcula_bin(aux, Lx, Ly)
  
  write(*, *)
  write(*, *) 'Best results with bin < ', aux

  call ending(back, finish, start, clock_rate) ! Finaliza programa e mostra tempo de processamento
  
end program suave_gauss
