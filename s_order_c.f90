program suave_order
  
  use types
  use variables
  use funcproc

  use gmxfort_trajectory
  use gmxfort_utils
  
  !============================
  type(Trajectory) :: trj
  real :: copy(3)
  !============================
  
  call startup(outer, bin, p_grid, coord, ind, ind2, rmsd, map, ind3, &
     l_coarse, begin, end, skip, lipid, rough, slices, inside, range, &
     n_grid, bin_out, fr_in, fr_end, n_skip, n_lipid, get, div, 's_order   ', version)

  !Leitura do PDB
12 format(a6, i5.1, a5, a5, a1, i4.1, a4, 3f8.3)
    
  !==Lendo os arquivos de index=========================
  
  if (.not.outer)then
     
     call abre_ndx(ind, in_num, n_index)
     
  end if

  !=================definindo frames para inicio e fim========
  
  call def_frame(frame, fr_in, fr_end, skip, n_skip, end, begin)

  !=================calculando o espaçamento==================
  
  call def_bin(outer, bin_coarse, n_index, bin, n_grid, bin_out)

  !==================== abrindo trajetórias =====================
  
  traj_type = coord(len(trim(coord))-2:len(trim(coord)))
  
  if (traj_type=='xtc') call trj%open(coord)
  if (traj_type=='pdb') call abre_trj(1, coord)
    
  !=========abrindo arquivo do parâmetro de ordem médio=========================

  call abre('order_aver', 2,'xvg', back)

  write(2, '(a7, a7)') "#SuAVE ", version
  write(2, '(a14)') '#Command Line:'
  write(2, '(a10)', advance='no') '#s_order  '
  write(2, *) (trim(get(i)),"  ", i=1, 30)
  write(2, *) '@    title "Average Curvature Order Parameter X Frame"'
  write(2, *) '@    xaxis  label "#Frame"'
  write(2, *) '@    yaxis  label "Average Curvature Order Parameter"'

  !=================== alocando memoria do XPM====================

  allocate (xpm(n_grid, n_grid), stat = ierr)
  if (ierr /= 0) stop 'Not enough memory to initialize xpm matrix'

  allocate (r_xpm1(n_grid, n_grid), stat = ierr)
  if (ierr /= 0) stop 'Not enough memory to initialize xpm matrix'

  do i=1, n_grid

     do j=1, n_grid

        r_xpm1(i,j) = 0

     end do

  end do

  do k=1, 100

     hist(k) = 0

  end do

  if (rmsd)then

     call abre('rmsd      ', 3,'xvg', back)     
     
     write(3, '(a7, a7)') "#SuAVE ", version
     write(3, '(a14)') '#Command Line:'
     write(3, '(a10)', advance='no') '#s_order  '
     write(3, *) (trim(get(i)),"  ", i=1, 30)
     write(3, *) '@    title "RMSD X Frame"'
     write(3, *) '@    xaxis  label "#Frame"'
     write(3, *) '@    yaxis  label "RMSD [nm]"'

  end if

  call system_clock(start, clock_rate, clock_max)

! Leitura do arquivo .gro
  
  i_atom = 0
  n_index = 1
  frame = 0
  ierr = 0
  x_min = 1000
  y_min = 1000
  x_max = 0
  y_max = 0
  tot_frame = 0
  frame = 0
  minv = 100000
  maxv = -100000
  
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
              
              if (i_atom == in_num(n_index)) then
                 
                 store(n_index) = buff
                 n_index = n_index + 1
                 
                 x_max = max(x_max, buff%x)
                 x_min = min(x_min, buff%x)
                 y_max = max(y_max, buff%y)
                 y_min = min(y_min, buff%y)
                 
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

           if (i_atom == in_num(n_index)) then

              copy = trj%x(1, i_atom)
              store(n_index)%x = copy(1)*10
              store(n_index)%y = copy(2)*10
              store(n_index)%z = copy(3)*10
              n_index = n_index + 1

              x_max = max(x_max, store(n_index-1)%x)
              x_min = min(x_min, store(n_index-1)%x)
              y_max = max(y_max, store(n_index-1)%y)
              y_min = min(y_min, store(n_index-1)%y)

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
                 
                 if (dist<r_fit*r_fit)then
                    
                    peso = exp(-(dist*al*al)/pi)
                    s_grid = s_grid + peso
                    grid(i,j)%z = grid(i,j)%z + peso*coarse(k)%z
                    
                 end if
                 
              end do
              
              grid(i,j)%z = grid(i,j)%z/s_grid
              
           end do
           
        end do
        !OMP end parallel do
        
        ! Fim da estruturação
        
        ! Cálculo do ângulo de inclinação
        
        call calc_order(n_grid, grid, r_xpm1, aver, aver2, hist)
        
        desv = sqrt(aver2 - aver*aver)
        minv = min(aver-2*desv, minv)
        maxv = max(aver+2*desv, maxv)
        
        ! esses termos extras são por conta da distribuição
        ! não ser simétrica o que causava valores maiores que 1
        ! e menores que -0.5 para os limites do intervalo
        
        minv = max(minv, -0.5)
        maxv = min(maxv, 1.0)
        
        !=============================================
        
        ! n_grid*n_grid é o número de células do grid

        write(2, *) frame, aver
        
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

        gx = x_min
        gy = y_min
        n_index = 1
        i_atom = 0
        x_min = 1000
        y_min = 1000
        x_max = 0
        y_max = 0
        
     end if !======((frame<fr_in-1).and.(frame>fr_end+1))
     
     !====garante que fr_end sempre seja maior que frame ===
     !====caso essa variável não tenha sido fixada==========
     if (.not.end)then
        
        fr_end = frame + 1
        
     end if

     !======================================================

     if ((frame>=trj%NFRAMES).and.(traj_type=='xtc')) ierr = -1 ! apenas quando lendo o tipo XTC
     
  end do
  
  if (traj_type=='xtc') call trj%close()
  if (traj_type=='pdb') close(1)
  
  close(2)
  
  if (rmsd) then
     
     close(3)
     
  end if
  
  call system_clock(finish, clock_rate, clock_max)
  
  ! escrita do arquivo XPM========================================
  if (range) then
     
     write(*, *)
     write(*, *) 'Calculated range = [', minv, ';', maxv, ']'
     write(*, *)
     write(*, '(a19)', advance='no') ' Inferior limit :  '
     read(*, *) minv
     write(*, *)
     write(*, '(a19)', advance='no') ' Superior limit :  '
     read(*, *) maxv

  end if

  call abre('order     ', 4, 'xpm', back)
  
  write(4, '(a7, a7)') "#SuAVE ", version
  write(4, '(a14)') '#Command Line:'
  write(4, '(a10)', advance='no') '#s_order  '
  write(4, *) (trim(get(i)),"  ", i=1, 30)
  write(4, *) '/* XPM */'
  write(4, *)'/* This matrix is generated by s_order.f90 */'
  write(4, *)'/* title:   "Order Parameter" */'
  write(4, *)'/* x-label: "x axis [nm]" */'
  write(4, *)'/* y-label: "y axis [nm]" */'
  write(4, *)'/* type:    "Continuous" */'
  write(4, *)'static char * gv_xpm[] = {'
  write(4, *) '"',n_grid,n_grid,' 7 1",'

  do i=1, n_grid
     
     do j=1, n_grid
        
        r_xpm1(i,j) = r_xpm1(i,j)/tot_frame
                
     end do
     
  end do
  
  del = (maxv - minv)/6

  call print_xpm(n_grid, del, dx, dy, gx, gy, xpm, r_xpm1, minv, 'r')
  
  close(4)
  
  ! Final da escrita do XPM

  !=========escrita do histograma===========================

  call abre('hist      ', 4, 'xvg', back)

  write(4, '(a7, a7)') "#SuAVE ", version
  write(4, '(a14)') '#Command Line:'
  write(4, '(a10)', advance='no') '#s_order  '
  write(4, *) (trim(get(i)),"  ", i=1, 30)
  write(4, *) '@    title "Curvature angles"'
  write(4, *) '@    xaxis  label "Angles [\So\N]"'
  write(4, *) '@    yaxis  label "%"'

  do i=1, 100
     
     graph = hist(i)
     graph = graph/(tot_frame*1*n_grid*n_grid)
     
     write(4, *) i-1, graph*100
     
  end do

  close(4)
  
  !== escrita do GRID========================================

  if (p_grid)then

     call print_grid_xpm(grid, r_xpm1, n_grid, n_grid, 'grid1', back)
     
  end if

  call ending(back, finish, start, clock_rate) ! Finaliza programa e mostra tempo de processamento
  
  deallocate(r_xpm1)
  deallocate(xpm)
  
end program suave_order
