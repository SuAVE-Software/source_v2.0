program suave_thick
  
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
       n_grid, bin_out, fr_in, fr_end, n_skip, n_lipid, get, div, 's_thick   ', version)
  
  !Leitura do PDB
12 format(a6, i5.1, a5, a5, a1, i4.1, a4, 3f8.3)

  if (.not.outer)then

     if (ind2=='miss')then

        write(*, *)
        write(*, *)'Second index file is missing'
        write(*, *)
        stop

     end if

  end if

  if (.not.lipid) then

     write(*, *)
     write(*, *) 'Number of lipids is missing'
     write(*, *)
     stop

  end if

  !==Lendo os arquivos de index=========================

  if (.not.outer)then

     call abre_ndx(ind, in_num, n_index)
     call abre_ndx(ind2, in_num2, n_index)

  end if

  !=================definindo frames para inicio e fim========
  
  call def_frame(frame, fr_in, fr_end, skip, n_skip, end, begin)

  !=================calculando o espaçamento==================
  
  call def_bin(outer, bin_coarse, n_index, bin, n_grid, bin_out)
  
  !==================== abrindo trajetórias =====================
  
  traj_type = coord(len(trim(coord))-2:len(trim(coord)))
  
  if (traj_type=='xtc') call trj%open(coord)
  if (traj_type=='pdb') call abre_trj(1, coord)
  
  !==================== abrindo outputs =========================

  call abre('thick     ', 2, 'xvg', back)
  
  write(2, '(a7, a7)') "#SuAVE ", version
  write(2, '(a14)') '#Command Line:'
  write(2, '(a10)', advance='no') '#s_thick  '
  write(2, *) (trim(get(i)),"  ", i=1, 30)
  write(2, *) '@    title "Thickness X Frame"'
  write(2, *) '@    xaxis  label "#Frame"'
  write(2, *) '@    yaxis  label "Average Thickness [nm]"'
  
  if (rmsd)then

     call abre('rmsd      ', 3, 'xvg', back)
     
     write(3, '(a7, a7)') "#SuAVE ", version
     write(3, '(a14)') '#Command Line:'
     write(3, '(a10)', advance='no') '#s_thick  '
     write(3, *) (trim(get(i)),"  ", i=1, 30)
     write(3, *) '@    title "RMSD X Frame"'
     write(3, *) '@    xaxis  label "#Frame"'
     write(3, *) '@    yaxis  label "RMSD [nm]"'

  end if

  call abre('volume    ', 4, 'xvg', back)
  
  write(4, '(a7, a7)') "#SuAVE ", version
  write(4, '(a14)') '#Command Line:'
  write(4, '(a10)', advance='no') '#s_thick  '
  write(4, *) (trim(get(i)),"  ", i=1, 30)
  write(4, *) '@    title "Volume X Frame"'
  write(4, *) '@    xaxis  label "#Frame"'
  write(4, *) '@    yaxis  label "Volume per Lipid [nm\S3\N]"'

  
  !Alocando XPM ======================================================
  !===================================================================

  allocate (xpm(n_grid, n_grid), stat = ierr)
  if (ierr /= 0) stop 'Not enough memory to initialize xpm matrix'

  allocate (r_xpm1(n_grid, n_grid), stat = ierr)
  if (ierr /= 0) stop 'Not enough memory to initialize xpm matrix'

  do i=1, n_grid

     do j=1, n_grid

        r_xpm1(i,j) = 0
        
     end do

  end do
  
  !===================================================================
  
  call system_clock(start, clock_rate, clock_max)

! Leitura do arquivo .pdb

  minv = 100000
  maxv = -100000
  i_atom = 0
  n_index = 1
  frame = 0
  ierr = 0
  num = 1
  num2 = 1
  x_min = 1000
  y_min = 1000
  x_max = 0
  y_max = 0
  tot_frame = 0

  do while (ierr >= 0)

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
              
              i_atom  = i_atom + 1
              
              buff%n_atom = i_atom
              
              if (i_atom == in_num(num)) then
                 
                 store(num) = buff
                 n_index = n_index + 1
                 num = num + 1
                 
                 x_max = max(x_max, buff%x)
                 x_min = min(x_min, buff%x)
                 y_max = max(y_max, buff%y)
                 y_min = min(y_min, buff%y)
                 
              end if
              
              if (i_atom == in_num2(num2)) then
                 
                 store2(num2) = buff
                 n_index = n_index + 1
                 num2 = num2 + 1
                 
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
           
           if (i_atom == in_num(num)) then
              
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

           if (i_atom == in_num2(num2)) then

              copy = trj%x(1, i_atom)
              store2(num2)%x = copy(1)*10
              store2(num2)%y = copy(2)*10
              store2(num2)%z = copy(3)*10
              n_index = n_index + 1
              num2 = num2 + 1

              x_max = max(x_max, store2(num2-1)%x)
              x_min = min(x_min, store2(num2-1)%x)
              y_max = max(y_max, store2(num2-1)%y)
              y_min = min(y_min, store2(num2-1)%y)

           end if

        end do

     end if

     ! Fitting process =============================================================
     
     eval_skip = ((mod((frame-fr_in),(n_skip+1))==0))

     if ((frame>fr_in-1).and.(frame<fr_end+1).and.(eval_skip).and.(ierr>=0)) then

        tot_frame = tot_frame + 1 ! para programas que calculam propriedades medias
        
        noi1 = num - 1
        noi2 = num2 - 1
        
        if (l_coarse)then
           
           ! Estruturação do primeiro coarse grid
           
           dx = (x_max - x_min)/bin_coarse
           dy = (y_max - y_min)/bin_coarse
           
           ! adaptação ao r_fit e ao alfa
           
           call param(x_max, x_min, y_max, y_min, num, r_fit, al, rough)
           
           !$OMP parallel do private(s_grid, dist, peso, k, aux, j)
           do i=1, bin_coarse + 1
              
              do j=1, bin_coarse + 1
                 
                 aux = (i-1)*(bin_coarse + 1) + j
                 coarse(aux)%x = (i-1)*dx + x_min
                 coarse(aux)%y = (j-1)*dy + y_min
                 
                 s_grid = 0
                 coarse(aux)%z = 0
                 
                 do k=1, num - 1
                    
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
           !OMP end parallel do
           
           num = (bin_coarse + 1)**2 + 1
           
           ! Estruturação do segundo coarse grid
           
           call param(x_max, x_min, y_max, y_min, num2, r_fit, al, rough)
           
           !$OMP parallel do private(s_grid, dist, peso, k, aux, j)
           do i=1, bin_coarse + 1
              
              do j=1, bin_coarse + 1
                 
                 aux = (i-1)*(bin_coarse + 1) + j
                 coarse2(aux)%x = (i-1)*dx + x_min
                 coarse2(aux)%y = (j-1)*dy + y_min
                 
                 s_grid = 0
                 coarse2(aux)%z = 0
                 
                 do k=1, num2 - 1
                    
                    dist = (coarse2(aux)%x - store2(k)%x)**2
                    dist = dist + (coarse2(aux)%y - store2(k)%y)**2
                    
                    if (dist<r_fit*r_fit)then
                       
                       peso = exp(-(dist*al*al)/pi)
                       coarse2(aux)%z = coarse2(aux)%z + peso*store2(k)%z
                       s_grid = s_grid + peso
                       
                    end if
                    
                 end do
                 
                 coarse2(aux)%z = coarse2(aux)%z/s_grid
                 
              end do
              
           end do
           !$OMP end parallel do
           
           num2 = (bin_coarse + 1)**2 + 1
           
        else
           
           coarse = store
           coarse2 = store2
           
        end if
        
        !Estruturação do primeiro grid de alta resolução
        
        dx = (x_max - x_min)/n_grid
        dy = (y_max - y_min)/n_grid
        
        call param(x_max, x_min, y_max, y_min, num, r_fit, al, rough)
        
        !$OMP parallel do private(s_grid, dist, peso, k, j)              
        do i=1, n_grid+1
           
           do j=1, n_grid+1
              
              grid(i,j)%x = (i-1)*dx + x_min
              grid(i,j)%y = (j-1)*dy + y_min
              
              s_grid = 0
              grid(i,j)%z = 0
              
              do k=1, num - 1
                 
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
        !$OMP end parallel do
        
        ! Estruturação do segundo grid de alta resolução
        
        call param(x_max, x_min, y_max, y_min, num2, r_fit, al, rough)
        
        !$OMP parallel do private(s_grid, dist, peso, k, j)
        do i=1, n_grid+1
           
           do j=1, n_grid+1
              
              grid2(i,j)%x = (i-1)*dx + x_min
              grid2(i,j)%y = (j-1)*dy + y_min
              
              s_grid = 0
              grid2(i,j)%z = 0
              
              do k=1, num2 - 1
                 
                 dist = (grid2(i,j)%x - coarse2(k)%x)**2
                 dist = dist + (grid2(i,j)%y - coarse2(k)%y)**2
                 
                 if (dist<r_fit*r_fit)then
                    
                    peso = exp(-(dist*al*al)/pi)
                    s_grid = s_grid + peso
                    grid2(i,j)%z = grid2(i,j)%z + peso*coarse2(k)%z
                    
                 end if
                 
              end do
              
              grid2(i,j)%z = grid2(i,j)%z/s_grid
              
           end do
           
        end do
        !$OMP end parallel do
        
        ! Calculando o volume e o thickness  =================
        
        do i=1, n_grid+1
           
           do j=1, n_grid+1
              
              grid3(i,j)%z = (grid(i,j)%z + grid2(i,j)%z)/2
              grid3(i,j)%x = (i-1)*dx + x_min
              grid3(i,j)%y = (j-1)*dy + y_min
              
           end do
           
        end do
        
        call calc_thick(n_grid, aver, aver2, s_v, dx, dy, grid, grid2, grid3, r_xpm1)
        
        desv = sqrt(aver2 - aver*aver)
        minv = min(aver-2*desv, minv)
        maxv = max(aver+2*desv, maxv)
        
        write(2, *) frame, aver
        write(4, *) frame, s_v/(1000*n_lipid*2)
        
        !calculando o RMSD
        
        if (rmsd)then
           
           write(3, *) calc_rmsd(store, store2, grid, grid2, x_min, y_min, dx, dy, noi1, noi2)/10   
           
        end if
        
        ! Fim do cálculo

        gx = x_min
        gy = y_min
        n_index = 1
        i_atom = 0
        num = 1
        num2 = 1
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
  
  if (rmsd)then
     
     close(3)
     
  end if
  
  close(4)
  close(2)
  
  call system_clock(finish, clock_rate, clock_max)
  
  ! Fim da leitura do arquivo .pdb

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

  call abre('thick     ', 4, 'xpm', back)
  
  write(4, '(a7, a7)') "#SuAVE ", version
  write(4, '(a14)') '#Command Line:'
  write(4, '(a10)', advance='no') '#s_thick  '
  write(4, *) (trim(get(i)),"  ", i=1, 30)
  write(4, *) '/* XPM */'
  write(4, *)'/* This matrix is generated by Thick.f90 */'
  write(4, *)'/* title:   "Thickness" */'
  write(4, *)'/* x-label: "x axis [nm]" */'
  write(4, *)'/* y-label: "y axis [nm]" */'
  write(4, *)'/* type:    "Continuous" */'
  write(4, *)'static char * gv_xpm[] = {'
  write(4, *) '"',n_grid,n_grid,' 7 1",'
  
  do i=1, n_grid
     
     do j=1, n_grid
        
        r_xpm1(i,j) = r_xpm1(i,j)/tot_frame/10

     end do
     
  end do

  del = (maxv - minv)/6

  call print_xpm(n_grid, del, dx, dy, gx, gy, xpm, r_xpm1, minv, 'b')

  ! Final da escrita do XPM
  
  if (p_grid)then
     
     call print_grid_xpm(grid, r_xpm1, n_grid, n_grid, 'grid1', back)
     call print_grid_xpm(grid2, r_xpm1, n_grid, n_grid, 'grid2', back)  
     
  end if
  
  call ending(back, finish, start, clock_rate) ! Finaliza programa e mostra tempo de processamento

  close(4)

  deallocate(r_xpm1)
  deallocate(xpm)
    
end program suave_thick
