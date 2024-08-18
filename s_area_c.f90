program suave_area
  
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
     n_grid, bin_out, fr_in, fr_end, n_skip, n_lipid, get, div, 's_area    ', version)

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
  
  !=================== alocando memoria do XPM====================

  allocate (r_xpm1(n_grid+1, n_grid+1), stat = ierr)
  if (ierr /= 0) stop 'Not enough memory to initialize xpm matrix'

  allocate (r_xpm2(n_grid+1, n_grid+1), stat = ierr)
  if (ierr /= 0) stop 'Not enough memory to initialize xpm matrix'

  do i=1, n_grid+1

     do j=1, n_grid+1

        r_xpm1(i,j) = 0
        r_xpm2(i,j) = 0

     end do

  end do

  !==================== abrindo trajetórias =====================

  traj_type = coord(len(trim(coord))-2:len(trim(coord)))

  if (traj_type=='xtc') call trj%open(coord)
  if (traj_type=='pdb') call abre_trj(1, coord)

  !==================== abrindo outputs =========================
  
  call abre('area-med  ', 2, 'xvg', back)

  write(2, '(a7, a7)') "#SuAVE ", version
  write(2, '(a14)') '#Command Line:'
  write(2, '(a9)', advance='no') '#s_area  '
  write(2, *) (trim(get(i)),"  ", i=1, 30)
  write(2, *) '@    title "Area per lipid X Frame"'
  write(2, *) '@    xaxis  label "#Frame"'
  write(2, *) '@    yaxis  label "Area per lipid [nm\S2\N]"'

  if (rmsd)then

     call abre('rmsd      ', 3, 'xvg', back)

     write(3, '(a7, a7)') "#SuAVE ", version
     write(3, '(a14)') '#Command Line:'
     write(3, '(a9)', advance='no') '#s_area  '
     write(3, *) (trim(get(i)),"  ", i=1, 30)
     write(3, *) '@    title "RMSD X Frame"'
     write(3, *) '@    xaxis  label "#Frame"'
     write(3, *) '@    yaxis  label "RMSD [nm]"'

  end if

  call system_clock (start, clock_rate, clock_max)

! Leitura da trajetória   =================================

  i_atom = 0
  n_index = 1
  ierr = 0
  num = 1
  num2 = 1
  x_min = 1000
  y_min = 1000
  x_max = 0
  y_max = 0
  tot_frame = 0
  frame = 0
  
  do while (ierr >= 0) ! até o final da trajetória 

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
              
              if (outer)then
                 
                 out(i_atom) = buff
                 n_index = n_index + 1
                 
                 x_max = max(x_max, buff%x)
                 x_min = min(x_min, buff%x)
                 y_max = max(y_max, buff%y)
                 y_min = min(y_min, buff%y)
                 
              else
                 
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
        
        if (l_coarse) then
           
           ! Estruturação do primeiro coarse grid
           
           dx = (x_max - x_min)/bin_coarse
           dy = (y_max - y_min)/bin_coarse
           
           !obtendo parametros 
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
           !$OMP end parallel do
           
           ! atualização do num (número de pontos do coarse grid)

           num = (bin_coarse + 1)**2 + 1
           
           ! Estruturação do segundo coarse grid
           ! obtendo parametros
           
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
        
        ! estruturação do prmeiro grid de alta resolução
        
        dx = (x_max - x_min)/n_grid
        dy = (y_max - y_min)/n_grid
        
        call param(x_max, x_min, y_max, y_min, num, r_fit, al, rough)
        
        !$OMP parallel do private(s_grid, dist, peso, k, j) reduction(min:minv)
        do i=1, n_grid+1
           
           do j=1, n_grid+1
              
              grid(i,j)%x = (i-1)*dx + x_min
              grid(i,j)%y = (j-1)*dy + y_min
              
              s_grid = 0
              grid(i,j)%z = 0
              minv = 1000
              
              do k=1, num - 1 
                 
                 dist = (grid(i,j)%x - coarse(k)%x)**2
                 dist = dist + (grid(i,j)%y - coarse(k)%y)**2
                 
                 minv = min(minv, dist)
                 
                 if (dist<r_fit*r_fit)then
                    
                    peso = exp(-(dist*al*al)/pi)
                    s_grid = s_grid + peso
                    grid(i,j)%z = grid(i,j)%z + peso*coarse(k)%z
                    
                 end if
                 
              end do
              
              grid(i,j)%z = grid(i,j)%z/s_grid
              r_xpm1(i,j) = r_xpm1(i,j) + minv
              
           end do
           
        end do
        !$OMP end parallel do
        
        ! estruturação do segundo grid de alta resolução
        
        call param(x_max, x_min, y_max, y_min, num2, r_fit, al, rough)
        
        !$OMP parallel do private(s_grid, dist, peso, k, j) reduction(min:minv)
        do i=1, n_grid+1
           
           do j=1, n_grid+1
              
              grid2(i,j)%x = (i-1)*dx + x_min
              grid2(i,j)%y = (j-1)*dy + y_min
              
              s_grid = 0
              grid2(i,j)%z = 0
              minv = 1000
              
              do k=1, num2 - 1
                 
                 dist = (grid2(i,j)%x - coarse2(k)%x)**2
                 dist = dist + (grid2(i,j)%y - coarse2(k)%y)**2
                 
                 minv = min(minv, dist)
                 
                 if (dist<r_fit*r_fit)then
                    
                    peso = exp(-(dist*al*al)/pi)
                    s_grid = s_grid + peso
                    grid2(i,j)%z = grid2(i,j)%z + peso*coarse2(k)%z
                    
                 end if
                 
              end do
              
              grid2(i,j)%z = grid2(i,j)%z/s_grid
              r_xpm2(i,j) = r_xpm2(i,j) + minv
              
           end do
           
        end do
        !$OMP end parallel do
        
        ! Fim da estruturação
        
        !$OMP parallel do private(j)
        do i=1, n_grid+1
           
           do j=1, n_grid+1
              
              grid3(i,j)%z = (grid(i,j)%z + grid2(i,j)%z)/2
              grid3(i,j)%x = grid(i,j)%x
              grid3(i,j)%y = grid(i,j)%y
              
           end do
           
        end do
        !$OMP end parallel do
        
        ! Calculando a área pela fórmula de Heron
        write(2, *) frame, calc_area(n_grid, grid3)/100/n_lipid
        
        !calculo do RMSD =======================================

        if (rmsd) then
           
           write(3, *) calc_rmsd(store, store2, grid, grid2, x_min, y_min, dx, dy, noi1, noi2)/10
           
        end if
        !======================================================================

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
          
  call system_clock (finish, clock_rate, clock_max)

  ! Fim da leitura do arquivo .pdb

  if (traj_type=='xtc') call trj%close()
  if (traj_type=='pdb') close(1)

  close(2)

  if (rmsd)then
     
     close(3)

  end if

  !creating the grid color map for area ====================

  do i=1, n_grid+1

     do j=1, n_grid+1

        r_xpm1(i,j) = r_xpm1(i,j)*pi/(100*tot_frame)
        r_xpm2(i,j) = r_xpm2(i,j)*pi/(100*tot_frame)
        
     end do

  end do

  !end of calculation ======================================
  
  if (p_grid)then

     call print_grid_xpm(grid, r_xpm1, n_grid, n_grid, 'grid1', back)
     call print_grid_xpm(grid2, r_xpm2, n_grid, n_grid, 'grid2', back)     
     call print_grid(grid3, n_grid+1, n_grid+1, 'grid3', back)
     
  end if

  call ending(back, finish, start, clock_rate) ! Finaliza programa e mostra tempo de processamento
  
  deallocate(r_xpm1)
  deallocate(r_xpm2)
  
end program suave_area
