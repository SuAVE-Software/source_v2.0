program suave_grid
  
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
     n_grid, bin_out, fr_in, fr_end, n_skip, n_lipid, get, div, 's_grid    ', version)

  !Leitura do PDB
12 format(a6, i5.1, a5, a5, a1, i4.1, a4, 3f8.3)

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
 
  !==================== abrindo outputs =========================

  if (rmsd)then

     call abre('rmsd      ', 3, 'xvg', back)

     write(3, '(a7, a7)') "#SuAVE ", version
     write(3, '(a14)') '#Command Line:'
     write(3, '(a9)', advance='no') '#s_grid  '
     write(3, *) (trim(get(i)),"  ", i=1, 20)
     write(3, *) '@    title "RMSD X Frame"'
     write(3, *) '@    xaxis  label "#Frame"'
     write(3, *) '@    yaxis  label "RMSD [nm]"'

  end if

  !==========================================================

  call abre('grid      ', 4, 'pdb', back)

  call system_clock(start, clock_rate, clock_max)
  
  ! Leitura do arquivo .pdb
  
  i_atom = 0
  frame = 0
  ierr = 0
  num = 1
  n_index = 1
  x_min = 1000
  y_min = 1000
  z_min = 1000
  x_max = 0
  y_max = 0
  z_max = 0
  
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
                 z_max = max(z_max, buff%z)
                 z_min = min(z_min, buff%z)
           
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
              store(num)%x = copy(1)*10
              store(num)%y = copy(2)*10
              store(num)%z = copy(3)*10
              n_index = n_index + 1
              num = num + 1

              x_max = max(x_max, store(num-1)%x)
              x_min = min(x_min, store(num-1)%x)
              y_max = max(y_max, store(num-1)%y)
              y_min = min(y_min, store(num-1)%y)
              z_max = max(z_max, store(num-1)%z)
              z_min = min(z_min, store(num-1)%z)

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
                 
           ! refinamento do r_fit e do alfa
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
              
        ! refinamento do r_fit e do alfa

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
        !$OMP end parallel do
              
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
              
        ! escreve cada frame da trajetoria do grid
        call print_pdb(grid, n_grid+1, n_grid+1, 4, x_max-x_min, y_max-y_min, z_max-z_min, frame, back)              

        !fim da escrita da trajetoria

        n_atom = n_index-1
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
     if (.not.end)then
              
        fr_end = frame + 1
              
     end if
     !======================================================
     
     if ((frame>=trj%NFRAMES).and.(traj_type=='xtc')) ierr = -1 ! apenas quando lendo o tipo XTC

  end do

  call system_clock (finish, clock_rate, clock_max)

  ! Fim da estruturação

  if (traj_type=='xtc') call trj%close()
  if (traj_type=='pdb') close(1)

  close(4)

  if (rmsd)then
     
     close(3)
     
  end if
    
  ! desenho dos pontos usados no ajuste
  
  call abre('adjust    ', 3, 'pdb', back)

  write(3, *) 'TITLE     Rectangular Regular GRID'
  write(3, *) 'REMARK    THIS IS A SIMULATION BOX'
  write(3, *) 'MODEL        1'

  do i=1, noi1
     
     write(3, 12, iostat=ierr) 'ATOM  ', i, '  DOT', &
          '  DOT',' ', 1,'    ', store(i)%x, &
          store(i)%y, store(i)%z

  end do
  
  write(3, '(a3)') 'TER'
  write(3, '(a6)') 'ENDMDL'
  
  close(3)
  
  if (l_coarse) then
     
     call abre('coarse    ', 3, 'pdb', back)
     
     write(3, *) 'TITLE     Rectangular Regular GRID'
     write(3, *) 'REMARK    THIS IS A SIMULATION BOX'
     write(3, *) 'MODEL        1'
     
     do i=1, (bin_coarse + 1)*(bin_coarse + 1)
        
        write(3, 12, iostat=ierr) 'ATOM  ', i, '  DOT', &
             '  DOT',' ', 1,'    ', coarse(i)%x, &
             coarse(i)%y, coarse(i)%z
        
     end do
     
     write(3, '(a3)') 'TER'
     write(3, '(a6)') 'ENDMDL'
     
     close(3)
     
  end if
  
  call ending(back, finish, start, clock_rate) ! Finaliza programa e mostra tempo de processamento
       
end program suave_grid