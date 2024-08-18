program suave_gridsph

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
     n_grid, bin_out, fr_in, fr_end, n_skip, n_lipid, get, div, 's_gridsph ', version)

  !Leitura do PDB
12 format(a6, i5.1, a5, a5, a1, i4.1, a4, 3f8.3)
  
 !==Lendo os arquivos de index=========================

  if (.not.outer)then
     
     call abre_ndx(ind, in_num, n_index)
     
  end if

  !=================definindo frames para inicio e fim========
  
  call def_frame(frame, fr_in, fr_end, skip, n_skip, end, begin)

  !=================calculando o espaçamento==================
  
  call def_bin_sph(bin_coarse, n_index, bin, n_grid)
  
  !==================== abrindo trajetórias =====================
 
  traj_type = coord(len(trim(coord))-2:len(trim(coord)))

  if (traj_type=='xtc') call trj%open(coord)
  if (traj_type=='pdb') call abre_trj(1, coord)
 
  !==================== abrindo outputs =========================

  if (rmsd)then

     call abre('rmsd      ', 3, 'xvg', back)
     
     write(3, '(a7, a7)') "#SuAVE ", version
     write(3, '(a14)') '#Command Line:'
     write(3, '(a12)', advance='no') '#s_gridsph  '
     write(3, *) (trim(get(i)),"  ", i=1, 20)
     write(3, *) '@    title "RMSD X Frame"'
     write(3, *) '@    xaxis  label "#Frame"'
     write(3, *) '@    yaxis  label "RMSD [nm]"'

  end if

  call abre('grid      ', 4, 'pdb', back) 

  call system_clock(start, clock_rate, clock_max)

! Leitura do arquivo .pdb

  i_atom = 0
  n_index = 1
  ierr = 0
  num = 1
  cent_x = 0
  cent_y = 0
  cent_z = 0
  x_min = 1000
  y_min = 1000
  z_min = 1000
  x_max = 0
  y_max = 0
  z_max = 0
  
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

              if ((i_atom == in_num(num)))then 

                 spher(n_index) = buff
                 cent_x = (cent_x*(n_index-1) + buff%x)/(n_index)
                 cent_y = (cent_y*(n_index-1) + buff%y)/(n_index)
                 cent_z = (cent_z*(n_index-1) + buff%z)/(n_index)
                 num = num  + 1
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

           if (i_atom == in_num(num)) then

              copy = trj%x(1, i_atom)
              spher(n_index)%x = copy(1)*10
              spher(n_index)%y = copy(2)*10
              spher(n_index)%z = copy(3)*10
              cent_x = (cent_x*(n_index-1) + copy(1)*10)/(n_index)
              cent_y = (cent_y*(n_index-1) + copy(2)*10)/(n_index)
              cent_z = (cent_z*(n_index-1) + copy(3)*10)/(n_index)
              num = num  + 1
              n_index = n_index + 1

              x_max = max(x_max, copy(1)*10)
              x_min = min(x_min, copy(1)*10)
              y_max = max(y_max, copy(2)*10)
              y_min = min(y_min, copy(2)*10)
              z_max = max(z_max, copy(3)*10)
              z_min = min(z_min, copy(3)*10)
              
           end if

        end do

     end if

     ! Fitting process =============================================================                                                                                                                        

     eval_skip = ((mod((frame-fr_in),(n_skip+1))==0))

     if ((frame>fr_in-1).and.(frame<fr_end+1).and.(eval_skip).and.(ierr>=0)) then

        tot_frame = tot_frame + 1 ! para programas que calculam propriedades medias         
          
        num = 1
        r_med1 = 0
              
        ! transformando em Coordenadas Esfericas ======================================================
              
        do i=1, n_index - 1

           call cart2sphe(num, r_med1, cent_x, cent_y, cent_z, store(num), spher(i))                 
              
        end do
              
        !! A partir desse ponto as coordenadas são esféricas e estão com o centro transladado para
        !! o centro da micela. Logo, quando retornarmos para as coordenadas cartesianas teremos que
        !! transladar novamente 
              
        noi1 = num - 1
              
        if (l_coarse) then
                 
           ! Estruturação do primeiro coarse grid=============================================
                 
           dph = (2*pi)/(2*(bin_coarse - int(bin_coarse/2)))
           dth = (2*pi)/bin_coarse

           call param_esf(r_med1, num, r_fit, al, rough)

           !$OMP parallel do private(s_grid, dist, j, k, aux)
           do i=1, (bin_coarse - int(bin_coarse/2)) + 1
                    
              do j=1, bin_coarse + 1
                       
                 aux = (j-1)*(bin_coarse - int(bin_coarse/2) + 1) + i
                 coarse(aux)%phi = (i-1)*dph
                 coarse(aux)%theta = (j-1)*dth 
                       
                 s_grid = 0
                 coarse(aux)%rho = 0
                       
                 do k=1, num - 1
                          
                    dist = sin(coarse(aux)%phi)*sin(store(k)%phi)*cos(coarse(aux)%theta-store(k)%theta)
                    dist = dist + cos(coarse(aux)%phi)*cos(store(k)%phi)
                    dist = store(k)%rho*acos(dist)
                    dist = al*dist
                          
                    if (dist/al<r_fit)then
                             
                       s_grid = s_grid + exp(-dist*dist/pi)
                       coarse(aux)%rho = coarse(aux)%rho + exp(-dist*dist/pi)*store(k)%rho
                             
                    end if
                          
                 end do
                       
                 coarse(aux)%rho = coarse(aux)%rho/s_grid
                       
              end do
                    
           end do
           !$OMP end parallel do
                 
           num = (bin_coarse + 1)*(bin_coarse - int(bin_coarse/2) + 1) + 1
                 
        else
                 
           coarse = store
                 
        end if
              
        ! estruturação do prmeiro grid de alta resolução========================
              
        dph = (2*pi)/(2*(n_grid - int(n_grid/2)))
        dth = (2*pi)/n_grid
              
        call param_esf(r_med1, num, r_fit, al, rough)

        !$OMP parallel do private(s_grid, dist, j, k)
        do i=1, (n_grid - int(n_grid/2)) + 1
                 
           do j=1, n_grid+1
                    
              grid(i,j)%phi = (i-1)*dph
              grid(i,j)%theta = (j-1)*dth
                    
              s_grid = 0
              grid(i,j)%rho = 0
                    
              do k=1, num - 1 
                       
                 dist = sin(grid(i,j)%phi)*sin(coarse(k)%phi)*cos(grid(i,j)%theta-coarse(k)%theta)
                 dist = dist + cos(grid(i,j)%phi)*cos(coarse(k)%phi)
                 dist = coarse(k)%rho*acos(dist)
                 dist = al*dist
                       
                 if (dist/al<r_fit)then
                          
                    s_grid = s_grid + exp(-dist*dist/pi)
                    grid(i,j)%rho = grid(i,j)%rho + exp(-dist*dist/pi)*coarse(k)%rho
                          
                 end if
                       
              end do
                    
              grid(i,j)%rho = grid(i,j)%rho/s_grid

           end do
                 
        end do
        !$OMP end parallel do
              
        ! Fim da estruturação========================================

        ! Calculando o RMSD========================================

        if (rmsd)then
                 
           noir = 0
                 
           do i=1, noi1
                    
              a = nint((store(i)%phi)/dph) + 1
              b = nint((store(i)%theta)/dth) + 1
              dist_z = store(i)%rho-grid(a,b)%rho
              noir = noir + dist_z*dist_z/noi1
                    
           end do
                                  
           noir = sqrt(noir)
                 
           write(3, *) noir/10
                 
        end if
        
        ! Fim do calculo do RMSD===================================

        ! Escreve cada frame da trajetoria do grid=================

        call sphe2cart((n_grid-int(n_grid/2))+1, n_grid+1, cent_x, cent_y, cent_z, grid, grid3)
        call print_pdb(grid3, (n_grid-int(n_grid/2))+1, n_grid+1, 4, x_max-x_min, y_max-y_min, z_max-z_min, frame, back)
                  
        !fim da escrita da trajetoria======================================

        n_index = 1
        i_atom = 0
        num = 1

        if (ierr<0) then

           cent_x = center%x
           cent_y = center%y
           cent_z = center%z

        end if
        
        center%x = cent_x
        center%y = cent_y
        center%z = cent_z
        cent_x = 0
        cent_y = 0
        cent_z = 0
        x_min = 1000
        y_min = 1000
        z_min = 1000
        x_max = 0
        y_max = 0
        z_max = 0
              
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

  if (rmsd)then
     
     close(3)

  end if
  
  close(4)

  ! Desenho dos pontos usados no ajuste

  call abre('adjust    ', 3, 'pdb', back)
  
  write(3, *) 'TITLE     Rectangular Regular GRID'
  write(3, *) 'REMARK    THIS IS A SIMULATION BOX'
  write(3, *) 'MODEL        1'

  do i=1, noi1
     
     write(3, 12, iostat=ierr) 'ATOM  ', i, '  DOT', &
          '  DOT',' ', 1,'    ', spher(i)%x, &
          spher(i)%y, spher(i)%z

  end do

  close(3)
  
  if (l_coarse) then

     call abre('coarse    ', 3, 'pdb', back)
     
     write(3, *) 'TITLE     Rectangular Regular GRID'
     write(3, *) 'REMARK    THIS IS A SIMULATION BOX'
     write(3, *) 'MODEL        1'
     
     do i=1, (bin_coarse - int(bin_coarse/2)+1)*(bin_coarse+1)

        buff%x = coarse(i)%rho*sin(coarse(i)%phi)*cos(coarse(i)%theta) + center%x
        buff%y = coarse(i)%rho*sin(coarse(i)%phi)*sin(coarse(i)%theta) + center%y
        buff%z = coarse(i)%rho*cos(coarse(i)%phi) + center%z
        
        write(3, 12, iostat=ierr) 'ATOM  ', i, '  DOT', &
             '  DOT',' ', 1,'    ', buff%x, &
             buff%y, buff%z

     end do

  end if
  
  ! Fim Arquivo escrito

  close(3)

  call ending(back, finish, start, clock_rate) ! Finaliza programa e mostra tempo de processamento
  
end program suave_gridsph