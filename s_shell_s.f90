program suave_shell

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
     n_grid, bin_out, fr_in, fr_end, n_skip, n_lipid, get, div, 's_shell   ', version)

  !Leitura do PDB
12 format(a6, i5.1, a5, a5, a1, i4.1, a4, 3f8.3)

  if (ind2=='miss')then

     write(*, *)
     write(*, *)'Second index file is missing'
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
  
  call def_bin_sph(bin_coarse, n_index, bin, n_grid)

  !==================== abrindo trajetórias =====================
 
  traj_type = coord(len(trim(coord))-2:len(trim(coord)))

  if (traj_type=='xtc') call trj%open(coord)
  if (traj_type=='pdb') call abre_trj(1, coord)
 
  !==================== abrindo outputs =========================
  
  call abre('thick     ', 2, 'xvg', back)
  write(2, '(a7, a7)') "#SuAVE ", version
  write(2, '(a14)') '#Command Line:'
  write(2, '(a10)', advance='no') '#s_shell  '
  write(2, *) (trim(get(i)),"  ", i=1, 20)
  write(2, *) '@    title "Thickness X Frame"'
  write(2, *) '@    xaxis  label "#Frame"'
  write(2, *) '@    yaxis  label "Average Thickness [nm]"'

  if (rmsd)then

     call abre('rmsd      ', 3, 'xvg', back)
     
     write(3, '(a7, a7)') "#SuAVE ", version
     write(3, '(a14)') '#Command Line:'
     write(3, '(a10)', advance='no') '#s_shell  '
     write(3, *) (trim(get(i)),"  ", i=1, 20)
     write(3, *) '@    title "RMSD X Frame"'
     write(3, *) '@    xaxis  label "#Frame"'
     write(3, *) '@    yaxis  label "RMSD [nm]"'

  end if

  call abre('volume    ', 4, 'xvg', back)

  write(4, '(a7, a7)') "#SuAVE ", version
  write(4, '(a14)') '#Command Line:'
  write(4, '(a10)', advance='no') '#s_shell  '
  write(4, *) (trim(get(i)),"  ", i=1, 20)
  write(4, *) '@    title "Volume X Frame"'
  write(4, *) '@    xaxis  label "#Frame"'
  write(4, *) '@    yaxis  label "Total Spherical Shell Volume[nm\S3\N]"'
  
  ! =====================================================================

  allocate (xpm1((n_grid - int(n_grid/2)), n_grid), stat = ierr)
  if (ierr /= 0) stop 'Not enough memory to initialize xpm matrix'

  allocate (r_xpm1((n_grid - int(n_grid/2)), n_grid), stat = ierr)
  if (ierr /= 0) stop 'Not enough memory to initialize xpm matrix'

  do i=1, n_grid - int(n_grid/2) 

     do j=1, n_grid 

        r_xpm1(i,j) = 0
        
     end do

  end do
  
  call system_clock(start, clock_rate, clock_max)

! Leitura do arquivo .pdb

  i_atom = 0
  n_index = 1
  ierr = 0
  num = 1
  num2 = 1
  cent_x = 0
  cent_y = 0
  cent_z = 0
  minv = 100000
  maxv = -100000
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

              if ((i_atom == in_num(num)))then 

                 spher(n_index) = buff
                 cent_x = (cent_x*(n_index-1) + buff%x)/(n_index)
                 cent_y = (cent_y*(n_index-1) + buff%y)/(n_index)
                 cent_z = (cent_z*(n_index-1) + buff%z)/(n_index)
                 num = num  + 1
                 n_index = n_index + 1

              end if

              if ((i_atom == in_num2(num2)))then

                 spher(n_index) = buff
                 cent_x = (cent_x*(n_index-1) + buff%x)/(n_index)
                 cent_y = (cent_y*(n_index-1) + buff%y)/(n_index)
                 cent_z = (cent_z*(n_index-1) + buff%z)/(n_index)
                 num2 = num2  + 1
                 n_index = n_index + 1

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
              spher(n_index)%n_atom = i_atom
              spher(n_index)%x = copy(1)*10
              spher(n_index)%y = copy(2)*10
              spher(n_index)%z = copy(3)*10
              cent_x = (cent_x*(n_index-1) + copy(1)*10)/(n_index)
              cent_y = (cent_y*(n_index-1) + copy(2)*10)/(n_index)
              cent_z = (cent_z*(n_index-1) + copy(3)*10)/(n_index)
              num = num  + 1
              n_index = n_index + 1
              
           end if

           if ((i_atom == in_num2(num2)))then

              copy = trj%x(1, i_atom)
              spher(n_index)%n_atom = i_atom
              spher(n_index)%x = copy(1)*10
              spher(n_index)%y = copy(2)*10
              spher(n_index)%z = copy(3)*10
              cent_x = (cent_x*(n_index-1) + copy(1)*10)/(n_index)
              cent_y = (cent_y*(n_index-1) + copy(2)*10)/(n_index)
              cent_z = (cent_z*(n_index-1) + copy(3)*10)/(n_index)
              num2 = num2  + 1
              n_index = n_index + 1
              
           end if

        end do

     end if

     ! Fitting process =============================================================                                                                                                                        

     eval_skip = ((mod((frame-fr_in),(n_skip+1))==0))

     if ((frame>fr_in-1).and.(frame<fr_end+1).and.(eval_skip).and.(ierr>=0)) then

        tot_frame = tot_frame + 1 ! para programas que calculam propriedades medias
         
        num = 1
        num2 = 1
        r_med1 = 0
        r_med2 = 0

        ! transformando em Coordenadas Esfericas ======================================================
              
        do i=1, n_index - 1
                 
           if ((spher(i)%n_atom == in_num(num))) then
                    
              call cart2sphe(num, r_med1, cent_x, cent_y, cent_z, store(num), spher(i))
                    
           else
                    
              call cart2sphe(num2, r_med2, cent_x, cent_y, cent_z, store2(num2), spher(i))
                    
           end if
                 
        end do
                     
        !! A partir desse ponto as coordenadas são esféricas e estão com o centro transladado para
        !! o centro da micela. Logo, quando retornarmos para as coordenadas cartesianas teremos que
        !! transladar novamente 
              
        noi1 = num - 1
        noi2 = num2 - 1
              
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
                 
           ! Estruturação do segundo coarse grid=======================================
                 
           call param_esf(r_med2, num2, r_fit, al, rough)

           !$OMP parallel do private(s_grid, dist, j, k, aux)
           do i=1, (bin_coarse - int(bin_coarse/2)) + 1
                    
              do j=1, bin_coarse + 1
                       
                 aux = (j-1)*(bin_coarse - int(bin_coarse/2) + 1) + i 
                 coarse2(aux)%phi = (i-1)*dph 
                 coarse2(aux)%theta = (j-1)*dth
                       
                 s_grid = 0
                 coarse2(aux)%rho = 0
                       
                 do k=1, num2 - 1
                          
                    dist = sin(coarse2(aux)%phi)*sin(store2(k)%phi)*cos(coarse2(aux)%theta-store2(k)%theta)
                    dist = dist + cos(coarse2(aux)%phi)*cos(store2(k)%phi)
                    dist = store2(k)%rho*acos(dist)
                    dist = al*dist
                          
                    if (dist/al<r_fit)then
                             
                       coarse2(aux)%rho = coarse2(aux)%rho + exp(-dist*dist/pi)*store2(k)%rho
                       s_grid = s_grid + exp(-dist*dist/pi)
                             
                    end if
                          
                 end do
                       
                 coarse2(aux)%rho = coarse2(aux)%rho/s_grid
                       
              end do
                    
           end do
           !$OMP end parallel do
                 
           num2 = (bin_coarse + 1)*(bin_coarse - int(bin_coarse/2) + 1) + 1
                 
        else
                 
           coarse = store
           coarse2 = store2
                 
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
              
        ! estruturação do segundo grid de alta resolução========================
              
        call param_esf(r_med2, num2, r_fit, al, rough)

        !$OMP parallel do private(s_grid, dist, j, k)
        do i=1, (n_grid - int(n_grid/2)) + 1
                 
           do j=1, n_grid + 1
                    
              grid2(i,j)%phi = (i-1)*dph
              grid2(i,j)%theta = (j-1)*dth
                    
              s_grid = 0
              grid2(i,j)%rho = 0
                    
              do k=1, num2 - 1
                       
                 dist = sin(grid2(i,j)%phi)*sin(coarse2(k)%phi)*cos(grid2(i,j)%theta-coarse2(k)%theta)
                 dist = dist + cos(grid2(i,j)%phi)*cos(coarse2(k)%phi)
                 dist = coarse2(k)%rho*acos(dist)
                 dist = al*dist
                       
                 if (dist/al<r_fit)then
                          
                    s_grid = s_grid + exp(-dist*dist/pi)
                    grid2(i,j)%rho = grid2(i,j)%rho + exp(-dist*dist/pi)*coarse2(k)%rho
                          
                 end if
                       
              end do
                    
              grid2(i,j)%rho = grid2(i,j)%rho/s_grid
                    
           end do
                 
        end do
        !$OMP end parallel do
              
        ! Fim da estruturação========================================
              
        ! Calculando o volume do grid==============================
              
        call sphe2cart((n_grid-int(n_grid/2))+1, n_grid+1, cent_x, cent_y, cent_z, grid, grid3)

        call calc_area_sph((n_grid-int(n_grid/2))+1, n_grid+1, s_area, s_vol, dph, dth, grid, grid3)

        ! Calculando o volume do grid2 ============================

        call sphe2cart((n_grid-int(n_grid/2))+1, n_grid+1, cent_x, cent_y, cent_z, grid2, grid3)
        
        call calc_area_sph((n_grid-int(n_grid/2))+1, n_grid+1, s_area2, s_vol2, dph, dth, grid2, grid3)

        write(4, *) frame, abs(s_vol - s_vol2)/(1000)
              
        ! Fim do calculo ==========================================

        ! Calculando o thickness médio=============================
              
        call calc_thick_sph((n_grid-int(n_grid/2))+1, n_grid+1, aver, aver2, r_xpm1, grid, grid2)
              
        desv = sqrt(aver2 - aver*aver)
        minv = min(aver-2*desv, minv)
        maxv = max(aver+2*desv, maxv)

        write(2, *) frame, aver
              
        ! Calculando o RMSD========================================

        if (rmsd)then
                 
           write(3, *) calc_rmsd_sph(noi1, noi2, store, store2, grid, grid2, dph, dth)/10
                 
        end if
        
        ! Fim do calculo do RMSD===================================

        n_index = 1
        i_atom = 0
        num = 1
        num2 = 1

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
  
  close(4)

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
  write(4, '(a10)', advance='no') '#s_shell  '
  write(4, *) (trim(get(i)),"  ", i=1, 20)
  write(4, *) '/* XPM */'
  write(4, *)'/* This matrix is generated by s_shell.f90 */'
  write(4, *)'/* title:   "Thickness" */'
  write(4, *)'/* x-label: "Theta angle " */'
  write(4, *)'/* y-label: "Phi angle " */'
  write(4, *)'/* type:    "Continuous" */'
  write(4, *)'static char * gv_xpm[] = {'
  write(4, *) '"',n_grid, (n_grid - int(n_grid/2)),' 7 1",'

  do i=1, n_grid - int(n_grid/2)
     
     do j=1, n_grid 
        
        r_xpm1(i,j) = r_xpm1(i,j)/tot_frame/10

     end do
     
  end do

  del = (maxv - minv)/6    

  call  print_xpm_sph(n_grid-int(n_grid/2), n_grid, del, dth, dph, xpm1, r_xpm1, minv, 'b')
  
  close(4)

  !Creating grid files ==============================================================

  if (p_grid)then

     call sphe2cart((n_grid-int(n_grid/2))+1, n_grid+1, center%x, center%y, center%z, grid, grid3)
     call print_grid_xpm(grid3, r_xpm1, (n_grid-int(n_grid/2)), n_grid, 'grid1', back)
     
     call sphe2cart((n_grid-int(n_grid/2))+1, n_grid+1, center%x, center%y, center%z, grid2, grid3)
     call print_grid_xpm(grid3, r_xpm1, (n_grid-int(n_grid/2)), n_grid, 'grid2', back)
     
  end if

  call ending(back, finish, start, clock_rate) ! Finaliza programa e mostra tempo de processamento

  deallocate(r_xpm1)
  deallocate(xpm1)
    
end program suave_shell