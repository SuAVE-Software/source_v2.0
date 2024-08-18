program suave_spher

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
     n_grid, bin_out, fr_in, fr_end, n_skip, n_lipid, get, div, 's_densph  ', version)

  !Leitura do PDB
12 format(a6, i5.1, a5, a5, a1, i4.1, a4, 3f8.3)

  if (ind2=='miss')then

     write(*, *)
     write(*, *)'Second index file is missing'
     write(*, *)
     stop

  end if

  if (ind3=='miss')then

     write(*, *)
     write(*, *) 'Density index file is missing'
     write(*, *)
     stop

  end if

  if (.not.slices) then

     div = 10

     write(*, *) ""
     write(*, *) "STD_SLICES = ", div

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

  !============================================================

  call abre_ndx(ind3, in_dens, n_index)
  
  !==================== abrindo trajetórias =====================
 
  traj_type = coord(len(trim(coord))-2:len(trim(coord)))

  if (traj_type=='xtc') call trj%open(coord)
  if (traj_type=='pdb') call abre_trj(1, coord)
 
  !==================== abrindo outputs =========================
  
  if (rmsd)then

     call abre('rmsd      ', 3, 'xvg', back)

     write(3, '(a7, a7)') "#SuAVE ", version
     write(3, '(a14)') '#Command Line:'
     write(3, '(a11)', advance='no') '#s_densph  '
     write(3, *) (trim(get(i)),"  ", i=1, 20)
     write(3, *) '@    title "RMSD X Frame"'
     write(3, *) '@    xaxis  label "#Frame"'
     write(3, *) '@    yaxis  label "RMSD [nm]"'

  end if
  
  ! Limpando o histograma==============================

  do i=1, 1000

     hist(i) = 0

  end do

  ! histograma limpo===================================

  call system_clock(start, clock_rate, clock_max)

! Leitura do arquivo .pdb

  i_atom = 0
  a_dens = 1
  n_index = 1
  ierr = 0
  num = 1
  num2 = 1
  cent_x = 0
  cent_y = 0
  cent_z = 0
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

              if (i_atom == in_dens(a_dens))then

                 dens_spher(a_dens) = buff
                 a_dens = a_dens + 1

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

           if (i_atom == in_dens(a_dens))then

              copy = trj%x(1, i_atom)
              dens_spher(a_dens)%x = copy(1)*10
              dens_spher(a_dens)%y = copy(2)*10
              dens_spher(a_dens)%z = copy(3)*10
              a_dens = a_dens + 1

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

        ! Transformando em Coordenadas Esfericas ======================================================

        do i=1, n_index - 1

           if ((spher(i)%n_atom == in_num(num))) then

              call cart2sphe(num, r_med1, cent_x, cent_y, cent_z, store(num), spher(i))

           else

              call cart2sphe(num2, r_med2, cent_x, cent_y, cent_z, store2(num2), spher(i))

           end if

        end do
              
        do i=1, a_dens - 1

           !aux e aux2 são apenas variáveis auxiliares e não vão guardar estes valores
           call cart2sphe(aux, aux2, cent_x, cent_y, cent_z, dens(i), dens_spher(i))
                 
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

        ! Contrução do grid médio para cálculo de densidade =========

        do i=1, (n_grid - int(n_grid/2))+1
                 
           do j=1, n_grid+1
                    
              grid3(i,j)%x = (grid(i,j)%rho + grid2(i,j)%rho)*sin(grid(i,j)%phi)*cos(grid(i,j)%theta)/2 + cent_x
              grid3(i,j)%y = (grid(i,j)%rho + grid2(i,j)%rho)*sin(grid(i,j)%phi)*sin(grid(i,j)%theta)/2 + cent_y
              grid3(i,j)%z = (grid(i,j)%rho + grid2(i,j)%rho)*cos(grid(i,j)%phi)/2 + cent_z

           end do
                 
        end do

        ! Calculando volume da superficie media ====================

        call calc_area_sph((n_grid-int(n_grid/2))+1, n_grid+1, s_area, s_vol, dph, dth, grid, grid3)
              
        ! cálculo da densidade esférica ============================
        ! o calculo é em funcao de K ===============================

        call calc_dens_sph(a_dens, hist, s_vol, dens, grid, grid2, div, del, dth, dph)
              
        ! Calculando o RMSD=========================================

        if (rmsd)then
                 
           write(3, *) calc_rmsd_sph(noi1, noi2, store, store2, grid, grid2, dph, dth)/10
                 
        end if

        ! Fim do calculo do RMSD====================================

        a_dens = 1
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
           
  if (rmsd)then
     
     close(3)

  end if
  
  ! Imprimindo arquivo de densidade ===========================================

  call abre('density   ', 5, 'xvg', back)  

  write(5, '(a7, a7)') "#SuAVE ", version
  write(5, '(a14)') '#Command Line:'
  write(5, '(a11)', advance='no') '#s_densph  '
  write(5, *) (trim(get(i)),"  ", i=1, 20)
  write(5, *) '@    title "System Density"'
  write(5, *) '@    xaxis  label "R/R\smed\N"'
  write(5, *) '@    yaxis  label "Density [nm\S-3\N]"'

  do i=500, 500+2*div

     write(5, *) (i-500)*del, hist(i)/tot_frame

  end do

  close(5)

  !Arquivo escrito ===========================================================

  if (p_grid)then

     call sphe2cart((n_grid-int(n_grid/2))+1, n_grid+1, center%x, center%y, center%z, grid, grid3)
     call print_grid(grid3, (n_grid-int(n_grid/2))+1, n_grid+1, 'grid1', back)
     
     call sphe2cart((n_grid-int(n_grid/2))+1, n_grid+1, center%x, center%y, center%z, grid2, grid3)
     call print_grid(grid3, (n_grid-int(n_grid/2))+1, n_grid+1, 'grid2', back)
     
  end if

  call ending(back, finish, start, clock_rate) ! Finaliza programa e mostra tempo de processamento
  
end program suave_spher