program s_inertia

  use types
  use variables
  use funcproc

  use gmxfort_trajectory
  use gmxfort_utils
 
  ! ==================================  
  type(Trajectory) :: trj
  real :: copy(3)
  double precision :: MI(3,3)
  double precision :: Q(3,3)
  double precision :: W(3), PCA(3), ORD(3)
  ! ==================================
  
  call startup(outer, bin, p_grid, coord, ind, ind2, rmsd, map, ind3, &
     l_coarse, begin, end, skip, lipid, rough, slices, inside, range, &
     n_grid, bin_out, fr_in, fr_end, n_skip, n_lipid, get, div, 's_inertia ', version)

  !Leitura do PDB
12 format(a6, i5.1, a5, a5, a1, i4.1, a4, 3f8.3)
  
  if (coord=='miss')then

     write(*, *)
     write(*, *)'PDB file is missing'
     write(*, *)
     stop

  end if

  if (ind=='miss')then
     
     write(*, *)
     write(*, *)'Index file is missing'
     write(*, *)
     stop
     
  end if
  
  !==Lendo os arquivos de index=========================

  call abre_ndx(ind, in_num, n_index)
  
  !=================definindo frames para inicio e fim========

  !=================definindo frames para inicio e fim========
  
  call def_frame(frame, fr_in, fr_end, skip, n_skip, end, begin)
  
  !=================calculando o espaçamento==================
  
  bin_coarse = nint(sqrt(real(2*(n_index-1))))
  
  if (.not.bin)then
     
     n_grid = nint(sqrt(real(2*(n_index-1))))
     write(*, *)
     write(*, *) 'STD_BIN = ', n_grid
     
  end if
  
  !==================== abrindo trajetórias =====================
 
  traj_type = coord(len(trim(coord))-2:len(trim(coord)))

  if (traj_type=='xtc') call trj%open(coord)
  if (traj_type=='pdb') call abre_trj(1, coord)
 
  !==================== abrindo outputs =========================

  call abre('inertia   ', 2, 'xvg', back)
  
  write(2, '(a7, a7)') "#SuAVE ", version
  write(2, '(a14)') '#Command Line:'
  write(2, '(a12)', advance='no') '#s_inertia  '
  write(2, *) (trim(get(i)),"  ", i=1, 20)
  write(2, *) '@    title "Principal Moments of Inertia X Frame"'
  write(2, *) '@    xaxis  label "#Frame"'
  write(2, *) '@    yaxis  label "Principal Moments of Inertia [nm\S2\N]"'

  call abre('roundness ', 4, 'xvg', back)
  
  write(4, '(a7, a7)') "#SuAVE ", version
  write(4, '(a14)') '#Command Line:'
  write(4, '(a12)', advance='no') '#s_inertia  '
  write(4, *) (trim(get(i)),"  ", i=1, 20)
  write(4, *) '@    title "Inertia Roundness X Frame"'
  write(4, *) '@    xaxis  label "#Frame"'
  write(4, *) '@    yaxis  label "Roundness"'
  write(4, *) '# Frame    Inert. Round.     Circle Round.'

  if (rmsd)then

     call abre('rmsd      ', 3, 'xvg', back)
     
     write(3, '(a7, a7)') "#SuAVE ", version
     write(3, '(a14)') '#Command Line:'
     write(3, '(a12)', advance='no') '#s_inertia  '
     write(3, *) (trim(get(i)),"  ", i=1, 20)
     write(3, *) '@    title "RMSD X Frame"'
     write(3, *) '@    xaxis  label "#Frame"'
     write(3, *) '@    yaxis  label "RMSD [nm]"'

  end if
  !====================================================================

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
                 
           dz = 2.0/(bin_coarse)
           dth = (2*pi)/bin_coarse

           call param_esf(r_med1, num, r_fit, al, rough)

           !$OMP parallel do private(s_grid, dist, j, k, aux)
           do i=1, (bin_coarse) 
                    
              do j=1, bin_coarse + 1
                       
                 aux = (j-1)*bin_coarse + i
                 coarse(aux)%phi = acos(i*dz-1-dz/2)
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
                 
           num = bin_coarse*(bin_coarse + 1)        
                 
        else
                 
           coarse = store
                 
        end if
              
        ! estruturação do prmeiro grid de alta resolução========================
              
        dz = 2.0/(n_grid)
        dth = (2*pi)/n_grid
              
        call param_esf(r_med1, num, r_fit, al, rough)
              
        minv = 100000
        maxv = -100000

        !$OMP parallel do private(s_grid, dist, j, k) reduction(min:minv) reduction(max:maxv)
        do i=1, (n_grid) 
                
           do j=1, n_grid+1

              grid(i,j)%phi = acos(i*dz-1-dz/2)
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
              maxv = max(maxv, grid(i,j)%rho)
              minv = min(minv, grid(i,j)%rho)
                    
           end do
                 
        end do
        !$OMP end parallel do
                            
        ! Fim da estruturação========================================
        ! Calculo dos momentos principais de inercia ================

        call calc_inertia(n_grid, grid, grid3, MI)
        call DSYEVJ3(MI, Q, PCA)

        !PCA possue os valores xx, yy e zz (variâncias)
        !W possue os valores de momento de inércia

        W(1) = PCA(1)
        W(2) = PCA(2)
        W(3) = PCA(3)

        call ordem(W, ORD)

        write(2, *) frame, ORD
        write(4, *) frame, ORD(3)/ORD(1), minv/maxv

        ! Fim do calculo ==========================================
              
        ! Calculando o RMSD========================================
        
        if (rmsd)then
                              
           write(3, *) frame, calc_rmsd_inert(noi1, store, dz, dth, grid)/10
                 
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
  close(4)

  if (rmsd)then
     
     close(3)

  end if
  
  !Creating grid files ==============================================================
  
  if (p_grid)then

     call sphe2cart(n_grid, n_grid, center%x, center%y, center%z, grid, grid3)
     call print_grid(grid3, n_grid, n_grid, 'grid1', back)

  end if

  call ending(back, finish, start, clock_rate) ! Finaliza programa e mostra tempo de processamento
  
end program s_inertia
