program s_count

  use types
  use variables
  use funcproc

  use gmxfort_trajectory
  use gmxfort_utils
 
  !============================
  type(Trajectory) :: trj
  real :: copy(3)
  logical :: consta, molec
  integer :: inner
  real :: r_max
  ! ==================================
  
  call startup(outer, bin, p_grid, coord, ind, ind2, rmsd, map, ind3, &
     l_coarse, begin, end, skip, lipid, rough, slices, inside, range, &
     n_grid, bin_out, fr_in, fr_end, n_skip, n_lipid, get, div, 's_count   ', version)

  !Leitura do PDB
12 format(a6, i5.1, a5, a5, a1, i4.1, a4, 3f8.3)

  consta = .false.
  molec = .false.

  do i=1, 20

     call getarg(i, get(i))

  end do

  do i=1, 20

     if (get(i)=='-molec')then

        molec = .true.

     end if

     if (get(i)=='-const')then

        consta = .true.
        read(get(i+1), *, iostat=ierr) const
        call erro(ierr, "CONST")

     end if

  end do
  
  if (ind3=='miss')then

     write(*, *)
     write(*, *) 'Density index file is missing'
     write(*, *)
     stop

  end if

  if (.not.consta)then

     write(*, *)
     write(*, *) 'Using CONST_DEF = 1.0'
     const = 1.0

  end if

  !==Lendo os arquivos de index=========================
  
  if (.not.outer)then
     
     call abre_ndx(ind, in_num, n_index)
     
  end if

  call abre_ndx(ind3, in_dens, n_index)
  
  !=================definindo frames para inicio e fim========
  
  call def_frame(frame, fr_in, fr_end, skip, n_skip, end, begin)

  !=================calculando o espaçamento==================
  
  call def_bin_sph(bin_coarse, n_index, bin, n_grid)
  
  !==================== abrindo trajetórias =====================
 
  traj_type = coord(len(trim(coord))-2:len(trim(coord)))

  if (traj_type=='xtc') call trj%open(coord)
  if (traj_type=='pdb') call abre_trj(1, coord)
 
  !==================== abrindo outputs =========================

  call abre('molecules ', 2, 'xvg', back)

  write(2, '(a7, a7)') "#SuAVE ", version
  write(2, '(a14)') '#Command Line:'
  write(2, '(a9)', advance='no') '#s_count  '
  write(2, *) (trim(get(i)),"  ", i=1, 20)
  write(2, *) '@    title " Molecules X Frame"'
  write(2, *) '@    xaxis  label "#Frame"'
  write(2, *) '@    yaxis  label "Molecules [nm\S3\N]"'
  
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
  a_dens = 1
  
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
              spher(n_index)%x = copy(1)*10
              spher(n_index)%y = copy(2)*10
              spher(n_index)%z = copy(3)*10
              cent_x = (cent_x*(n_index-1) + copy(1)*10)/(n_index)
              cent_y = (cent_y*(n_index-1) + copy(2)*10)/(n_index)
              cent_z = (cent_z*(n_index-1) + copy(3)*10)/(n_index)
              num = num  + 1
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
        r_max = 0
              
        ! transformando em Coordenadas Esfericas ======================================================
              
        do i=1, n_index - 1

           call cart2sphe(num, r_med1, cent_x, cent_y, cent_z, store(num), spher(i))                 
           r_max = max(store(num-1)%rho, r_max)
                 
        end do

        do i=1, a_dens - 1
                 
           call cart2sphe(aux, r_med1, cent_x, cent_y, cent_z, dens(i), dens_spher(i))
                 
        end do                    

        !! A partir desse ponto as coordenadas são esféricas e estão com o centro transladado para
        !! o centro da micela. Logo, quando retornarmos para as coordenadas cartesianas teremos que
        !! transladar novamente 
              
        noi1 = num - 1
              
        if (l_coarse) then
                 
           ! Estruturação do primeiro coarse grid=============================================
                 
           dph = (2*pi)/(2*(bin_coarse - int(bin_coarse/2)))
           dth = (2*pi)/bin_coarse

           call param_esf(r_max, num, r_fit, al, rough)
                 
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
                             
                       s_grid = s_grid + exp(-dist*dist/pi)*exp(store(k)%rho)
                       coarse(aux)%rho = coarse(aux)%rho + exp(-dist*dist/pi)*exp(store(k)%rho)*store(k)%rho
                             
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
              
        call param_esf(r_max, num, r_fit, al, rough)

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
                          
                    s_grid = s_grid + exp(-dist*dist/pi)*exp(coarse(k)%rho)
                    grid(i,j)%rho = grid(i,j)%rho + exp(-dist*dist/pi)*exp(coarse(k)%rho)*coarse(k)%rho
                          
                 end if

              end do

              grid(i,j)%rho = grid(i,j)%rho/s_grid
                    
           end do
                 
        end do
        !$OMP end parallel do

        !========== MOLECULAS NO PORO =================================

        inner = 0

        do i=1, a_dens-1
                 
           a = nint(dens(i)%phi/dph) + 1
           b = nint(dens(i)%theta/dth) + 1
                 
           if (dens(i)%rho<=grid(a,b)%rho*const) then

              inner = inner + 1

           end if

        end do

        write(2, *) frame, inner

        ! Fim do calculo ===========================================

        n_index = 1
        aux = a_dens
        a_dens = 1
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
  
  if (p_grid)then

     call sphe2cart((n_grid-int(n_grid/2))+1, n_grid+1, center%x, center%y, center%z, grid, grid3)
     call print_grid(grid3, (n_grid-int(n_grid/2))+1, n_grid+1, 'grid1', back)
    
  end if

  !escreve os átomos dentro da estrutura para o último frame

  if (molec) then

     call abre('molecules ', 3, 'pdb', back)
         
     do i=1, aux-1
        
        a = nint(dens(i)%phi/dph) + 1
        b = nint(dens(i)%theta/dth) + 1
        
        if (dens(i)%rho<=grid(a,b)%rho*const) then

           dens_spher(i)%x = dens(i)%rho*sin(dens(i)%phi)*cos(dens(i)%theta) + center%x
           dens_spher(i)%y = dens(i)%rho*sin(dens(i)%phi)*sin(dens(i)%theta) + center%y
           dens_spher(i)%z = dens(i)%rho*cos(dens(i)%phi) + center%z

           write(3, 12, iostat=ierr) 'ATOM  ', i , dens_spher(i)%atom, &
                dens_spher(i)%resid,' ', 1,'    ', dens_spher(i)%x, &
                dens_spher(i)%y, dens_spher(i)%z
           
        end if
        
     end do
     
  end if
  
  close(3)

  ! finaliza escrita ================================
  
  call ending(back, finish, start, clock_rate) ! Finaliza programa e mostra tempo de processamento  
  
end program s_count