module grd

  implicit none

  real, parameter :: pi = 3.141592654
  character(len=5), parameter :: version = "2.0.0"

  type vet1

     character(len=50) :: atom, resid, ident, code
     real ::  x, y, z
     integer :: n_atom, n_resid

  end type vet1

  type vet2

     real :: x, y, z

  end type vet2


end module grd

program grid_project

  use grd

  logical :: ex, bin, outer, down, rmsd, l_coarse
  logical :: begin, end, skip, eval_skip, help, back
  
  integer :: i, j, k, ierr, n_grid, frame, aux, num, noi1
  integer :: n_index, n_atom, i_atom, bin_out, bin_coarse, a, b
  integer :: fr_in, fr_end, n_skip, start, finish, clock_rate, clock_max
  integer, dimension(50000) :: in_num

  real :: dx, dy, dist, s_grid, r_fit, hour, minu, sec
  real :: x_max, x_min, y_max, y_min, noir, al
  real :: z_max, z_min, peso, rough
  real :: total_lines, lines, progress

  character(len=30) :: coord,  index, atom
  character(len=30), dimension(20) :: get

  type(vet1) :: buff
  type(vet1), dimension(50000) :: store, coarse
  type(vet1), dimension(500000) :: out
  type(vet2), dimension(:,:), allocatable :: mat
  type(vet2), dimension(1001,1001) :: grid

  write(*, *) "  .--.--.                   ,---,                      ,---,. "
  write(*, *) " /  /    '.                '  .' \            ,---.  ,'  .' | "
  write(*, *) "|  :  /`. /          ,--, /  ;    '.         /__./|,---.'   | "
  write(*, *) ";  |  |--`         ,'_ /|:  :       \   ,---.;  ; ||   |   .' "
  write(*, *) "|  :  ;_      .--. |  | ::  |   /\   \ /___/ \  | |:   :  |-, "
  write(*, *) " \  \    `. ,'_ /| :  . ||  :  ' ;.   :\   ;  \ ' |:   |  ;/| "
  write(*, *) "  `----.   \|  ' | |  . .|  |  ;/  \   \\   \  \: ||   :   .' "
  write(*, *) "  __ \  \  ||  | ' |  | |'  :  | \  \ ,' ;   \  ' .|   |  |-, "
  write(*, *) " /  /`--'  /:  | : ;  ; ||  |  '  '--'    \   \   ''   :  ;/| "
  write(*, *) "'--'.     / '  :  `--'   \  :  :           \   `  ;|   |    \ "
  write(*, *) "  `--'---'  :  ,      .-./  | ,'            :   \ ||   :   .' "
  write(*, *) "             `--`----'   `--''               '---: |   | ,'   "
  write(*, *) "                                                   `----'     "
  write(*, *) ""
  write(*, *) ""  
  write(*, *) "                       ** s_grid **"
  write(*, *) ""
  write(*, *) "                   ** VERSION ", version, " **"
  write(*, *) ""
  write(*, *) "             Santos, D. E. S.; Soares, T. A."
  write(*, *) ""
  write(*, *) "Please cite "
  write(*, *)
  write(*, *) "Santos, D. E. S.; Coutinho, K.; Soares, T. A. (2022) Surface "
  write(*, *) "Assessment Grid Evaluation (SuAVE) for Every Surface Curvature"
  write(*, *) "and Cavity Shape. Journal of Chemical Information and Modeling,"
  write(*, *) "v. 62, p. 4690–4701"
  write(*, *)
  write(*, *) "Santos, D. E. S.; Pontes, J. F. S.; Lins, R. D.; Coutinho, K.; "
  write(*, *) "Soares, T. A. (2020) SuAVE: A Tool for Analyzing Curvature-Dependent"
  write(*, *) "Properties in Chemical Interfaces. Journal of Chemical Information "
  write(*, *) "and Modeling, v. 60, p. 473-484."

!
! pegando os arquivos de entrada
!

  back = .false.
  outer = .false.
  bin = .false.
  down = .false.
  coord = 'miss'
  index = 'miss'
  rmsd = .false.
  l_coarse = .false.
  begin = .false.
  end = .false.
  skip = .false.
  help = .false.
  rough = 1.0
  
  do i=1, 20

     call getarg(i, get(i))

  end do

  do i=1, 20

     if (get(i)=='-in')then

        coord = get(i+1)

     end if

     if (get(i)=='-ind')then

        index = get(i+1)

     end if

     if (get(i)=='-bin')then

        bin= .true.
        read(get(i+1), *, iostat=ierr) n_grid
        call erro(ierr, " BIN ")
        
     end if

     if (get(i)=='-outer')then

        outer= .true.
        read(get(i+1), *, iostat=ierr) bin_out
        call erro(ierr, "OUTER")

     end if

     if (get(i)=='-down')then

        down= .true.

     end if

     if (get(i)=='-rmsd')then

        rmsd = .true.

     end if

     if (get(i)=='-coarse')then

        l_coarse = .true.

     end if

     if (get(i)=='-begin')then

        begin = .true.
        read(get(i+1), *, iostat=ierr) fr_in
        call erro(ierr, "BEGIN")
                
     end if

     if (get(i)=='-end')then

        end = .true.
        read(get(i+1), *, iostat=ierr) fr_end
        call erro(ierr, " END ")

     end if

     if (get(i)=='-skip')then

        skip = .true.
        read(get(i+1), *, iostat=ierr) n_skip
        call erro(ierr, " SKIP")

     end if

     if (get(i)=='-rough')then

        read(get(i+1), *, iostat=ierr) rough
        call erro(ierr, "REDUC")

     end if
     
     if (get(i)=='-help')then

        help = .true.

     end if

  end do

  !HELP begins==================================================================
  if (help) then

     write(*, *) ""
     write(*, *) ""
     write(*, *) "s_grid builds a grid per frame throughout a trajectory file." 
     write(*, *) "Its output is useful to verify how accurate is the fitting of" 
     write(*, *) "the calculated grid on the chemical surface.  "
     write(*, *) ""
     write(*, *) "Usage: s_grid -in file.pdb -ind file.ndx"
     write(*, *) ""
     write(*, *) "file.pdb ---- atomic coordinates in PDB format"
     write(*, *) "file.ndx ---- index file containing user-selected atoms used" 
     write(*, *) "to fit the grid points to the chemical surface."
     write(*, *) ""
     write(*, *) "Options:"
     write(*, *) ""
     write(*, *) "-bin             defines the number of rectangular partition "
     write(*, *) "bins along the x- and y-axes"
     write(*, *) ""
     write(*, *) "-outer           automatically selects the "
     write(*, *) "surface/interface outmost atoms to fit the grid points. This" 
     write(*, *) "option overwrites the user-selected index files."
     write(*, *) ""
     write(*, *) "-rmsd            calculates the RMSD between the fitted" 
     write(*, *) "grid and the selected atoms in the index files. This" 
     write(*, *) "estimates how precisely is the grid surface fitted to the" 
     write(*, *) "chemical surface throughout the trajectory file."
     write(*, *) ""
     write(*, *) "-coarse          generates a coarse grid over the surface" 
     write(*, *) "index atoms from which a finer grid will be generated. This" 
     write(*, *) "is recommended for surfaces defined by atoms which greatly" 
     write(*, *) "fluctuate throughout the trajectory. "
     write(*, *) ""
     write(*, *) "-begin           first frame to use in the calculations"
     write(*, *) ""
     write(*, *) "-end             last frame to use in the calculations"
     write(*, *) ""
     write(*, *) "-skip            number of trajectory frames to be skipped" 
     write(*, *) "during the analysis "
     write(*, *) ""
     Write(*, *) "-down            used to generate a PDB file with the last "
     Write(*, *) "external lowest fitting grid. This option is only used when" 
     Write(*, *) "followed by option –outer."  
     Write(*, *) ""
     write(*, *) "-rough           percentage of the original surface roughness"
     write(*, *) "the user wants to keep (1.0 by default, meaning 100%)"
     write(*, *) ""
     Write(*, *) "-help            prints HELP information and quit."
     stop
  end if
  !HELP end=====================================================================

  if (coord=='miss')then

     write(*, *)
     write(*, *)'PDB file is missing'
     write(*, *)
     stop

  end if

  if (.not.outer)then

     if (index=='miss')then

        write(*, *)
        write(*, *)'Index file is missing'
        write(*, *)
        stop

     end if

  end if
  
  !==Lendo os arquivos de index=========================
  if (.not.outer)then
     
     open(2, file=index, status='old', iostat=ierr)
     
     if(ierr /= 0) then

        write(*, *)
        write(*, *) 'Unable to open file', index
        write(*, *)
        stop
        
     end if

     ! Inicio da leitura do index ============================
     
     n_index = 1
     
     do while (ierr==0)
        
        read(2, *, iostat=ierr) aux
        
        if (ierr == 0)then
           
           in_num(n_index) = aux
           n_index = n_index + 1
           
        end if
        
     end do
     
     n_atom  = n_index - 1
     
     if(ierr>0) then

        write(*, *)
        write(*, *) 'Problem by reading ', index
        write(*, *)
        stop
        
     end if
     
     do i=n_index, 50000
        
        in_num(i)=-1
        
     end do

     write(*, *)
     write(*, *) index, " input file is OK"
     
     close(2)

  end if
  
  !=================definindo frames para inicio e fim========
  !=================definindo skip============================
  
  frame = 0
  
  if (.not.begin) then
     
     fr_in = 1

  end if
  
  if (.not.end) then

     fr_end = fr_in + 2 !==para que fr_end seja maior que frame

  end if

  if (.not.skip) then

     n_skip = 0
     
  end if
  !============================================================
    
  !=======================================================
  ! calculando o espaçamento

  if (.not.outer)then

     bin_coarse = nint(sqrt(real(n_index-1)) - 1)

     if (.not.bin)then
        
        n_grid = nint(sqrt(real(n_index-1)) - 1)
        write(*, *)
        write(*, *) 'STD_BIN = ', n_grid
        
     end if

  else

     bin_coarse = bin_out

     if (.not.bin)then

        n_grid = bin_out
        write(*, *)
        write(*, *) 'STD_BIN = ', n_grid

     end if

  end if
  !==================================================
  
  !==================================================
  ! definição da matriz de pontos externos
  
  if (outer)then

     allocate (mat(bin_out+1, bin_out+1), stat = ierr)
     if (ierr /= 0) stop 'Not enough memory to initialize mat matrix'
     
     do i=1, bin_out + 1
        
        do j=1, bin_out + 1

           if (down)then
              
              mat(i,j)%z = 1000

           else

              mat(i,j)%z = -1000

           end if

        end do
        
     end do

  end if

  !==================================================

  open(1, file=coord, status='old', iostat=ierr)

  if(ierr /=0) then

     write(*, *)
     write(*, *) 'Unable to open file ', coord
     write(*, *)
     stop

  endif

  ! adquirindo numero total de linhas =====
  
  total_lines = 0

  do while (ierr >= 0)
     
     read(1, *, iostat=ierr)
     total_lines = total_lines + 1
     
  end do

  !========================================
  
  rewind(1) ! reinicializando o arquivo 1

  if (rmsd)then

     inquire(file='rmsd.xvg', exist=ex)

     if (ex) then

        call execute_command_line("mv rmsd.xvg rmsd_p.xvg")
        back = .true.

     end if

     open(3, file='rmsd.xvg', status='new', iostat=ierr)

     if(ierr /= 0) then

        write(*, *)
        write(*, *) 'Unable to open file rmsd.xvg'
        write(*, *)
        stop

     endif

     write(3, '(a7, a5)') "#SuAVE ", version
     write(3, '(a14)') '#Command Line:'
     write(3, '(a9)', advance='no') '#s_grid  '
     write(3, *) (trim(get(i)),"  ", i=1, 20)
     write(3, *) '@    title "RMSD X Frame"'
     write(3, *) '@    xaxis  label "#Frame"'
     write(3, *) '@    yaxis  label "RMSD [nm]"'

  end if

  !==========================================================
  inquire(file='grid.pdb', exist=ex)

  if (ex) then

     call execute_command_line("mv grid.pdb grid_p.pdb")
     back = .true.
     
  end if

  open(4, file='grid.pdb', status='new', iostat=ierr)

  if(ierr /= 0) then

     write(*, *)
     write(*, *) 'Unable to open file grid.pdb'
     write(*, *)
     stop

  end if

  call system_clock(start, clock_rate, clock_max)
  
  ! Leitura do arquivo .pdb
  
  i_atom = 0
  frame = 0
  ierr = 0
  n_index = 1
  x_min = 1000
  y_min = 1000
  z_min = 1000
  x_max = 0
  y_max = 0
  z_max = 0
  line = 0
  progress = 0.01
  
  write(*, *)
  write(*, '(a17)', advance='no') ' Progress : 00.0%'
 
  
  do while (ierr>=0)
  
     read(1, 12, iostat=ierr) atom, buff%n_atom, buff%atom, &
          buff%resid, buff%ident, buff%n_resid, buff%code, buff%x, &
          buff%y, buff%z

12   format(a6, i5.1, a5, a5, a1, i4.1, a4, 3f8.3)

     lines = lines + 1
     
     if (atom.eq.'ATOM  ') then

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
           z_max = max(z_max, buff%z)
           z_min = min(z_min, buff%z)

        else
           
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
        
     else

        if (n_index > 1)then

           frame  = frame + 1

           eval_skip = ((mod((frame-fr_in),(n_skip+1))==0))
           
           if ((frame>fr_in-1).and.(frame<fr_end+1).and.(eval_skip)) then
              
              if (outer)then
                 
                 dx = (x_max - x_min)/bin_out
                 dy = (y_max - y_min)/bin_out
                 
                 do k=1, n_index-1
                    
                    i = nint((out(k)%x-x_min)/dx + 1)
                    
                    j = nint((out(k)%y-y_min)/dy + 1)
                    
                    if (down) then
                       
                       if (out(k)%z<mat(i,j)%z)then
                          
                          mat(i,j)%z = out(k)%z
                          mat(i,j)%x = out(k)%x
                          mat(i,j)%y = out(k)%y
                          
                       end if
                       
                    else 
                       
                       if (out(k)%z>mat(i,j)%z)then
                          
                          mat(i,j)%z = out(k)%z
                          mat(i,j)%x = out(k)%x
                          mat(i,j)%y = out(k)%y
                          
                       end if
                       
                    end if
                    
                 end do
                 
                 num = n_index
                 n_index = 1
                 
                 do i=1,bin_out+1
                    
                    do j=1,bin_out+1
                       
                       if (down)then
                          
                          if (mat(i,j)%z<1000)then
                             
                             store(n_index)%x = mat(i,j)%x
                             store(n_index)%y = mat(i,j)%y
                             store(n_index)%z = mat(i,j)%z
                             
                             n_index = n_index + 1
                             
                          end if
                          
                       else
                          
                          if (mat(i,j)%z>-1000)then
                             
                             store(n_index)%x = mat(i,j)%x
                             store(n_index)%y = mat(i,j)%y
                             store(n_index)%z = mat(i,j)%z
                             
                             n_index = n_index + 1
                             
                          end if
                          
                       end if
                       
                    end do
                    
                 end do
                 
              end if! if outer
              
              noi1 = n_index - 1
              
              if (l_coarse) then
                 
                 !
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

3             format(a8, f7.2, a2, f7.2, a2, f7.2, a35)
              write(4, '(a26)') 'REMARK Generated by s_grid'
              write(4, '(a34)') 'TITLE     Rectangular Regular GRID'
              write(4, '(a22)') 'REMARK trajectory file'
              write(4, 3) 'CRYST1  ', x_max-x_min, '  ', y_max-y_min, '  ', z_max-z_min, '  90.00  90.00  90.00 P  1       1 '
              write(4, '(a6, i4)') 'MODEL ', frame

              do i=1, n_grid+1
                 
                 do j=1, n_grid+1
                    
                    write(4, 12, iostat=ierr) 'ATOM  ', (i-1)*(n_grid+1) + j, '  DOT', &
                         '  DOT',' ', 1,'    ', grid(i,j)%x, &
                         grid(i,j)%y, grid(i,j)%z
                    
                 end do
                 
              end do
              
              write(4, '(a3)') 'TER'
              write(4, '(a6)') 'ENDMDL'
              
           !fim da escrita da trajetoria

           end if !======((frame<fr_in-1).and.(frame>fr_end+1))
              
           n_atom = n_index-1
           n_index = 1
           i_atom = 0
           x_min = 1000
           y_min = 1000
           x_max = 0
           y_max = 0

           !====garante que fr_end sempre seja maior que frame ===
           !====caso essa variável não tenha sido fixada==========
           if (.not.end)then
              
              fr_end = frame + 1
              
           end if
           !======================================================

           
           if (outer)then
              
              do i=1, bin_out + 1
                 
                 do j=1, bin_out + 1

                    if (down)then

                       mat(i,j)%z = 1000

                    else

                       mat(i,j)%z = -1000

                    end if
                    
                 end do
                 
              end do

           end if
           
        end if
        
     end if 

     if ((lines/total_lines) > progress) then ! imprimindo barra de avanco
        
        write(*,'(1a1, a12, f5.1, a1, $)') char(13), ' Progress : ', progress*100, '%' 
        progress = progress + 0.01
        
     end if

     
  end do
  
  ! Fim da estruturação

  close(1)
  close(4)

  if (rmsd)then
     
     close(3)
     
  end if
  
  call system_clock(finish, clock_rate, clock_max)

  
  ! desenho dos pontos usados no ajuste
  
  inquire(file='adjust.pdb', exist=ex)
  
  if (ex) then
     
     call execute_command_line("mv adjust.pdb adjust_p.pdb")
     back = .true.
     
  end if
  
  open(3, file='adjust.pdb', status='new', iostat=ierr)
  
  if(ierr /= 0) then

     write(*, *)
     write(*, *) 'Unable to open file ajuste.pdb'
     write(*, *)
     stop

  end if

  write(3, *) 'TITLE     Rectangular Regular GRID'
  write(3, *) 'REMARK    THIS IS A SIMULATION BOX'
  write(3, *) 'MODEL        1'

  do i=1, noi1
     
     write(3, 12, iostat=ierr) 'ATOM  ', i, '  DOT', &
          '  DOT',' ', 1,'    ', store(i)%x, &
          store(i)%y, store(i)%z

  end do

  close(3)
  
  if (l_coarse) then
     
     inquire(file='coarse.pdb', exist=ex)
     
     if (ex) then
        
        call execute_command_line("mv coarse.pdb coarse_p.pdb")
        back = .true.
        
     end if
     
     open(3, file='coarse.pdb', status='new', iostat=ierr)
     
     if(ierr /= 0) then
        
        write(*, *)
        write(*, *) 'Unable to open file coarse.pdb'
        write(*, *)
        stop
        
     end if
     
     write(3, *) 'TITLE     Rectangular Regular GRID'
     write(3, *) 'REMARK    THIS IS A SIMULATION BOX'
     write(3, *) 'MODEL        1'
     
     do i=1, (bin_coarse + 1)*(bin_coarse + 1)
        
        write(3, 12, iostat=ierr) 'ATOM  ', i, '  DOT', &
             '  DOT',' ', 1,'    ', coarse(i)%x, &
             coarse(i)%y, coarse(i)%z
        
     end do

  end if
  
  write(*, *)
  write(*, *) 
  write(*, *) "Finished"
  write(*, *)
  if (back) write(*, '(a31)') ' Previous files were backed up!'
  
! Fim. Arquivo escrito

  close(3)

  hour = real(finish-start)/(3600*real(clock_rate))
  minu = (hour - int(hour))*60
  sec = (minu-int(minu))*60

  if (back)  write(*, *)
  write(*, 1) " Processing time : ", int(hour), " h", int(minu), " min", int(sec), " sec"
1 format(a19, i6.3, a2, i6.3, a4, i6.3, a4)
  write(*, *)
  
  if (outer)then
     
     deallocate(mat)

  end if
     
end program grid_project


subroutine erro(i, name)

  integer, intent(in) :: i
  character(len=5), intent(in) :: name
  character(len=60) :: aux

  if (i/=0) then

     aux = '    Problem by passing '//name//' parameter'
     write(*, *)
     write(*, *)"==================ERROR==================="
     write(*, *) aux
     write(*, *)"==================ERROR==================="
     write(*, *)
     stop

  end if

end subroutine erro

subroutine param(x_max, x_min, y_max, y_min, num, r_fit, al, rough)

  real, intent(in) :: x_max, x_min, y_max, y_min, rough
  real :: r_fit, al
  integer, intent(in) :: num

  r_fit = 3*sqrt((x_max - x_min)**2 + (y_max - y_min)**2)/sqrt(real(num)-1)

  al = (num-1)*100/((x_max - x_min)*(y_max - y_min))
  al = exp(0.4247*rough*log(al)-1.3501/rough)

end subroutine param
