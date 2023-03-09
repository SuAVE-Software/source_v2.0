module top
 
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


end module top

program topog
  
  use top

  logical :: ex, bin, outer, p_grid, rmsd, l_coarse
  logical :: begin, end, skip, eval_skip, range, help, back
  
  integer :: i, j, k, ierr, n_grid, frame, aux, bin_out, noi1, noi2
  integer :: n_index, num, num2, i_atom, bin_coarse, a, b
  integer :: fr_in, fr_end, n_skip, tot_frame
  integer :: start, finish, clock_rate, clock_max
  integer, dimension(50000) :: in_num, in_num2

  real :: dist_z, hour, minu, sec, r_fit, noir, desv
  real :: x_max, x_min, y_max, y_min, dx, dy, dist, gx, gy, s_grid, peso
  real :: la, lb, lc, minv, maxv, del, al, aver, aver2, rough
  real :: total_lines, lines, progress
  real, dimension(:, :), allocatable :: r_xpm

  character(len=30) :: coord,  ind, atom, ind2
  character(len=1), dimension(:,:), allocatable :: xpm
  character(len=30), dimension(30) :: get

  type(vet1), dimension(50000) :: store, store2, coarse, coarse2
  type(vet1) :: buff
  type(vet1), dimension(500000) :: out
  type(vet2), dimension(:,:), allocatable :: mat1, mat2
  type(vet2), dimension(1001,1001) :: grid, grid2

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
  write(*, *) "                       ** s_topog **"
  write(*, *) ""
  write(*, *) "                   ** VERSION ", version, " **"
  write(*, *) ""
  write(*, *) "              Santos, D. E. S.; Soares, T. A."
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
  range = .false.
  outer = .false.
  bin = .false.
  p_grid = .false.
  coord = 'miss'
  ind = 'miss'
  ind2 = 'miss'
  rmsd = .false.
  l_coarse = .false.
  begin = .false.
  end = .false.
  skip = .false.
  help = .false.
  rough = 1.0
  
  do i=1, 30

     call getarg(i, get(i))

  end do

  do i=1, 30

     if (get(i)=='-in')then

        coord = get(i+1)

     end if

     if (get(i)=='-ind1')then

        ind = get(i+1)

     end if

     if (get(i)=='-ind2')then

        ind2 = get(i+1)

     end if
     
     if (get(i)=='-bin')then

        bin = .true.
        read(get(i+1), *, iostat=ierr) n_grid
        call erro(ierr, " BIN ")

     end if

     if (get(i)=='-outer')then

        outer = .true.
        read(get(i+1), *, iostat=ierr) bin_out
        call erro(ierr, "OUTER")

     end if

     if(get(i)=='-grid')then

        p_grid = .true.

     end if

     if (get(i)=='-rmsd')then

        rmsd = .true.

     end if

     if (get(i)=='-coarse')then

        l_coarse = .true.

     end if

     if (get(i)=='-range')then

        range = .true.

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

  !HELP begins ===============================================================
  if (help) then

     write(*, *) ""
     write(*, *) ""
     write(*, *) "s_topog calculates the level surface of the chemical interface" 
     write(*, *) ""
     write(*, *) "Usage: s_topog -in file.pdb -ind1 file1.ndx -ind2 file2.ndx"
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
     write(*, *) "-grid            generates a PDB file containing the grid" 
     write(*, *) "points used in the fitting for the last frame in the "
     write(*, *) "trajectory file."
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
     write(*, *) "-range           defines the range specified in the XPM file"
     write(*, *) ""
     write(*, *) "-rough           percentage of the original surface roughness"
     write(*, *) "the user wants to keep (1.0 by default, meaning 100%)"
     write(*, *) ""     
     write(*, *) "-help            prints HELP information and quits"
     stop
     
  end if
  !HELP ends ==================================================================

  if (coord=='miss')then

     write(*, *)
     write(*, *)'PDB file is missing'
     write(*, *)
     stop

  end if

  if (.not.outer)then
  
     if (ind=='miss')then
        
        write(*, *)
        write(*, *)'First index file is missing'
        write(*, *)
        stop
        
     end if
     
     if (ind2=='miss')then
        
        write(*, *)
        write(*, *)'Second index file is missing'
        write(*, *)
        stop
        
     end if

  end if

  !==Lendo os arquivos de index=========================
  if (.not.outer)then
     
     open(2, file=ind, status='old', iostat=ierr)
     
     if(ierr /= 0) then

        write(*, *)
        write(*, *) 'Unable to open file ', ind
        write(*, *)
        stop
        
     endif
     
     open(3, file=ind2, status='old', iostat=ierr)
     
     if(ierr/=0) then

        write(*, *)
        write(*, *) 'Unable to open file', ind2
        write(*, *)
        stop
        
     end if

     ! Inicio da leitura do index======================

     ierr = 0
     n_index = 1
     
     do while (ierr==0)
        
        read(2, *, iostat=ierr) aux
        
        if (ierr==0) then
           
           in_num(n_index) = aux
           n_index = n_index + 1
           
        end if
        
     end do
     
     if(ierr>0) then

        write(*, *)
        write(*, *) 'Problem by reading ', ind
        write(*, *)
        stop
        
     end if
     
     do i=n_index, 50000
        
        in_num(i) = -1
        
     end do

     write(*, *)
     write(*, *) ind, " input file is OK"
     
     ! Fim da leitura do index ==========================

     close(2)
     
     ! Inicio da leitura do ind2 ========================

     ierr = 0
     n_index = 1
     
     do while (ierr==0)
        
        read(3, *, iostat=ierr) aux
        
        if (ierr==0)then
           
           in_num2(n_index) = aux
           n_index = n_index + 1
           
        end if
        
     end do
     
     if(ierr>0) then

        write(*, *)
        write(*, *) 'Problem by reading ', ind2
        write(*, *)
        stop
        
     end if
     
     do i=n_index, 50000
        
        in_num2(i) = -1
        
     end do

     write(*, *)
     write(*, *) ind2, " input file is OK"
     
     ! Fim da leitura do ind2 ================================

     close(3)
     
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

  !============================================================
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
  !===================================================

  !===================================================
  ! definição da matriz de pontos externos

  if (outer)then

     allocate (mat1(bin_out+1, bin_out+1), stat = ierr)
     if (ierr /= 0) stop 'Not enough memory to initialize mat matrix'

     allocate (mat2(bin_out+1, bin_out+1), stat = ierr)
     if (ierr /= 0) stop 'Not enough memory to initialize mat matrix'

     do i=1, bin_out + 1

        do j=1, bin_out + 1

           mat1(i,j)%z = -1000
           mat2(i,j)%z = 1000

        end do

     end do

  end if

  !=============================================================

  open(1, file=coord, status='old', iostat=ierr)

  if(ierr /=0) then

     write(*, *)
     write(*, *) 'Unable to open file',coord
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
     write(3, '(a10)', advance='no') '#s_topog  '
     write(3, *) (trim(get(i)),"  ", i=1, 30)
     write(3, *) '@    title "RMSD X Frame"'
     write(3, *) '@    xaxis  label "#Frame"'
     write(3, *) '@    yaxis  label "RMSD [nm]"'

  end if
  
  allocate (xpm(n_grid, n_grid), stat = ierr)
  if (ierr /= 0) stop 'Not enough memory to initialize xpm matrix'

  allocate (r_xpm(n_grid, n_grid), stat = ierr)
  if (ierr /= 0) stop 'Not enough memory to initialize xpm matrix'

  do i=1, n_grid

     do j=1, n_grid

        r_xpm(i,j) = 0
        
     end do

  end do

  call system_clock(start, clock_rate, clock_max)

! Leitura do arquivo .pdb

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
  aver = 0
  aver2 = 0
  minv = 100000
  maxv = -100000
  lines = 0
  progress = 0.01
  
  write(*, *)
  write(*, '(a17)', advance='no') ' Progress : 00.0%'
  
  do while (ierr >= 0)

     read(1, 12, iostat=ierr) atom, buff%n_atom, buff%atom, &
          buff%resid, buff%ident, buff%n_resid, buff%code, buff%x, &
          buff%y, buff%z
     
12   format(a6, i5.1, a5, a5, a1, i4.1, a4, 3f8.3)

     lines = lines + 1

     if ((atom.eq.'TER  ')) cycle
     
     if (atom.eq.'ATOM  ') then
        
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
        
     else

        if (n_index > 1)then

           frame  = frame + 1

           eval_skip = ((mod((frame-fr_in),(n_skip+1))==0))
           
           if ((frame>fr_in-1).and.(frame<fr_end+1).and.(eval_skip)) then

              tot_frame = tot_frame + 1
           
              if (outer)then
                 
                 dx = (x_max - x_min)/bin_out
                 dy = (y_max - y_min)/bin_out
                 
                 do k=1, n_index-1
                    
                    i = nint((out(k)%x-x_min)/dx + 1)
                    
                    j = nint((out(k)%y-y_min)/dy + 1)
                    
                    if (out(k)%z>mat1(i,j)%z)then
                       
                       mat1(i,j)%z = out(k)%z
                       mat1(i,j)%x = out(k)%x
                       mat1(i,j)%y = out(k)%y
                       
                    end if
                    
                    if (out(k)%z<mat2(i,j)%z)then
                       
                       mat2(i,j)%z = out(k)%z
                       mat2(i,j)%x = out(k)%x
                       mat2(i,j)%y = out(k)%y
                       
                    end if
                    
                 end do
                 
                 n_index = 1
                 num = 1
                 num2 = 1
                 
                 do i=1,bin_out+1
                    
                    do j=1,bin_out+1
                       
                       if (mat1(i,j)%z>-1000)then
                          
                          store(num)%x = mat1(i,j)%x
                          store(num)%y = mat1(i,j)%y
                          store(num)%z = mat1(i,j)%z
                          
                          num = num + 1
                          n_index = n_index + 1
                          
                       end if
                       
                       if (mat2(i,j)%z<1000)then
                          
                          store2(num2)%x = mat2(i,j)%x
                          store2(num2)%y = mat2(i,j)%y
                          store2(num2)%z = mat2(i,j)%z
                          
                          num2 = num2 + 1
                          n_index = n_index + 1
                          
                       end if
                       
                    end do
                    
                 end do
                 
              end if ! if outer
              
              noi1 = num - 1
              noi2 = num2 - 1
              
              if (l_coarse) then
                 
                 ! Estruturação do primeiro grid coarse
                 
                 dx = (x_max - x_min)/bin_coarse
                 dy = (y_max - y_min)/bin_coarse
                 
                 !adaptação ao r_fit e ao alfa

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
                 
                 num = (bin_coarse + 1)**2 + 1
                 
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
              
              ! estruturação do primeiro grid de alta resolução
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
                       
                       if (dist<r_fit*r_fit) then
                          
                          peso = exp(-(dist*al*al)/pi)
                          s_grid = s_grid + peso
                          grid(i,j)%z = grid(i,j)%z + peso*coarse(k)%z
                          
                       end if
                       
                    end do
                    
                    grid(i,j)%z = grid(i,j)%z/s_grid
                    
                 end do
                 
              end do
              !$OMP end parallel do
              
              ! estruturação do segundo grid de alta resolução
              
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
           
              ! calcula topografia ===============================================
              !Topografia da superfície média
              lc = 0
              aver = 0
              aver2 = 0

              do i=2, n_grid+1
                 
                 do j=2, n_grid+1
                    
                    la = (grid(i-1,j-1)%z +grid(i-1,j)%z + grid(i,j-1)%z + grid(i,j)%z)/4
                    lb = (grid2(i-1,j-1)%z + grid2(i-1,j)%z + grid2(i,j-1)%z + grid2(i,j)%z)/4 
                    lc = (la+lb)
                    
                    r_xpm(i-1,j-1) = r_xpm(i-1, j-1) +lc/10

                    aver = aver + lc/10/(n_grid*n_grid)
                    aver2 = aver2 + lc*lc/100/(n_grid*n_grid)
         
                 end do
                 
              end do
              
              desv = sqrt(aver2 - aver*aver)
              minv = min(aver-2*desv, minv)
              maxv = max(aver+2*desv, maxv)
              ! calcula topografia ===============================================

              if (rmsd)then
                 
                 noir = 0
                 
                 do i=1, noi1
                    
                    a = nint((store(i)%x-x_min)/dx) + 1
                    b = nint((store(i)%y-y_min)/dy) + 1
                    dist_z = store(i)%z-grid(a,b)%z
                    noir = noir + dist_z*dist_z/(noi1+noi2)
                    
                 end do
                 
                 do i=1, noi2
                    
                    a = nint((store2(i)%x-x_min)/dx) + 1
                    b = nint((store2(i)%y-y_min)/dy) + 1
                    dist_z = store2(i)%z-grid2(a,b)%z
                    noir = noir + dist_z*dist_z/(noi1+noi2)
                    
                 end do
                 
                 noir = sqrt(noir)
                 
                 write(3, *) noir/10
                 
              end if
              
              ! Fim do cálculo 

           end if !======((frame<fr_in-1).and.(frame>fr_end+1))
           
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

           !====garante que fr_end sempre seja maior que frame ===
           !====caso essa variável não tenha sido fixada==========
           if (.not.end)then
              
              fr_end = frame + 1
              
           end if
           !======================================================
           
           if (outer)then

              do i=1, bin_out + 1

                 do j=1, bin_out + 1

                    mat1(i,j)%z = -1000
                    mat2(i,j)%z = 1000

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
  
  call system_clock(finish, clock_rate, clock_max)

  close(1)
  close(3)
  
  ! Fim da leitura do arquivo .pdb
  
  ! escrita do arquivo XPM

  r_xpm = r_xpm/tot_frame

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
  
  !=============================================================

  inquire(file='topog.xpm', exist=ex)

  if (ex) then

     call execute_command_line("mv topog.xpm topog_p.xpm")
     back = .true.

  end if

  open(4, file='topog.xpm', status='new', iostat=ierr)

  if(ierr /= 0) then

     write(*, *)
     write(*, *) 'Unable to open file topog.xpm'
     write(*, *)
     stop

  endif

  write(4, '(a7, a5)') "#SuAVE ", version
  write(4, '(a14)') '#Command Line:'
  write(4, '(a10)', advance='no') '#s_thick  '
  write(4, *) (trim(get(i)),"  ", i=1, 30)
  write(4, *) '/* XPM */'
  write(4, *)'/* This matrix was generated by s_topog.f90 */'
  write(4, *)'/* title:   "Topog" */'
  write(4, *)'/* x-label: "x axis [nm]" */'
  write(4, *)'/* y-label: "y axis [nm]" */'
  write(4, *)'/* type:    "Continuous" */'
  write(4, *)'static char * gv_xpm[] = {'
  write(4, *) '"',n_grid,n_grid,' 7 1",'
  
  do i=1, n_grid

     do j=1, n_grid

       r_xpm(i,j) = r_xpm(i,j) - minv

     end do

  end do

  del = (maxv - minv)/6
7 format(a20, f4.1, a5)
  write(4, 7) '"A  c #ffffff " /* "', minv+000.0,'" */,'
  write(4, 7) '"B  c #87cefa " /* "',minv+del,'" */,'
  write(4, 7) '"C  c #00bfff " /* "',minv+2*del,'" */,'
  write(4, 7) '"D  c #1e90ff " /* "',minv+3*del,'" */,'
  write(4, 7) '"E  c #4169e1 " /* "',minv+4*del,'" */,'
  write(4, 7) '"F  c #0000ff " /* "',minv+5*del,'" */,'
  write(4, 7) '"G  c #00008b " /* "',minv+6*del,'" */,'

  do i=1, n_grid

     do j=1, n_grid

        if (r_xpm(i,j)<del)then
           xpm(i,j)= 'A'
        else
           if (r_xpm(i,j)<2*del)then
              xpm(i,j) = 'B'
           else
              if (r_xpm(i,j)<3*del)then
                 xpm(i,j) = 'C'
              else
                 if (r_xpm(i,j)<4*del) then
                    xpm(i,j) = 'D'
                 else
                    if (r_xpm(i,j)<5*del) then
                       xpm(i,j) = 'E'
                    else
                       if (r_xpm(i,j)<6*del) then
                          xpm(i,j) = 'F'
                       else
                          xpm(i,j) = 'G'

                       end if
                       
                    end if
                    
                 end if
                 
              end if
              
           end if
           
        end if

     end do

  end do

  k = 1

3 format(a10)
4 format(f6.2, a1)

  do i=1, int(n_grid/20)+1

     write(4, 3, advance='no') '/* y-axis:'
     j = 1

     do while ((j<=20).and.(k<=n_grid))
        
        write(4, 4, advance='no') ((k-1/2)*dy + gy)/10, ' '
        k = k+1
        j = j+1

     end do
     write(4,*) '*/'

  end do

  k = 1
  do i=1, int(n_grid/20)+1

     write(4, 3, advance='no') '/* x-axis:'
     j = 1

     do while ((j<=20).and.(k<=n_grid))

        write(4, 4, advance='no') ((k-1/2)*dx + gx)/10, ' '
        k = k+1
        j = j+1

     end do
     write(4,*) '*/'

  end do

5 format(a1)
  do j=n_grid, 1, -1
    
     write(4, 5, advance='no') '"'
     
     do i=1, n_grid
        
        write(4, 5, advance='no') xpm(i,j)

     end do

     write(4, 5, advance='no')'"'
     write(4, 5) ','
     
  end do

  close(4)
  
! Final da escrita do XPM
  
  if (p_grid)then

     inquire(file='grid1.pdb', exist=ex)

     if (ex) then

        call execute_command_line("mv grid1.pdb grid1_p.pdb")
        back = .true.

     end if

     open(3, file='grid1.pdb', status='new', iostat=ierr)

     if(ierr /= 0) then

        write(*, *)
        write(*, *) 'Unable to open file grid1.pdb'
        write(*, *)
        stop

     end if

     write(3, *) 'TITLE     Rectangular Regular GRID'
     write(3, *) 'REMARK    THIS IS A SIMULATION BOX'
     write(3, *) 'MODEL        1'

13   format(a6, i5.1, a5, a5, a1, i4.1, a4, 3f8.3, f5.1)
     
     do i=1, n_grid

        do j=1, n_grid

           write(3, 13, iostat=ierr) 'ATOM  ', (i-1)*(n_grid+1) + j, '  DOT', &
                '  DOT',' ', 1,'    ', grid(i,j)%x, &
                grid(i,j)%y, grid(i,j)%z, r_xpm(i,j)

        end do

     end do

     close(3)
     ! first grid created
     
     inquire(file='grid2.pdb', exist=ex)

     if (ex) then

        call execute_command_line("mv grid2.pdb grid2_p.pdb")
        back = .true.

     end if

     open(3, file='grid2.pdb', status='new', iostat=ierr)

     if(ierr /= 0) then

        write(*, *)
        write(*, *) 'Unable to open file grid2.pdb'
        write(*, *)
        stop

     end if

     write(3, *) 'TITLE     Rectangular Regular GRID'
     write(3, *) 'REMARK    THIS IS A SIMULATION BOX'
     write(3, *) 'MODEL        1'

     do i=1, n_grid

        do j=1, n_grid

           write(3, 13, iostat=ierr) 'ATOM  ', (i-1)*(n_grid+1) + j, '  DOT', &
                '  DOT',' ', 1,'    ', grid2(i,j)%x, &
                grid2(i,j)%y, grid2(i,j)%z, r_xpm(i,j)

        end do

     end do

     close(3)
     ! second grid created
     
  end if
  
  deallocate(r_xpm)
  deallocate(xpm)

  if (outer)then

     deallocate(mat1)
     deallocate(mat2)

  end if

  write(*, *)
  write(*, *) 
  write(*, *) "Finished"
  write(*, *)
  if (back) write(*, '(a31)') ' Previous files were backed up!'

  hour = real(finish-start)/(3600*real(clock_rate))
  minu = (hour - int(hour))*60
  sec = (minu-int(minu))*60

  if (back) write(*, *)
  write(*, 1) " Processing time : ", int(hour), " h", int(minu), " min", int(sec), " sec"
1 format(a19, i6.3, a2, i6.3, a4, i6.3, a4)
  write(*, *)
  
end program topog


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
