module area

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


end module area

program ar_grid

  use area

  logical :: ex, bin, outer, p_grid, rmsd, l_coarse, begin, end, skip
  logical :: eval_skip, lipid, help, back

  integer :: i, j, k, l, ierr, n_grid, frame, aux, bin_out, noi1, noi2
  integer :: n_index, n_lipid, num, num2, i_atom, bin_coarse, a, b
  integer :: fr_in, fr_end, n_skip, tot_frame
  integer :: start, finish, clock_rate, clock_max
  integer, dimension(50000) :: in_num, in_num2

  real :: dist_x, dist_y, dist_z, s_soma, hour, minu, sec, noir, minv
  real :: x_max, x_min, y_max, y_min, dx, dy, dist, s_grid, al, peso
  real :: s_area, la, lb, lc, r_fit, aux2, gridx, gridy, tempo, rough
  real :: total_lines, lines, progress
  real, dimension(: , :), allocatable :: r_xpm1, r_xpm2
  
  character(len=30) :: coord,  ind, atom, ind2
  character(len=30), dimension(30) :: get

  type(vet1), dimension(50000) :: store, store2, coarse, coarse2
  type(vet1) :: buff
  type(vet1), dimension(500000) :: out
  type(vet2), dimension(:,:), allocatable :: mat1, mat2
  type(vet2), dimension(1001,1001) :: grid, grid2, grid3

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
  write(*, *) "                       ** s_area **"
  write(*, *) ""
  write(*, *) "                   ** VERSION ", version, " **"
  write(*, *) ""
  write(*, *) "             Santos, D. E. S.; Soares, T. A."
  write(*, *) ""
  write(*, *) "Please cite SuAVE: A Tool for Analyzing Curvature-Dependent"
  write(*, *) "Properties in Chemical Interfaces (2020) Denys E. S. Santos,"
  write(*, *) "Frederico J. S. Pontes, Roberto D. Lins, Kaline Coutinho,"
  write(*, *) "Thereza A. Soares. J. Chem. Inf. Model., v. 60(2), p. 473-484."

!
! pegando os arquivos de entrada =======================================
!

  back = .false.
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
  lipid = .false.
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

     if (get(i)=='-grid')then

        p_grid = .true.

     end if

     if (get(i)=='-rmsd')then

        rmsd = .true.

     end if
     
     if (get(i)=='-coarse')then

        l_coarse = .true.

     end if

     if (get(i)=='-begin')then

        begin = .true.
        read(get(i+1),*, iostat=ierr) fr_in
        call erro(ierr, "BEGIN")
        
     end if

     if (get(i)=='-end')then

        end = .true.
        read(get(i+1), *, iostat=ierr) fr_end
        call erro(ierr, " END ")

     end if

     if (get(i)=='-skip')then

        skip = .true.
        read(get(i+1), *, iostat=ierr)n_skip
        call erro(ierr, " SKIP")

     end if
     
     if (get(i)=='-lipid')then
        
        lipid = .true.
        read(get(i+1), *, iostat=ierr) n_lipid
        call erro(ierr, "LIPID")

     end if

     if (get(i)=='-rough')then

        read(get(i+1), *, iostat=ierr) rough
        call erro(ierr, "REDUC")

     end if
     
     if (get(i)=='-help')then

        help = .true.

     end if

  end do
  
  !HELP begins ==================================================================
    if (help)then
     
     write(*, *) ""
     write(*, *) ""
     write(*, *) "s_area calculates the time-dependent area per molecule of the" 
     write(*, *) "interface."
     write(*, *) ""
     write(*, *) "Usage: s_area -in file.pdb -ind1 file1.ndx -ind2 file2.ndx"
     write(*, *) "              -lipid N_lipids"
     write(*, *) ""
     write(*, *) "file.pdb ---- atomic coordinates in PDB format"
     write(*, *) "file.ndx ---- index file containing user-selected atoms used" 
     write(*, *) "to fit the grid points to the chemical surface."
     write(*, *) "N_lipids ---- number of lipids composing the leaflet"
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
     write(*, *) "-rough           percentage of the original surface roughness"
     write(*, *) "the user wants to keep (1.0 by default, meaning 100%)"
     write(*, *) ""
     write(*, *) "-help            prints HELP information and quits"
     stop
  end if
  !HELP ends=====================================================================
  
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

  if (.not.lipid) then

     write(*, *)
     write(*, *) 'Number of lipids is missing'
     write(*, *)
     stop

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

     ! Inicio da leitura do index1====================

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
     
     ! Fim da leitura do ind2===================

     close(3)

  end if
  
2 format(a10)
  
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
        
     end if
          
  else

     bin_coarse = bin_out
     
     if (.not.bin)then
        
        n_grid = bin_out
        
     end if

  end if  
  !===================================================

  !===================================================

 ! write(*, *) ""
 ! if (.not. bin) write(*, '(a15, i3)') " DEF_BIN = ", n_grid 
 ! if (rough == 1.0) write(*, '(a15, f3.1)') " DEF_ROUGH = ", rough
 ! if (.not. skip) write(*, '(a15, i3)') " DEF_SKIP = ", n_skip
 ! if (l_coarse) write(*, '(a15, i3)') " DEF_COARSE = ", bin_coarse


  !===================================================
  
  !===================================================
  ! definição da matriz de pontos externos
  
  if (outer)then
     
     allocate (mat1(bin_out+1, bin_out+1), stat = ierr)
     if (ierr /= 0) stop 'Not enough memory to initialize mat matrix'
     
     allocate (mat2(bin_out+1, bin_out+1), stat = ierr)
     if (ierr /= 0) stop 'Not enough memory to initialize mat matrix'

     mat1(1:bin_out+1,1:bin_out+1)%z = -1000
     mat2(1:bin_out+1,1:bin_out+1)%z = 1000
          
  end if
  
  !=================================================

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

  !===============================================================

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
  
  inquire(file='area-med.xvg', exist=ex) ! avaliando existencia do arquivo area-med.xvg
  
  if (ex) then
     
     call execute_command_line("mv area-med.xvg area-med_p.xvg")
     back = .true.
     
  end if

  open(2, file='area-med.xvg', status='new', iostat=ierr)

  if(ierr /= 0) then

     write(*, *)
     write(*, *) 'Unable to open file area-med.xvg'
     write(*, *)
     stop

  endif

  write(2, '(a7, a5)') "#SuAVE ", version
  write(2, '(a14)') '#Command Line:'
  write(2, '(a9)', advance='no') '#s_area  '
  write(2, *) (trim(get(i)),"  ", i=1, 30)
  write(2, *) '@    title "Area per lipid X Frame"'
  write(2, *) '@    xaxis  label "#Frame"'
  write(2, *) '@    yaxis  label "Area per lipid [nm\S2\N]"'

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
     write(3, '(a9)', advance='no') '#s_area  '
     write(3, *) (trim(get(i)),"  ", i=1, 30)
     write(3, *) '@    title "RMSD X Frame"'
     write(3, *) '@    xaxis  label "#Frame"'
     write(3, *) '@    yaxis  label "RMSD [nm]"'

  end if

  call system_clock (start, clock_rate, clock_max)

! Leitura do arquivo .pdb =================================

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
        
     else
        
        if (n_index > 1)then

           frame  = frame + 1

           eval_skip = ((mod((frame-fr_in),(n_skip+1))==0))
           
           if ((frame>fr_in-1).and.(frame<fr_end+1).and.(eval_skip)) then

              tot_frame = tot_frame + 1 ! para programas que calculam propriedades medias
              
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
                 
              end if! if outer   
              
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
              s_area = 0

              do i=2, n_grid+1

                 do j=2, n_grid+1
                    
                    dist_x = abs(grid3(i-1,j-1)%x - grid3(i-1,j)%x)
                    dist_y = abs(grid3(i-1,j-1)%y - grid3(i-1,j)%y)
                    dist_z = abs(grid3(i-1,j-1)%z - grid3(i-1,j)%z)
                    
                    la = sqrt(dist_x**2 + dist_y**2 + dist_z**2)
                    
                    dist_x = abs(grid3(i-1,j-1)%x - grid3(i,j-1)%x)
                    dist_y = abs(grid3(i-1,j-1)%y - grid3(i,j-1)%y)
                    dist_z = abs(grid3(i-1,j-1)%z - grid3(i,j-1)%z)
                    
                    lb = sqrt(dist_x**2 + dist_y**2 + dist_z**2)
                    
                    dist_x = abs(grid3(i-1,j)%x - grid3(i,j-1)%x)
                    dist_y = abs(grid3(i-1,j)%y - grid3(i,j-1)%y)
                    dist_z = abs(grid3(i-1,j)%z - grid3(i,j-1)%z)
                    
                    lc = sqrt(dist_x**2 + dist_y**2 + dist_z**2)
                    
                    s_soma = (la + lb + lc)/2
                    
                    s_area = s_area + sqrt(s_soma*(s_soma-la)*(s_soma-lb)*(s_soma-lc))
                    
                    ! segundo triângulo
                    
                    dist_x = abs(grid3(i-1,j)%x - grid3(i,j)%x)
                    dist_y = abs(grid3(i-1,j)%y - grid3(i,j)%y)
                    dist_z = abs(grid3(i-1,j)%z - grid3(i,j)%z)
                    
                    la = sqrt(dist_x**2 + dist_y**2 + dist_z**2)
                    
                    dist_x = abs(grid3(i,j-1)%x - grid3(i,j)%x)
                    dist_y = abs(grid3(i,j-1)%y - grid3(i,j)%y)
                    dist_z = abs(grid3(i,j-1)%z - grid3(i,j)%z)
                    
                    lb = sqrt(dist_x**2 + dist_y**2 + dist_z**2)
                    
                    s_soma = (la + lb + lc)/2
                    
                    s_area = s_area + sqrt(s_soma*(s_soma-la)*(s_soma-lb)*(s_soma-lc))
                    
                 end do
                 
              end do
                            
              write(2, *) frame, s_area/100/n_lipid
              !FIM do calculo ================================================
              
              !calculo do RMSD =======================================
              if (rmsd) then
                 
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
              !======================================================================
              
           end if !======((frame<fr_in-1).and.(frame>fr_end+1))
           
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

  call system_clock (finish, clock_rate, clock_max)

  ! Fim da leitura do arquivo .pdb

  close(1)
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
     
     do i=1, n_grid+1

        do j=1, n_grid+1

           write(3, 13, iostat=ierr) 'ATOM  ', (i-1)*(n_grid+1) + j, '  DOT', &
                '  DOT',' ', 1,'    ', grid(i,j)%x, &
                grid(i,j)%y, grid(i,j)%z, r_xpm1(i,j)

        end do

     end do

     close(3)
     !first grid created

     inquire(file='grid2.pdb', exist=ex)

     if (ex) then

        call execute_command_line("mv grid2.pdb grid2_p.pdb")
        back  = .true.
        
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

     do i=1, n_grid+1

        do j=1, n_grid+1

           write(3, 13, iostat=ierr) 'ATOM  ', (i-1)*(n_grid+1) + j, '  DOT', &
                '  DOT',' ', 1,'    ', grid2(i,j)%x, &
                grid2(i,j)%y, grid2(i,j)%z, r_xpm2(i,j)

        end do

     end do

     close(3)
     !second grid created

     inquire(file='grid3.pdb', exist=ex)

     if (ex) then

        call execute_command_line("mv grid3.pdb grid3_p.pdb")
        back = .true.
        
     end if

     open(3, file='grid3.pdb', status='new', iostat=ierr)

     if(ierr /= 0) then

        write(*, *)
        write(*, *) 'Unable to open file grid3.pdb'
        write(*, *)
        stop

     end if

     write(3, *) 'TITLE     Rectangular Regular GRID'
     write(3, *) 'REMARK    THIS IS A SIMULATION BOX'
     write(3, *) 'MODEL        1'

     do i=1, n_grid+1

        do j=1, n_grid+1

           write(3, 12, iostat=ierr) 'ATOM  ', (i-1)*(n_grid+1) + j, '  DOT', &
                '  DOT',' ', 1,'    ', grid3(i,j)%x, &
                grid3(i,j)%y, grid3(i,j)%z

        end do

     end do

     close(3)
     !third grid created
     
  end if

  write(*, *)
  write(*, *) 
  write(*, *) "Finished"
  write(*, *)
  if (back) write(*, '(a31)') ' Previous files were backed up!'
  
  hour = real(finish-start)/(3600*real(clock_rate))
  minu = (hour - int(hour))*60
  sec = (minu-int(minu))*60

  if (back)  write(*, *)
  write(*, 1) " Processing time : ", int(hour), " h", int(minu), " min", int(sec), " sec"
1 format(a19, i6.3, a2, i6.3, a4, i6.3, a4)
  write(*, *)
  
  if (outer)then

     deallocate(mat1)
     deallocate(mat2)

  end if

  deallocate(r_xpm1)
  deallocate(r_xpm2)
  
end program ar_grid


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
