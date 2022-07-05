module histo

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


end module histo

program histograma

  use histo

  logical :: ex, bin, outer, p_grid, down, rmsd, l_coarse
  logical :: begin, end, skip, eval_skip, range, help, back
  
  integer :: i, j, ierr, n_grid, frame, bini, a, b, noi1, noi2
  integer :: n_index, i_atom, aux, bin_out, bin_coarse
  integer :: fr_in, fr_end, n_skip, tot_frame
  integer :: start, finish, clock_rate, clock_max
  integer, dimension(50000) :: in_num
  integer, dimension(100) :: hist

  real :: dist_x, dist_y, dist_z, dx, dy, dist, noir
  real :: x_max, x_min, y_max, y_min, hour, minu, sec, r_fit, al
  real :: la, graph, minv, maxv, del, gx, gy, aux2, rough
  real :: total_lines, lines, progress
  real :: aver, aver2, desv, peso
  real, dimension(: , :), allocatable :: r_xpm

  character(len=30) :: coord, index, atom
  character, dimension(: , :), allocatable :: xpm
  character(len=30), dimension(30) :: get
  
  type(vet1) :: buff
  type(vet1), dimension(500000) :: out
  type(vet2), dimension(:,:), allocatable :: mat
  type(vet1), dimension(50000) :: store, coarse
  type(vet2), dimension(1001,1001) :: grid
  type(vet2) :: v1, v2, v3

  
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
  write(*, *) "                       ** s_order **"
  write(*, *) ""
  write(*, *) "                   ** VERSION ", version, " **"
  write(*, *) ""
  write(*, *) "              Santos, D. E. S.; Soares, T. A."
  write(*, *) ""
  write(*, *) "Please cite SuAVE: A Tool for Analyzing Curvature-Dependent"
  write(*, *) "Properties in Chemical Interfaces (2020) Denys E. S. Santos,"
  write(*, *) "Frederico J. S. Pontes, Roberto D. Lins, Kaline Coutinho,"
  write(*, *) "Thereza A. Soares. J. Chem. Inf. Model., v. 60(2), p. 473-484."

!
! pegando os arquivos de entrada
!

  back = .true.
  outer = .false.
  bin = .false.
  down = .false.
  p_grid = .false.
  coord = 'miss'
  index = 'miss'
  rmsd = .false.
  l_coarse = .false.
  begin = .false.
  end = .false.
  skip = .false.
  range = .false.
  help = .false.
  rough = 1.0
  
  do i=1, 30

     call getarg(i, get(i))

  end do

  do i=1, 30

     if (get(i)=='-in')then

        coord = get(i+1)

     end if

     if (get(i)=='-ind')then

        index = get(i+1)

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

     if (get(i)=='-down')then

        down= .true.

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
     
     if (get(i)=='-help') then

        help = .true.

     end if
     
  end do

  !HELP begins =========================================================
  if (help) then

     write(*, *) ""
     write(*, *) ""
     Write(*, *) "s_order calculates the distribution of angles between the "
     Write(*, *) "vectors normal to the rectangular partition surfaces and the z" 
     Write(*, *) "axis of the system. s_order also calculates the curvature" 
     Write(*, *) "order parameter P(theta) or SC."
     write(*, *) ""
     write(*, *) "Usage: s_order -in file.pdb -ind file.ndx"
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
     write(*, *) "-down            used to generate a PDB file with the last "
     write(*, *) "external lowest fitting grid. This option is only used when" 
     write(*, *) "followed by option –outer."
     write(*, *) ""
     write(*, *) "-range           defines the range specified in the XPM file"
     write(*, *) ""
     write(*, *) "-rough           percentage of the original surface roughness"
     write(*, *) "the user wants to keep (1.0 by default, meaning 100%)"
     write(*, *) ""
     write(*, *) "-help            prints HELP information and quits"
     stop

  end if
  !HELP ends ===========================================================

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
     
     if(ierr/=0) then

        write(*, *)
        write(*, *) 'Unable to open file ', index
        write(*, *)
        stop
        
     end if

     ! Inicio da leitura do index

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
        write(*, *) 'Problem by reading ', index
        write(*, *)
        stop
        
     end if
     
     do i=n_index, 50000
        
        in_num(i) = -1
        
     end do
     
     close(2)

     write(*, *)
     write(*, *) index, " input file is OK"
     
     ! Fim da leitura do index   

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
  
  !============================================================

  
  !============================================================
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

  !=========abrindo arquivo do parâmetro de ordem médio=========================
  inquire(file='order_aver.xvg', exist=ex)

  if (ex) then

     call execute_command_line("mv order_aver.xvg order_aver_p.xvg")
     back = .true.

  end if

  open(2, file='order_aver.xvg', status='new', iostat=ierr)

  if(ierr /= 0) then

     write(*, *)
     write(*, *) 'Unable to open file order_aver.xvg'
     write(*, *)
     stop

  endif

  write(2, '(a7, a5)') "#SuAVE ", version
  write(2, '(a14)') '#Command Line:'
  write(2, '(a10)', advance='no') '#s_order  '
  write(2, *) (trim(get(i)),"  ", i=1, 30)
  write(2, *) '@    title "Average Curvature Order Parameter X Frame"'
  write(2, *) '@    xaxis  label "#Frame"'
  write(2, *) '@    yaxis  label "Average Curvature Order Parameter"'

  !=================== alocando memoria do XPM====================
  allocate (xpm(n_grid, n_grid), stat = ierr)
  if (ierr /= 0) stop 'Not enough memory to initialize xpm matrix'

  allocate (r_xpm(n_grid, n_grid), stat = ierr)
  if (ierr /= 0) stop 'Not enough memory to initialize xpm matrix'

  do i=1, n_grid

     do j=1, n_grid

        r_xpm(i,j) = 0

     end do

  end do

  do k=1, 100

     hist(k) = 0

  end do

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
     write(3, '(a10)', advance='no') '#s_order  '
     write(3, *) (trim(get(i)),"  ", i=1, 30)
     write(3, *) '@    title "RMSD X Frame"'
     write(3, *) '@    xaxis  label "#Frame"'
     write(3, *) '@    yaxis  label "RMSD [nm]"'

  end if

  call system_clock(start, clock_rate, clock_max)

! Leitura do arquivo .gro
  
  i_atom = 0
  n_index = 1
  frame = 0
  ierr = 0
  x_min = 1000
  y_min = 1000
  x_max = 0
  y_max = 0
  tot_frame = 0
  minv = 100000
  maxv = -100000
  lines = 0
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
           
        else

           if (i_atom == in_num(n_index)) then

              store(n_index) = buff
              n_index = n_index + 1

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
                 
                 ! Estruturação do grid
                 
                 dx = (x_max - x_min)/bin_coarse
                 dy = (y_max - y_min)/bin_coarse
                 
                 ! adaptação dos valores de r_fit e alfa

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
              !OMP end parallel do
              
              ! Fim da estruturação
              
              ! Cálculo do ângulo de inclinação

              aver = 0
              aver2 = 0

              do i=2, n_grid+1
                 
                 do j=2, n_grid+1
                    
                    v1%x = grid(i,j)%x-grid(i-1,j-1)%x
                    v1%y = grid(i,j)%y-grid(i-1,j-1)%y
                    v1%z = grid(i,j)%z-grid(i-1,j-1)%z
                    
                    v2%x = grid(i-1,j)%x-grid(i,j-1)%x
                    v2%y = grid(i-1,j)%y-grid(i,j-1)%y
                    v2%z = grid(i-1,j)%z-grid(i,j-1)%z
                    
                    v3%x = - v2%y*v1%z + v2%z*v1%y
                    v3%y = - v2%z*v1%x + v2%x*v1%z
                    v3%z = - v2%x*v1%y + v2%y*v1%x
                    
                    la = acos(v3%z/sqrt(v3%x**2 + v3%y**2 + v3%z**2))*180/pi
 
                    aux2 = 0.5*(3*cos(la*pi/180)**2 - 1)
                   
                    r_xpm(i-1, j-1) = r_xpm(i-1, j-1) + aux2
                    
                    aver = aver + aux2/(n_grid*n_grid)
                    aver2 = aver2 + aux2*aux2/(n_grid*n_grid)

                    if ((la<0).or.(la>90)) then
                       
                       write(*, *) "problem ... la = ", la
                       stop
                       
                    end if

                    bini = nint(la) + 1
                                        
                    hist(bini) = hist(bini) + 1
                    
                 end do
                 
              end do

              desv = sqrt(aver2 - aver*aver)

              minv = min(aver-2*desv, minv)
              maxv = max(aver+2*desv, maxv)
              
              ! esses termos extras são por conta da distribuição
              ! não ser simétrica o que causava valores maiores que 1
              ! e menores que -0.5 para os limites do intervalo
              
              minv = max(minv, -0.5)
              maxv = min(maxv, 1.0)

              !=============================================

              ! n_grid*n_grid é o número de células do grid
              write(2, *) frame, aver
              
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

           !! ******* Fim do cálculo

           end if !======((frame<fr_in-1).and.(frame>fr_end+1))
           
           gx = x_min
           gy = y_min
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

  close(1)
  close(2)
  
  if (rmsd) then
     
     close(3)

  end if
  
  call system_clock(finish, clock_rate, clock_max)
  
  ! escrita do arquivo XPM========================================

  
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

  inquire(file='order.xpm', exist=ex)

  if (ex) then

     call execute_command_line("mv order.xpm order_p.xpm")
     back = .true.

  end if

  open(5, file='order.xpm', status='new', iostat=ierr)

  if(ierr /= 0) then

     write(*, *)
     write(*, *) 'Unable to open file order.xpm'
     write(*, *)
     stop

  endif

  write(5, '(a7, a5)') "#SuAVE ", version
  write(5, '(a14)') '#Command Line:'
  write(5, '(a10)', advance='no') '#s_order  '
  write(5, *) (trim(get(i)),"  ", i=1, 30)
  write(5, *) '/* XPM */'
  write(5, *)'/* This matrix is generated by s_order.f90 */'
  write(5, *)'/* title:   "Order Parameter" */'
  write(5, *)'/* x-label: "x axis [nm]" */'
  write(5, *)'/* y-label: "y axis [nm]" */'
  write(5, *)'/* type:    "Continuous" */'
  write(5, *)'static char * gv_xpm[] = {'
  write(5, *) '"',n_grid,n_grid,' 7 1",'

  do i=1, n_grid
     
     do j=1, n_grid
        
        r_xpm(i,j) = r_xpm(i,j)/tot_frame
                
     end do
     
  end do
  
  del = (maxv - minv)/6
7 format(a20, f4.1, a5)
  write(5, 7) '"A  c #993333" /* "',minv,'" */,'
  write(5, 7) '"B  c #cc0000" /* "',(minv+del),'" */,'
  write(5, 7) '"C  c #ff5050" /* "',(minv+2*del),'" */,'
  write(5, 7) '"D  c #ff6666" /* "',(minv+3*del),'" */,'
  write(5, 7) '"E  c #ff9999" /* "',(minv+4*del),'" */,'
  write(5, 7) '"F  c #ffcccc" /* "',(minv+5*del),'" */,'
  write(5, 7) '"G  c #ffffff" /* "',(minv+6*del),'" */,'

  do i=1, n_grid
     
     do j=1, n_grid
        
        if (r_xpm(i,j)<minv+del)then
           xpm(i,j)= 'A'
        else
           if (r_xpm(i,j)<minv+2*del)then
              xpm(i,j) = 'B'
           else
              if (r_xpm(i,j)<minv+3*del)then
                 xpm(i,j) = 'C'
              else
                 if (r_xpm(i,j)<minv+4*del) then
                    xpm(i,j) = 'D'
                 else
                    if (r_xpm(i,j)<minv+5*del) then
                       xpm(i,j) = 'E'
                    else
                       if (r_xpm(i,j)<minv+6*del) then
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

     write(5, 3, advance='no') '/* y-axis:'
     j = 1

     do while ((j<=20).and.(k<=n_grid))
        
        write(5, 4, advance='no') ((k-1/2)*dy + gy)/10, ' '
        k = k+1
        j = j+1

     end do
     write(5,*) '*/'

  end do

  k = 1
  do i=1, int(n_grid/20)+1

     write(5, 3, advance='no') '/* x-axis:'
     j = 1

     do while ((j<=20).and.(k<=n_grid))

        write(5, 4, advance='no') ((k-1/2)*dx + gx)/10, ' '
        k = k+1
        j = j+1

     end do
     write(5,*) '*/'

  end do

5 format(a1)
  do j=n_grid, 1, -1
    
     write(5, 5, advance='no') '"'
     
     do i=1, n_grid
        
        write(5, 5, advance='no') xpm(i,j)

     end do

     write(5, 5, advance='no')'"'
     write(5, 5) ','
     
  end do

  close(5)
  
  ! Final da escrita do XPM

  !=========escrita do histograma===========================
  
  inquire(file='hist.xvg', exist=ex)

  if (ex) then

     call execute_command_line("mv hist.xvg hist_p.xvg")
     back = .true.

  end if

  open(4, file='hist.xvg', status='new', iostat=ierr)

  if(ierr/=0) then

     write(*, *)
     write(*, *) 'Unable to open file hist.xvg'
     write(*, *)
     stop

  end if

  write(4, '(a7, a5)') "#SuAVE ", version
  write(4, '(a14)') '#Command Line:'
  write(4, '(a10)', advance='no') '#s_order  '
  write(4, *) (trim(get(i)),"  ", i=1, 30)
  write(4, *) '@    title "Curvature angles"'
  write(4, *) '@    xaxis  label "Angles [\So\N]"'
  write(4, *) '@    yaxis  label "%"'

  do i=1, 100
     
     graph = hist(i)
     graph = graph/(tot_frame*1*n_grid*n_grid)
     
     write(4, *) i-1, graph*100
     
  end do

  close(4)
  
  !== escrita do GRID========================================
  if (p_grid)then

     inquire(file='grid.pdb', exist=ex)

     if (ex) then

        call execute_command_line("mv grid.pdb grid_p.pdb")
        back = .true.

     end if

     open(3, file='grid.pdb', status='new', iostat=ierr)

     if(ierr /= 0) then

        write(*, *)
        write(*, *) 'Unable to open file grid.pdb'
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
     ! grid created
     
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

! Fim da leitura do arquivo .pdb

  deallocate(r_xpm)
  deallocate(xpm)

  if (outer)then

     deallocate(mat)

  end if
  
end program histograma
  

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
