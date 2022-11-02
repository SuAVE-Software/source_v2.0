module gauss

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


end module gauss

program curvature

  use gauss

  logical :: ex, bin, outer, p_grid, down, rmsd, l_coarse
  logical :: begin, end, skip, eval_skip, help, back
  
  integer :: i, j, ierr, n_grid, frame, a, b, noi1
  integer :: n_index, i_atom, aux, bin_out, bin_coarse
  integer :: fr_in, fr_end, n_skip, tot_frame, bini_g, bini_h
  integer :: start, finish, clock_rate, clock_max
  integer, dimension(50000) :: in_num

  real :: dist_z, dx, dy, dist, noir
  real :: x_max, x_min, y_max, y_min, hour, minu, sec, r_fit, al
  real :: del, gx, gy, peso, rough, Lx, Ly
  real :: total_lines, lines, progress
  real :: maxg, ming, maxh, minh  
  real, dimension(1000, 1000) :: r_xpmg, r_xpmh
  real :: aveg, ave2g, aveh, ave2h, desvg, desvh
  real, dimension(1000) :: hist_g, hist_h
  
  double precision :: fx, fy, fxx, fyy, fxy, kg, hm
  
  character(len=30) :: coord, index, atom
  character, dimension(1000 , 1000):: xpmg, xpmh
  character(len=30), dimension(30) :: get
  
  type(vet1) :: buff
  type(vet1), dimension(500000) :: out
  type(vet2), dimension(:,:), allocatable :: mat
  type(vet1), dimension(50000) :: store, coarse
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
  write(*, *) "                       ** s_gauss **"
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
! pegando os arquivos de entrada =======================================
!

  back = .false.
  outer = .false.
  down = .false.
  p_grid = .false.
  coord = 'miss'
  index = 'miss'
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

     if (get(i)=='-ind')then

        index = get(i+1)

     end if

     if (get(i)=='-outer')then

        outer = .true.
        read(get(i+1), *, iostat=ierr) bin_out
        call erro(ierr, "OUTER")

     end if

     if (get(i)=='-down')then

        down= .true.

     end if

     if (get(i)=='-bin')then

        bin = .true.
        read(get(i+1), *, iostat=ierr) n_grid
        call erro(ierr, " BIN ")

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
     Write(*, *) "s_gauss calculates the Gaussian and Mean Curvatures for"
     write(*, *) "the interface" 
     write(*, *) ""
     write(*, *) "Usage: s_gauss -in file.pdb -ind file.ndx"
     write(*, *) ""
     write(*, *) "file.pdb ---- atomic coordinates in PDB format"
     write(*, *) "file.ndx ---- index file containing user-selected atoms used" 
     write(*, *) "to fit the grid points to the chemical surface."
     write(*, *) ""
     write(*, *) "Options:"
     write(*, *) ""
     write(*, *) "-bin             defines the number of rectangular partition"
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
  
  !=================== zerando memoria do XPM====================
  !============== e outras variáveis ============================

  do i=1, 1000

     do j=1, 1000

        r_xpmg(i,j) = 0
        r_xpmh(i,j) = 0

     end do
     
     hist_g(i) = 0
     hist_h(i) = 0
     
  end do

  maxg = -10000000
  ming = 10000000
  maxh = -100000000
  minh = 10000000

  !====================================================================
  
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
     write(3, '(a10)', advance='no') '#s_gauss  '
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
              !$OMP end parallel do
              !=================== Fim da estruturação ==========================

              
              !==================================================================
              ! Cálculo das curvaturas gaussiana e média
              ! Pela metodologia implementada n_grid deve ser PAR!!!!
              ! Pois estamos pulando de 2 em 2

              aveg = 0
              ave2g = 0
              aveh = 0
              ave2h = 0

              do i=2, n_grid, 2
                 
                 do j=2, n_grid, 2
                    
	            fx = (grid(i+1,j)%z-grid(i-1,j)%z)/(2*dx)
                    fy = (grid(i,j+1)%z-grid(i,j-1)%z)/(2*dy)
                    
                    fxx = (grid(i+1,j)%z-2*grid(i,j)%z+grid(i-1,j)%z)/(dx*dx)
                    fyy = (grid(i,j+1)%z-2*grid(i,j)%z+grid(i,j-1)%z)/(dy*dy)
                    fxy = (grid(i+1,j+1)%z-grid(i-1,j+1)%z-grid(i+1,j-1)%z+grid(i-1,j-1)%z)/(4*dy*dx)
                    
                    kg = (fxx*fyy-fxy*fxy)/((1+fx*fx+fy*fy)**2)
                    hm = ((1+fy*fy)*fxx-fx*fy*fxy+(1+fx*fx)*fyy)/(2*(1+fx*fx+fy*fy)*sqrt(1+fx*fx+fy*fy))
                    
                    r_xpmg(i, j) = r_xpmg(i, j) + kg
                    r_xpmh(i, j) = r_xpmh(i, j) + hm

                    bini_g = nint(kg/0.0005) + 500 ! parametros ad hoc
                    hist_g(bini_g) = hist_g(bini_g) + 1
                    bini_h = nint(hm/0.0005) + 500 
                    hist_h(bini_h) = hist_h(bini_h) + 1
                    
                    aveg = aveg + kg/(n_grid*n_grid/4)
                    ave2g = ave2g + kg*kg/(n_grid*n_grid/4)
                    aveh = aveh + hm/(n_grid*n_grid/4)
                    ave2h = ave2h + hm*hm/(n_grid*n_grid/4)

                 end do
                 
              end do
             
              desvg = sqrt(ave2g - aveg*aveg)
              desvh = sqrt(ave2h - aveh*aveh)
              
              maxg = aveg + 2*desvg
              ming = aveg - 2*desvg
              maxh = aveh + 2*desvh
              minh = aveh - 2*desvh

              !=================================================================

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

           Lx = x_max - x_min
           Ly = y_max - y_min
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
  close(4)
  
  if (rmsd) then
     
     close(3)

  end if
  
  call system_clock(finish, clock_rate, clock_max)
  
  ! escrita do arquivo XPMG========================================

  inquire(file='gaussian.xpm', exist=ex)

  if (ex) then

     call execute_command_line("mv gaussian.xpm gaussian_p.xpm")
     back = .true.

  end if

  open(5, file='gaussian.xpm', status='new', iostat=ierr)

  if(ierr /= 0) then

     write(*, *)
     write(*, *) 'Unable to open file gaussian.xpm'
     write(*, *)
     stop

  endif

  write(5, '(a7, a5)') "#SuAVE ", version
  write(5, '(a14)') '#Command Line:'
  write(5, '(a9)', advance='no') '#s_gauss  '
  write(5, *) (trim(get(i)),"  ", i=1, 30)
  write(5, *) '/* XPM */'
  write(5, *)'/* This matrix is generated by s_gauss.f90 */'
  write(5, *)'/* title:   "Gaussian Curvature" */'
  write(5, *)'/* x-label: "x axis [nm]" */'
  write(5, *)'/* y-label: "y axis [nm]" */'
  write(5, *)'/* type:    "Continuous" */'
  write(5, *)'static char * gv_xpm[] = {'
  write(5, *) '"',n_grid/2,n_grid/2,' 7 1",'

  do i=2, n_grid, 2
     
     do j=2, n_grid, 2
        
        r_xpmg(i,j) = r_xpmg(i,j)/tot_frame
        
     end do
     
  end do
  
  del = (maxg - ming)/6
7 format(a20, f6.3, a5)
  write(5, 7) '"A  c #993333" /* "',ming,'" */,'
  write(5, 7) '"B  c #cc0000" /* "',(ming+del),'" */,'
  write(5, 7) '"C  c #ff5050" /* "',(ming+2*del),'" */,'
  write(5, 7) '"D  c #ff6666" /* "',(ming+3*del),'" */,'
  write(5, 7) '"E  c #ff9999" /* "',(ming+4*del),'" */,'
  write(5, 7) '"F  c #ffcccc" /* "',(ming+5*del),'" */,'
  write(5, 7) '"G  c #ffffff" /* "',(ming+6*del),'" */,'

  do i=2, n_grid, 2
     
     do j=2, n_grid, 2
        
        if (r_xpmg(i,j)<ming+del)then
           xpmg(i,j)= 'A'
        else
           if (r_xpmg(i,j)<ming+2*del)then
              xpmg(i,j) = 'B'
           else
              if (r_xpmg(i,j)<ming+3*del)then
                 xpmg(i,j) = 'C'
              else
                 if (r_xpmg(i,j)<ming+4*del) then
                    xpmg(i,j) = 'D'
                 else
                    if (r_xpmg(i,j)<ming+5*del) then
                       xpmg(i,j) = 'E'
                    else
                       if (r_xpmg(i,j)<ming+6*del) then
                          xpmg(i,j) = 'F'
                       else
                          xpmg(i,j) = 'G'

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

  do i=1, int(n_grid/40)+1

     write(5, 3, advance='no') '/* y-axis:'
     j = 1

     do while ((j<=20).and.(k<=n_grid/2))
        
        write(5, 4, advance='no') ((k*2-1)*dy + gy)/10, ' '
        k = k+1
        j = j+1

     end do
     write(5,*) '*/'

  end do
  write(5,*) '*/'

  k = 1
  do i=1, int(n_grid/40)+1

     write(5, 3, advance='no') '/* x-axis:'
     j = 1

     do while ((j<=20).and.(k<=n_grid/2))

        write(5, 4, advance='no') ((k*2-1)*dx + gx)/10, ' '
        k = k+1
        j = j+1

     end do
     write(5,*) '*/'

  end do

5 format(a1)
  do j=n_grid, 2, -2
    
     write(5, 5, advance='no') '"'
     
     do i=2, n_grid, 2
        
        write(5, 5, advance='no') xpmg(i,j)

     end do

     write(5, 5, advance='no')'"'
     write(5, 5) ','
     
  end do

  close(5)
  
  ! Final da escrita do XPMG


  ! escrita do arquivo XPMH========================================

  inquire(file='mean.xpm', exist=ex)

  if (ex) then

     call execute_command_line("mv mean.xpm mean_p.xpm")
     back = .true.

  end if

  open(5, file='mean.xpm', status='new', iostat=ierr)

  if(ierr /= 0) then

     write(*, *)
     write(*, *) 'Unable to open file mean.xpm'
     write(*, *)
     stop

  endif

  write(5, '(a7, a5)') "#SuAVE ", version
  write(5, '(a14)') '#Command Line:'
  write(5, '(a9)', advance='no') '#s_gauss  '
  write(5, *) (trim(get(i)),"  ", i=1, 30)
  write(5, *) '/* XPM */'
  write(5, *)'/* This matrix is generated by s_gauss.f90 */'
  write(5, *)'/* title:   "Mean Curvature" */'
  write(5, *)'/* x-label: "x axis [nm]" */'
  write(5, *)'/* y-label: "y axis [nm]" */'
  write(5, *)'/* type:    "Continuous" */'
  write(5, *)'static char * gv_xpm[] = {'
  write(5, *) '"',n_grid/2,n_grid/2,' 7 1",'

  do i=2, n_grid, 2
     
     do j=2, n_grid, 2
        
        r_xpmh(i,j) = r_xpmh(i,j)/tot_frame
                
     end do
     
  end do
  
  del = (maxh - minh)/6

  write(5, 7) '"A  c #993333" /* "',minh,'" */,'
  write(5, 7) '"B  c #cc0000" /* "',(minh+del),'" */,'
  write(5, 7) '"C  c #ff5050" /* "',(minh+2*del),'" */,'
  write(5, 7) '"D  c #ff6666" /* "',(minh+3*del),'" */,'
  write(5, 7) '"E  c #ff9999" /* "',(minh+4*del),'" */,'
  write(5, 7) '"F  c #ffcccc" /* "',(minh+5*del),'" */,'
  write(5, 7) '"G  c #ffffff" /* "',(minh+6*del),'" */,'

  do i=2, n_grid, 2
     
     do j=2, n_grid, 2
        
        if (r_xpmh(i,j)<minh+del)then
           xpmh(i,j)= 'A'
        else
           if (r_xpmh(i,j)<minh+2*del)then
              xpmh(i,j) = 'B'
           else
              if (r_xpmh(i,j)<minh+3*del)then
                 xpmh(i,j) = 'C'
              else
                 if (r_xpmh(i,j)<minh+4*del) then
                    xpmh(i,j) = 'D'
                 else
                    if (r_xpmh(i,j)<minh+5*del) then
                       xpmh(i,j) = 'E'
                    else
                       if (r_xpmh(i,j)<minh+6*del) then
                          xpmh(i,j) = 'F'
                       else
                          xpmh(i,j) = 'G'

                       end if
                       
                    end if
                    
                 end if
                 
              end if
              
           end if
           
        end if

     end do

  end do

  k = 1

  do i=1, int(n_grid/40)+1

     write(5, 3, advance='no') '/* y-axis:'
     j = 1

     do while ((j<=20).and.(k<=n_grid/2))
        
        write(5, 4, advance='no') ((k*2-1)*dy + gy)/10, ' '
        k = k+1
        j = j+1

     end do
     write(5,*) '*/'

  end do

  k = 1
  do i=1, int(n_grid/40)+1

     write(5, 3, advance='no') '/* x-axis:'
     j = 1

     do while ((j<=20).and.(k<=n_grid/2))

        write(5, 4, advance='no') ((k*2-1)*dx + gx)/10, ' '
        k = k+1
        j = j+1

     end do
     write(5,*) '*/'

  end do

  do j=n_grid, 2, -2
    
     write(5, 5, advance='no') '"'
     
     do i=2, n_grid, 2
        
        write(5, 5, advance='no') xpmh(i,j)

     end do

     write(5, 5, advance='no')'"'
     write(5, 5) ','
     
  end do

  close(5)
  
  ! Final da escrita do XPMH ================================
  
  ! Escrita das distribuições de probabilidade ==============
  
  inquire(file='hist_gaussian.xvg', exist=ex)

  if (ex) then

     call execute_command_line("mv hist_gaussian.xvg hist_gaussian_p.xvg")
     back = .true.

  end if

  open(5, file='hist_gaussian.xvg', status='new', iostat=ierr)

  if(ierr /=0) then

     write(*, *)
     write(*, *) 'Unable to open file hist_gaussian.xvg'
     write(*, *)
     stop

  endif

  write(5, '(a7, a5)') "#SuAVE ", version
  write(5, '(a14)') '#Command Line:'
  write(5, '(a9)', advance='no') '#s_gauss  '
  write(5, *) (trim(get(i)),"  ", i=1, 30)
  write(5, *) '@    title "Gaussian Curvature Distribution"'
  write(5, *) '@    xaxis  label "Gaussian Curvature"'
  write(5, *) '@    yaxis  label "PDF"'

  do i=200, 800

     write(5, *) ((i-500)*0.0005), hist_g(i)/(tot_frame*0.0005*(n_grid*n_grid/4))

  end do

  close(5)
  
  ! Mean Curvature ==========================================
  inquire(file='hist_mean.xvg', exist=ex)

  if (ex) then

     call execute_command_line("mv hist_mean.xvg hist_mean_p.xvg")
     back = .true.

  end if

  open(5, file='hist_mean.xvg', status='new', iostat=ierr)

  if(ierr /=0) then

     write(*, *)
     write(*, *) 'Unable to open file hist_mean.xvg'
     write(*, *)
     stop

  endif

  write(5, '(a7, a5)') "#SuAVE ", version
  write(5, '(a14)') '#Command Line:'
  write(5, '(a9)', advance='no') '#s_gauss  '
  write(5, *) (trim(get(i)),"  ", i=1, 30)
  write(5, *) '@    title "Mean Curvature Distribution"'
  write(5, *) '@    xaxis  label "Mean Curvature"'
  write(5, *) '@    yaxis  label "PDF"'

  do i=1, 1000

     write(5, *) ((i-500)*0.0005), hist_h(i)/(tot_frame*0.0005*(n_grid*n_grid/4))

  end do
  
  ! Fim da escrita ==========================================
  
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

13   format(a6, i5.1, a5, a5, a1, i4.1, a4, 3f8.3, 2f5.1)
     
     do i=2, n_grid, 2

        do j=2, n_grid, 2

           write(3, 13, iostat=ierr) 'ATOM  ', (i-1)*(n_grid+1) + j, '  DOT', &
                '  DOT',' ', 1,'    ', grid(i,j)%x, &
                grid(i,j)%y, grid(i,j)%z, (r_xpmg(i,j)-ming)/(maxg-ming), (r_xpmh(i,j)-minh)/(maxh-minh)

        end do

     end do

     close(3)
     ! grid created
     
  end if

  call calcula_bin(aux, Lx, Ly)
  
  write(*, *)
  write(*, *)
  write(*, *) 'Best results with bin < ', aux
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

! Fim da leitura do arquivo .pdb


  if (outer)then

     deallocate(mat)

  end if
  
end program curvature
  

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

subroutine calcula_bin(n_grid, X, Y)

  integer :: n_grid
  real :: X, Y

  !o valor de 0.08 veio de um estudo com diferentes funções
  !onde verificou-se que valores de dx ou dy menores que 0.08
  !poderiam levar a erros nos calculos das curvaturas

  if (X>Y) then

     n_grid = int(Y/0.08)

  else

     n_grid = int(X/0.08)

  end if

  if (mod(n_grid,2) /= 0.0)  n_grid = n_grid + 1

  if (n_grid>500) n_grid = 500
  
end subroutine calcula_bin
 
