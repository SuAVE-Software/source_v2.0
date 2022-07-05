module shellmod

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

  type vet3

     real ::  rho, phi, theta

  end type vet3

end module shellmod

program shell
  
  use shellmod

  logical :: ex, bin, p_grid, rmsd, l_coarse, begin, end, skip
  logical :: eval_skip, range, help, back

  integer :: i, j, ierr, n_grid, frame, aux, noi1, noi2
  integer :: n_index, num, num2, i_atom, bin_coarse, a, b
  integer :: fr_in, fr_end, n_skip, tot_frame
  integer :: start, finish, clock_rate, clock_max
  integer, dimension(50000) :: in_num, in_num2

  real :: dist_x, dist_y, dist_z, s_soma, hour, minu, sec, noir
  real :: dist, s_grid, r_med1, r_med2, aux2, aver, aver2
  real :: s_area, s_area2, la, lb, lc, r_fit, al
  real :: cent_x, cent_y, cent_z, dph, dth, s_vol, s_vol2
  real :: minv, maxv, desv, c_angle, rough
  real :: total_lines, lines, progress
  real, dimension(:, :), allocatable :: r_xpm

  character(len=30) :: coord,  ind, atom, ind2
  character(len=1), dimension(:,:), allocatable :: xpm
  character(len=30), dimension(20) :: get

  type(vet3), dimension(50000) :: store, store2, coarse, coarse2
  type(vet1) :: buff
  type(vet1), dimension(500000) :: spher
  type(vet3), dimension(1001,1001) :: grid, grid2
  type(vet2), dimension(1001,1001) :: grid3
  type(vet2) :: center 
  
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
  write(*, *) "                       ** s_shell **"
  write(*, *) ""
  write(*, *) "                   ** VERSION ", version, " **"
  write(*, *) ""
  write(*, *) "               Santos, D. E. S.; Soares, T. A."
  write(*, *) ""
  write(*, *) "Please cite SuAVE: A Tool for Analyzing Curvature-Dependent"
  write(*, *) "Properties in Chemical Interfaces (2020) Denys E. S. Santos,"
  write(*, *) "Frederico J. S. Pontes, Roberto D. Lins, Kaline Coutinho,"
  write(*, *) "Thereza A. Soares. J. Chem. Inf. Model., v. 60(2), p. 473-484."

!
! pegando os arquivos de entrada
!

  back = .false.
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
  range=.false.
  help = .false.
  rough = 1.0
  
  do i=1, 20

     call getarg(i, get(i))

  end do

  do i=1, 20

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

        bin= .true.
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
        read(get(i+1), *, iostat=ierr) n_skip
        call erro(ierr, " SKIP")

     end if
     
     if (get(i)=='-range')then

        range = .true.

     end if

     if (get(i)=='-rough')then

        read(get(i+1), *, iostat=ierr) rough
        call erro(ierr, "REDUC")

     end if
     
     if (get(i)=='-help')then

        help = .true.

     end if

  end do

  !HELP begins ============================================================
  if (help) then

       write(*, *) ""
       write(*, *) ""
       write(*, *) "s_shell calculates the average thickness and volume "
       write(*, *) "for spherical shells along the simulation time."
       write(*, *) ""
       write(*, *) "Usage: s_shell -in file.pdb -ind1 file1.ndx -ind2 file2.ndx"
       write(*, *) ""
       write(*, *) "file.pdb ---- atomic coordinates in PDB format"
       write(*, *) "file.ndx ---- index file containing user-selected atoms used"
       write(*, *) "to fit the grid points to the chemical surface."
       write(*, *) ""
       write(*, *) "Options:"
       write(*, *) ""
       write(*, *) "-bin             defines the number of rectangular partition "
       write(*, *) "bins along the phi- and psi-angles"
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
       write(*, *)"-help             prints HELP information and quits"
       stop
       
  end if
  !HELP ends ========================================================================
  
  if (coord=='miss')then

     write(*, *)
     write(*, *)'PDB file is missing'
     write(*, *)
     stop

  end if

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

 !==Lendo os arquivos de index=========================

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
  
  ! Fim da leitura do index=====================

  close(2)
  
  ! Inicio da leitura do ind2===================

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
  
  ! Fim da leitura do ind2=======================

  close(3)

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
  
  bin_coarse = nint(sqrt(real(2*(n_index-1))))
  
  if (.not.bin)then
     
     n_grid = nint(sqrt(real(2*(n_index-1))))
     write(*, *)
     write(*, *) 'STD_BIN = ', n_grid
     
  end if
     
  !============================================================
     
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
  
  inquire(file='thick.xvg', exist=ex)

  if (ex) then

     call execute_command_line("mv thick.xvg thick_p.xvg")
     back = .true.
     
  end if
  
  open(2, file='thick.xvg', status='new', iostat=ierr)

  if(ierr /= 0) then

     write(*, *)
     write(*, *) 'Unable to open file thick.xvg'
     write(*, *)
     stop

  endif

  write(2, '(a7, a5)') "#SuAVE ", version
  write(2, '(a14)') '#Command Line:'
  write(2, '(a10)', advance='no') '#s_shell  '
  write(2, *) (trim(get(i)),"  ", i=1, 20)
  write(2, *) '@    title "Thickness X Frame"'
  write(2, *) '@    xaxis  label "#Frame"'
  write(2, *) '@    yaxis  label "Average Thickness [nm]"'

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
     write(3, '(a10)', advance='no') '#s_shell  '
     write(3, *) (trim(get(i)),"  ", i=1, 20)
     write(3, *) '@    title "RMSD X Frame"'
     write(3, *) '@    xaxis  label "#Frame"'
     write(3, *) '@    yaxis  label "RMSD [nm]"'

  end if

  allocate (xpm((n_grid - int(n_grid/2)+1), n_grid+1), stat = ierr)
  if (ierr /= 0) stop 'Not enough memory to initialize xpm matrix'

  allocate (r_xpm((n_grid - int(n_grid/2)+1), n_grid+1), stat = ierr)
  if (ierr /= 0) stop 'Not enough memory to initialize xpm matrix'

  do i=1, n_grid - int(n_grid/2) + 1

     do j=1, n_grid + 1

        r_xpm(i,j) = 0
        
     end do

  end do

  inquire(file='volume.xvg', exist=ex)

  if (ex) then

     call execute_command_line("mv volume.xvg volume_p.xvg")
     back = .true.

  end if

  open(4, file='volume.xvg', status='new', iostat=ierr)

  if(ierr /= 0) then

     write(*, *)
     write(*, *) 'Unable to open file volume.xvg'
     write(*, *)
     stop

  endif

  write(4, '(a7, a5)') "#SuAVE ", version
  write(4, '(a14)') '#Command Line:'
  write(4, '(a10)', advance='no') '#s_shell  '
  write(4, *) (trim(get(i)),"  ", i=1, 20)
  write(4, *) '@    title "Volume X Frame"'
  write(4, *) '@    xaxis  label "#Frame"'
  write(4, *) '@    yaxis  label "Total Spherical Shell Volume[nm\S3\N]"'
  
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
     
     if (atom.eq.'ATOM  ') then
        
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
        
     else

        if (n_index > 1)then
           
           frame  = frame + 1

           eval_skip = ((mod((frame-fr_in),(n_skip+1))==0))
           
           if ((frame>fr_in-1).and.(frame<fr_end+1).and.(eval_skip)) then

              tot_frame = tot_frame + 1
              
              num = 1
              num2 = 1
              r_med1 = 0
              r_med2 = 0
              
              ! transformando em Coordenadas Esfericas ======================================================
              
              do i=1, n_index - 1
                 
                 if ((spher(i)%n_atom == in_num(num))) then
                    
                    store(num)%rho = sqrt((spher(i)%x-cent_x)**2 + (spher(i)%y-cent_y)**2 + (spher(i)%z-cent_z)**2)
                    store(num)%phi = acos((spher(i)%z-cent_z)/store(num)%rho)
                    
                    if ((spher(i)%y-cent_y)==0)then
                       
                       store(num)%theta = 0.000

                    else
                       
                       store(num)%theta = acos((spher(i)%x-cent_x)/sqrt((spher(i)%x-cent_x)**2 + (spher(i)%y-cent_y)**2))
                       
                       if ((spher(i)%y-cent_y)<=0)then
                          
                          store(num)%theta = 2*pi - store(num)%theta
                          
                       end if

                    end if

                    r_med1 = (r_med1*(num-1) + store(num)%rho)/(num)

                    num = num + 1
                    
                 else
                    
                    store2(num2)%rho = sqrt((spher(i)%x-cent_x)**2 + (spher(i)%y-cent_y)**2 + (spher(i)%z-cent_z)**2)
                    store2(num2)%phi = acos((spher(i)%z-cent_z)/store2(num2)%rho)

                    if ((spher(i)%y-cent_y)==0) then

                       store2(num2)%theta = 0.000

                    else

                       store2(num2)%theta = acos((spher(i)%x-cent_x)/sqrt((spher(i)%x-cent_x)**2 + (spher(i)%y-cent_y)**2))

                       if ((spher(i)%y-cent_y)<0)then
                          
                          store2(num2)%theta = 2*pi - store2(num2)%theta
                          
                       end if

                    end if

                    r_med2 = (r_med2*(num2-1) + store2(num2)%rho)/(num2)
                    
                    num2 = num2 + 1
                    
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

              do i=1, (n_grid - int(n_grid/2))+1
                 
                 do j=1, n_grid+1
                    
                    grid3(i,j)%x = grid(i,j)%rho*sin(grid(i,j)%phi)*cos(grid(i,j)%theta) + cent_x
                    grid3(i,j)%y = grid(i,j)%rho*sin(grid(i,j)%phi)*sin(grid(i,j)%theta) + cent_y
                    grid3(i,j)%z = grid(i,j)%rho*cos(grid(i,j)%phi) + cent_z
                    
                 end do
                 
              end do
              
              ! Calculando a primeira area do grid==============================
              
              s_vol  = 0
              
              do i=2, (n_grid - int(n_grid/2))+1
                 
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
                    
                    aux2 = sqrt(s_soma*abs(s_soma-la)*abs(s_soma-lb)*abs(s_soma-lc))
                    c_angle = ang(grid3(i-1,j-1), grid3(i-1,j), grid3(i,j-1), (i-1-0.5)*dph, (j-1-0.5)*dth)
                    
                    s_vol = s_vol + c_angle*aux2*grid(i-1,j-1)%rho/3

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
                    
                    aux2 = sqrt(s_soma*abs(s_soma-la)*abs(s_soma-lb)*abs(s_soma-lc))
                    c_angle = ang(grid3(i-1,j), grid3(i,j-1), grid3(i,j), (i-1-0.5)*dph, (j-1-0.5)*dth)

                    s_vol = s_vol + c_angle*aux2*grid(i,j)%rho/3

                 end do
                 
              end do
              
              ! Calculando para a segunda area, do grid2=========================
              do i=1, (n_grid - int(n_grid/2))+1
                 
                 do j=1, n_grid+1
                    
                    grid3(i,j)%x = grid2(i,j)%rho*sin(grid2(i,j)%phi)*cos(grid2(i,j)%theta) + cent_x
                    grid3(i,j)%y = grid2(i,j)%rho*sin(grid2(i,j)%phi)*sin(grid2(i,j)%theta) + cent_y
                    grid3(i,j)%z = grid2(i,j)%rho*cos(grid2(i,j)%phi) + cent_z
                    
                 end do
                 
              end do

              s_vol2  = 0

              do i=2, (n_grid - int(n_grid/2))+1
                 
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
                    
                    aux2 = sqrt(s_soma*abs(s_soma-la)*abs(s_soma-lb)*abs(s_soma-lc))
                    c_angle = ang(grid3(i-1,j-1), grid3(i-1,j), grid3(i,j-1), (i-1-0.5)*dph, (j-1-0.5)*dth)

                    s_vol2 = s_vol2 + c_angle*aux2*grid2(i-1,j-1)%rho/3
                    
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
                    
                    aux2 = sqrt(s_soma*abs(s_soma-la)*abs(s_soma-lb)*abs(s_soma-lc))
                    c_angle = ang(grid3(i-1,j), grid3(i,j-1), grid3(i,j), (i-1-0.5)*dph, (j-1-0.5)*dth)

                    s_vol2 = s_vol2 + c_angle*aux2*grid2(i,j)%rho/3
                    
                 end do
                 
              end do

              write(4, *) frame, abs(s_vol - s_vol2)/(1000)
              
              ! Fim do calculo ==========================================

              ! Calculando o thickness médio=============================
              
              aver = 0
              aver2 = 0
              
              do i=1, (n_grid - int(n_grid/2))+1
                 
                 do j=1, n_grid+1
                    
                    aux2 = abs((grid(i,j)%rho-grid2(i,j)%rho))/10
                    aver = aver + aux2/((n_grid+1)*(n_grid - int(n_grid/2)+1))
                    aver2 = aver2 + aux2*aux2/((n_grid+1)*(n_grid - int(n_grid/2)+1))

                    r_xpm(i,j) = r_xpm(i,j) + abs(grid(i,j)%rho-grid2(i,j)%rho)
                    
                 end do
                 
              end do
              
              desv = sqrt(aver2 - aver*aver)

              minv = min(aver-2*desv, minv)
              maxv = max(aver+2*desv, maxv)

              write(2, *) frame, aver
              
              ! Calculando o RMSD========================================
              if (rmsd)then
                 
                 noir = 0
                 
                 do i=1, noi1
                    
                    a = nint((store(i)%phi)/dph) + 1
                    b = nint((store(i)%theta)/dth) + 1
                    dist_z = store(i)%rho-grid(a,b)%rho
                    noir = noir + dist_z*dist_z/(noi1+noi2)
                    
                 end do
                 
                 do i=1, noi2
                    
                    a = nint((store2(i)%phi)/dph) + 1
                    b = nint((store2(i)%theta)/dth) + 1
                    dist_z = store2(i)%rho-grid2(a,b)%rho
                    noir = noir + dist_z*dist_z/(noi1+noi2)
                    
                 end do
                 
                 noir = sqrt(noir)
                 
                 write(3, *) noir/10
                 
              end if
              ! Fim do calculo do RMSD===================================
              
           end if !======((frame<fr_in-1).and.(frame>fr_end+1)) 
              
           n_index = 1
           i_atom = 0
           num = 1
           num2 = 1
           center%x = cent_x
           center%y = cent_y
           center%z = cent_z
           cent_x = 0
           cent_y = 0
           cent_z = 0

           !====garante que fr_end sempre seja maior que frame ===
           !====caso essa variável não tenha sido fixada==========
	   if (.not.end)then

              fr_end = frame + 1

           end if
           !======================================================
           
        end if
        
     end if

     if ((lines/total_lines) > progress) then ! imprimindo barra de avanco
        
        write(*,'(1a1, a12, f5.1, a1, $)') char(13), ' Progress : ', progress*100, '%' 
        progress = progress + 0.01
        
     end if


  end do
  
  call system_clock(finish, clock_rate, clock_max)
  
  ! Fim da leitura do arquivo .pdb

  close(1)
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

  inquire(file='thick.xpm', exist=ex)
  
  if (ex) then
     
     call execute_command_line("mv thick.xpm thick_p.xpm")
     back = .true.

  end if
  
  open(4, file='thick.xpm', status='new', iostat=ierr)

  if(ierr /= 0) then

     write(*, *)
     write(*, *) 'Unable to open file thick.xpm'
     write(*, *)
     stop

  endif

  write(4, '(a7, a5)') "#SuAVE ", version
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
  write(4, *) '"',n_grid+1, (n_grid - int(n_grid/2)+1),' 7 1",'

  do i=1, n_grid - int(n_grid/2)+1
     
     do j=1, n_grid + 1
        
        r_xpm(i,j) = r_xpm(i,j)/tot_frame/10

     end do
     
  end do

  del = (maxv - minv)/6
7 format(a20, f4.1, a5)
  write(4, 7) '"A  c #ffffff " /* "',minv,'" */,'
  write(4, 7) '"B  c #87cefa " /* "',(minv+del),'" */,'
  write(4, 7) '"C  c #00bfff " /* "',(minv+2*del),'" */,'
  write(4, 7) '"D  c #1e90ff " /* "',(minv+3*del),'" */,'
  write(4, 7) '"E  c #4169e1 " /* "',(minv+4*del),'" */,'
  write(4, 7) '"F  c #0000ff " /* "',(minv+5*del),'" */,'
  write(4, 7) '"G  c #00008b " /* "',(minv+6*del),'" */,'

  do i=1, n_grid - int(n_grid/2) + 1

     do j=1, n_grid + 1

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

  k = n_grid - int(n_grid/2)
3 format(a10)
4 format(f6.2, a1)

  do i=1, int((n_grid - int(n_grid/2))/20)+1

     write(4, 3, advance='no') '/* y-axis:'
     j = 1

     do while ((j<=20).and.(k>=1))
        
        write(4, 4, advance='no') ((k-1/2)*dph), ' '
        k = k-1
        j = j+1

     end do
     write(4,*) '*/'

  end do

  k = 1
  do i=1, int(n_grid/20)+1

     write(4, 3, advance='no') '/* x-axis:'
     j = 1

     do while ((j<=20).and.(k<=n_grid))

        write(4, 4, advance='no') ((k-1/2)*dth), ' '
        k = k+1
        j = j+1

     end do
     write(4,*) '*/'

  end do

5 format(a1)
  do i=1, n_grid - int(n_grid/2)+1
    
     write(4, 5, advance='no') '"'
     
     do j=1, n_grid +1
        
        write(4, 5, advance='no') xpm(i,j)

     end do

     write(4, 5, advance='no')'"'
     write(4, 5) ','
     
  end do

  close(4)

  !creating grid files ==============================================================
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
     
     do i=1, (n_grid - int(n_grid/2))+1

        do j=1, n_grid+1

           grid3(i,j)%x = grid(i,j)%rho*sin(grid(i,j)%phi)*cos(grid(i,j)%theta) + center%x
           grid3(i,j)%y = grid(i,j)%rho*sin(grid(i,j)%phi)*sin(grid(i,j)%theta) + center%y
           grid3(i,j)%z = grid(i,j)%rho*cos(grid(i,j)%phi) + center%z

           write(3, 13, iostat=ierr) 'ATOM  ', (i-1)*(n_grid +1) + j, '  DOT', &
                '  DOT',' ', 1,'    ', grid3(i,j)%x, &
                grid3(i,j)%y, grid3(i,j)%z, r_xpm(i,j)

        end do

     end do

     close(3)
     !first grid created========================================================

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

     do i=1, (n_grid - int(n_grid/2))+1

        do j=1, n_grid+1

           grid3(i,j)%x = grid2(i,j)%rho*sin(grid2(i,j)%phi)*cos(grid2(i,j)%theta) + center%x
           grid3(i,j)%y = grid2(i,j)%rho*sin(grid2(i,j)%phi)*sin(grid2(i,j)%theta) + center%y
           grid3(i,j)%z = grid2(i,j)%rho*cos(grid2(i,j)%phi) + center%z

           write(3, 13, iostat=ierr) 'ATOM  ', (i-1)*(n_grid +1) + j, '  DOT', &
                '  DOT',' ', 1,'    ', grid3(i,j)%x, &
                grid3(i,j)%y, grid3(i,j)%z, r_xpm(i,j)

        end do

     end do

     close(3)
     !second grid created =======================================================
     
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

  deallocate(r_xpm)
  deallocate(xpm)
    
end program shell


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

real function ang(p1, p2, p3, phi, theta)

  use shellmod

  type(vet2), intent(in) :: p1, p2, p3
  type(vet2) :: v1, v2, v3

  real :: la, phi, theta

  v1%x = p2%x-p1%x
  v1%y = p2%y-p1%y
  v1%z = p2%z-p1%z
  
  v2%x = p3%x-p2%x
  v2%y = p3%y-p2%y
  v2%z = p3%z-p2%z
  
  v3%x = - v2%y*v1%z + v2%z*v1%y
  v3%y = - v2%z*v1%x + v2%x*v1%z
  v3%z = - v2%x*v1%y + v2%y*v1%x
  
  !reutilizando a variável v1 para facilitar o cálculo ====== 
  
  v1%x = sin(phi)*cos(theta)
  v1%y = sin(phi)*sin(theta)
  v1%z = cos(phi)
  
  !========================================================== 

  if (v3%x**2 + v3%y**2 + v3%z**2<0.000001) then

     !serve para os casos das extremidades de phi, onde a área
     !do triangulo é nula. Neste caso não importa o ângulo.
     !Por isso utilizamos o valor 1, pois a área já é nula

     ang = 1

  else

     la = v1%x*v3%x + v1%y*v3%y + v1%z*v3%z
     la = la/sqrt(v3%x**2 + v3%y**2 + v3%z**2)
     la = la/sqrt(v1%x**2 + v1%y**2 + v1%z**2)

     ang = abs(la)

  end if

end function ang

subroutine param_esf(r_med, num, r_fit, al, rough)

  real, parameter :: pi = 3.141592654
  real, intent(in) :: r_med, rough
  real :: r_fit, al
  integer :: num

  r_fit = 6*r_med*pi/sqrt(real(num-1))

  al = (num-1)*100/(4*pi*r_med*r_med)
  al = exp(0.4984*rough*log(al) - 1.06016110229/rough)

end subroutine param_esf
