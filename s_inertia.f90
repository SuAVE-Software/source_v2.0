module inert

  implicit none

  real, parameter :: pi = 3.14159265358979323846
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

end module inert

program bend
  
  use inert

  logical :: ex, bin, p_grid, rmsd, l_coarse, begin, end, skip
  logical :: eval_skip, help, back

  integer :: i, j, ierr, n_grid, frame, aux, noi1
  integer :: n_index, num, num2, i_atom, bin_coarse, a, b
  integer :: fr_in, fr_end, n_skip, tot_frame
  integer :: start, finish, clock_rate, clock_max
  integer, dimension(50000) :: in_num

  double precision :: MI(3,3)
  double precision :: Q(3,3)
  double precision :: W(3), PCA(3), ORD(3) 

  real :: dist_z, hour, minu, sec, noir
  real :: dist, s_grid, r_med1, averx, avery, averz
  real :: r_fit, al, rough
  real :: cent_x, cent_y, cent_z, dth, dz
  real :: total_lines, lines, progress
  real :: minv, maxv

  character(len=30) :: coord,  ind, atom, ind2
  character(len=30), dimension(20) :: get

  type(vet3), dimension(50000) :: store, coarse
  type(vet1) :: buff
  type(vet1), dimension(500000) :: spher
  type(vet3), dimension(1001,1001) :: grid
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
  write(*, *) "                       ** s_inertia **"
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
  
  do i=1, 20

     call getarg(i, get(i))

  end do

  do i=1, 20

     if (get(i)=='-in')then

        coord = get(i+1)

     end if

     if (get(i)=='-ind')then

        ind = get(i+1)

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

     if (get(i)=='-rough')then

        read(get(i+1), *, iostat=ierr) rough
        call erro(ierr, "REDUC")

     end if
     
     if (get(i)=='-help')then

        help = .true.

     end if

  end do
  
  !HELP begins =================================================================
  if (help) then

     write(*, *) ""
     write(*, *) ""
     Write(*, *) "s_inertia calculates the Principal Moments of Inertia of the"
     write(*, *) "closed surface"
     write(*, *) ""
     write(*, *) "Usage: s_inertia -in file.pdb -ind file.ndx "
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

  if (ind=='miss')then
     
     write(*, *)
     write(*, *)'Index file is missing'
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
  
  !=========abrindo arquivos de saída =====================================
  inquire(file='inertia.xvg', exist=ex)

  if (ex) then

     call execute_command_line("mv inertia.xvg inertia_p.xvg")
     back = .true.

  end if

  open(2, file='inertia.xvg', status='new', iostat=ierr)

  if(ierr /= 0) then

     write(*, *)
     write(*, *) 'Unable to open file inertia.xvg'
     write(*, *)
     stop

  endif

  write(2, '(a7, a5)') "#SuAVE ", version
  write(2, '(a14)') '#Command Line:'
  write(2, '(a12)', advance='no') '#s_inertia  '
  write(2, *) (trim(get(i)),"  ", i=1, 20)
  write(2, *) '@    title "Principal Moments of Inertia X Frame"'
  write(2, *) '@    xaxis  label "#Frame"'
  write(2, *) '@    yaxis  label "Principal Moments of Inertia [nm\S2\N]"'

        
  inquire(file='roundness.xvg', exist=ex)

  if (ex) then

     call execute_command_line("mv roundness.xvg roundness_p.xvg")
     back = .true.

  end if

  open(4, file='roundness.xvg', status='new', iostat=ierr)

  if(ierr /= 0) then

     write(*, *)
     write(*, *) 'Unable to open file roundness.xvg'
     write(*, *)
     stop

  endif

  write(4, '(a7, a5)') "#SuAVE ", version
  write(4, '(a14)') '#Command Line:'
  write(4, '(a12)', advance='no') '#s_inertia  '
  write(4, *) (trim(get(i)),"  ", i=1, 20)
  write(4, *) '@    title "Inertia Roundness X Frame"'
  write(4, *) '@    xaxis  label "#Frame"'
  write(4, *) '@    yaxis  label "Roundness"'
  write(4, *) '# Frame    Inert. Round.     Circle Round.'

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
        
     else

        if (n_index > 1)then
           
           frame  = frame + 1

           eval_skip = ((mod((frame-fr_in),(n_skip+1))==0))
           
           if ((frame>fr_in-1).and.(frame<fr_end+1).and.(eval_skip)) then

              tot_frame = tot_frame + 1
              
              num = 1
              r_med1 = 0

              ! transformando em Coordenadas Esfericas ======================================================
              
              do i=1, n_index - 1
                 
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

               data MI/0, 0, 0, 0, 0, 0, 0, 0, 0/

               averx = 0
               avery = 0
               averz = 0

               do i=1, (n_grid) 
                  
                  do j=1, n_grid 
                     
                    grid3(i,j)%x = grid(i,j)%rho*sin(grid(i,j)%phi)*cos(grid(i,j)%theta)/10
                    grid3(i,j)%y = grid(i,j)%rho*sin(grid(i,j)%phi)*sin(grid(i,j)%theta)/10
                    grid3(i,j)%z = grid(i,j)%rho*cos(grid(i,j)%phi)/10

                    averx = averx + grid3(i,j)%x
                    avery = avery + grid3(i,j)%y
                    averz = averz + grid3(i,j)%z
                    
                 end do
                 
              end do
              
              averx = averx/(n_grid*n_grid)
              avery = avery/(n_grid*n_grid)
              averz = averz/(n_grid*n_grid)
              
              do i=1, (n_grid) 
                 
                 do j=1, n_grid 
                                        
                    MI(1,1) = MI(1,1) + (grid3(i,j)%y-avery)**2 + (grid3(i,j)%z-averz)**2
                    MI(2,2) = MI(2,2) + (grid3(i,j)%x-averx)**2 + (grid3(i,j)%z-averz)**2
                    MI(3,3) = MI(3,3) + (grid3(i,j)%x-averx)**2 + (grid3(i,j)%y-avery)**2
                    MI(1,2) = MI(1,2) - (grid3(i,j)%x-averx)*(grid3(i,j)%y-avery)
                    MI(1,3) = MI(1,3) - (grid3(i,j)%x-averx)*(grid3(i,j)%z-averz)
                    MI(2,3) = MI(2,3) - (grid3(i,j)%y-avery)*(grid3(i,j)%z-averz)
                    
                 end do
                 
              end do
              
              MI(1,1) = MI(1,1)/(n_grid*n_grid)
              MI(2,2) = MI(2,2)/(n_grid*n_grid)
              MI(3,3) = MI(3,3)/(n_grid*n_grid)
              MI(1,2) = MI(1,2)/(n_grid*n_grid)
              MI(1,3) = MI(1,3)/(n_grid*n_grid)
              MI(2,3) = MI(2,3)/(n_grid*n_grid)

              MI(2,1) = MI(1,2)
              MI(3,1) = MI(1,3)
              MI(3,2) = MI(2,3)

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
                 
                 noir = 0
                 
                 do i=1, noi1

                    !esse calculo é único para este programa
                    a = nint((cos(store(i)%phi)+1)/dz + 1/2)
                    b = nint((store(i)%theta)/dth) + 1

                    dist_z = store(i)%rho-grid(a,b)%rho
                    
                    noir = noir + dist_z*dist_z/(noi1)
                    
                 end do
                                  
                 noir = sqrt(noir)
                 
                 write(3, *) frame, noir/10
                 
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
  close(4)

  if (rmsd)then
     
     close(3)

  end if
  
  !creating grid files ==============================================================
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
     
     do i=1, (n_grid )

        do j=1, n_grid

           grid3(i,j)%x = grid(i,j)%rho*sin(grid(i,j)%phi)*cos(grid(i,j)%theta) + center%x
           grid3(i,j)%y = grid(i,j)%rho*sin(grid(i,j)%phi)*sin(grid(i,j)%theta) + center%y
           grid3(i,j)%z = grid(i,j)%rho*cos(grid(i,j)%phi) + center%z

           write(3, 13, iostat=ierr) 'ATOM  ', (i-1)*(n_grid +1) + j, '  DOT', &
                '  DOT',' ', 1,'    ', grid3(i,j)%x, &
                grid3(i,j)%y, grid3(i,j)%z

        end do

     end do

     close(3)
     
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


end program bend


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

subroutine param_esf(r_med, num, r_fit, al, rough)

  real, parameter :: pi = 3.141592654
  real, intent(in) :: r_med
  real :: r_fit, al
  integer :: num

  r_fit = 6*r_med*pi/sqrt(real(num-1))

  al = (num-1)*100/(4*pi*r_med*r_med)
  al = exp(0.4984*rough*log(al) - 1.06016110229/rough)

end subroutine param_esf


subroutine ordem(W, ORD)

  integer :: i
  double precision :: W(3), ORD(3)

  ORD(1) = maxval(W)
  ORD(3) = minval(W)

  do i=1, 3

     if (W(i) < ORD(1)) then

        if (W(i) > ORD(3)) then
           
           ORD(2) = W(i)

        end if
     
     end if

  end do

end subroutine ordem
