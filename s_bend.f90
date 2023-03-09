module bendmod

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

end module bendmod

program bend
  
  use bendmod

  logical :: ex, bin, p_grid, rmsd, l_coarse, begin, end, skip
  logical :: eval_skip, range, help, back

  integer :: i, j, ierr, n_grid, frame, aux, noi1
  integer :: n_index, num, num2, i_atom, bin_coarse, a, b
  integer :: fr_in, fr_end, n_skip, tot_frame, bini
  integer :: start, finish, clock_rate, clock_max
  integer, dimension(50000) :: in_num
  integer, dimension(200) :: hist

  real :: dist_z, hour, minu, sec, noir
  real :: dist, s_grid, r_med1, aux2, aver, aver2
  real :: la, r_fit, al
  real :: cent_x, cent_y, cent_z, dph, dth
  real :: minv, maxv, desv, graph, rough
  real, dimension(:, :), allocatable :: r_xpm

  character(len=30) :: coord,  ind, atom, ind2
  character(len=1), dimension(:,:), allocatable :: xpm
  character(len=30), dimension(20) :: get

  type(vet3), dimension(50000) :: store, coarse
  type(vet1) :: buff
  type(vet1), dimension(500000) :: spher
  type(vet3), dimension(1001,1001) :: grid
  type(vet2), dimension(1001,1001) :: grid3
  type(vet2) :: center, v1, v2, v3 
  
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
  write(*, *) "                       ** s_bend **"
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
  
  !HELP begins =================================================================
  if (help) then

     write(*, *) ""
     write(*, *) ""
     Write(*, *) "s_bend calculates the distribution of angles between the "
     Write(*, *) "vectors normal to the partition surfaces and the gradient "
     Write(*, *) "vector of the sphere circumscribing the compact surface. "
     write(*, *) "s_bend also calculates the angular deflection among these "
     write(*, *) "two vectors as an order parameter P(theta)"
     Write(*, *) "or SC."
     write(*, *) ""
     write(*, *) "Usage: s_bend -in file.pdb -ind file.ndx "
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
  write(2, '(a9)', advance='no') '#s_bend  '
  write(2, *) (trim(get(i)),"  ", i=1, 20)
  write(2, *) '@    title "Average Curvature Order Parameter X Frame"'
  write(2, *) '@    xaxis  label "#Frame"'
  write(2, *) '@    yaxis  label "Average Curvature Order Parameter"'

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
     write(3, '(a9)', advance='no') '#s_bend  '
     write(3, *) (trim(get(i)),"  ", i=1, 20)
     write(3, *) '@    title "RMSD X Frame"'
     write(3, *) '@    xaxis  label "#Frame"'
     write(3, *) '@    yaxis  label "RMSD [nm]"'

  end if

  do k=1, 200

     hist(k) = 0

  end do

  allocate (xpm((n_grid - int(n_grid/2)), n_grid), stat = ierr)
  if (ierr /= 0) stop 'Not enough memory to initialize xpm matrix'

  allocate (r_xpm((n_grid - int(n_grid/2)), n_grid), stat = ierr)
  if (ierr /= 0) stop 'Not enough memory to initialize xpm matrix'

  do i=1, n_grid - int(n_grid/2) 

     do j=1, n_grid 

        r_xpm(i,j) = 0
        
     end do

  end do
  
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
                 
              else
                 
                 coarse = store
                 
              end if
              
              ! estruturação do prmeiro grid de alta resolução========================
              
              dph = (2*pi)/(2*(n_grid - int(n_grid/2)))
              dth = (2*pi)/n_grid
              
              call param_esf(r_med1, num, r_fit, al, rough)

	      !$OMP parallel do private(s_grid, dist, k, j)
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
                            
              ! Fim da estruturação========================================

              !$OMP parallel do private(j)
              do i=1, (n_grid - int(n_grid/2))+1
                 
                 do j=1, n_grid+1
                    
                    grid3(i,j)%x = grid(i,j)%rho*sin(grid(i,j)%phi)*cos(grid(i,j)%theta)
                    grid3(i,j)%y = grid(i,j)%rho*sin(grid(i,j)%phi)*sin(grid(i,j)%theta)
                    grid3(i,j)%z = grid(i,j)%rho*cos(grid(i,j)%phi)
                    
                 end do
                 
              end do
              !$OMP end parallel do

              ! Cálculo do ângulo de inclinação =================================

              aver = 0
              aver2 = 0

              do i=2, (n_grid - int(n_grid/2))+1

                 do j=2, n_grid+1

                    v1%x = grid3(i,j)%x-grid3(i-1,j-1)%x
                    v1%y = grid3(i,j)%y-grid3(i-1,j-1)%y
                    v1%z = grid3(i,j)%z-grid3(i-1,j-1)%z

                    v2%x = grid3(i-1,j)%x-grid3(i,j-1)%x
                    v2%y = grid3(i-1,j)%y-grid3(i,j-1)%y
                    v2%z = grid3(i-1,j)%z-grid3(i,j-1)%z

                    v3%x = - v2%y*v1%z + v2%z*v1%y
                    v3%y = - v2%z*v1%x + v2%x*v1%z
                    v3%z = - v2%x*v1%y + v2%y*v1%x

                    !reutilizando a variável v1 para facilitar o cálculo ======
                    ! v1 é o gradiente definido entre os quatro pontos avaliados

                    v1%x = sin((i-1-0.5)*dph)*cos((j-1-0.5)*dth)
                    v1%y = sin((i-1-0.5)*dph)*sin((j-1-0.5)*dth)
                    v1%z = cos((i-1-0.5)*dph)

                    !==========================================================

                    la = v1%x*v3%x + v1%y*v3%y + v1%z*v3%z
                    la = la/sqrt(v3%x**2 + v3%y**2 + v3%z**2)
                    la = la/sqrt(v1%x**2 + v1%y**2 + v1%z**2)

                    !erros de arredondamento levam a valores 
                    !levemente maiores que 1.00 (1.00000012)
                    la = min(la, 1.00)
                    
                    aux2 = 0.5*(3*(la)**2 - 1)

                    r_xpm(i-1, j-1) = r_xpm(i-1, j-1) + aux2

                    aver = aver + aux2/(n_grid*(n_grid - int(n_grid/2)))
                    aver2 = aver2 + aux2*aux2/(n_grid*(n_grid - int(n_grid/2)))

                    la = acos(la)*180/pi !colocando em graus

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

	      write(2, *) frame, aver
              
              ! esses termos extras são por conta da distribuição
              ! não ser simétrica o que causava valores maiores que 1
              ! e menores que -0.5 para os limites do intervalo

              minv = max(minv, -0.5)
              maxv = min(maxv, 1.0)

              !=============================================

              ! Fim do calculo ==========================================
              
              ! Calculando o RMSD========================================
              if (rmsd)then
                 
                 noir = 0
                 
                 do i=1, noi1
                    
                    a = nint((store(i)%phi)/dph) + 1                    
                    b = nint((store(i)%theta)/dth) + 1
                    dist_z = store(i)%rho-grid(a,b)%rho
                    noir = noir + dist_z*dist_z/(noi1)
                    
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
  write(5, '(a9)', advance='no') '#s_bend  '
  write(5, *) (trim(get(i)),"  ", i=1, 20)
  write(5, *) '/* XPM */'
  write(5, *)'/* This matrix is generated by s_shell.f90 */'
  write(5, *)'/* title:   "Order Parameter" */'
  write(5, *)'/* x-label: "Theta angle " */'
  write(5, *)'/* y-label: "Phi angle " */'
  write(5, *)'/* type:    "Continuous" */'
  write(5, *)'static char * gv_xpm[] = {'
  write(5, *) '"',n_grid, (n_grid - int(n_grid/2)),' 7 1",'

  do i=1, n_grid - int(n_grid/2)
     
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

  do i=1, n_grid - int(n_grid/2) 

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

  k = n_grid - int(n_grid/2)
3 format(a10)
4 format(f6.2, a1)

  do i=1, int((n_grid - int(n_grid/2))/20)+1

     write(5, 3, advance='no') '/* y-axis:'
     j = 1

     do while ((j<=20).and.(k>=1))
        
        write(5, 4, advance='no') ((k-1/2)*dph), ' '
        k = k-1
        j = j+1

     end do
     write(5,*) '*/'

  end do

  k = 1
  do i=1, int(n_grid/20)+1

     write(5, 3, advance='no') '/* x-axis:'
     j = 1

     do while ((j<=20).and.(k<=n_grid))

        write(5, 4, advance='no') ((k-1/2)*dth), ' '
        k = k+1
        j = j+1

     end do
     write(5,*) '*/'

  end do

5 format(a1)
  do i=1, n_grid - int(n_grid/2)
    
     write(5, 5, advance='no') '"'
     
     do j=1, n_grid 
        
        write(5, 5, advance='no') xpm(i,j)

     end do

     write(5, 5, advance='no')'"'
     write(5, 5) ','
     
  end do
  
  close(5)

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
  write(4, '(a9)', advance='no') '#s_bend  '
  write(4, *) (trim(get(i)),"  ", i=1, 20)
  write(4, *) '@    title "Curvature angles"'
  write(4, *) '@    xaxis  label "Angles [\So\N]"'
  write(4, *) '@    yaxis  label "%"'

  do i=1, 100

     graph = hist(i)
     graph = graph/(tot_frame*1*n_grid*(n_grid - int(n_grid/2)))

     write(4, *) i-1, graph*100

  end do

  close(4)

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
     
     do i=1, (n_grid - int(n_grid/2))

        do j=1, n_grid

           grid3(i,j)%x = grid(i,j)%rho*sin(grid(i,j)%phi)*cos(grid(i,j)%theta) + center%x
           grid3(i,j)%y = grid(i,j)%rho*sin(grid(i,j)%phi)*sin(grid(i,j)%theta) + center%y
           grid3(i,j)%z = grid(i,j)%rho*cos(grid(i,j)%phi) + center%z

           write(3, 13, iostat=ierr) 'ATOM  ', (i-1)*(n_grid +1) + j, '  DOT', &
                '  DOT',' ', 1,'    ', grid3(i,j)%x, &
                grid3(i,j)%y, grid3(i,j)%z, r_xpm(i,j)

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

  if(back) write(*, *)
  write(*, 1) " Processing time : ", int(hour), " h", int(minu), " min", int(sec), " sec"
1 format(a19, i6.3, a2, i6.3, a4, i6.3, a4)
  write(*, *)

  deallocate(r_xpm)
  deallocate(xpm)
    
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
  real, intent(in) :: r_med, rough
  real :: r_fit, al
  integer :: num

  r_fit = 6*r_med*pi/sqrt(real(num-1))

  al = (num-1)*100/(4*pi*r_med*r_med)
  al = exp(0.4984*rough*log(al) - 1.06016110229/rough)

end subroutine param_esf
