module spherical

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

end module spherical

program esferico
  
  use spherical

  logical :: ex, bin, p_grid, rmsd, l_coarse, begin, end, skip
  logical :: eval_skip, slices, help, back

  integer :: i, j, ierr, n_grid, frame, aux, noi1, noi2
  integer :: n_index, num, num2, i_atom, bin_coarse, a, b
  integer :: fr_in, fr_end, n_skip, a_dens, div, bini, tot_frame
  integer :: start, finish, clock_rate, clock_max
  integer, dimension(1000000) :: in_num, in_num2, in_dens

  real :: dist_x, dist_y, dist_z, s_soma, hour, minu, sec, noir
  real :: dist, s_grid, r_med1, r_med2, aux2
  real :: la, lb, lc, r_fit, al, del, rad_aver, rough, c_angle
  real :: cent_x, cent_y, cent_z, dph, dth, s_vol, k_inf, k_sup
  real :: total_lines, lines, progress  
  real, dimension(1000):: hist

  character(len=30) :: coord,  ind, atom, ind2, ind3
  character(len=30), dimension(30) :: get

  type(vet3), dimension(1000000) :: store, store2, coarse, coarse2
  type(vet1) :: buff
  type(vet3), dimension(1000000) :: dens
  type(vet1), dimension(1000000) :: spher, dens_spher
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
  write(*, *) "                       ** s_densph **"
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
  slices = .false.
  bin = .false.
  p_grid = .false.
  coord = 'miss'
  ind = 'miss'
  ind2 = 'miss'
  ind3 = 'miss'
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

        bin= .true.
        read(get(i+1), *, iostat=ierr) n_grid
        call erro(ierr, " BIN ")
        
     end if
     
     if (get(i)=='-grid')then

        p_grid = .true.

     end if

     if (get(i)=='-slices')then

	slices = .true.
	read(get(i+1),*, iostat=ierr) div
        call erro(ierr, "SLICE")
 
     end if
     
     if (get(i)=='-dens')then
        
	ind3 = get(i+1)
 
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
     write(*, *) "s_densph calculates the density profile for compact surfaces"
     write(*, *) "taking into account their curvature"
     write(*, *) ""
     write(*, *) "Usage: s_densph -in file.pdb -ind1 file1.ndx -ind2 file2.ndx"
     write(*, *) "-dens dens.ndx"
     write(*, *) ""
     write(*, *) "file.pdb ---- atomic coordinates in PDB format"
     write(*, *) "file.ndx ---- index file containing user-selected atoms used" 
     write(*, *) "to fit the grid points to the chemical surface."
     Write(*, *) "dens.ndx ---- index file containing user-selected atoms used "
     Write(*, *) "to calculate density profile."
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
     write(*, *) "-slices          defines the number of slices along the axis"
     write(*, *) "normal to the system used to calculate the density profile"
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

  if (ind3=='miss')then

     write(*, *)
     write(*, *) 'Density index file is missing'
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
  
  do i=n_index, 1000000
     
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
  
  do i=n_index, 1000000
     
     in_num2(i) = -1
     
  end do
  
  write(*, *)
  write(*, *) ind2, " input file is OK"
  
  ! Fim da leitura do ind2=======================

  close(3)

  ! inicio da leitura do index da densidade===============
  open(4, file=ind3, status='old', iostat=ierr)

  if(ierr /= 0) then

     write(*, *)
     write(*, *) 'Unable to open file ', ind3
     write(*, *)
     stop

  end if

  ierr = 0
  n_index = 1

  do while (ierr==0)

     read(4, *, iostat=ierr) aux

     if (ierr==0) then

        in_dens(n_index) = aux
        n_index = n_index + 1

     end if

  end do

  if(ierr>0) then

     write(*, *)
     write(*, *) 'Problem by reading ', ind3
     write(*, *)
     stop

  end if

  do i=n_index, 1000000

     in_dens(i) = -1

  end do

  write(*, *)
  write(*, *) ind3, 'input file is OK'
    
  ! Fim da leitura do index3
  
  close(4)

  if (.not.slices) then

     div = 10

     write(*, *) ""
     write(*, *) "STD_SLICES = ", div
     
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
     write(3, '(a11)', advance='no') '#s_densph  '
     write(3, *) (trim(get(i)),"  ", i=1, 30)
     write(3, *) '@    title "RMSD X Frame"'
     write(3, *) '@    xaxis  label "#Frame"'
     write(3, *) '@    yaxis  label "RMSD [nm]"'

  end if
  
  ! Limpando o histograma==============================

  do i=1, 1000

     hist(i) = 0

  end do

  ! histograma limpo===================================

  call system_clock(start, clock_rate, clock_max)

! Leitura do arquivo .pdb

  i_atom = 0
  a_dens = 1
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

        if ((i_atom == in_num2(num2)))then

           spher(n_index) = buff
           cent_x = (cent_x*(n_index-1) + buff%x)/(n_index)
           cent_y = (cent_y*(n_index-1) + buff%y)/(n_index)
           cent_z = (cent_z*(n_index-1) + buff%z)/(n_index)
           num2 = num2  + 1
           n_index = n_index + 1

        end if

        if (i_atom == in_dens(a_dens))then

           dens_spher(a_dens) = buff
           a_dens = a_dens + 1

        end if
        
     else

        if ((n_index > 1).or.(a_dens > 1))then
           
           frame  = frame + 1

           eval_skip = ((mod((frame-fr_in),(n_skip+1))==0))
           
           if ((frame>fr_in-1).and.(frame<fr_end+1).and.(eval_skip)) then

              tot_frame = tot_frame + 1 ! para programas que calculam propriedades medias

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

              do i=1, a_dens - 1

                 dens(i)%rho = sqrt((dens_spher(i)%x-cent_x)**2 + (dens_spher(i)%y-cent_y)**2 + (dens_spher(i)%z-cent_z)**2)
                 dens(i)%phi = acos((dens_spher(i)%z-cent_z)/dens(i)%rho)

                 if ((dens_spher(i)%y-cent_y)==0)then

                    dens(i)%theta = 0.000

                 else

                    dens(i)%theta = acos((dens_spher(i)%x-cent_x)/sqrt((dens_spher(i)%x-cent_x)**2 + (dens_spher(i)%y-cent_y)**2))
                    
                    if((dens_spher(i)%y-cent_y)<0)then
                       
                       dens(i)%theta = 2*pi - dens(i)%theta
                       
                    end if

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
              ! Contrução do grid médio para cálculo de densidade =========

              !$OMP parallel do
              do i=1, (n_grid - int(n_grid/2))+1
                 
                 do j=1, n_grid+1
                    
                    grid3(i,j)%x = (grid(i,j)%rho + grid2(i,j)%rho)*sin(grid(i,j)%phi)*cos(grid(i,j)%theta)/2 + cent_x
                    grid3(i,j)%y = (grid(i,j)%rho + grid2(i,j)%rho)*sin(grid(i,j)%phi)*sin(grid(i,j)%theta)/2 + cent_y
                    grid3(i,j)%z = (grid(i,j)%rho + grid2(i,j)%rho)*cos(grid(i,j)%phi)/2 + cent_z

                 end do
                 
              end do
              !$OMP end parallel do

              ! Calculando volume da superficie media ====================
              
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
              
              ! Fim do calculo ==========================================

              ! cálculo da densidade esférica ===========================
              ! o calculo eh em funcao de K =============================
              
              do i=1, a_dens - 1

                 a = nint((dens(i)%phi)/dph) + 1
                 b = nint((dens(i)%theta)/dth) + 1

                 rad_aver = (grid(a,b)%rho + grid2(a,b)%rho)/2
                 del = 2.0/div

                 bini = int((dens(i)%rho/rad_aver)/del)
                 k_inf = bini*del
                 k_sup = bini*del + del 
                 
                 hist(bini+500) = hist(bini+500) + 1*1000/(s_vol*(k_sup**3 - k_inf**3))
                 
              end do

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
              
           a_dens = 1
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

  if (rmsd)then
     
     close(3)

  end if
  
  ! Imprimindo arquivo de densidade ===========================================
  inquire(file='density.xvg', exist=ex)

  if (ex) then

     call execute_command_line("mv density.xvg density_p.xvg")
     back = .true.

  end if

  open(5, file='density.xvg', status='new', iostat=ierr)

  if(ierr /=0) then

     write(*, *)
     write(*, *) 'Unable to open file density.xvg'
     write(*, *)
     stop

  endif

  write(5, '(a7, a5)') "#SuAVE ", version
  write(5, '(a14)') '#Command Line:'
  write(5, '(a11)', advance='no') '#s_densph  '
  write(5, *) (trim(get(i)),"  ", i=1, 30)
  write(5, *) '@    title "System Density"'
  write(5, *) '@    xaxis  label "R/R\smed\N"'
  write(5, *) '@    yaxis  label "Density [nm\S-3\N]"'

  do i=500, 500+2*div

     write(5, *) (i-500)*del, hist(i)/tot_frame

  end do

  close(5)

  !Arquivo escrito ===========================================================

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

     do i=1, (n_grid - int(n_grid/2))+1

        do j=1, n_grid+1

           grid3(i,j)%x = grid(i,j)%rho*sin(grid(i,j)%phi)*cos(grid(i,j)%theta) + center%x
           grid3(i,j)%y = grid(i,j)%rho*sin(grid(i,j)%phi)*sin(grid(i,j)%theta) + center%y
           grid3(i,j)%z = grid(i,j)%rho*cos(grid(i,j)%phi) + center%z

           write(3, 12, iostat=ierr) 'ATOM  ', (i-1)*(n_grid +1) + j, '  DOT', &
                '  DOT',' ', 1,'    ', grid3(i,j)%x, &
                grid3(i,j)%y, grid3(i,j)%z

        end do

     end do

     close(3)
     !first grid created

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

           write(3, 12, iostat=ierr) 'ATOM  ', (i-1)*(n_grid +1) + j, '  DOT', &
                '  DOT',' ', 1,'    ', grid3(i,j)%x, &
                grid3(i,j)%y, grid3(i,j)%z

        end do

     end do

     close(3)
     !second grid created
     
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
  
  
end program esferico


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

  use spherical

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
