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
  logical :: eval_skip, help, consta, molec, back

  integer :: i, j, ierr, n_grid, frame, aux, noi1
  integer :: n_index, num, num2, i_atom, bin_coarse, a, b
  integer :: fr_in, fr_end, n_skip, inner, a_dens
  integer :: start, finish, clock_rate, clock_max
  integer, dimension(500000) :: in_num, in_dens

  real :: hour, minu, sec
  real :: s_grid, r_max, dist
  real :: r_fit, al, rough
  real :: total_lines, lines, progress  
  real :: cent_x, cent_y, cent_z, dph, dth

  character(len=30) :: coord,  ind, atom, ind3
  character(len=30), dimension(20) :: get

  type(vet3), dimension(500000) :: store, coarse
  type(vet1) :: buff
  type(vet1), dimension(500000) :: spher
  type(vet3), dimension(1001,1001) :: grid
  type(vet2), dimension(1001,1001) :: grid3
  type(vet2) :: center
  type(vet3) :: dens(500000)
  type(vet1), dimension(500000) :: dens_spher  

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
  write(*, *) "                       ** s_count **"
  write(*, *) ""
  write(*, *) "                    ** VERSION ", version, " **"
  write(*, *) ""
  write(*, *) "               Santos, D. E. S.; Soares, T. A."
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
  ind3 = 'miss'
  rmsd = .false.
  l_coarse = .false.
  begin = .false.
  end = .false.
  skip = .false.
  help = .false.
  consta = .false.
  molec = .false.
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

     if (get(i)=='-dens')then

        ind3 = get(i+1)

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

     if (get(i)=='-molec')then

        molec = .true.

     end if

     if (get(i)=='-const')then

        consta = .true.
        read(get(i+1), *, iostat=ierr) const
        call erro(ierr, "CONST")

     end if

  end do

  !HELP begins =================================================================
  if (help) then

     write(*, *) ""
     write(*, *) ""
     write(*, *) "s_count calculates the number of molecules inside a specified"
     write(*, *) "structure surrounded by SuAVE grid."
     write(*, *) ""
     write(*, *) "Usage: s_count -in file.pdb -ind file1.ndx  -dens dens.ndx "
     write(*, *) ""
     write(*, *) "file.pdb ---- atomic coordinates in PDB format"
     write(*, *) "file.ndx ---- index file containing the atoms composing the"
     write(*, *) "the structure "
     Write(*, *) "dens.ndx ---- index file indicating the molecules that have"
     Write(*, *) "to be counted"
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
     write(*, *) "-const           defines the constant multiplying the radius"
     write(*, *) ""
     write(*, *) "-molec           generates a PDB file containing the molecules"
     write(*, *) "observed inside the structure in each frame of the trajectory"
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

  ! inicio da leitura do index das moleculas===============
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

  do i=n_index, 500000

     in_dens(i) = -1

  end do

  write(*, *)
  write(*, *) ind3, 'input file is OK'
    
  ! Fim da leitura do index3
  
  close(4)
  
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
  
  inquire(file='molecules.xvg', exist=ex)
  
  if (ex) then
     
     call execute_command_line("mv molecules.xvg molecules_p.xvg")
     back = .true.

  end if

  open(2, file='molecules.xvg', status='new', iostat=ierr)

  if(ierr /= 0) then

     write(*, *)
     write(*, *) 'Unable to open file molecules.xvg'
     write(*, *)
     stop

  endif

  write(2, '(a7, a5)') "#SuAVE ", version
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
        
        if (i_atom == in_dens(a_dens))then
           
           dens_spher(a_dens) = buff
           a_dens = a_dens + 1
           
        end if

     else

        if (n_index > 1)then
           
           frame  = frame + 1

           eval_skip = ((mod((frame-fr_in),(n_skip+1))==0))
           
           if ((frame>fr_in-1).and.(frame<fr_end+1).and.(eval_skip)) then

              num = 1
              r_max = 0
              
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
                 
                 r_max = max(store(num)%rho, r_max)
                 
                 num = num + 1                    
                 
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
              
              
           end if !======((frame<fr_in-1).and.(frame>fr_end+1)) 

           n_index = 1
           aux = a_dens
           a_dens = 1
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
  close(5)

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
    
  end if

  !escreve os átomos dentro da estrutura para o último frame
  if (molec) then

     inquire(file='molecules.pdb', exist=ex)
     
     if (ex) then
        
        call execute_command_line("mv molecules.pdb molecules_p.pdb")
        back = .true.
        
     end if
     
     open(3, file='molecules.pdb', status='new', iostat=ierr)
     
     if(ierr /= 0) then
        
        write(*, *)
        write(*, *) 'Unable to open file molecules.pdb'
        write(*, *)
        stop
        
     endif
    
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


subroutine param_esf(r_med, num, r_fit, al, rough)

  real, parameter :: pi = 3.141592654
  real, intent(in) :: r_med, rough
  real :: r_fit, al
  integer :: num

  r_fit = 6*r_med*pi/sqrt(real(num-1))

  al = (num-1)*100/(4*pi*r_med*r_med)
  al = exp(0.4984*rough*log(al) - 1.06016110229/rough)

end subroutine param_esf
