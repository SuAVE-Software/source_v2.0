module separa

  implicit none

  character(len=5), parameter :: version = "2.0.0"
  
  type vet1

     character(len=30) :: atom, resid, ident, code
     real ::  x, y, z
     integer :: n_atom, n_resid
     
  end type vet1

end module separa


program separa_indice

  use separa

  logical :: ex, res, sph, igual, join, gro, bound, help, back
  
  integer :: i, j, ierr, res_i, res_j, id(500)
  integer :: n_index, i_atom, n_sph, aux2, n_index2
  integer, dimension(15)  :: ax
  
  real :: cent_x, cent_y, cent_z, rho, dist_z, aux3
  
  character(len=30) :: coord, cond, index, sai, atom, resid, gro_file
  character(len=30), dimension(30) :: get
  character(len=30) :: atomid(500,500), residue(500), aux
  
  type(vet1) :: buff 
  type(vet1), dimension(1000000) :: store, store2, sphstore

                                                              
                                                              
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
  write(*, *) "                       ** s_index ** " 
  write(*, *) "" 
  write(*, *) "                   ** VERSION ", version, " **"
  write(*, *) ""
  write(*, *) "               Santos, D. E. S; Soares, T. A."
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
  coord = 'miss'
  cond = 'miss'
  res = .false.
  sph = .false.
  join = .false.
  gro = .false.
  bound = .false.
  help = .false.
  
  do i=1, 30

     call getarg(i, get(i))

  end do

  do i=1, 30

     if (get(i)=='-in')then

        coord = get(i+1)

     end if

     if (get(i)=='-residue')then

        res = .true.

     end if

     if (get(i)=='-sphere') then

        sph = .true.

     end if

     if (get(i)=='-join') then

        join = .true.

     end if
     
     if (get(i)=='-gromacs')then

        gro = .true.
        gro_file = get(i+1)
        
     end if
     
     if (get(i)=='-bound') then

        bound = .true.

     end if
     
     if (get(i)=='-help') then
        
        help = .true.
        
     end if

  end do

  !HELP begins ===============================================================
  if (help) then

     write(*, *) ""
     write(*, *) "" 
     write(*, *) "s_index creates the input file containing indexes assigned to" 
     write(*, *) "the atoms, which will be used to define the surface/interface." 
     write(*, *) "The user should select atoms distributed along the "
     write(*, *) "full surface, and preferentially, with low atomic "
     write(*, *) "fluctuation along time."
     write(*, *) " "
     write(*, *) "Usage: s_index -in file.pdb"
     write(*, *) ""
     write(*, *) "file.pdb ---- atomic coordinates in PDB format"
     write(*, *) ""
     write(*, *) "Options: "
     write(*, *) ""
     write(*, *) "-residue         selects a full residue to fit the grid or" 
     write(*, *) "to be used for the calculation of the density profile."  
     write(*, *) "Otherwise the s_index will ask for specific atoms."
     write(*, *) ""
     write(*, *) "-sphere          defines compact systems (micelles," 
     write(*, *) "vesicles, etc.)"
     write(*, *) ""
     write(*, *) "-join            joins all indexes into a single file."
     write(*, *) "Creates the index files for the calculation of the" 
     write(*, *) "density profile. This option must be used to calculate" 
     write(*, *) "monolayers as opposed to bilayers."
     write(*, *) ""
     write(*, *) "-gromacs         converts density index files from GROMACS" 
     write(*, *) "format to SuAVE format"
     write(*, *) ""
     write(*, *) "-bound           defines the position of placement of a" 
     write(*, *) "divisor plane between the two leaflets of a bilayer. If the" 
     write(*, *) "flag is not used, s_index will assume the boundary at the" 
     write(*, *) "average distance between the two surfaces defined by the" 
     write(*, *) "atoms listed in the index file."
     write(*, *) ""
     write(*, *) "-help            prints HELP information and quits."
     stop

  end if
  !HELP end =================================================================

  if ((coord=='miss').and.(.not.gro))then

     write(*, *)
     write(*, *)'PDB file is missing'
     write(*, *)
     stop

  end if

  if (.not.gro)then 
     
     open(1, file=coord, status='old', iostat=ierr)
     
     if(ierr /=0) then
        
        write(*, *)
        write(*, *) 'Unable to open file ', coord
        write(*, *)
        stop
        
     endif

  end if
  
  !===============GROMACS to SuAVE================================================

  if (gro) then
     
     write(*, *)
     write(*, *) "=================================================="
     write(*, *) " This action will convert all complete lines with"
     write(*, *) " just 15 columns. If necessary, complete the last"
     write(*, *) " one with zeros!"
     write(*, *) "=================================================="
     
     open(1, file=gro_file, status='old', iostat=ierr)
     
     if(ierr /=0) then
        
        write(*, *)
        write(*, *) 'Unable to open file ', gro_file
        write(*, *)
        stop
        
     endif
     
     inquire(file='single_ind.ndx', exist=ex)
     
     if (ex) then

        call execute_command_line("mv single_ind.ndx single_ind_p.ndx")
        back = .true.

     end if

     open(2, file='single_ind.ndx', status='new', iostat=ierr)

     if(ierr /= 0) then

        write(*, *)
        write(*, *) 'Unable to open file single_ind.ndx'
        write(*, *)
        stop

     end if

4    format(i10)
     
     ierr = 0
     
     do while (ierr==0)

        read(1, *, iostat=ierr) ax

        if(ierr>0) then
           
           write(*, *)
           write(*, *) 'Problem by reading ', gro_file
           write(*, *)
           stop
           
        end if

        if (ierr==0) then
           
           do i=1, 15
              
              if (ax(i)>0)then
                 
                 write(2, 4) ax(i)
                 
              end if
              
           end do
           
        end if

     end do
        
     write(*, *)
     write(*, *) 'Finished'
     write(*, *)
     write(*, *) 'Output files created'
     write(*, *) 
     
     close(1)
     close(2)
     stop
     
  end if
  
!!==========================inicio do bloco de listagem============================================

  ierr = 0
  res_i = 0

  do i=1, 500

     id(i) = 0
     
  end do
  
! Leitura do arquivo .pdb

  do while (ierr >= 0)

     read(1, 12, iostat=ierr) atom, buff%n_atom, buff%atom, &
          buff%resid, buff%ident, buff%n_resid, buff%code, buff%x, &
          buff%y, buff%z

12   format(a4, i7.1, a5, a5, a1, i4.1, a4, 3f8.3)

     if (atom.eq.'ATOM  ') then
        
        if(ierr > 0) then

           write(*, *)
           write(*, *) 'Problem reading atomic position in block 1!'
           write(*, *)
           stop

        endif

        i = 1
        igual  = .false.
        
        do while ((.not.igual).and.(i<=res_i)) ! procurando os residuos anotados

           if (buff%resid .eq. residue(i)) then
              
              igual = .true.
              
           end if

           i = i + 1
              
        end do

        if (igual)then ! se sim, entao o residuo do atomo em questao ja foi anotado

           j = 1
           igual = .false.

           do while ((.not.igual).and.(j<=id(i-1))) ! procurando os atomos anotados
              
              if (buff%atom .eq. atomid(i-1, j)) then
                 
                 igual = .true.
                 
              else
              
                 j = j + 1
                 
              end if
              
           end do
           
           if (.not.igual)then !adicionando um novo atomo a lista do resíduo em analise

              id(i-1) = id(i-1) + 1
              atomid(i-1, j) = buff%atom

           end if

        else

           res_i = i
           residue(res_i) = buff%resid !adicionado novo residuo na lista

           id(i) = id(i) + 1
           atomid(i, 1) = buff%atom

        end if

     end if

  end do

!!==========================fim do bloco de listagem=================================================

  write(*, *)
  write(*, *) " Residue groups : "
  write(*, *)

  call imprime(residue, res_i)

1 format(a28)
  write(*, *)
  write(*, 1, advance='no') " Choose the residue group : "
  read(*, *) i
  
  resid = residue(i)

  if (.not.res)then
     
     write(*, *)
     write(*, *) " Atoms inside : "
     write(*, *)
     
     call imprime(atomid(i,:), id(i))
     
2    format(a19)
     write(*, *)
     write(*, 2, advance='no') " Choose the atom : "
     read(*, *) j
     
     index = atomid(i,j)

  end if


  rewind(1)

  i_atom = 0
  n_sph = 1
  n_index = 1
  n_index2 = 1
  ierr = 0
  cent_x = 0
  cent_y = 0
  cent_z = 0
  dist_z = 0
  
! Leitura do arquivo .pdb==============================================================

  write(*, *)
  write(*, *)'Processing ..........'

  do while (ierr >= 0)

     read(1, 12, iostat=ierr) atom, buff%n_atom, buff%atom, &
          buff%resid, buff%ident, buff%n_resid, buff%code, buff%x, &
          buff%y, buff%z


     if (atom.eq.'ATOM  ') then

        if(ierr > 0) then

           write(*, *)
           write(*, *) 'Problem reading atomic position!'
           write(*, *)
           stop
           
        endif
        
        i_atom = i_atom + 1
        buff%n_atom = i_atom
        
        if (buff%resid.eq.resid)then

           if (sph) then

              sphstore(n_sph) = buff
              cent_x = (cent_x*(n_sph-1) + buff%x)/(n_sph)
              cent_y = (cent_y*(n_sph-1) + buff%y)/(n_sph)
              cent_z = (cent_z*(n_sph-1) + buff%z)/(n_sph)
              n_sph = n_sph + 1
              
           else
              
              if (res) then
                 
                 sphstore(n_sph) = buff
                 dist_z = (dist_z*(n_sph-1) + buff%z)/(n_sph)
                 n_sph = n_sph + 1

              else 
                 
                 if ((buff%atom.eq.index)) then
                    
                    sphstore(n_sph) = buff
                    dist_z = (dist_z*(n_sph-1) + buff%z)/(n_sph)
                    n_sph = n_sph + 1

                 end if
                 
              end if
              
           end if
           
        end if

     end if
     
  end do

  ! separando as lamelas ==========================
  if (.not.sph) then

     if ((.not.join)) then

       if ((bound)) then
          
3         format(a33)
          write(*, *)
          write(*, *) "Separating the index files "
          write(*, *)
          write(*, 3, advance='no') 'Write the boundary coordinate :  '
          read(*, *) dist_z
          
       else
          
          write(*, *)
          write(*, *) "Separating the index files "
          write(*, *)
          write(*, *) "Calculated boundary coordinate = ", dist_z/10, " nm"
          
       end if

    end if
     
     if (join) then
        
        do i=1, n_sph-1
           
           store(n_index) = sphstore(i)
           n_index = n_index + 1
           
        end do
        
     else
        
        do i=1, n_sph-1
           
           if (sphstore(i)%z>dist_z) then
              
              store(n_index) = sphstore(i)
              n_index = n_index + 1
              
           end if
           
           if (sphstore(i)%z<dist_z) then
              
              store2(n_index2) = sphstore(i)
              n_index2 = n_index2 + 1
              
           end if
           
        end do
        
     end if

  end if
  
  !caso esferico=======================================================================
  if (sph) then
     
     n_index = 1
     n_index2 = 1
     dist_z = 0
     
     do i=1, n_sph-1
        
        rho = sqrt((sphstore(i)%x-cent_x)**2 + (sphstore(i)%y-cent_y)**2 + (sphstore(i)%z-cent_z)**2)
        dist_z = (dist_z*(i-1) + rho)/(i)
        
     end do

     if ((.not.join)) then

       if ((bound)) then
          
          write(*, *)
          write(*, *) "Separating the index files "
          write(*, *)
          write(*, 3, advance='no') 'Write the boundary coordinate :  '
          read(*, *) dist_z
          
       else
          
          write(*, *)
          write(*, *) "Separating the index files "
          write(*, *)
          write(*, *) "Calculated boundary coordinate = ", dist_z/10, " nm"
          
       end if

    end if
     
     do i=1, n_sph-1

	rho = sqrt((sphstore(i)%x-cent_x)**2 + (sphstore(i)%y-cent_y)**2 + (sphstore(i)%z-cent_z)**2)
        
        if (res) then

           if (join) then

              store(n_index) = sphstore(i)
              n_index = n_index + 1
              
           else
              
              if (rho>dist_z) then
                 
                 store(n_index) = sphstore(i)
                 n_index = n_index + 1
                 
              end if
              
              if (rho<dist_z) then
                 
                 store2(n_index2) = sphstore(i)
                 n_index2 = n_index2 + 1
                 
              end if

           end if
           
        else

           if ((join).and.(sphstore(i)%atom.eq.index)) then

              store(n_index) = sphstore(i)
              n_index = n_index + 1

           else
              
              if ((sphstore(i)%atom.eq.index).and.(rho>dist_z)) then
                 
                 store(n_index) = sphstore(i)
                 n_index = n_index + 1
                 
              end if
              
              if ((sphstore(i)%atom.eq.index).and.(rho<dist_z)) then
                 
                 store2(n_index2) = sphstore(i)
                 n_index2 = n_index2 + 1
                 
              end if

           end if
           
        end if
        
     end do
     
  end if

! Fim da leitura inicial

  if (join) then

     inquire(file='single_ind.ndx', exist=ex)
     
     if (ex) then
        
        call execute_command_line("mv single_ind.ndx single_ind_p.ndx")
        back = .true.
        
     end if
     
     open(2, file='single_ind.ndx', status='new', iostat=ierr)
     
     if(ierr /= 0) then
        
        write(*, *)
        write(*, *) 'Unable to open file single_ind.ndx'
        write(*, *)
        stop
        
     end if

     do i=1, n_index-1
        
        write(2, 4) store(i)%n_atom
        
     end do
     
  else
     
     inquire(file='ind2.ndx', exist=ex)
     
     if (ex) then
        
        call execute_command_line("mv ind2.ndx ind2_p.ndx")
        back = .true.
        
     end if
  
     open(2, file='ind2.ndx', status='new', iostat=ierr)
     
     if(ierr /= 0) then
        
        write(*, *)
        write(*, *) 'Unable to open file ind2.ndx'
        write(*, *)
        stop
        
     end if

     inquire(file='ind1.ndx', exist=ex)
     
     if (ex) then
        
        call execute_command_line("mv ind1.ndx ind1_p.ndx")
        back = .true.
        
     end if
     
     open(3, file='ind1.ndx', status='new', iostat=ierr)
     
     if(ierr /= 0) then
        
        write(*, *)
        write(*, *) 'Unable to open file ind1.ndx'
        write(*, *)
        stop
        
     end if
     
     do i=1, n_index-1
        
        write(3, 4) store(i)%n_atom
        
     end do
     
     do i=1, n_index2-1
        
        write(2, 4) store2(i)%n_atom
        
     end do

  end if
  
! arquivo escrito                       

  write(*, *) 
  write(*, *) "Finished"
  write(*, *)
  if (back) write(*, '(a31)') ' Previous files were backed up!'

  close(1)
  close(2)

  if (.not.join)then
  
     close(3)
     
  end if
     
end program separa_indice

subroutine imprime(vetor, indice)

    character(len=30), intent(in) :: vetor(500)
    integer, intent(in) :: indice

5   format(i4, a10)
    
    do i=1, int((indice+4)/5)
       
       do j=1, 5
          
          if ((i-1)*5+j<=indice)then
             
             write(*, 5, advance='no') (i-1)*5+j, vetor((i-1)*5+j)
             
          end if
          
       end do
       
       write(*,*)
       
    end do
    
end subroutine imprime
