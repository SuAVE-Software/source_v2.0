subroutine startup(outer, bin, p_grid, coord, ind, ind2, rmsd, map, ind3, &
     l_coarse, begin, end, skip, lipid, rough, slices, inside, range, &
     n_grid, bin_out, fr_in, fr_end, n_skip, n_lipid, get, div, name, version)
  
  logical :: bin, outer, p_grid, rmsd, l_coarse, begin, end, skip
  logical :: lipid, map, slices, inside, range
  
  integer :: i, ierr, div
  integer :: n_grid, bin_out, n_lipid, fr_in, fr_end, n_skip
  
  real :: rough

  character(len=7) :: version
  character(len=10) :: name
  character(len=30) :: coord,  ind, ind2, ind3 
  character(len=30), dimension(30) :: get    
  
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
  write(*, *) "                    ** "//trim(name)//" **"
  write(*, *) ""
  write(*, *) "               ** VERSION ", version, " **"
  write(*, *) ""
  write(*, *) "            Santos, D. E. S.; Soares, T. A."
  write(*, *) ""
  write(*, *) "Colaborators: "
  write(*, *) "Coutinho, K.        Caetano, D. L. Z"
  write(*, *) ""
  write(*, *) "Please cite "
  write(*, *)
  write(*, *) "Santos, D. E. S.; Coutinho, K.; Soares, T. A. (2022) Surface "
  write(*, *) "Assessment via Grid Evaluation (SuAVE) for Every Surface Curvature"
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

  slices = .false.
  outer = .false.
  bin = .false.
  inside = .false.
  p_grid = .false.
  coord = 'miss'
  ind = 'miss'
  ind2 = 'miss'
  ind3 = 'miss'
  rmsd = .false.
  l_coarse = .false.
  map = .false.
  begin = .false.
  end = .false.
  skip = .false.
  range = .false.
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

     if (get(i)=='-dens')then

        ind3 = get(i+1)

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

     if (get(i)=='-inside')then

        inside = .true.
        
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

     if (get(i)=='-slices')then

        slices = .true.
        read(get(i+1), *, iostat=ierr) div
        call erro(ierr, "SLICE")

     end if

     if (get(i)=='-map')then

        map = .true.

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
        call erro(ierr, "ROUGH")

     end if
          
     if (get(i)=='-help')then
        
        call write_help(name(3:10))
        
     end if

     if (get(i)=='-lipid')then
        
        lipid = .true.
        read(get(i+1), *, iostat=ierr) n_lipid
        call erro(ierr, "LIPID")
        
     end if

     if (get(i)=='-range')then
        
        range = .true.
        
     end if
     
  end do
  
  
  ! Realizando consideracoes ==============================
  
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
          
  end if
    
end subroutine startup

subroutine startup_filter(back, signal, hfilter, lfilter, sfilter, paramet, &
     hfreq, lfreq, low, get, name, version)
  
  logical :: help, hfilter, lfilter, sfilter, back
  
  character(len=30) :: get(20), signal
  character(len=10) :: name
  character(len=7) :: version
  
  integer :: i, ierr
  
  real :: low, lfreq, hfreq, paramet


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
  write(*, *) "                    ** "//trim(name)//" **"
  write(*, *) ""
  write(*, *) "               ** VERSION ", version, " **"
  write(*, *) ""
  write(*, *) "            Santos, D. E. S.; Soares, T. A."
  write(*, *) ""
  write(*, *) "  Colaborators: "
  write(*, *) "  Caetano, D. L. Z"
  write(*, *) ""
  write(*, *) "Please cite "
  write(*, *)
  write(*, *) "Santos, D. E. S.; Coutinho, K.; Soares, T. A. (2022) Surface "
  write(*, *) "Assessment via Grid Evaluation (SuAVE) for Every Surface Curvature"
  write(*, *) "and Cavity Shape. Journal of Chemical Information and Modeling,"
  write(*, *) "v. 62, p. 4690–4701"
  write(*, *)
  write(*, *) "Santos, D. E. S.; Pontes, J. F. S.; Lins, R. D.; Coutinho, K.; "
  write(*, *) "Soares, T. A. (2020) SuAVE: A Tool for Analyzing Curvature-Dependent"
  write(*, *) "Properties in Chemical Interfaces. Journal of Chemical Information "
  write(*, *) "and Modeling, v. 60, p. 473-484."

1 format(a19, i6.3, a2, i6.3, a4, i6.3, a4)

  back = .false.
  signal = 'miss'
  help = .false.
  hfilter = .false.
  lfilter = .false.
  sfilter = .false.
  hfreq = 1000000000000000.0
  lfreq = -1000000000000000.0
  low = -1.0

  do i=1, 20
     
     call getarg(i, get(i))
     
  end do

  do i=1, 20
     
     if (get(i)=='-in')then
        
        signal = get(i+1)

     end if
     
     if (get(i)=='-hfreq')then
        
        hfilter = .true.
        read(get(i+1), *, iostat=ierr) hfreq
        call erro(ierr, "HFREQ ")

     end if

     if (get(i)=='-lfreq')then

        lfilter = .true.
        read(get(i+1), *, iostat=ierr) lfreq
        call erro(ierr, "LFREQ ")

     end if

     if (get(i)=='-lowint')then

        read(get(i+1), *, iostat=ierr) low
        call erro(ierr, "LOWINT")
        
     end if

     if (get(i)=='-param')then
        
        sfilter = .true.
        read(get(i+1), *, iostat=ierr) paramet
        call erro(ierr, "PARAM ")
        
     end if

     if (get(i)=='-help')then

        call write_help(name(3:10))

     end if
     
  end do

  if (signal=='miss')then

     write(*, *)
     write(*, *)'Input file is missing'
     write(*, *)
     stop

  end if

end subroutine startup_filter
  
subroutine erro(i, name)
  
  integer, intent(in) :: i
  character(len=5), intent(in) :: name
  character(len=60) :: aux
  
  if (i/=0) then
     
     aux = '    Problem by passing '//name//' parameter'
     write(*, *)
     write(*, *) "        ==================ERROR==================="
     write(*, *) "        "//aux
     write(*, *) "        ==================ERROR==================="
     write(*, *)
     stop
     
  end if
  
end subroutine erro
