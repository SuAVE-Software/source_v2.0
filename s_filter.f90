program fourier

  real, parameter :: pi = 3.141592654
  character(len=5), parameter :: version = "2.0.0"
  
  logical :: ex, help, hfilter, lfilter, sfilter, back
  
  character(len=30) :: get(20), signal

  integer :: n_index, i, j, ierr
  integer :: start, finish, clock_rate, clock_max

  real :: aux, aux2, low, lfreq, hfreq, inicial 
  real :: hour, minu, sec, param, dt, dw, t, w
  real, dimension(3000000) ::func, re, im, ref, imf


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
  write(*, *) "                       ** s_filter **"
  write(*, *) ""
  write(*, *) "                   ** VERSION ", version, " **"
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
        read(get(i+1), *, iostat=ierr) param
        call erro(ierr, "PARAM ")
        
     end if

     if (get(i)=='-help')then

        help = .true.

     end if

  end do

  !HELP begins ============================================================
  if (help) then

     write(*, *) ""
     write(*, *) ""
     write(*, *) "s_filter performs a Discrete Fourier Transform (DFT)"
     write(*, *) "on the signal described in the input file and performs"
     write(*, *) "filtering treatments"
     write(*, *) ""
     write(*, *) "Usage: s_filter -in file.xvg"
     write(*, *) ""
     write(*, *) "file.xvg ---- contains a SuAVE output"
     write(*, *) ""
     write(*, *) "Options: "
     write(*, *) ""
     write(*, *) "-hfreq           defines the highest angular frequency to be"
     write(*, *) "used in the Inverse Fourier Transform (default: no filter)"
     write(*, *) ""
     write(*, *) "-lfreq           defines the lowest angular frequency to be"
     write(*, *) "used in the Inverse Fourier Transform (default: no filter)"
     write(*, *) ""
     write(*, *) "-lowint          defines the lowest intensity to be used in the"
     write(*, *) "Inverse Fourier Transform (default: -1)"
     write(*, *) ""
     write(*, *) "-param           defines the parameter used on the unidimensional" 
     write(*, *) "SuAVE fitting process. In this case the DFT is not applied"
     write(*, *) ""
     write(*, *) "-help            prints HELP information and quits"
     stop

  end if
  !HELP ends ==============================================================

  if (signal=='miss')then

     write(*, *)
     write(*, *)'Input file is missing'
     write(*, *)
     stop

  end if

  open(1, file=signal, status='old', iostat=ierr)
  
  if(ierr /= 0) then
     
     write(*, *)
     write(*, *) ' Unable to open file ', signal
     write(*, *)
     stop
     
  endif

  if (.not.sfilter) then
     
     inquire(file='Ftransf.xvg', exist=ex)
     
     if (ex) then
        
        call execute_command_line("mv Ftransf.xvg Ftransf_p.xvg")
        back = .true.
        
     end if
  
     open(2, file='Ftransf.xvg', status='new', iostat=ierr)
     
     if(ierr /= 0) then
        
        write(*, *)
        write(*, *) ' Unable to open file Ftransf.xvg'
        write(*, *)
        stop
        
     endif
     
     write(2, '(a7, a5)') "#SuAVE ", version
     write(2, '(a14)') '#Command Line:'
     write(2, '(a11)', advance='no') '#s_filter  '
     write(2, *) (trim(get(i)),"  ", i=1, 20)
     write(2, *) '@    title "Discrete Fourier Transform"'
     write(2, *) '@    xaxis  label "w [rad/(time unity)]"'
     write(2, *) '@    yaxis  label "|F(w)|"'
     
     
     inquire(file='Ffiltered.xvg', exist=ex)
     
     if (ex) then
        
        call execute_command_line("mv Ffiltered.xvg Ffiltered_p.xvg")
        back = .true.
        
     end if
     
     open(3, file='Ffiltered.xvg', status='new', iostat=ierr)
     
     if(ierr /= 0) then
        
        write(*, *)
        write(*, *) ' Unable to open file Ffiltered.xvg'
        write(*, *)
        stop
        
     endif
     
     write(3, '(a7, a5)') "#SuAVE ", version
     write(3, '(a14)') '#Command Line:'
     write(3, '(a11)', advance='no') '#s_filter  '
     write(3, *) (trim(get(i)),"  ", i=1, 20)
     write(3, *) '@    title "Inverse Fourier Transform"'
     write(3, *) '@    xaxis  label "t"'
     write(3, *) '@    yaxis  label "f(t)"'
     
  end if

  do i=1, 3000000
     
     func(i) = 0
     re(i) = 0
     im(i) = 0

  end do
  
  if (sfilter) then

     inquire(file='Sfiltered.xvg', exist=ex)
     
     if (ex) then
        
        call execute_command_line("mv Sfiltered.xvg Sfiltered_p.xvg")
        back = .true.
        
     end if
     
     open(4, file='Sfiltered.xvg', status='new', iostat=ierr)
     
     if(ierr /= 0) then
        
        write(*, *)
        write(*, *) ' Unable to open file Sfiltered.xvg'
        write(*, *)
        stop
        
     endif
     
     write(4, '(a7, a5)') "#SuAVE ", version
     write(4, '(a14)') '#Command Line:'
     write(4, '(a11)', advance='no') '#s_filter  '
     write(4, *) (trim(get(i)),"  ", i=1, 20)
     write(4, *) '@    title "SuAVE Unidimensional Fitting Curve"'
     write(4, *) '@    xaxis  label "t"'
     write(4, *) '@    yaxis  label "f(t)"'
     
  end if

  !===================================

  n_index = 0
  
  do while (ierr>=0)

     read(1, *, iostat=ierr) aux, aux2
     
     if (ierr == 0)then

        if(n_index == 0) inicial = aux

        n_index = n_index + 1
        func(n_index) = aux2

     end if

  end do

  dt = (aux-inicial)/(n_index-1) ! dividindo pelo numero de intervalos

  if(n_index>1000000) then

     write(*, *)
     write(*, *) ' Too many points (>1000000)'
     write(*, *)
     stop

  end if

  call system_clock(start, clock_rate, clock_max)
  
  !SuAVE fitting Curve ===============

  if (sfilter)  then

     do i=1, n_index

        aux = 0
        sum = 0
        
        do j=1, n_index
           
           aux = aux + func(j)*exp(-param*param*(j-i)*(j-i))
           sum = sum + exp(-param*param*(j-i)*(j-i))
           
        end do
        
        aux = aux/sum
        write(4, *) (i-1)*dt+inicial, aux
        
     end do

     close(4)
     close(1)
          
  else
  
     !===================================
     ! calculando a transformada
     
     dw = (2*pi/dt)/(n_index) 

     !$OMP parallel do 
     do i=1, n_index 

        w = -pi/dt + (i-1)*dw ! com isso w vai de -pi/dt até pi/dt
        
        re(i) = 0
        im(i) = 0
        
        do j=1, n_index 
           
           t = inicial + (j-1)*dt

           re(i) = re(i) + func(j)*cos(w*t)
           im(i) = im(i) - func(j)*sin(w*t)

        end do
        
     end do
     !$OMP end parallel do 
     
     do i=1, n_index
        
        w = -pi/dt + (i-1)*dw

        write(2, *, iostat=ierr) w, sqrt(re(i)*re(i) + im(i)*im(i))
        
        if(ierr>0) then
           
           write(*, *)
           write(*, *) ' Problem by writing output '
           write(*, *)
           stop
           
        end if
        
     end do
     
     close(2)
     
     !==inserindo FILTROS =========================
     !$OMP parallel do
     do i=1, n_index
        
        if (sqrt(re(i)*re(i) + im(i)*im(i))<=low) then
           
           re(i) = 0
           im(i) = 0
           
        end if
        
        if ((hfilter).and.(-pi/dt + (i-1)*dw>hfreq)) then
           
           re(i) = 0
           im(i) = 0
           
        end if
        
        if ((lfilter).and.(-pi/dt + (i-1)*dw<lfreq)) then
           
           re(i) = 0
           im(i) = 0
           
        end if
        
     end do
     !$OMP end parallel do
     
     !=============================================
     ! calculando a funcao atraves da transformada
     !$OMP parallel do
     do j=1, n_index 

        t = inicial + (j-1)*dt

        ref(j) = 0
        imf(j) = 0
        
        do i=1, n_index 
           
           w = -pi/dt + (i-1)*dw ! com isso w vai de -pi/dt até pi/dt

           ref(j) = ref(j) + re(i)*cos(w*t) - im(i)*sin(w*t)
           imf(j) = imf(j) + re(i)*sin(w*t) + im(i)*cos(w*t)
           
        end do
        
     end do
     !$OMP end parallel do
     
     do j=1, n_index 
        
        t = inicial + (j-1)*dt
        write(3, *, iostat=ierr) t, sinal(ref(j))*sqrt(ref(j)*ref(j) + imf(j)*imf(j))/(n_index) 
        
        if(ierr>0) then
           
           write(*, *)
           write(*, *) 'Problem by writing output '
           write(*, *)
           stop
           
        end if
        
     end do
     
     close(1)
     close(3)

  end if ! if (s_filter) then

  call system_clock(finish, clock_rate, clock_max)

  write(*, *) 
  write(*, *) "Finished"
  write(*, *)
  if (back) write(*, '(a31)') ' Previous files were backed up!'

  hour = real(finish-start)/(3600*real(clock_rate))
  minu = (hour - int(hour))*60
  sec = (minu-int(minu))*60

  if (back) write(*, *)
  write(*, 1) " Processing time : ", int(hour), " h", int(minu), " min", int(sec), " sec"
  write(*, *)
  
end program fourier


subroutine erro(i, name)

  integer, intent(in) :: i
  character(len=6), intent(in) :: name
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


real function sinal(re)

  real :: re

  sinal = 1

  if (re<0) then
     
     sinal = -1

  end if

end function sinal
