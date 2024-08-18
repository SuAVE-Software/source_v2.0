program suave_filter

  use funcproc
  
  ! Permanecendo com as variaveis para facilitar =======
  real, parameter :: pi = 3.141592654
  character(len=5), parameter :: version = "2.0.0"
  
  logical :: ex, help, hfilter, lfilter, sfilter, back
  
  character(len=30) :: get(20), signal

  integer :: n_index, i, j, ierr
  integer :: start, finish, clock_rate, clock_max

  real :: aux, aux2, low, lfreq, hfreq, inicial 
  real :: hour, minu, sec, paramet, dt, dw, t, w
  real, dimension(3000000) ::func, re, im, ref, imf
  ! ====================================================

  call startup_filter(back, signal, hfilter, lfilter, sfilter, paramet, &
       hfreq, lfreq, low, get, "s_filter  ", version)
  
  open(1, file=signal, status='old', iostat=ierr)
  
  if(ierr /= 0) then
     
     write(*, *)
     write(*, *) ' Unable to open file ', signal
     write(*, *)
     stop
     
  endif
  
  if (.not.sfilter) then

     call abre('Ftransf   ', 2, 'xvg', back)
          
     write(2, '(a7, a7)') "#SuAVE ", version
     write(2, '(a14)') '#Command Line:'
     write(2, '(a11)', advance='no') '#s_filter  '
     write(2, *) (trim(get(i)),"  ", i=1, 20)
     write(2, *) '@    title "Discrete Fourier Transform"'
     write(2, *) '@    xaxis  label "w [rad/(time unity)]"'
     write(2, *) '@    yaxis  label "|F(w)|"'
     

     call abre('Ffiltered ', 3, 'xvg', back)
          
     write(3, '(a7, a7)') "#SuAVE ", version
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

     call abre('Sfiltered ', 4, 'xvg', back)
          
     write(4, '(a7, a7)') "#SuAVE ", version
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

     call filter_suave(n_index, paramet, func, dt, inicial)
     
     close(4)
     close(1)
          
  else
  
     !===================================
     ! calculando a transformada
     
     call filter_ft(n_index, dt, dw, inicial, func, re, im)
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

     call filter_ift(n_index, dt, inicial, func, re, im, ref, imf)

     close(1)
     close(3)

  end if ! if (s_filter) then

  call system_clock(finish, clock_rate, clock_max)

  call ending(back, finish, start, clock_rate) ! Finaliza programa e mostra tempo de processamento

end program suave_filter

