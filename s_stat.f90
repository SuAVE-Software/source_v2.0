program suave_stat

  use types
  use variables
  use funcproc

  call startup(outer, bin, p_grid, coord, ind, ind2, rmsd, map, ind3, &
     l_coarse, begin, end, skip, lipid, rough, slices, inside, range, &
     n_grid, bin_out, fr_in, fr_end, n_skip, n_lipid, get, div, 's_stat    ', version)

  back = .false.
  signal = 'miss'
  
  do i=1, 20
     
     call getarg(i, get(i))
     
  end do

  do i=1, 20
     
     if (get(i)=='-in')then
        
        signal = get(i+1)

     end if
     
  end do

  if (signal=='miss')then

     write(*, *)
     write(*, *)'Input file is missing'
     write(*, *)
     stop

  end if

  call abre_trj(1, signal)

  call abre('pdf       ', 2, 'xvg', back)
  
  write(2, '(a7, a7)') "#SuAVE ", version
  write(2, '(a14)') '#Command Line:'
  write(2, '(a9)', advance='no') '#s_stat  '
  write(2, *) (trim(get(i)),"  ", i=1, 20)
  write(2, *) '@    title "Probability Density Function"'
  write(2, *) '@    xaxis  label "x"'
  write(2, *) '@    yaxis  label "PDF(x)"'
  write(2, *) '@    s0 comment "pdf.xvg"'
  write(2, *) '@    s0 legend  "Original Data"'
  write(2, *) '@    s1 comment "pdf.xvg"'
  write(2, *) '@    s1 legend  "Gaussian Model"'
  write(2, *) '@    s1 line linewidth 2.0'
  write(2, *) '@    s2 comment "pdf.xvg"'
  write(2, *) '@    s2 legend  "Lower Quartile"'
  write(2, *) '@    s2 symbol 1'
  write(2, *) '@    s2 symbol color 3'
  write(2, *) '@    s2 symbol fill pattern 1'
  write(2, *) '@    s3 comment "pdf.xvg"'
  write(2, *) '@    s3 legend  "Upper Quartile"'
  write(2, *) '@    s3 symbol 1'
  write(2, *) '@    s3 symbol color 4'
  write(2, *) '@    s3 symbol fill pattern 1'
  write(2, *) '@    s4 comment "pdf.xvg"'
  write(2, *) '@    s4 legend  "Lower Decile"'
  write(2, *) '@    s4 symbol 1'
  write(2, *) '@    s4 symbol color 5'
  write(2, *) '@    s4 symbol fill pattern 1'
  write(2, *) '@    s5 comment "pdf.xvg"'
  write(2, *) '@    s5 legend  "Upper Decile"'
  write(2, *) '@    s5 symbol 1'
  write(2, *) '@    s5 symbol color 6'
  write(2, *) '@    s5 symbol fill pattern 1'

  call abre('acf       ', 3, 'xvg', back)  

  write(3, '(a7, a7)') "#SuAVE ", version
  write(3, '(a14)') '#Command Line:'
  write(3, '(a9)', advance='no') '#s_stat  '
  write(3, *) (trim(get(i)),"  ", i=1, 20)
  write(3, *) '@    title "Autocorrelation Function X Frame"'
  write(3, *) '@    xaxis  label "Frame"'
  write(3, *) '@    yaxis  label "ACF(frame)"'

  !===================================

  do i=1, 10000000

     func(i) = 0

  end do
  
  do i=1, 1000

     hist(i) = 0

  end do

  !===================================

  n_index = 1
  maxf = -1000000
  minf = 1000000
  
  do while (ierr>=0)
     
     read(1, *, iostat=ierr) aux, aux2

     if (ierr == 0)then

        func(n_index) = aux2
        n_index = n_index + 1
        maxf = max(aux2, maxf)
        minf = min(aux2, minf)

     end if
     
  end do

  n_index = n_index - 1
  
  if(n_index>10000000) then

     write(*, *)
     write(*, *) ' Too many points'
     write(*, *)
     stop

  end if
  
  call system_clock(start, clock_rate, clock_max)
  
  !===========================================!
  ! calculando a medidas de tendência central !
  ! e dispersão                               !
  !===========================================!

  del = (maxf-minf)/800
  call calc_stat_aver(aver, aver2, desv, skew, kurt, st_mom, n_index, func)
  call do_histogram(n_index, hist, minf, desv, aver, aux, del, func, 2)
  call calc_stat_all(n_index, hist, del, minf, skew, kurt, aver, desv)  
  call calc_acf(n_index, aver, desv, func)

  call system_clock(finish, clock_rate, clock_max)
  
  close(1)
  close(2)

  call ending(back, finish, start, clock_rate) ! Finaliza programa e mostra tempo de processamento
  
end program suave_stat