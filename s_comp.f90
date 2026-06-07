program comp

  use types
  use variables
  use funcproc

  integer :: N_sample, n_rand, locate, n_Ka, n_kc
  double precision, dimension(10000000) :: Ka_true, Ka_proj, kc, kc2
  
  call startup(outer, bin, p_grid, coord, ind, ind2, rmsd, map, ind3, &
     l_coarse, begin, end, skip, lipid, rough, slices, inside, range, &
     n_grid, bin_out, fr_in, fr_end, n_skip, n_lipid, get, div, 's_comp    ', version)

  back = .false.
  signal2 = 'miss'
  signal = 'miss'
  ac_time = 0
  n = 1
  temp = 300 ! DEFAULT
  
  do i=1, 20
     
     call getarg(i, get(i))
     
  end do

  do i=1, 20
     
     if (get(i)=='-true')then
        
        signal = get(i+1)

     end if

     if (get(i)=='-proj') then

        signal2 = get(i+1)
        
     end if
     
     if (get(i)=='-size')then

        read(get(i+1), *, iostat=ierr) size
                        
     end if
     
     if (get(i)=='-temp') then

        read(get(i+1), *, iostat=ierr) temp
        
     end if

     if (get(i)=='-ac')then
        
	read(get(i+1), *, iostat=ierr) ac_time
 
     end if
     
  end do

  !==============================================================
  if (signal=='miss')then

     write(*, *)
     write(*, *)'True area input file is missing'
     write(*, *)
     stop

  end if

  if (ac_time==0)then
     
     write(*, *)
     write(*, *)'You must provide the autocorr. time [in frames]'
     write(*, *)
     stop
     
  end if
  
  call abre_trj(1, signal)

  if (signal2=='miss')then

     write(*, *)
     write(*, *)'Projected area input file is missing'
     write(*, *)
     stop

  end if

  call abre_trj(3, signal2)
  !=============================================================

  call abre('Ka        ', 2, 'xvg', back)
    
  write(2, '(a7, a5)') "#SuAVE ", version
  write(2, '(a14)') '#Command Line:'
  write(2, '(a9)', advance='no') '#s_stat  '
  write(2, *) (trim(get(i)),"  ", i=1, 20)
  write(2, *) '@    title "Area Compressibility"'
  write(2, *) '@    xaxis  label "Moving Block Bootstrap Sample"'
  write(2, *) '@    yaxis  label "K\sA\N [N/m]"'


  call abre('kc        ', 4, 'xvg', back)

  write(4, '(a7, a5)') "#SuAVE ", version
  write(4, '(a14)') '#Command Line:'
  write(4, '(a9)', advance='no') '#s_comp  '
  write(4, *) (trim(get(i)),"  ", i=1, 20)
  write(4, *) '@    title "Bending Modulus"'
  write(4, *) '@    xaxis  label "Moving Block Bootstrap Sample"'
  write(4, *) '@    yaxis  label "kc [J]"'
  
  !=============================================================

  do i=1, 10000000

     func(i) = 0
     func2(i) = 0
     
  end do

  do i=1, 1000

     hist(i) = 0

  end do
     
  !===================================

  n_index = 1
  
  do while (ierr>=0)
     
     read(1, *, iostat=ierr) aux, aux_true
     read(3, *, iostat=ierr) aux, aux_proj
     
     if (ierr == 0)then

        func(n_index) = aux_true
        func2(n_index) = aux_proj
        n_index = n_index + 1

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
  
  !===================================
  ! calculando parametro de compressibilidade
  ! e modulo de elasticidade

  j = 0
  maxf = -100000
  minf = 100000
  n_Ka = 0
  n_kc = 0
  
  do i=size, n_index

     n_Ka = n_Ka + 1
     call calc_running_aver(aver, aver2, desv2, func, size, i)
     Ka_true(n_Ka) = aver*kb*temp/(desv2*1.0e-18) ! valor em N/m

     call calc_running_aver(aver, aver2, desv2, func2, size, i)
     Ka_proj(n_Ka) = aver*kb*temp/(desv2*1.0e-18) ! valor em N/m

     A0 = aver*1E-18          ! corrigindo para m2
     
     aux = 1/Ka_true(n_Ka) - 1/Ka_proj(n_Ka)

     !guardando os valores para avaliar o filtro e remoção de outliers

     if ((aux<0).and.(aux>-10000)) then

        n_kc = n_kc + 1
        kc(n_kc) = sqrt(-(A0*kb*temp)/(32*pi*pi*pi*aux)) ! valor em J
        
        maxf = max(maxf, kc(n_kc))
        minf = min(minf, kc(n_kc))
           
     end if
     
  end do

  !criando histograma ===========================
  !aqui já temos o perfil de kc, e agora vamos extrair
  !informações estatísticas para termos o IQR como
  !filtro de outliers ===========================

  del = (maxf - minf)/100

  do i=1, n_kc
     
     bini = nint((kc(i)-minf)/del) + 100
     hist(bini) = hist(bini) + 1/(n_kc*del)
     
  end do
  
  !Calculando o IQR
  call calc_stat_IQR(n_kc, hist, del, minf, quart1, quart3)
  IQR = quart3 - quart1
  
  k = 0
  
  do i=1, n_kc

     if ((kc(i) > quart1 - 1.5*IQR).and.(kc(i) < quart3 + 1.5*IQR)) then

        k = k + 1
        kc2(k) = kc(i)
        
     end if

  end do

  n_kc = k ! atualizando o número de pontos em kc

  call mbb(Ka_true, n_Ka, ac_time, int(1000*start/clock_rate), 2)
  call mbb(kc2, n_kc, ac_time, int(1000*start/clock_rate), 4)
  
  !==================================

  call system_clock(finish, clock_rate, clock_max)
  
  close(1)
  close(2)
  close(3)
  close(4)

  call ending(back, finish, start, clock_rate) ! Finaliza programa e mostra tempo de processamento
    
end program comp

