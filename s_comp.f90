program comp

  use types
  use variables
  use funcproc

  call startup(outer, bin, p_grid, coord, ind, ind2, rmsd, map, ind3, &
     l_coarse, begin, end, skip, lipid, rough, slices, inside, range, &
     n_grid, bin_out, fr_in, fr_end, n_skip, n_lipid, get, div, 's_comp    ', version)

  back = .false.
  signal2 = 'miss'
  signal = 'miss'
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
     
  end do

  !==============================================================
  if (signal=='miss')then

     write(*, *)
     write(*, *)'True area input file is missing'
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
  write(2, *) '@    xaxis  label "Frame"'
  write(2, *) '@    yaxis  label "KA [N/m]"'


  call abre('kc        ', 4, 'xvg', back)

  write(4, '(a7, a5)') "#SuAVE ", version
  write(4, '(a14)') '#Command Line:'
  write(4, '(a9)', advance='no') '#s_comp  '
  write(4, *) (trim(get(i)),"  ", i=1, 20)
  write(4, *) '@    title "Area Compressibility"'
  write(4, *) '@    xaxis  label "Sample"'
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

  do i=size, n_index

     call calc_running_aver(aver, aver2, desv2, func, size, i)
     Ka_true = aver*kb*temp/(desv2*1.0e-18) ! valor em N/m

     call calc_running_aver(aver, aver2, desv2, func2, size, i)
     Ka_proj = aver*kb*temp/(desv2*1.0e-18) ! valor em N/m
          
     write(2, *) i, Ka_true, Ka_proj


     A0 = aver*1E-18          ! corrigindo para m2
     
     aux = 1/Ka_true - 1/Ka_proj

     !guardando os valores para avaliar o filtro e remoção de outliers

     if ((aux<0).and.(aux>-10000)) then

        j = j + 1
        kc = sqrt(-(A0*kb*temp)/(16.6*pi*pi*pi*aux)) ! valor em J
        kc_v1(j) = kc

        maxf = max(maxf, kc)
        minf = min(minf, kc)
           
     end if
     
  end do

  !criando histograma ===========================
  !aqui já temos o perfil de kc, e agora vamos extrair
  !informações estatísticas para termos o IQR como
  !filtro de outliers ===========================

  del = (maxf - minf)/100

  do i=1, j
     
     bini = nint((kc_v1(i)-minf)/del) + 100
     hist(bini) = hist(bini) + 1/(j*del)
     
  end do
  
  !Calculando o IQR
  call calc_stat_IQR(j, hist, del, minf, quart1, quart3)
  IQR = quart3 - quart1
  
  k = 0
  
  do i=1, j

     if ((kc_v1(i) > quart1 - 1.5*IQR).and.(kc_v1(i) < quart3 + 1.5*IQR)) then

        k = k + 1
        kc_v2(k) = kc_v1(i)
        write(4, *) k, kc_v2(k)
        
     end if

  end do
  
  !==================================

  call system_clock(finish, clock_rate, clock_max)
  
  close(1)
  close(2)
  close(3)
  close(4)

  call ending(back, finish, start, clock_rate) ! Finaliza programa e mostra tempo de processamento
    
end program comp

