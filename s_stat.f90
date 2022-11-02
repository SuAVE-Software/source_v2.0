program fourier

  real, parameter :: pi = 3.141592654
  character(len=5), parameter :: version = "2.0.0"
  
  logical :: ex, help, back
  
  character(len=30) :: get(20), signal

  integer :: n_index, i, j, ierr, bini, class_modal
  integer :: class_med, class_q1, class_q3, class_d1, class_d9

  double precision :: aux, aux2, aver, aver2, desv, kurt, delta1, delta2
  double precision :: del, hist(1000), sum, acum, moda, mediana, st_mom
  double precision :: acum25, acum75, quart1, quart3, acum10, acum90
  double precision :: decil1, decil9, skew, acf, maxf, minf

  real, dimension(10000000) ::func


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
  write(*, *) "                       ** s_stat **"
  write(*, *) ""
  write(*, *) "                   ** VERSION ", version, " **"
  write(*, *) ""
  write(*, *) "              Santos, D. E. S.; Soares, T. A."
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


  back = .false.
  signal = 'miss'
  help = .false.

  do i=1, 20
     
     call getarg(i, get(i))
     
  end do

  do i=1, 20
     
     if (get(i)=='-in')then
        
        signal = get(i+1)

     end if
     
     if (get(i)=='-help')then

        help = .true.

     end if

  end do

  !HELP begins ============================================================
  if (help) then

     write(*, *) ""
     write(*, *) ""
     write(*, *) "s_stat performs central tendency and dispersion measurements,"
     write(*, *) "as well as assimetry and kurtosis upon the Probability "
     write(*, *) "Density Function describing the input data."
     write(*, *) ""
     write(*, *) "Usage: s_stat -in file.xvg"
     write(*, *) ""
     write(*, *) "file.xvg ---- contains a SuAVE output"
     write(*, *) ""
     write(*, *) "Options: "
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

  inquire(file='pdf.xvg', exist=ex)

  if (ex) then
     
     call execute_command_line("mv pdf.xvg pdf_p.xvg")
     back = .true.
     
  end if
  
  open(2, file='pdf.xvg', status='new', iostat=ierr)
  
  if(ierr /= 0) then
     
     write(*, *)
     write(*, *) ' Unable to open file pdf.xvg'
     write(*, *)
     stop
     
  endif
  
  write(2, '(a7, a5)') "#SuAVE ", version
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

  inquire(file='acf.xvg', exist=ex)

  if (ex) then

     call execute_command_line("mv acf.xvg acf_p.xvg")
     back = .true.

  end if

  open(3, file='acf.xvg', status='new', iostat=ierr)

  if(ierr /= 0) then

     write(*, *)
     write(*, *) ' Unable to open file acf.xvg'
     write(*, *)
     stop

  endif

  write(3, '(a7, a5)') "#SuAVE ", version
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
  
  
  !===================================
  ! calculando a medidas de tendência central
  ! e dispersão

  aver = 0
  aver2 = 0
  desv = 0

  do i=1, n_index

     aver = aver + func(i)/n_index
     aver2 = aver2 + func(i)*func(i)/n_index

  end do
  
  desv = sqrt((aver2 - aver*aver))
  
  write(*, *)
  write(*, '(a45)')"****************** SUMMARY ******************"
  write(*, '(a45)')"*                                           *"
  write(*, '(a13, f10.3, a22)')"* Average  = ", aver, "                   *"
  write(*, '(a13, f10.3, a22)')"* S. dev.  = ", desv, "                   *"

  skew = 0
  kurt = 0
  do i=1, n_index

     st_mom = (func(i) - aver)/desv
     kurt = kurt + st_mom*st_mom*st_mom*st_mom/n_index
     skew = skew + st_mom*st_mom*st_mom/n_index

  end do

  del = (maxf-minf)/800

  !criando histograma ===========================
  do i=1, n_index

     bini = nint((func(i)-minf)/del)+100
     hist(bini) = hist(bini) + 1/(n_index*del)
     
  end do

       !Original =========
  do i=1, 1000

     write(2, *) (i-100)*del+minf, hist(i)

  end do

  write(2, *)'&'

       !Modelo Gaussiano ========
  do i=1, 1000

     aux = (i-100)*del+minf
     aux = exp(-(aux-aver)*(aux-aver)/(2*desv*desv))
     aux = aux/(desv*sqrt(2*pi))
     write(2, *) (i-100)*del+minf, aux

  end do
  !==============================================

  !Calculando mediana, moda, quartil e decil=====
  sum = 0
  class_modal = 1
  class_med = 1
  class_q1 = 1
  class_q3 = 1
  class_d1 = 1
  class_d9 = 1

  do i=1, 1000

     if (hist(i)>hist(class_modal)) class_modal = i
     
     sum = sum + hist(i)*del*n_index

     if (sum<n_index/2)then

        class_med = i+1
        acum = sum

     end if

     if (sum<n_index/4)then

        class_q1 = i+1
        acum25 = sum

     end if

     if (sum<3*n_index/4)then

        class_q3 = i+1
        acum75 = sum

     end if

     if (sum<n_index/10)then

        class_d1 = i+1
        acum10 = sum

     end if

     if (sum<9*n_index/10)then

        class_d9 = i+1
        acum90 = sum

     end if

  end do

  !moda de Czuber para dados tabulados====
  delta1 = hist(class_modal) - hist(class_modal-1)
  delta2 = hist(class_modal) - hist(class_modal+1)
  moda = (((class_modal-100)*del+minf)-del/2)
  moda = moda + del*(delta1/(delta1+delta2))

  !mediana para dados tabulados===========
  mediana = (((class_med-100)*del+minf)-del/2)
  mediana = mediana + (n_index/2 - acum)/(hist(class_med)*n_index)

  !primeiro Quartil para dados tabulados =
  quart1 = (((class_q1-100)*del+minf)-del/2)
  quart1 = quart1 + (n_index/4 - acum25)/(hist(class_q1)*n_index)

  !terceiro Quartil para dados tabulados =
  quart3 = (((class_q3-100)*del+minf)-del/2)
  quart3 = quart3 + (3*n_index/4 - acum75)/(hist(class_q3)*n_index)

  !primeiro Decil para dados tabulados====
  decil1 = (((class_d1-100)*del+minf)-del/2)
  decil1 = decil1 + (n_index/10 - acum10)/(hist(class_d1)*n_index)

  !nono Decil para dados tabulados========
  decil9 = (((class_d9-100)*del+minf)-del/2)
  decil9 = decil9 + (9*n_index/10 - acum90)/(hist(class_d9)*n_index)

  write(2, *) '&'
  write(2, *) quart1, hist(class_q1)
  write(2, *) '&'
  write(2, *) quart3, hist(class_q3)
  write(2, *) '&'
  write(2, *) decil1, hist(class_d1)
  write(2, *) '&'
  write(2, *) decil9, hist(class_d9)


  write(*, '(a13, f10.3, a22)')"* Mode     = ", moda, "                   *"
  write(*, '(a13, f10.3, a22)')"* Median   = ", mediana, "                   *"
  write(*, '(a45)') "*                                           *"
  write(*, '(a20, f10.3, a15)')"* Lower Quartile  = ", quart1, "          *" 
  write(*, '(a20, f10.3, a15)')"* Upper Quartile  = ", quart3, "          *"  
  write(*, '(a20, f10.3, a15)')"* Lower Decile    = ", decil1, "          *"
  write(*, '(a20, f10.3, a15)')"* Upper Decile    = ", decil9, "          *"
  write(*, '(a45)') "*                                           *"
  write(*, '(a30, f10.3, a5)')"* Coefficient of Variation  = ", desv/aver, "    *"
  write(*, '(a45)')"*===================================        *"
  
  if (desv*desv/aver<1) write(*, '(a45)') "* Under-Dispersed Distribution              *"
  if (desv*desv/aver>1) write(*, '(a45)') "* Over-Dispersed Distribution               *"

  write(*, '(a45)') "*===================================        *"
  write(*, '(a45)') "*                                           *"
  write(*, '(a14, f10.3, a21)') "* Skewness  = ", skew, "                    *" 
  write(*, '(a45)') "*===================================        *"

  if (skew<0) write(*, '(a45)') "* Negative Skew                             *"
  if (skew==0) write(*, '(a45)') "* Symmetrical Distribution                  *"
  if (skew>0) write(*, '(a45)') "* Positive Skew                             *"
  
  write(*, '(a45)') "*===================================        *"
  write(*, '(a45)') "*                                           *"
  write(*, '(a21, f10.3, a14)')"* Excess Kurtosis  = ", kurt-3, "             *"
  write(*, '(a45)') "*===================================        *"

  if (kurt-3<0) write(*, '(a45)') "* Platykurtic Distribution                  *"
  if (kurt-3==0) write(*, '(a45)') "* Mesokurtic Distribution                   *"
  if (kurt-3>0) write(*, '(a45)') "* Leptokurtic Distribution                  *"

  write(*, '(a45)') "*===================================        *"
  write(*, '(a45)') "*                                           *"
  write(*, '(a45)') "****************** SUMMARY ******************"

  !==============================================
  !==============================================

  ! Calculating autocorrelation function ========

  write(*, *)
  write(*, *) "Calculating Autocorrelation Function ......"

  do j=1,  int(n_index/2)
  
     acf = 0
     do i=1, int(n_index/2)

        acf = acf + (func(i)-aver)*(func(i+j)-aver)/(int(n_index/2)*desv*desv)

     end do

     write(3, *) j, acf

  end do

  !==============================================

  close(1)
  close(2)

  write(*, *) 
  write(*, *) "Finished"
  write(*, *)
  if (back) write(*, '(a31)') ' Previous files were backed up!'
  if (back) write(*, *)
    
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
