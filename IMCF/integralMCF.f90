program integralMCF
    use omp_lib
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer :: i, j, k, batches, histories, sum_hist
    real(dp) :: x, a, b, IMC, mean, hist_real, lit_val, start_time, end_time, sum_fx
    real(dp), allocatable :: IMC_val(:), variance(:), stdv(:), rel_err(:), batch_time(:)
    integer(dp), allocatable :: history_count(:)
    call random_seed()

    a = 0.0_dp ! lower range of integration
    b = 5.0_dp ! upper range of integration
    batches = 800 ! how many integrations to compute
    histories = 1000
    lit_val = 0.52791728116532241384461568_dp
! HEY YOU SHOULD ADD A CONFIDENCE INTERVAL FEATURE!!! *&()*@()#*()@&)!(*#$&#*)(!@#$*&()#@*&)($*&()^*)B^&(#)B&*(^@#)*($^*@#)$^@)#@#$*)(@#(*$*@#()$))
    allocate(IMC_val(batches), history_count(batches), variance(batches), stdv(batches), rel_err(batches), batch_time(batches))

    ! Creates an array of histories that each batch takes its respective index of
    do k = 1, batches
        history_count(k) = histories
        histories = histories * 1.01
    end do

    sum_hist = sum(history_count)

    !$OMP PARALLEL DO PRIVATE(j, i, histories, x, sum_fx, IMC, start_time, end_time) &
    !$OMP& SHARED(a, b, IMC_val, batches, batch_time, lit_val)
    do j = 1, batches
        start_time = omp_get_wtime()
        ! BEGIN CALCULATION ***********************************************************************************
        histories = history_count(j)
        sum_fx = 0.0_dp

        do i = 1, histories
            call random_number(x)
            x = a + x * (b - a) ! scales x to be within a range from b to a
            sum_fx = sum_fx + f(x) ! adds rand f(x) value to sum
        end do

        IMC = (b - a) * (sum_fx / histories) ! CALCULATES INTEGRAL
        ! END CALCULATION *************************************************************************************
        end_time = omp_get_wtime()
        batch_time(j) = end_time - start_time

        IMC_val(j) = IMC

        sum_hist = sum_hist - histories
        print *, 'Remaining Histories:',sum_hist
    end do
    !$OMP END PARALLEL DO

    do j = 1, batches
        mean = sum(IMC_val(1:j)) / j  ! Use all calculated IMC values up to the current batch
        hist_real = real(history_count(j))
        variance(j) = sum((IMC_val(1:j) - mean) ** 2) / (hist_real - 1)  ! Calculate variance
        stdv(j) = sqrt(variance(j))
        rel_err(j) = abs((IMC_val(j) - lit_val)) / lit_val
    end do

    ! Write arrays to CSV
    open(unit=1, file='results.csv', status='replace')
    write(1,*)'batch,IMC_val,history_count,variance,stdv,rel_err,batch_time'
    do j = 1, batches
        write(1, '(I0,",",ES15.7,",",I0,",",ES15.7,",",ES15.7,",",ES15.7,",",ES15.7)') &
              j, IMC_val(j), history_count(j), variance(j), stdv(j), rel_err(j), batch_time(j)
    end do
    close(1)
    print *, 'CSV file written successfully!'


    deallocate(IMC_val, history_count, variance, stdv, rel_err)
contains

    function f(x) result(y)
        real(dp) :: x, y
        y = x
    end function f

end program integralMCF
