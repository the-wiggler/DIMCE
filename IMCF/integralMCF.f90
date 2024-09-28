program integralMCF
    use omp_lib
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer(dp) :: i, j, k, batches, histories, hist_sum
    real(dp) :: x, a, b, IMC, mean, hist_real, lit_val
    real(dp), allocatable :: IMC_val(:), r_val(:), variance(:), stdv(:), rel_err(:)
    integer(dp), allocatable :: history_count(:)
    call random_seed()

    a = 0.0_dp ! lower range of integration
    b = 5.0_dp ! upper range of integration
    batches = 800 ! how many integrations to compute
    histories = 1000
    lit_val = 0.52791728116532241384461568_dp

    allocate(IMC_val(batches), history_count(batches), variance(batches), stdv(batches), rel_err(batches))

    do k = 1, batches
        history_count(k) = histories
        histories = histories * 1.01
    end do

    hist_sum = sum(history_count)

    do j = 1, batches
        ! BEGIN CALCULATION ***********************************************************************************
        histories = history_count(j)
        allocate(r_val(histories))
        !$OMP PARALLEL DO
        do i = 1, histories
            call random_number(x)
            x = a + x * (b - a) ! scales x to be within a range from b to a
            r_val(i) = f(x) ! appends rand f(x) val to array
        end do
        !$OMP END PARALLEL DO

        IMC = (b - a) * (sum(r_val) / histories) ! CALCULATES INTEGRAL
        ! END CALCULATION *************************************************************************************
        IMC_val(j) = IMC
        deallocate(r_val)
        write(*, '(A, I3, A, I0, A, I0)') 'Batch:', j, ' | Histories:', history_count(j), ' | Remaining Histories:', hist_sum
        hist_sum = hist_sum - histories
    end do

    do j = 1, batches
        mean = sum(IMC_val(1:j)) / j  ! Use all calculated IMC values up to the current batch
        hist_real = real(j)
        variance(j) = sum((IMC_val(1:j) - mean) ** 2) / hist_real  ! Calculate variance
        stdv(j) = sqrt(variance(j))
        rel_err(j) = abs((IMC - lit_val)) / lit_val
    end do

    ! Write arrays to CSV
    open(unit=1, file='results.csv', status='replace')
    write(1,*) 'batch,IMC_val,history_count,variance,stdv,rel_err'
    do j = 1, batches
        write(1, '(I0,",",ES15.7,",",I0,",",ES15.7,",",ES15.7,",",ES15.7)') &
              j, IMC_val(j), history_count(j), variance(j), stdv(j), rel_err(j)
    end do
    close(1)
    print *, 'CSV file written successfully!'


    deallocate(IMC_val, history_count, variance, stdv, rel_err)
contains

    function f(x) result(y)
        real(dp) :: x, y
        y = sin(x**2)
    end function f

end program integralMCF
