program integralMCF
    use omp_lib
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer(dp) :: i, j, k, total_samples, batches, histories
    real(dp) :: x, a, b, IMC, mean, hist_real
    real(dp), allocatable :: IMC_val(:), r_val(:), variance(:), stdv(:)
    integer(dp), allocatable :: history_count(:)
    call random_seed()

    a = 0.0_dp ! lower range of integration
    b = 5.0_dp ! upper range of integration
    batches = 50000
    histories = 1000

    allocate(IMC_val(batches), history_count(batches), variance(batches), stdv(batches))

    do j = 1, batches
        ! BEGIN CALCULATION ***********************************************************************************
        allocate(r_val(histories))
        do i = 1, histories
            call random_number(x)
            x = a + x * (b - a) ! scales x to be within a range from b to a
            r_val(i) = f(x) ! appends rand f(x) val to array
        end do

        IMC = (b - a) * (sum(r_val) / histories) ! CALCULATES INTEGRAL
        ! END CALCULATION *************************************************************************************
        IMC_val(j) = IMC

        mean = sum(IMC_val(1:j)) / j  ! Use all calculated IMC values up to the current batch
        hist_real = real(j)
        variance(j) = sum((IMC_val(1:j) - mean) ** 2) / hist_real  ! Calculate variance
        stdv(j) = sqrt(variance(j))

        history_count(j) = histories

        print *, 'Batch:', j, 'Histories:', histories, 'IMC:', IMC, 'Variance:', variance(j)

        histories = histories + 1.001 ! the growth factor of histories between batches
        deallocate(r_val)
    end do
    

    ! Write arrays to CSV
    open(unit=1, file='results.csv', status='replace')
    write(1,*) 'batch,IMC_val,history_count,variance,stdv'
    do j = 1, batches
        write(1, '(I0,",",ES15.7,",",I0,",",ES15.7,",",ES15.7)') &
              j, IMC_val(j), history_count(j), variance(j), stdv(j)  ! Corrected references
    end do
    close(1)
    print *, 'CSV file written successfully!'


    deallocate(IMC_val, history_count, variance, stdv)
contains

    function f(x) result(y)
        real(dp) :: x, y
        y = sin(x**2)
        ! 0.527917
    end function f

end program integralMCF
