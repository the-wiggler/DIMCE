program integralMCF
    use omp_lib
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer(dp) :: i, j, k, total_samples, batches, histories
    real(dp) :: x, a, b, IMC, mean
    real(dp), allocatable :: IMC_val(:), r_val(:)
    call random_seed()

    a = 0.0_dp ! lower range of integration
    b = 5.0_dp ! upper range of integration
    batches = 500
    histories = 100

    allocate(IMC_val(batches))

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

        mean = sum(IMC_val) / histories
        print *, mean

        print *, 'Batch:', j, 'Histories:', histories, 'IMC:', IMC, 'Variance:'

        histories = histories * 1.1 ! the growth factor of histories between batches
        deallocate(r_val)
    end do
    deallocate(IMC_val)
contains

    function f(x) result(y)
        real(dp) :: x, y
        y = sin(x**2)
        ! 0.527917
    end function f

end program integralMCF
