program monte_carlo_integration
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer :: points_per_batch = 10000
    real(dp) :: x, sum_curve, final_integral_estimate
    integer :: i, j, k
    real(dp) :: a, b
    integer :: total_checks, batches, bat_countdown, iter_pb
    real(dp), allocatable :: batch_results(:)
    real(dp), allocatable :: calc_int(:)
   
    call random_seed()
    a = 0.0_dp
    b = 6.0_dp
    batches = 3
    iter_pb = 5

    allocate(calc_int(batches))
    allocate(batch_results(iter_pb))
    bat_countdown = batches
    do j = 1, batches
        do k = 1, iter_pb
        sum_curve = 0.0_dp
        total_checks = 0
            do i = 1, points_per_batch
                call random_number(x)
                x = a + x * (b - a)
                sum_curve = sum_curve + f(x)
                total_checks = total_checks + 1
            end do
        batch_results(k) = (b - a) * (sum_curve / total_checks)
        end do
        ! print *, 'Time remaining:', bat_countdown, '| Estimates:', points_per_batch
        bat_countdown = bat_countdown - 1
        points_per_batch = points_per_batch + 100
        print *, batch_results
        calc_int = sum(batch_results) / iter_pb
    end do
    print *, 'Estimated integral value:', calc_int
    

    deallocate(batch_results)
    deallocate(calc_int)

contains
    function f(x) result(y)
        real(dp) :: x, y
        y = x ** 2
    end function f
end program monte_carlo_integration