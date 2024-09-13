program monte_carlo_integration
    implicit none

    integer(8), parameter :: points_per_batch = 1000000000
    real(8) :: x, integral_estimate, sum_curve
    integer(8) :: i, total_checks
    real(8) :: a, b

    call random_seed()

    a = 0.0_8
    b = 6.0_8

    sum_curve = 0.0_8
    total_checks = 0

    do i = 1, points_per_batch
        call random_number(x)
        x = a + x * (b - a)
        sum_curve = sum_curve + f(x)
        total_checks = total_checks + 1
    end do

    integral_estimate = (b - a) * (sum_curve / total_checks)

    print *, "Estimated integral value:", integral_estimate

contains
    function f(x) result(y)
        real(8) :: x, y
        y = x ** 2
    end function f

end program monte_carlo_integration
