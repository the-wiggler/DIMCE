program test
    use openacc
    use omp_lib
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer :: histories = 10000
    real(dp) :: x, sum_curve, mean, variance, stddev
    real(dp) :: a, b, start_time, end_time
    integer :: i, j, k, total_checks, batches, bat_countdown, iter_pb, x_track
    real(dp), allocatable :: batch_results(:), calc_int(:), calc_stddev(:), batch_times(:), random_x(:)
    integer, allocatable :: history_count(:)

    call random_seed()

    a = 0.0_dp ! lower range of integration
    b = 6.0_dp ! upper range of integration
    batches = 1000 ! how many times to perform an integral estimation
    iter_pb = 3 ! how many iterations should be performed in each batch (to be averaged together)

    allocate(calc_int(batches), calc_stddev(batches), history_count(batches), batch_times(batches), batch_results(iter_pb))

    allocate(random_x(batches*iter_pb*histories))
    !$OMP PARALLEL DO
    do j = 1, batches * iter_pb * histories
        call random_number(x)
        x = a + x * (b - a)
        random_x(j) = x
    end do
    !$OMP END PARALLEL DO
    print *, random_x
end program test