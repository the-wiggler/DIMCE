program integralMCF
    use openacc
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer :: histories = 1000000
    real(dp) :: x, sum_curve, mean, variance, stddev
    real(dp) :: a, b, start_time, end_time
    integer(dp) :: i, j, k, total_checks, batches, bat_countdown, iter_pb, x_track
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

    bat_countdown = batches
    x_track = batches * iter_pb * histories
    !$acc kernels
    do j = 1, batches ! a loop that recursively creates new estimations with an increasing history rate for data analysis
        ! start_time = omp_get_wtime() ! batch time start
        do k = 1, iter_pb
            sum_curve = 0.0_dp ! the sum of output values
            total_checks = 0 ! the total points taken in this iteration
            do i = 1, histories
                sum_curve = sum_curve + f(random_x(x_track)) ! calculates an output value based on a random x and adds it to the total
                total_checks = total_checks + 1
            end do
            batch_results(k) = (b - a) * (sum_curve / total_checks) ! calculates the integral by averaging the range of y values and calculating the area based on width x avg height
            x_track = x_track - 1
        end do

        mean = sum(batch_results) / iter_pb ! Calculate mean of the batch
        calc_int(j) = mean
        ! end_time = omp_get_wtime() ! batch time end
        ! batch_times(j) = end_time - start_time

        variance = sum((batch_results - mean) ** 2) / (iter_pb - 1) ! calculates variance of batch_results set over {iter_pb} trials
        stddev = sqrt(variance) ! calculates standard deviation
        calc_stddev(j) = stddev

        print *, 'Time remaining:', bat_countdown, '| Estimates:', histories
        bat_countdown = bat_countdown - 1
        history_count(j) = histories
        histories = histories + 100
    end do
    !$acc end kernels
    

    print *, 'Estimated integral values and standard deviations:'
    do j = 1, batches
        print *, 'Batch ', j, ': Integral =', calc_int(j), 'Std Dev =', calc_stddev(j), &
                 'Histories =', history_count(j)!, 'Time =', batch_times(j)
    end do


    ! write arrays
    open(unit=1, file='results.csv', status='replace')
    write(1,*) 'batch,history,calc_int,stddev,batch_time'
    do j = 1, batches
        write(1,'(I0,",",I0,",",ES15.7,",",ES15.7,",",ES15.7)') &
              j, history_count(j), calc_int(j), calc_stddev(j), batch_times(j)
    end do
    close(1)
    print *, 'CSV file written successfully!'


    deallocate(calc_int, calc_stddev, history_count, batch_times)

contains
    function f(x) result(y)
        real(dp) :: x, y
        y = x **2 + 2 * x + 2
    end function f
end program integralMCF
