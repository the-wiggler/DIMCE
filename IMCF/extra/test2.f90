program integralMCF
    use omp_lib
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer(dp) :: histories = 10
    real(dp) :: x, sum_curve, mean, variance, stddev
    real(dp) :: a, b, start_time, end_time
    integer(dp) :: i, j, k, total_checks, batches, bat_countdown, iter_pb, start, finish, rate
    real(dp), allocatable :: batch_results(:), calc_int(:), calc_stddev(:), batch_times(:)
    integer(dp), allocatable :: history_count(:), thread_num(:)

    call random_seed()
    call system_clock(count_rate=rate)

    a = 0.0_dp ! lower range of integration
    b = 1.0_dp ! upper range of integration
    batches = 180 ! how many times to perform an integral estimation
    iter_pb = 3 ! how many iterations should be performed in each batch (to be averaged together)

    allocate(calc_int(batches), calc_stddev(batches), history_count(batches), batch_times(batches), batch_results(iter_pb) &
        , thread_num(batches))
    bat_countdown = batches

    call system_clock(start)
    !$OMP PARALLEL DO PRIVATE(j, k, i, sum_curve, total_checks, x) &
    !$OMP& PRIVATE(variance, stddev, mean) &
    !$OMP& SHARED(a, b, histories, calc_int, calc_stddev, history_count, batch_times)
    do j = 1, batches ! a loop that recursively creates new estimations with an increasing history rate for data analysis
        do k = 1, iter_pb
            sum_curve = 0.0_dp ! the sum of output values
            total_checks = 0 ! the total points taken in this iteration
            do i = 1, histories
                call random_number(x)
                x = a + x * (b - a) ! scales x to be within a range from b to a
                sum_curve = sum_curve + f(x) ! calculates an output value based on a random x and adds it to the total
                total_checks = total_checks + 1
            end do
            batch_results(k) = (b - a) * (sum_curve / total_checks) ! calculates the integral by averaging the range of y values and calculating the area based on width x avg height
        end do

        mean = sum(batch_results) / iter_pb ! Calculate mean of the batch

        calc_int(j) = mean

        variance = sum((batch_results - mean) ** 2) / (iter_pb - 1) ! calculates variance of batch_results set over {iter_pb} trials
        stddev = sqrt(variance) ! calculates standard deviation 
        calc_stddev(j) = stddev

        print *, 'Time remaining:', bat_countdown, '| Estimates:', histories
        bat_countdown = bat_countdown - 1
        history_count(j) = histories
        histories = histories * 1.1
    end do
    !$OMP END PARALLEL DO
    call system_clock(finish)

    print *, 'Estimated integral values and standard deviations:'
    do j = 1, batches
        print *, 'Relative Batch ', j, ': Integral =', calc_int(j), 'Std Dev =', calc_stddev(j), &
                 'Histories =', history_count(j), 'Time =', batch_times(j), 'Thread =', thread_num(j)
    end do

    ! write arrays
    open(unit=1, file='resultsOMP.csv', status='replace')
    write(1,*) 'rel_batch,history,calc_int,stddev,batch_time,thread_num'
    do j = 1, batches
        write(1,'(I0,",",I0,",",ES15.7,",",ES15.7,",",ES15.7,",",I0)') &
              j, history_count(j), calc_int(j), calc_stddev(j), batch_times(j), thread_num(j)
    end do
    close(1)
    print *, 'CSV file written successfully!'


    deallocate(calc_int, calc_stddev, history_count, batch_times)

    print *, "Runtime: ", real(finish - start) / real(rate)
contains
    function f(x) result(y)
        real(dp) :: x, y
        y = log(1 + x) / x
    end function f
end program integralMCF
