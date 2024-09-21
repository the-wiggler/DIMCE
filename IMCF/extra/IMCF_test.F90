program integralMCF
    use omp_lib
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer(dp) :: histories = 10000000
    real(dp) :: x, sum_curve, mean, variance, stddev
    real(dp) :: a, b, start_time, end_time
    integer(dp) :: i, j, k, l, n, total_checks, batches, bat_countdown, iter_pb, start, finish, rate, rep_fac
    real(dp), allocatable :: batch_results(:), calc_int(:), calc_stddev(:), batch_times(:), time_list(:), compute_time(:)
    integer(dp), allocatable :: history_count(:), thread_num(:)

    call random_seed()
    call system_clock(count_rate=rate)

    a = 0.0_dp ! lower range of integration
    b = 1.0_dp ! upper range of integration
    batches = 1 ! how many times to perform an integral estimation
    iter_pb = 1 ! how many iterations should be performed in each batch (to be averaged together)
    rep_fac = 25

    allocate(compute_time(rep_fac))
    do l = 1, rep_fac
        allocate(calc_int(batches), calc_stddev(batches), history_count(batches), batch_times(batches), batch_results(iter_pb) &
            , thread_num(batches))
        bat_countdown = batches
        call system_clock(start)
        !$OMP PARALLEL DO PRIVATE(j, k, i, sum_curve, total_checks, x) &
        !$OMP& PRIVATE(variance, stddev, mean) &
        !$OMP& SHARED(a, b, histories, calc_int, calc_stddev, history_count, batch_times)
        do j = 1, batches ! a loop that recursively creates new estimations with an increasing history rate for data analysis
        thread_num(j) = OMP_GET_THREAD_NUM()
            start_time = omp_get_wtime() ! batch time start
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

            ! calc_int(j) = mean
            end_time = omp_get_wtime() ! batch time end
            ! batch_times(j) = end_time - start_time

            variance = sum((batch_results - mean) ** 2) / (iter_pb - 1) ! calculates variance of batch_results set over {iter_pb} trials
            stddev = sqrt(variance) ! calculates standard deviation 
            ! calc_stddev(j) = stddev

            print *, 'rep_num', batches
            bat_countdown = bat_countdown - 1
            ! history_count(j) = histories
        end do
        !$OMP END PARALLEL DO
        call system_clock(finish)
        compute_time(l) = real(finish - start) / real(rate)
        deallocate(calc_int, calc_stddev, history_count, batch_times, batch_results, thread_num)
        batches = batches + 1
    end do

    print *, compute_time
    open(unit=1, file='varbatchSERIAL.csv', status='replace')
    write(1,*) 'compute_time'
    do n = 1, rep_fac
        write(1,'(ES15.7)') compute_time(n)
    end do


contains
    function f(x) result(y)
        real(dp) :: x, y
        y = log(1 + x) / x
    end function f
end program integralMCF
