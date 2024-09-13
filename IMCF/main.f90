program integralMCF
    use omp_lib
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer :: histories = 50000
    real(dp) :: x, sum_curve, final_integral_estimate, mean, variance, stddev
    real(dp) :: a, b
    integer :: i, j, k, total_checks, batches, bat_countdown, iter_pb
    integer :: start, finish, rate
    real(dp), allocatable :: batch_results(:), calc_int(:), calc_stddev(:), batch_times(:)
    integer, allocatable :: history_count(:)
    call random_seed()
    call system_clock(count_rate=rate)

    a = 0.0_dp ! lower range of integration
    b = 6.0_dp ! upper range of integration
    batches = 100 ! how many times to perform an integral estimation
    iter_pb = 5 ! how many iterations should be performed in each batch (to be averaged together)

    allocate(calc_int(batches))
    allocate(calc_stddev(batches))
    allocate(history_count(batches))
    allocate(batch_times(batches))
    allocate(batch_results(iter_pb))
    bat_countdown = batches

    call omp_set_num_threads(16)
    !$OMP PARALLEL DO PRIVATE(j, k, i, sum_curve, total_checks, x, variance, stddev, mean) SHARED(a, b, histories, calc_int, calc_stddev, history_count, batch_times)
    do j = 1, batches ! a loop that recursively creates new estimations with an increasing history rate for data analysis
        call system_clock(start) ! batch time start
        ! allocate(batch_results(iter_pb))
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
        call system_clock(finish) ! batch time end
        batch_times(j) = real(finish - start) / real(rate)

        variance = sum((batch_results - mean) ** 2) / (iter_pb - 1) ! calculates variance of batch_results set over {iter_pb} trials
        stddev = sqrt(variance) ! calculates standard deviation 
        calc_stddev(j) = stddev

        print *, 'Time remaining:', bat_countdown, '| Estimates:', histories
        bat_countdown = bat_countdown - 1
        history_count(j) = histories
        histories = histories + 10000
        ! deallocate(batch_results)
    end do
    !$OMP END PARALLEL DO
    

    print *, 'Estimated integral value:'
    do j = 1, batches
        print *, 'Batch ', j, ': ', calc_int(j)
    end do
    print *, 'Standard deviation for each batch:'
    do j = 1, batches
        print *, 'Batch ', j, ': ', calc_stddev(j)
    end do
    print *, 'History Values:'
    do j = 1, batches
        print *, 'Batch ', j, ': ', history_count(j)
    end do
        print *, 'Batch Times:'
    do j = 1, batches
        print *, 'Batch ', j, ': ', batch_times(j)
    end do

    ! write arrays
    open(unit=1, file='IMCF_out/calc_int.bin', form='unformatted', status='replace')
    write(1) calc_int
    close(1)

    open(unit=2, file='IMCF_out/calc_stddev.bin', form='unformatted', status='replace')
    write(2) calc_stddev
    close(2)

    open(unit=3, file='IMCF_out/history_count.bin', form='unformatted', status='replace')
    write(3) history_count
    close(3)

    open(unit=4, file='IMCF_out/batch_times.bin', form='unformatted', status='replace')
    write(4) batch_times
    close(4)

    print *, 'files written to dir successfully!'

    deallocate(calc_int)
    deallocate(calc_stddev)
    deallocate(history_count)
    deallocate(batch_results)

contains
    function f(x) result(y)
        real(dp) :: x, y
        y = x **2 + 2 * x + 2
    end function f
end program integralMCF
