program monte_carlo_pi
  implicit none
  integer :: num_samples, i, inside_circle
  real :: x, y, pi_estimate

  integer :: seed
  call random_seed()

  num_samples = 10000000
  inside_circle = 0

  do i = 1, num_samples
      call random_number(x)
      call random_number(y)
      if (x**2 + y**2 <= 1.0) then
          inside_circle = inside_circle + 1
      end if
  end do

  pi_estimate = 4.0 * inside_circle / num_samples
  print *, "Estimated Pi:", pi_estimate
end program monte_carlo_pi
