from DIMCE_GPU import IntegralMC

integral_calc = IntegralMC(
    batches=200,  # Number of times that sets of integral should be estimated ( y in calc_int array )
    int_per_batch=2,  # How many estimates of integral should be output per batch ( x in calc_int array )
    chunk_size=10 ** 7,  # for VRAM management, chunks history count into "manageable bites" for your VRAM. decrease for less VRAM (current system uses RTX 4070) -- untested at higher values
    a=0,  # Lower bound of the integral
    b=5,  # Upper bound of the integral
    histories=10000000,  # number of histories to take per permutation
    hist_factor=100,  # factor in which histories are multiplied in order to create a gradient for data analysis,
    repetition_factor=20,
    f=lambda x: x ** 2  # the function to integrate
)

print("Starting integral calculation...")
integral_calc.calc()


from intplot import plot_arrays

plot_arrays.rungraphs()
