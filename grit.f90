
program grit
    use utils
    use iso_fortran_env, only: output_unit
    !$ use omp_lib
    implicit none

    ! ======================================== !
    ! Preparing Things for the Parallelization !
    ! ======================================== !

    call init_mpi()
    call set_structure_parameters()

    if (mpi_rank == 0) then
        call print_intro()

        call print_with_border("Information about Input-Parameters")

        call print_with_border2("Modes")
        write(*,'(A,L1)') " - Fine k-grid near Fermi surface    = ", Fermi_centered_k
        write(*,'(A,L1)') " - Fine q-grid near the Gamma-point  = ", Gamma_centered_q
        write(*,'(A,L1)') " - Interacting Electrons and phonons = ", interacting
        write(*,'(A,L1)') " - Debug mode enabled                = ", debug_mode
        write(*,*)

        call print_with_border2("Tasks")
        write(*,'(A,L1)') " - Create binary files                    = ", create_binary
        write(*,'(A,L1)') " - Calc N_states_per_UC (mu = E_F, T = 0) = ", calc_N_states_per_UC_mu_eq_E_F_T_eq_zero
        write(*,'(A,L1)') " - Calc N_per_UC (mu = E_F, T = 0)        = ", calc_N_per_UC_mu_eq_E_F_T_eq_zero
        write(*,'(A,L1)') " - Calc mu vs Temperature                 = ", calc_mu_vs_T
        write(*,'(A,L1)') " - Calc specific heat per UC              = ", calc_c_V_per_UC
        write(*,'(A,L1)') " - Calc Sommerfeld gamma from linear fit  = ", calc_gamma
        write(*,*)

        call print_with_border2("Temperature parameters")
        write(*,'(A,I0)') " - Number of temperature points = ", n_temp
        write(*,'(A,F12.4)') " - Minimum temperature (K)   = ", T_min
        write(*,'(A,F12.4)') " - Maximum temperature (K)   = ", T_max
        write(*,*)

        call print_with_border2("Threshold parameters")
        write(*,'(A,ES12.4)') " - Bisection threshold   = ", threshold_bisection
        write(*,'(A,ES12.4)') " - Omega threshold       = ", threshold_omegas
        write(*,'(A,ES12.4)') " - Coupling threshold    = ", threshold_couplings
        write(*,'(A,ES12.4)') " - Denominator threshold = ", threshold_denominator
        write(*,'(A,ES12.4)') " - Xi threshold          = ", threshold_xis
        write(*,*)

        call print_with_border2("GRIT convergence parameters")
        write(*,'(A,F12.4)') " - fsthick                           = ", fsthick
        write(*,'(A,ES12.4)') " - Eta (electronic broadening)       = ", eta
        write(*,'(A,ES12.4)') " - eta_q (phononic broadening)       = ", eta_q
        write(*,'(A,ES12.4)') " - Energy tails (DO(O)S integration) = ", E_tails
        write(*,'(A,ES12.4)') " - Energy radius (SH integration)    = ", E_Length_SH
        write(*,'(A,ES12.4)') " - Delta mu search window            = ", delta_mu
        write(*,'(A,I0)')     " - Number of energy points           = ", n_energy
        write(*,*)
        write(*,*)


        call print_with_border("Information about Parallelizations")

        write(*,'(A,I0,A)') " - MPI initialized with ", mpi_size, " processes"
        !$ write(*,'(A,I0,A)') " - OpenMP is enabled with ", omp_get_max_threads(), " threads available"
        
        write(*,*)
        write(*,*)
    end if 



    ! ================= !
    ! Setting the paths !
    ! ================= !

    call set_filepaths()



    ! ==================== !
    ! Task = create_binary !
    ! ==================== !

    ! Only one process, reads data from textfile and creates binary file, only in the case of interactions. 
    if (mpi_rank == 0 .and. interacting .and. create_binary) then

        call print_with_border("Creating binary-versions of the big Self-Energy-array")

        write(*,*) "- We have to read the Self-Energy k-grid and the Gamma-centered q-grid data to know the size of the arrays"
        write(*,*) 

        call read_SE_epsilons_data()
        call read_omegas_data("Self-Energy")
        call read_n_write_g_data("Self-Energy")
        call read_n_write_SE_epsilon_kplusq_data()

    end if

    ! Only one process, reads data from textfile and creates binary file, only in the case of interactions. 
    if (mpi_rank == 0 .and. interacting .and. create_binary .and. calc_c_V_per_UC) then

        call print_with_border("Creating binary-versions of the uniform coupling-arrays")

        call read_k_data()
        call read_n_write_g_data("Specific-Heat")

    end if

    ! Make sure that all processes wait until binary file is created
    call MPI_Barrier(MPI_COMM_WORLD, mpi_ierr)

    flush(output_unit)



    ! =========================================== !
    ! Loading the Data for subsequent other tasks !
    ! =========================================== !

    if (mpi_rank == 0) then
        call print_with_border("Loading the Data for subsequent other tasks")
    end if


    ! ------------------------ !
    ! Load Electronic Energies !
    ! ------------------------ !

    call read_epsilons_data("Non-Interacting")
    call read_omegas_data("Non-Interacting")
    
    if (interacting .or. Fermi_centered_k) call read_SE_epsilons_data()

    if (interacting) then 
        call read_omegas_data("Self-Energy")
        call read_coupling_array("Self-Energy" )
    end if 

    if (interacting .and. calc_c_V_per_UC) then 
        ! We need routines here: 
        call read_k_data()
        call read_epsilons_data("Uniform")
        call read_omegas_data("Uniform")
        call read_coupling_array("Uniform")
    end if 

    flush(output_unit)



    ! ====================================== !
    ! Print Summary of Structure-Information !
    ! ====================================== !

    if (mpi_rank == 0) then
        write(*,*)
        call print_with_border("Structure-Information")
        write(*,*) "- Structure = ", trim(structure_name)
        write(*,*)
        write(*,'(A,I0,A)') " - Number of branches           =  n_branches = ", n_branches
        write(*,'(A,I0,A)') " - Number of acoustic branches  =          dx = ", dx
        write(*,'(A,I0,A)') " - Number of optical branches   = r * dx - dx = ", r * dx - dx
        write(*,'(A,I0,A)') " - Number of relevant W90-bands =       N_bnd = ", N_bnd
        write(*,'(A,I0,A)') " - Number of states per UC      =   2 * N_bnd = ", (N_bnd * 2)
        write(*,'(A,I0,A)') " - Number of electrons per UC   =        n_UC = ", n_UC
        write(*,*)
        write(*,*)


        call print_with_border("EPW-Parameters")
        
        write(*,'(A,I0,A)') " - Number of k-points in each direction = Nk = ", Nk
        write(*,'(A,I0,A)') " - Number of q-points in each direction for Non-Interacting-contribution = Nq = ", Nq
        
        if (interacting) then 
            write(*,'(A,I0,A)') " - Number of q-points in each direction for SE-calculation = Nq_SE = ", Nq_SE
        end if

        if (calc_c_V_per_UC) then 
            write(*,'(A,I0,A)') " - Number of q-points in each direction for SH-calculation = Ni_uniform = ", Ni_uniform
        end if 
        
        write(*,*)
        write(*,*)
    end if 

    ! Find min and max energies (all processes need these)
    E_min = minval(epsilon_arr)
    E_max = maxval(epsilon_arr)

    E_min_rel = E_min - E_F
    E_max_rel = E_max - E_F
    
    if (mpi_rank == 0) then
        call print_with_border2("Information about Energies")
        write(*,'(A,F8.3,A)') " - The Fermi energy is (from nscf-calculation in QE)    E_F = ", E_F, " eV."
        write(*,'(A,F8.3,A,F8.3,A)') " - Extremal energies of the electronic system: E_min = ", E_min, " eV, E_max = ", E_max, " eV"
        write(*,*)
        write(*,'(A)') " - Boundaries for the Integration-Interval:"
        write(*,'(A,F8.3,A)') " - E_min - E_F = E_min_rel = ", E_min_rel, " eV."
        write(*,'(A,F8.3,A)') " - E_max - E_F = E_max_rel = ", E_max_rel, " eV."
        write(*,*)
        write(*,*)
    end if

    flush(output_unit)



    ! ======================================== !
    ! Task = calc_N_states_mu_eq_E_F_T_eq_zero !
    ! ======================================== !

    if (calc_N_states_per_UC_mu_eq_E_F_T_eq_zero) then
        
        if (mpi_rank == 0) then 
            call print_with_border("Task = calc_N_states_per_UC_mu_eq_E_F_T_eq_zero")
            write(*,'(A)') " - Computing the number of states per UC at T = 0 K for mu = E_F."
        end if 

        ! Time the states calculation
        call system_clock(start_time, clock_rate)

        call calculate_no_of_States_per_UC(0.0_dp, E_F)

        call system_clock(end_time)
        elapsed_time = real(end_time - start_time) / real(clock_rate)

        if (mpi_rank == 0) then
            ! Assessing the normalization of our Setup for mu = E_F
            write(*,*)
            write(*,'(A,F12.3,A)') &
                 " - The computed number of states per UC at T = 0 K for mu = E_F is " // &
                 "<N_states_per_UC> = ", N_states_per_UC

            write(*,'(A,F12.3,A)') &
                 " - The expected number of states per UC is " // &
                 "                                N_bnd * 2 = ", real(N_bnd * 2)

            write(*,'(A,F12.3,A)') &
                 " - The ratio is " // &
                 "                       <N_states_per_UC> / (N_bnd * 2) = ", &
                                 N_states_per_UC / real(N_bnd * 2)

            call print_elapsed_time("States calculation", elapsed_time)
        end if

        flush(output_unit)
    
    end if


    ! ======================================== !
    ! Task = calc_N_per_UC_mu_eq_E_F_T_eq_zero !
    ! ======================================== !

    if (calc_N_per_UC_mu_eq_E_F_T_eq_zero) then 

        if (mpi_rank == 0) then 
            call print_with_border("Task = calc_N_per_UC_mu_eq_E_F_T_eq_zero")
            write(*,'(A)') " - Computing the number of electrons per UC at T = 0 K for mu = E_F."
        end if 

        N_per_UC_mu_eq_E_F_T_eq_zero = calculate_particle_number_per_UC(0.0_dp, E_F, .true.) 

        if (mpi_rank == 0) then
            ! Printing the particle number at T = 0 K for mu = E_F
            write(*,*)

            write(*,'(A,F8.5,A)') " - The computed particle-number per UC is n_UC = ", &
    N_per_UC_mu_eq_E_F_T_eq_zero, " (electrons)."
            write(*,'(A,I0,A)') " - The expected particle-number per UC is n_UC = ", n_UC
            write(*,*)
            write(*,*)
            
        end if

    end if

    flush(output_unit)



    ! =================== !
    ! Task = calc_mu_vs_T !
    ! =================== !

    if (calc_mu_vs_T) then

        if (mpi_rank == 0) then 
            call print_with_border("Task = calc_mu_vs_T")
            write(*,'(A)') " - Computing the chemical potential for different temperatures."
            write(*,*)
        end if 

        call setup_temp_arr()
        
        ! Set up chemical potential array
        allocate(mu_arr(n_temp))

        ! Calculate chemical potential for each temperature
        do i_T = 1, n_temp
            mu_arr(i_T) = calculate_mu_T(temp_arr(i_T), E_F)
            flush(output_unit)
        end do

        ! Write chemical potential vs temperature to file (only mu values, one per line)
        if (mpi_rank == 0) then
            if (interacting) then 
                mu_vs_T_filename = trim(file_tree)//"GRIT_chemical_potential_vs_T_interacting.txt"
            else if (.not. interacting) then
                mu_vs_T_filename = trim(file_tree)//"GRIT_chemical_potential_vs_T_noninteracting.txt"
            end if
            open(newunit=unit, file=mu_vs_T_filename, status='replace', action='write', iostat=ios)
            if (ios /= 0) then
                write(*,*) "Warning: Could not open file for writing: ", trim(mu_vs_T_filename)
            else
                do i_T = 1, n_temp
                    write(unit, '(F16.10)') mu_arr(i_T)
                end do
                close(unit)
                write(*,*) "- Chemical potential values written to: ", trim(mu_vs_T_filename)
            end if
        end if

        deallocate(mu_arr)
        deallocate(temp_arr)

    end if 

    flush(output_unit)




    ! ====================== !
    ! Task = calc_c_V_per_UC !
    ! ====================== !

    if (calc_c_V_per_UC) then

        if (allocated(E_arr)) deallocate(E_arr)
        allocate(E_arr(n_energy))

        if (mpi_rank == 0) then 
            write(*,*)
            write(*,*)
            call print_with_border("Task = calc_c_V_per_UC")
            write(*,'(A)') " - Computing the specific heat capacity for different temperatures."
        end if

        call setup_temp_arr()

        ! Set up chemical potential array and read values from file
        ! All processes need to allocate the array before broadcast
        allocate(mu_arr(n_temp))
        allocate(c_V_per_uc(n_temp))
        
        ! Read chemical potential values from the corresponding file (only master process)
        if (mpi_rank == 0) then
            if (interacting) then 
                mu_vs_T_filename = trim(file_tree)//"GRIT_chemical_potential_vs_T.txt"
            else if (.not. interacting) then
                mu_vs_T_filename = trim(file_tree)//"GRIT_chemical_potential_vs_T_noninteracting.txt"
            end if
            
            open(newunit=unit, file=mu_vs_T_filename, status='old', action='read', iostat=ios)
            if (ios /= 0) then
                write(*,*) " - Could not open file for reading: ", trim(mu_vs_T_filename)
                write(*,*) " - We will use mu(T) = E_F = ", E_F, " eV for all temperatures."
                mu_arr = E_F
                write(*,*) " - Make sure that this is really what you want, or " , & 
                           "run the task calc_mu_vs_T first to generate the chemical potential data."
                write(*,*)
            else
                do i_T = 1, n_temp
                    read(unit, '(F16.10)', iostat=ios) mu_arr(i_T)
                    if (ios /= 0) then
                        write(*,*) "Error reading chemical potential data at temperature index: ", i_T
                        exit
                    end if
                end do
                close(unit)
                write(*,*) "- Chemical potential values for different temperatures read from: ", trim(mu_vs_T_filename)
                write(*,*)
            end if
        end if
        
        ! Broadcast mu_arr to all MPI processes (fix datatype)
        if (mpi_size > 1) then
            call MPI_Bcast(mu_arr, n_temp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
        end if



        ! Calculate specific heat capacity for each temperature
        do i_T = 1, n_temp
            ! Set up Energy-array for integration of DOOS
            c_V_per_uc(i_T) = calculate_c_V_per_uc(temp_arr(i_T), mu_arr(i_T))
        end do



        ! Write specific heat capacity and temperature values to file
        if (mpi_rank == 0) then
            
            if (interacting) then 
                c_V_per_uc_filename = trim(file_tree)//"GRIT_c_V_per_uc_vs_T_interacting.txt"
            else if (.not. interacting) then
                c_V_per_uc_filename = trim(file_tree)//"GRIT_c_V_per_uc_vs_T_noninteracting.txt"
            end if

            open(newunit=unit, file=c_V_per_uc_filename, status='replace', action='write', iostat=ios)
            if (ios /= 0) then
                write(*,*) "Warning: Could not open file for writing: ", trim(c_V_per_uc_filename)
            else
                do i_T = 1, n_temp
                    write(unit, '(E14.6)') c_V_per_uc(i_T)
                end do
                close(unit)
                write(*,*)
                write(*,*) "- Specific heat values written to: ", trim(c_V_per_uc_filename)
            end if
            
            T_filename = trim(file_tree)//"GRIT_temperature_values.txt"
            open(newunit=unit, file=T_filename, status='replace', action='write', iostat=ios)
            if (ios /= 0) then
                write(*,*) "Warning: Could not open file for writing: ", trim(T_filename)
            else
                do i_T = 1, n_temp
                    write(unit, '(F16.10)') temp_arr(i_T)
                end do
                close(unit)
                write(*,*) "- Temperature values written to: ", trim(T_filename)
                write(*,*)
                write(*,*)
            end if
        end if
        
    end if 

    flush(output_unit)



    ! ===================================================================== !
    ! Make a linear fit (in temperature) to c_V vs. T and extract the slope !
    ! ===================================================================== !

    if (calc_gamma) then 
        ! Perform a linear fit. (TBD) 

        ! Only the master process does this.

        if (mpi_rank == 0) then 
            call print_with_border("Task = calc_gamma")
            write(*,'(A)') " - Computing the electronic specific heat coefficient."
            write(*,*)

            ! Set up heat capacity potential array and read values from file
            ! All processes need to allocate the array before broadcast
            if (.not. allocated(temp_arr)) allocate(temp_arr(n_temp))
            if (.not. allocated(c_V_per_uc)) allocate(c_V_per_uc(n_temp))


            ! ======================================== !
            ! Setting up filepaths and reading in Data !
            ! ======================================== !

            ! Temperature values
            T_filename = trim(file_tree)//"GRIT_temperature_values.txt"

            open(newunit=unit, file=T_filename, status='old', action='read', iostat=ios)
            if (ios /= 0) then
                write(*,*) " - Could not open file for reading: ", trim(T_filename)
                write(*,*)
            else
                do i_T = 1, n_temp
                    read(unit, '(F16.10)', iostat=ios) temp_arr(i_T)
                    if (ios /= 0) then
                        write(*,*) "Error reading temperature-data at index: ", i_T
                        exit
                    end if
                end do
                close(unit)
                write(*,*) "- Temperature values read from: ", trim(T_filename)
            end if

            ! Specific heat values

            if (interacting) then 
                c_V_per_uc_filename = trim(file_tree)//"GRIT_c_V_per_uc_vs_T.txt"
            else if (.not. interacting) then
                c_V_per_uc_filename = trim(file_tree)//"GRIT_c_V_per_uc_vs_T_noninteracting.txt"
            end if
            
            open(newunit=unit, file=c_V_per_uc_filename, status='old', action='read', iostat=ios)
            if (ios /= 0) then
                write(*,*) " - Could not open file for reading: ", trim(c_V_per_uc_filename)
                write(*,*)
            else
                do i_T = 1, n_temp
                    read(unit, '(F16.10)', iostat=ios) c_V_per_uc(i_T)
                    if (ios /= 0) then
                        write(*,*) "Error reading specific heat data at temperature index: ", i_T
                        exit
                    end if
                end do
                close(unit)
                write(*,*) "- Specific heat values for different temperatures read from: ", trim(c_V_per_uc_filename)
                write(*,*)
            end if

            ! Now we will make a linear fit with the ansatz 
            !   c_V_per_uc(T) = gamma * T + const. * T^3
            ! Which is equivalent to fitting 
            !   c_V_per_uc(T) / T = gamma + const. * T^2

            ! For this we call linear_fit from utils.f90
            
            call linear_fit(temp_arr**2, c_V_per_uc / temp_arr, n_temp, const, gamma)

            ! Now we print out the value of gamma and const.
            write(*,'(A)') " - Linear fit results:"
            write(*,'(A,F14.6,A)') "   gamma (electronic specific heat coefficient) = ", gamma * 1.0e3, " mJ / K^2"
            write(*,'(A,F14.6,A)') "   const (phononic specific heat coefficient)    = ", const * 1.0e3, " mJ / K^4"

        end if

    end if



    ! =========== !
    ! Cleaning up !
    ! =========== !

    if (allocated(E_arr))          deallocate(E_arr)
    if (allocated(temp_arr))       deallocate(temp_arr)
    if (allocated(mu_arr))         deallocate(mu_arr)
    if (allocated(epsilon_arr))    deallocate(epsilon_arr)
    if (allocated(omega_arr))      deallocate(omega_arr)
    if (allocated(c_V_per_uc))     deallocate(c_V_per_uc)
    
    if (mpi_rank == 0) then
        call print_with_border("Program finished successfully.")
    end if
    
    call finalize_mpi()

end program grit
