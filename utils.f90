module utils
    use, intrinsic :: ieee_arithmetic
    use modmain 
    use mpi
    !$ use omp_lib
    implicit none

contains

    ! ============================================ !
    ! MPI initialization and finalization routines !
    ! ============================================ !

    subroutine init_mpi()
        call mpi_init(mpi_ierr)
        call mpi_comm_rank(MPI_COMM_WORLD, mpi_rank, mpi_ierr)
        call mpi_comm_size(MPI_COMM_WORLD, mpi_size, mpi_ierr)
    end subroutine init_mpi


    subroutine setup_mpi_kbounds(mode)
        character(len=*), intent(in) :: mode

        select case (trim(mode))
            
            case ("Non-Interacting")
                
                i_k_start =  mpi_rank      * Nk_cube / mpi_size + 1
                i_k_end   = (mpi_rank + 1) * Nk_cube / mpi_size

                if (mpi_rank == mpi_size - 1) i_k_end = Nk_cube ! Handle remainder

                Nk_local_MPI = i_k_end - i_k_start + 1

            case ("Self-Energy")
                i_k_start =  mpi_rank      * Nk_SE / mpi_size + 1
                i_k_end   = (mpi_rank + 1) * Nk_SE / mpi_size

                if (mpi_rank == mpi_size - 1) i_k_end = Nk_SE ! Handle remainder

                Nk_local_MPI = i_k_end - i_k_start + 1

            case ("Uniform")
                i_k_start =  mpi_rank      * Ni_uniform_cube / mpi_size + 1
                i_k_end   = (mpi_rank + 1) * Ni_uniform_cube / mpi_size

                if (mpi_rank == mpi_size - 1) i_k_end = Ni_uniform_cube ! Handle remainder

                Nk_local_MPI = i_k_end - i_k_start + 1

        end select 
        
        ! Print information for debug-info
        if (mpi_rank == 0 .and. debug_mode) then
            write(*,*)
            write(*,'(A)') repeat('-', 110)
            write(*,*) "- MPI distribution for fine coupling data:"
            write(*,*) "  Number of MPI processes: ", mpi_size
            write(*,*) "  Total number of fine k-points: ", Nk_SE
            write(*,*) "  Each MPI process handles approx.: ", Nk_SE / mpi_size, " k-points."
            write(*,*) "  Current MPI rank: ", mpi_rank
            write(*,*) "  Assigned k-point indices from i_k = ", i_k_start, " to i_k = ", i_k_end
            write(*,'(A)') repeat('-', 110)
            write(*,*)
        end if


    end subroutine setup_mpi_kbounds


    subroutine read_coupling_array(mode)
        character(len=*), intent(in) :: mode

        select case (trim(mode))
            
            case ("Self-Energy")
                call setup_mpi_kbounds("Self-Energy")

                if (allocated(SE_g_arr)) deallocate(SE_g_arr)

                ! Allocate distributed coupling array - only local slice
                allocate(SE_g_arr(N_bnd, N_bnd, n_branches, Nk_local_MPI, Nq_tot_SE))

                ! Read distributed coupling data
                call read_g_data_binary_distributed(i_k_start, i_k_end, mode)

            case ("Uniform")
                call setup_mpi_kbounds("Uniform")

                if (allocated(uniform_g_arr)) deallocate(uniform_g_arr)

                ! Allocate distributed coupling array - only local slice
                allocate(uniform_g_arr(N_bnd, N_bnd, n_branches, Nk_local_MPI, Ni_uniform_cube))

                ! Read distributed coupling data
                call read_g_data_binary_distributed(i_k_start, i_k_end, mode)

            case default
                if (mpi_rank == 0) then
                    write(*,*) "Error: unknown mode = ", trim(mode)
                end if
                call MPI_Abort(MPI_COMM_WORLD, 1, mpi_ierr)
                stop
        end select



        ! 
        select case (trim(mode))
            
            case ("Self-Energy")
                if (allocated(SE_epsilon_kplusq_arr)) deallocate(SE_epsilon_kplusq_arr)

                ! Allocate distributed epsilon_{ k + q } array - only local slice
                allocate(SE_epsilon_kplusq_arr(N_bnd, Nk_local_MPI, Nq_tot_SE))

                ! Read distributed epsilon_{ k + q } data
                call read_epsilon_kplusq_data_binary_distributed(i_k_start, i_k_end)

        end select        
        


    end subroutine read_coupling_array


    subroutine finalize_mpi()
        call mpi_finalize(mpi_ierr)
    end subroutine finalize_mpi
           

    subroutine read_n_write_g_data(mode)
        character(len=*), intent(in) :: mode

        character(250) :: binary_filename
        integer        :: i_k, i_q, ios, m, n, nu, N_failures, N_tot_couplings, unit, &
                          Nq_tot_local, Nk_tot_local
        logical        :: binary_exists
        real(dp)       :: temp_val

        select case (trim(mode))
        
            case ("Self-Energy")
                Nk_tot_local = Nk_SE
                Nq_tot_local = Nq_tot_SE
                allocate(SE_g_arr(N_bnd, N_bnd, n_branches, Nk_tot_local, Nq_tot_local))
            case ("Specific-Heat")
                Nk_tot_local = Ni_uniform_cube
                Nq_tot_local = Ni_uniform_cube
                allocate(uniform_g_arr(N_bnd, N_bnd, n_branches, Nk_tot_local, Nq_tot_local))
            case default

                if (mpi_rank == 0) then
                    write(*,*) "Error: unknown omega mode = ", trim(mode)
                end if
                call MPI_Abort(MPI_COMM_WORLD, 1, mpi_ierr)
                stop

        end select


        if (mpi_rank == 0) then
            ! Write array-size info
            write(*,*) "  - Coupling-Data-Information for mode: ", trim(mode)
            write(*,*) "  - Total number of q-points for coupling data: ", Nq_tot_local
            write(*,*) "  - Total number of k-points for coupling data: ", Nk_tot_local
            write(*,*) "  - Number of bands: ", N_bnd
            write(*,*) "  - Number of phonon branches: ", n_branches
            write(*,*)
        end if

        N_tot_couplings = N_bnd * N_bnd * n_branches * Nk_tot_local * Nq_tot_local
    
        !-------------------------- !
        ! Construct binary filename !
        !-------------------------- !

        ! Temporarily store the text-filename as the binary filename (and then change it)

        select case (trim(mode))
            case ("Self-Energy")
                binary_filename = SE_couplings_file
            case ("Specific-Heat")
                binary_filename = uniform_couplings_file
        end select
        
    
        if (binary_filename(len_trim(binary_filename)-3:) == '.txt') then
            binary_filename(len_trim(binary_filename)-3:) = '.bin'
        else
            binary_filename = trim(binary_filename) // '.bin'
        end if
    
    
        !---------------------- !
        ! Check for binary file !
        !-----------------------! 

        inquire(file=binary_filename, exist=binary_exists)
    
        if (binary_exists) then
            
            !----------------- !
            ! Load from binary !
            !----------------- !

            if (mpi_rank == 0) then
                write(*,*) "- Loading the EPW-g-data from binary file..."
            end if
    
            open(newunit=unit, file=binary_filename, status='old', action='read', &
                 form='unformatted', access='stream', iostat=ios)
    
            if (ios /= 0) then
                if (mpi_rank == 0) then
                    write(*,*) "Error opening binary file: ", trim(binary_filename)
                end if
                call mpi_finalize(mpi_ierr)
                call MPI_Abort(MPI_COMM_WORLD, 1, mpi_ierr)
            end if

            select case (trim(mode))
                case ("Self-Energy")
                    read(unit, iostat=ios) SE_g_arr
                case ("Specific-Heat")
                    read(unit, iostat=ios) uniform_g_arr
            end select

            if (ios /= 0) then
                if (mpi_rank == 0) then
                    write(*,*) "Error reading binary file."
                end if
                close(unit)
                call mpi_finalize(mpi_ierr)
                call MPI_Abort(MPI_COMM_WORLD, 1, mpi_ierr)
            end if
    
            close(unit)
    
        else

            !---------------------------------------- !
            ! Load from text and create binary backup !
            !---------------------------------------- !

            if (mpi_rank == 0) then
                write(*,*) "- Loading the EPW-g-data from text file (this may take a while)..."
            end if

            select case (trim(mode))
                case ("Self-Energy")
                    open(newunit=unit, file=SE_couplings_file, status='old', action='read', iostat=ios)
                case ("Specific-Heat")
                    open(newunit=unit, file=uniform_couplings_file, status='old', action='read', iostat=ios)
            end select
    
            if (ios /= 0) then
                if (mpi_rank == 0) then
                    write(*,*) "Error opening Coupling-file."
                end if
                call MPI_Abort(MPI_COMM_WORLD, 1, mpi_ierr)
            end if

            N_failures = 0
    
            do i_q = 1, Nq_tot_local
                do i_k = 1, Nk_tot_local
                    do m = 1, N_bnd
                        do n = 1, N_bnd
                            do nu = 1, n_branches
                                
                                read(unit, *, iostat=ios) temp_val
                                
                                if (ios /= 0) then
                                    if (mpi_rank == 0) then
                                        write(*,*) "Error reading coupling data"
                                    end if
                                    call MPI_Abort(MPI_COMM_WORLD, 1, mpi_ierr)
                                end if

                                if (temp_val > 500.0_dp .or. i_q .eq. 1) then
                                    N_failures = N_failures + 1
                                    temp_val = 0.0_dp
                                end if

                                select case (trim(mode))
        
                                    case ("Self-Energy")
                                        SE_g_arr(m,n,nu,i_k,i_q) = temp_val
                                    case ("Specific-Heat")
                                        uniform_g_arr(m,n,nu,i_k,i_q) = temp_val
                                end select


                            end do
                        end do
                    end do
                end do
            end do
    
            close(unit)


    
            if (mpi_rank == 0) then
                write(*,*) "   - Text coupling data loaded successfully."
                write(*,*) "   - Creating binary file."
                write(*,*)
                write(*,*) "   - Number of unusually large coupling values that were set to zero by hand = ", N_failures
                write(*,*) "   - Total number of coupling values read = ", N_tot_couplings
                ! Print ratio in percent for debug-info
                write(*,'(A,F8.3,A)') "   - Ratio of unusually large values = N_failures / N_tot_couplings * 100 =  ", &
                      real(N_failures) / real(N_tot_couplings) * 100.0_dp, " %"
                write(*,*)
            end if
    
            open(newunit=unit, file=binary_filename, status='replace', action='write', &
                 form='unformatted', access='stream', iostat=ios)
    
            if (ios == 0) then

                select case (trim(mode))
                    case ("Self-Energy")
                         write(unit, iostat=ios) SE_g_arr
                    case ("Specific-Heat")
                        if (symmetrize) then
                            call symmetrize_n_write_g_data()  ! This subroutine also creates the binary file, so we can skip the write here.
                            if (mpi_rank == 0) then
                                write(*,*) "   - The coupling data was symmetrized by hand to fix small numerical asymmetries."
                                write(*,*)
                            end if
                        else 
                            write(unit, iostat=ios) uniform_g_arr
                        end if 
                end select


                close(unit)
    
                if (ios == 0 .and. mpi_rank == 0) then
                    write(*,*) "- Binary backup created: ", trim(binary_filename)
                    write(*,*)
                end if
            else
                if (mpi_rank == 0) then
                    write(*,*) " - Warning: Could not create binary file"
                end if
            end if
        end if
    
        !-------------- !
        ! Final message !
        !-------------- !

        if (mpi_rank == 0) then
            write(*,*) "- The electron-phonon coupling-array was successfully loaded for the mode = ", trim(mode)
            write(*,*) "  Note: g_arr[m,n,nu,i_k,i_q] = g_{m,n,nu}(k_i,q_i)"
            write(*,*)
        end if

        if (allocated(SE_epsilon_arr)) deallocate(SE_epsilon_arr)
        if (allocated(SE_g_arr)) deallocate(SE_g_arr)
        if (allocated(uniform_g_arr)) deallocate(uniform_g_arr)
    
    end subroutine read_n_write_g_data


    subroutine read_k_data()
        integer :: unit, i, ios

        if (.not. allocated(uniform_wavevector_arr)) allocate(uniform_wavevector_arr(3, Ni_uniform_cube))

        open(newunit=unit, file=uniform_k_points_file, status='old', action='read', iostat=ios)

        if (ios /= 0) then
            if (mpi_rank == 0) then
                write(*,*) "Error opening file: ", trim(uniform_k_points_file)
            end if
            stop
        end if

        do i = 1, Ni_uniform_cube
            read(unit, *, iostat=ios) uniform_wavevector_arr(:, i)
            if (ios /= 0) then
                if (mpi_rank == 0) then
                    write(*,*) "Error reading k-point data at line ", i
                end if
                stop
            end if
        end do

        close(unit)

        if (mpi_rank == 0) then
            write(*,*) "- The wavevector-array was successfully loaded."
            write(*,*) "  Note: wavevector_arr[i] = k_i = q_i in reciprocal coordinates."
            write(*,*)
        end if
    end subroutine read_k_data


    function get_i_kplusq(i_k, i_q) result(i_kplusq)

        integer, intent(in) :: i_k, i_q
        integer             :: i_kplusq
        integer             :: i_tmp
        real(dp)            :: k(3), q(3), k_plus_q(3)
        real(dp)            :: shifted_vector(3), distance, min_distance

        ! Get k and q vectors
        k = uniform_wavevector_arr(:, i_k)
        q = uniform_wavevector_arr(:, i_q)

        ! Calculate k + q with modulo operation
        k_plus_q = modulo(k + q, 1.0)

        ! Find the closest wavevector to k_plus_q
        min_distance = huge(min_distance)
        i_kplusq = 1

        do i_tmp = 1, Ni_uniform_cube
            shifted_vector = uniform_wavevector_arr(:, i_tmp) - k_plus_q
            distance       = sqrt(sum(shifted_vector**2))  ! L2 norm

            if (distance < min_distance) then
                min_distance = distance
                i_kplusq     = i_tmp
            end if

        end do

    end function get_i_kplusq


    function get_i_minusq(i_q) result(i_minusq)

        integer, intent(in) :: i_q
        integer             :: i_minusq
        integer             :: i_tmp
        real(dp)            :: q(3), minusq(3)
        real(dp)            :: shifted_vector(3), distance, min_distance

        ! Get q vector
        q = uniform_wavevector_arr(:, i_q)

        ! Calculate -q with modulo operation
        minusq = modulo(-q, 1.0)

        ! Find the closest wavevector to -q
        min_distance = huge(min_distance)
        i_minusq = 1

        do i_tmp = 1, Ni_uniform_cube
            shifted_vector = uniform_wavevector_arr(:, i_tmp) - minusq
            distance       = sqrt(sum(shifted_vector**2))  ! L2 norm

            if (distance < min_distance) then
                min_distance = distance
                i_minusq     = i_tmp
            end if

        end do

    end function get_i_minusq


    subroutine symmetrize_n_write_g_data()
        ! Enforcing inversion and hermiticity "by hand".
        ! Remark: 
        ! If we treat the couplings as complex-numbers, we must take the hermitian conjugate.
        integer        :: unit, ios
        integer        :: m, n, nu, i_k, i_q
        integer        :: i_minusq, i_minusk, i_kplusq
        real(dp)       :: temp_val
        character(250) :: binary_filename

        do i_q = 1, Ni_uniform_cube
            i_minusq = get_i_minusq(i_q)
            do i_k = 1, Ni_uniform_cube
                i_minusk = get_i_minusq(i_k)
                i_kplusq = get_i_kplusq(i_k, i_q)
                do m = 1, N_bnd
                    do n = 1, N_bnd
                        do nu = 1, n_branches
                            ! Fixing inversion-symmetry "by hand":
                            temp_val = 0.5 * (uniform_g_arr(m, n, nu, i_k, i_q) + uniform_g_arr(m, n, nu, i_minusk, i_minusq))
                            uniform_g_arr(m, n, nu, i_k, i_q) = temp_val
                            uniform_g_arr(m, n, nu, i_minusk, i_minusq) = temp_val

                            ! Fixing hermiticity "by hand":
                            temp_val = 0.5 * (uniform_g_arr(m, n, nu, i_k, i_q) + uniform_g_arr(n, m, nu, i_kplusq, i_minusq))
                            uniform_g_arr(m, n, nu, i_k, i_q) = temp_val
                            uniform_g_arr(n, m, nu, i_kplusq, i_minusq) = temp_val
                        end do
                    end do
                end do
            end do
        end do
        
        ! Create binary filename by replacing .txt with .bin
        binary_filename = uniform_couplings_file

        if (binary_filename(len_trim(binary_filename)-3:) == '.txt') then
            binary_filename(len_trim(binary_filename)-3:) = '.bin'
        else
            binary_filename = trim(binary_filename) // '.bin'
        end if

        open(newunit=unit, file=binary_filename, status='replace', action='write', &
                 form='unformatted', access='stream', iostat=ios)
            
        if (ios == 0) then
            write(unit, iostat=ios) uniform_g_arr
            close(unit)
            if (ios == 0) then
                if (mpi_rank == 0) then
                    write(*,*)
                end if
            else
                if (mpi_rank == 0) then
                    write(*,*) " - Warning: Error writing binary file"
                end if
            end if
        else
            if (mpi_rank == 0) then
                write(*,*) " - Warning: Could not create binary file"
            end if
        end if
         
    end subroutine symmetrize_n_write_g_data



    ! ========================================================================= !
    ! Subroutines that are shared for non- and -interacting cases (coarse grid) !
    ! ========================================================================= !

    subroutine read_epsilons_data(mode)
        character(len=*), intent(in) :: mode
        integer  :: unit, i_k, m, Nk_tot_local, ios
        real(dp) :: temp_val



        ! Decide which file and which array to use
        select case (trim(mode))
        case ("Non-Interacting")
            Nk_tot_local = Nk_cube
            allocate(epsilon_arr(N_bnd, Nk_cube))
            open(newunit=unit, file=epsilons_file, status='old', action='read', iostat=ios)
        case ("Uniform")
            Nk_tot_local = Ni_uniform_cube
            allocate(uniform_epsilon_arr(N_bnd, Ni_uniform_cube))
            open(newunit=unit, file=uniform_epsilons_file, status='old', action='read', iostat=ios) 
        case default

            if (mpi_rank == 0) then
                write(*,*) "Error: unknown mode = ", trim(mode)
            end if
            call MPI_Abort(MPI_COMM_WORLD, 1, mpi_ierr)
            stop

        end select

        

        if (ios /= 0) then
            if (mpi_rank == 0) then
                write(*,*) "Error opening file: ", trim(epsilons_file)
            end if
            call MPI_Abort(MPI_COMM_WORLD, 1, mpi_ierr)
        end if

        do i_k = 1, Nk_tot_local
            do m = 1, N_bnd
                read(unit, *, iostat=ios) temp_val
                if (ios /= 0) then
                    if (mpi_rank == 0) then
                        write(*,*) "Error reading epsilon data"
                    end if
                    call MPI_Abort(MPI_COMM_WORLD, 1, mpi_ierr)
                end if


                select case (trim(mode))
                    case ("Non-Interacting")
                        epsilon_arr(m, i_k) = temp_val
                    case ("Uniform")
                        uniform_epsilon_arr(m, i_k) = temp_val
                end select 
                
            end do
        end do

        close(unit)

        if (mpi_rank == 0) then
            write(*,*) "- The electronic-energy-array was successfully loaded for mode = ", trim(mode)
            write(*,*) "  Note: epsilon_arr[m, i_k] = epsilon_m ( k_i )"
            write(*,*)
        end if

    end subroutine read_epsilons_data


    subroutine read_SE_epsilons_data()

        integer  :: i_k, ios, m, n_values, unit
        real(dp) :: temp_val

        n_values = 0

        open(newunit=unit, file=SE_epsilons_file, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            if (mpi_rank == 0) then
                write(*,*) "Error opening file: ", trim(SE_epsilons_file)
            end if
            call MPI_Abort(MPI_COMM_WORLD, 1, mpi_ierr)
        end if

        !----------------------------------------- ! 
        ! First pass: count total number of values !
        !----------------------------------------- ! 
        do
            read(unit, *, iostat=ios) temp_val
            if (ios < 0) exit        ! End of file
            if (ios > 0) then
                write(*,*) "Error while scanning epsilon file."
                call MPI_Abort(MPI_COMM_WORLD, 1, mpi_ierr)
            end if
            n_values = n_values + 1
        end do

        if (mod(n_values, N_bnd) /= 0) then
            write(*,*) "Error: epsilon file size not divisible by N_bnd."
            call MPI_Abort(MPI_COMM_WORLD, 1, mpi_ierr)
        end if

        Nk_SE = n_values / N_bnd

        !------------------------------- !
        ! Allocate (or reallocate) array !
        !------------------------------- !

        if (allocated(SE_epsilon_arr)) deallocate(SE_epsilon_arr)
        allocate(SE_epsilon_arr(N_bnd, Nk_SE))

        rewind(unit)

        !------------------------- !
        ! Second pass: actual read !
        !------------------------- !

        do i_k = 1, Nk_SE
            do m = 1, N_bnd
                read(unit, *, iostat=ios) SE_epsilon_arr(m, i_k)
                if (ios /= 0) then
                    write(*,*) "Error reading epsilon data."
                    call MPI_Abort(MPI_COMM_WORLD, 1, mpi_ierr)
                end if
            end do
        end do

        close(unit)

        if (mpi_rank == 0) then
            write(*,*) "- The Fermi-electronic-energy-array was successfully loaded."
            write(*,*) "  Detected number of k-points: ", Nk_SE
            write(*,*) "  Note: SE_epsilon_arr(m, i_k) = epsilon_m(k_i)"
            write(*,*)
        end if

    end subroutine read_SE_epsilons_data


    subroutine read_omegas_data(mode)
        character(len=*), intent(in) :: mode

        integer  :: i_q, ios, nu, n_values, Nq_tot_local, unit
        real(dp) :: temp_val

        n_values = 0

        ! Decide which file and which array to use
        select case (trim(mode))
        case ("Non-Interacting")
            open(newunit=unit, file=omegas_file, status='old', action='read', iostat=ios)
        case ("Self-Energy")
            open(newunit=unit, file=SE_omegas_file, status='old', action='read', iostat=ios)
        case ("Uniform")
            open(newunit=unit, file=uniform_omegas_file, status='old', action='read', iostat=ios)
        case default

            if (mpi_rank == 0) then
                write(*,*) "Error: unknown omega mode = ", trim(mode)
            end if
            call MPI_Abort(MPI_COMM_WORLD, 1, mpi_ierr)
            stop

        end select

        if (ios /= 0) then
            if (mpi_rank == 0) then
                write(*,*) "Error opening Omegas-file in read_omegas_data with mode: ", trim(mode)
            end if
            call MPI_Abort(MPI_COMM_WORLD, 1, mpi_ierr)
        end if

        !----------------------------------------- ! 
        ! First pass: count total number of values !
        !----------------------------------------- ! 
        do
            read(unit, *, iostat=ios) temp_val
            if (ios < 0) exit        ! End of file
            if (ios > 0) then
                write(*,*) "Error while scanning omega file."
                call MPI_Abort(MPI_COMM_WORLD, 1, mpi_ierr)
            end if
            n_values = n_values + 1
        end do

        if (mod(n_values, n_branches) /= 0) then
            write(*,*) "Error: omega file size not divisible by n_branches."
            call MPI_Abort(MPI_COMM_WORLD, 1, mpi_ierr)
        end if


        !------------------------------- !
        ! Allocate (or reallocate) array !
        !------------------------------- !

        select case (trim(mode))
        case ("Non-Interacting")
            Nq_tot_0    = n_values / n_branches
            Nq_tot_local = Nq_tot_0

            allocate(omega_arr(n_branches, Nq_tot_0))

        case ("Self-Energy")
            Nq_tot_SE  = n_values / n_branches
            Nq_tot_local = Nq_tot_SE
            
            if (.not. allocated(SE_omega_arr)) allocate(SE_omega_arr(n_branches, Nq_tot_SE))
        case ("Uniform")
                Nq_tot_local = Ni_uniform_cube
            if (.not. allocated(uniform_omega_arr)) allocate(uniform_omega_arr(n_branches, Ni_uniform_cube))
        case default
            if (mpi_rank == 0) then
                write(*,*) "Error: unknown omega mode = ", trim(mode)
            end if
            call MPI_Abort(MPI_COMM_WORLD, 1, mpi_ierr)
            stop
        end select

        rewind(unit)

        !------------------------- !
        ! Second pass: actual read !
        !------------------------- !

        do i_q = 1, Nq_tot_local
            do nu = 1, n_branches

                select case (trim(mode))
                    case ("Non-Interacting")
                        read(unit, *, iostat=ios) omega_arr(nu, i_q)
                    case ("Self-Energy")
                        read(unit, *, iostat=ios) SE_omega_arr(nu, i_q)
                    case ("Uniform")
                        read(unit, *, iostat=ios) uniform_omega_arr(nu, i_q)
                end select

                if (ios /= 0) then
                    write(*,*) "Error reading omega data."
                    call MPI_Abort(MPI_COMM_WORLD, 1, mpi_ierr)
                end if
            end do
        end do

        close(unit)

        if (mpi_rank == 0) then
            write(*,*) "- The omega-energy-array was successfully loaded for the mode : ", trim(mode)
            write(*,*) "  Detected number of q-points: ", Nq_tot_local

            select case (trim(mode))
            case ("Non-Interacting")
                write(*,*) "  Expected number of q-points: Nq_cube = ", Nq_cube
                write(*,*) "  Note: omega_arr(nu, i_q) = omega_nu ( q_i )"
            case ("Self-Energy")
                if (.not. Gamma_centered_q) then 
                    write(*,*) "  Expected number of q-points: Nq_SE_cube = ", Nq_SE_cube
                else if (Gamma_centered_q) then
                    write(*,*) "  Expected number of q-points can only be concluded by looking at the Q_Path.dat file."
                end if 
                write(*,*) "  Note: SE_omega_arr(nu, i_q) = omega_nu ( q_i )"
            case ("Uniform")
                write(*,*) "  Expected number of q-points: Ni_uniform_cube = ", Ni_uniform_cube
                write(*,*) "  Note: uniform_omega_arr(nu, i_q) = omega_nu ( q_i )"
            end select

            write(*,*)
        end if

    end subroutine read_omegas_data


    subroutine read_coupling_data_binary_slice(binary_filename, i_k_start, i_k_end, mode)
        character(len=*), intent(in) :: mode
        character(len=*), intent(in) :: binary_filename
        integer,          intent(in) :: i_k_end, i_k_start
        integer(kind=8)              :: block_size, elem_size, offset
        integer                      :: i_k, ios, local_i_k, unit
        

        select case (trim(mode))
            
            case ("Self-Energy")
                elem_size = storage_size(SE_g_arr(1,1,1,1,1)) / 8
                block_size = N_bnd * N_bnd * n_branches * Nq_tot_SE * elem_size
            case ("Uniform")
                elem_size = storage_size(uniform_g_arr(1,1,1,1,1)) / 8
                block_size = N_bnd * N_bnd * n_branches * Ni_uniform_cube * elem_size
            case default
                if (mpi_rank == 0) then
                    write(*,*) "Error: unknown mode = ", trim(mode)
                end if
                call MPI_Abort(MPI_COMM_WORLD, 1, mpi_ierr)
                stop

        end select
        

        open(newunit=unit, file=binary_filename, status='old', access='stream', &
             form='unformatted', action='read', iostat=ios)
        if (ios /= 0) stop "Error opening binary file"

        do local_i_k = 1, i_k_end - i_k_start + 1
            i_k    = i_k_start + local_i_k - 1
            offset = (i_k-1) * block_size

            select case (trim(mode))
            
                case ("Self-Energy")
                    read(unit, pos=offset+1, iostat=ios) SE_g_arr(:,:,:,local_i_k,:)
                case ("Uniform")
                    read(unit, pos=offset+1, iostat=ios) uniform_g_arr(:,:,:,local_i_k,:)

            end select

            if (ios /= 0) stop "Error reading slice"
        end do

        close(unit)
    end subroutine read_coupling_data_binary_slice


    subroutine read_g_data_binary_distributed(i_k_start, i_k_end, mode)
        character(len=*), intent(in) :: mode
        integer, intent(in)          :: i_k_start, i_k_end
    
        character(len=250) :: binary_filename
        logical            :: binary_exists
    

        ! Temporarily store the text-filename as the binary filename (and then change it)
        select case (trim(mode))
            
            case ("Self-Energy")
                binary_filename = SE_couplings_file
            case ("Uniform")
                binary_filename = uniform_couplings_file

        end select
        

        if (binary_filename(len_trim(binary_filename)-3:) == '.txt') then
            binary_filename(len_trim(binary_filename)-3:) = '.bin'
        else
            binary_filename = trim(binary_filename) // '.bin'
        end if
    


        !----------------------------------- !
        ! Check for existence of binary file !
        !----------------------------------- !

        inquire(file=binary_filename, exist=binary_exists)
    
        if (.not. binary_exists) then
            if (mpi_rank == 0) then
                call print_with_border("Error!")
                write(*,*) "  Error: MPI is non trivial and no binary file is available."
                write(*,*) "         Rerun the program with one MPI-process to create it."
                write(*,*)
            end if
            call MPI_Abort(MPI_COMM_WORLD, 1, mpi_ierr)
        end if
    


        !--------------------------- !
        ! Load the appropriate slice !
        !--------------------------- !
    
        call read_coupling_data_binary_slice(binary_filename, i_k_start, i_k_end, mode)

    
    end subroutine read_g_data_binary_distributed


    subroutine read_n_write_SE_epsilon_kplusq_data()
    
        character(250) :: binary_filename
        integer        :: i_k, i_q, ios, m, unit
        logical        :: binary_exists
        real(dp)       :: temp_val


        allocate(SE_epsilon_kplusq_arr(N_bnd, Nk_SE, Nq_tot_SE))
    
        !-------------------------- !
        ! Construct binary filename !
        !-------------------------- !

        ! Temporarily store the text-filename as the binary filename (and then change it)
        binary_filename = SE_epsilons_kplusq_file
    
        if (binary_filename(len_trim(binary_filename)-3:) == '.txt') then
            binary_filename(len_trim(binary_filename)-3:) = '.bin'
        else
            binary_filename = trim(binary_filename) // '.bin'
        end if
    
    
        !---------------------- !
        ! Check for binary file !
        !-----------------------! 

        inquire(file=binary_filename, exist=binary_exists)

        if (binary_exists) then
            
            !----------------- !
            ! Load from binary !
            !----------------- !

            if (mpi_rank == 0) then
                write(*,*) "- Loading the fine epsilon_{k + q}-data from binary file..."
            end if
    
            open(newunit=unit, file=binary_filename, status='old', action='read', &
                 form='unformatted', access='stream', iostat=ios)
    
            if (ios /= 0) then
                if (mpi_rank == 0) then
                    write(*,*) "Error opening binary file: ", trim(binary_filename)
                end if
                call mpi_finalize(mpi_ierr)
                call MPI_Abort(MPI_COMM_WORLD, 1, mpi_ierr)
            end if
    
            read(unit, iostat=ios) SE_epsilon_kplusq_arr
    
            if (ios /= 0) then
                if (mpi_rank == 0) then
                    write(*,*) "Error reading epsilon binary file."
                end if
                close(unit)
                call mpi_finalize(mpi_ierr)
                call MPI_Abort(MPI_COMM_WORLD, 1, mpi_ierr)
            end if
    
            close(unit)
    
        else

            !---------------------------------------- !
            ! Load from text and create binary backup !
            !---------------------------------------- !

            if (mpi_rank == 0) then
                write(*,*) "- Loading the fine epsilon_{ k + q } data from text file (this may take a while)..."
            end if
    
            open(newunit=unit, file=SE_epsilons_kplusq_file, status='old', action='read', iostat=ios)
    
            if (ios /= 0) then
                if (mpi_rank == 0) then
                    write(*,*) "Error opening file: ", trim(SE_epsilons_kplusq_file)
                end if
                call MPI_Abort(MPI_COMM_WORLD, 1, mpi_ierr)
            end if
    
            do i_q = 1, Nq_tot_SE
                do i_k = 1, Nk_SE
                    do m = 1, N_bnd
                            
                        read(unit, *, iostat=ios) temp_val
                        if (ios /= 0) then
                            if (mpi_rank == 0) then
                                write(*,*) "Error reading epsilon_{k + q} data"
                            end if
                            call MPI_Abort(MPI_COMM_WORLD, 1, mpi_ierr)
                        end if
                        
                        SE_epsilon_kplusq_arr(m,i_k,i_q) = temp_val
                            
                    end do
                end do
            end do
    
            close(unit)
    
            if (mpi_rank == 0) then
                write(*,*) "   - Text epsilon_{ k + q } data loaded successfully."
                write(*,*) "   - Creating binary file."
            end if
    
            open(newunit=unit, file=binary_filename, status='replace', action='write', &
                 form='unformatted', access='stream', iostat=ios)
    
            if (ios == 0) then
                write(unit, iostat=ios) SE_epsilon_kplusq_arr

                close(unit)
    
                if (ios == 0 .and. mpi_rank == 0) then
                    write(*,*) "- Binary backup created: ", trim(binary_filename)
                    write(*,*)
                end if
            else
                if (mpi_rank == 0) then
                    write(*,*) " - Warning: Could not create binary file"
                end if
            end if
        end if
    
        !-------------- !
        ! Final message !
        !-------------- !

        if (mpi_rank == 0) then
            write(*,*) "- The epsilon_{ k + q } data was successfully loaded."
            write(*,*) "  Note: SE_epsilon_kplusq_arr[m,i_k,i_q] = epsilon_{m} (k_i + q_i)"
            write(*,*)
        end if

        deallocate(SE_epsilon_kplusq_arr)
    
    end subroutine read_n_write_SE_epsilon_kplusq_data


    subroutine read_epsilon_kplusq_data_binary_slice(binary_filename, i_k_start, i_k_end)
        character(len=*), intent(in) :: binary_filename
        integer,          intent(in) :: i_k_start, i_k_end

        integer           :: i_k, ios, local_i_k, unit
        integer(kind=8)   :: offset, elem_size, block_size


        elem_size = storage_size(SE_epsilon_kplusq_arr(1,1,1)) / 8
        block_size = N_bnd * Nq_tot_SE * elem_size

        open(newunit=unit, file=binary_filename, status='old', access='stream', &
             form='unformatted', action='read', iostat=ios)
        if (ios /= 0) stop "Error opening binary file"

        do local_i_k = 1, i_k_end - i_k_start + 1
            i_k    = i_k_start + local_i_k - 1
            offset = (i_k-1) * block_size
            read(unit, pos=offset+1, iostat=ios) SE_epsilon_kplusq_arr(:,local_i_k,:)
            if (ios /= 0) stop "Error reading slice"
        end do

        close(unit)
    end subroutine read_epsilon_kplusq_data_binary_slice


    subroutine read_epsilon_kplusq_data_binary_distributed(i_k_start, i_k_end)
        integer, intent(in) :: i_k_end, i_k_start
    
        character(len=250) :: binary_filename
        logical            :: binary_exists
        

        binary_filename = SE_epsilons_kplusq_file
    
        !-------------------------------------- !
        ! Replace .txt by .bin (or append .bin) !
        !-------------------------------------- !

        if (binary_filename(len_trim(binary_filename)-3:) == '.txt') then
            binary_filename(len_trim(binary_filename)-3:) = '.bin'
        else
            binary_filename = trim(binary_filename) // '.bin'
        end if
    
        !----------------------------------- !
        ! Check for existence of binary file !
        !----------------------------------- !

        inquire(file=binary_filename, exist=binary_exists)
    
        if (.not. binary_exists) then
            if (mpi_rank == 0) then
                call print_with_border("Error!")
                write(*,*) "  Error: MPI is non trivial and no binary file is available."
                write(*,*) "         Rerun the program with one MPI-process to create it."
                write(*,*)
            end if
            call MPI_Abort(MPI_COMM_WORLD, 1, mpi_ierr)
        end if
    
        !--------------------------- !
        ! Load the appropriate slice !
        !--------------------------- !

        if (mpi_rank == 0 .and. debug_mode) then
            write(*,*) "- Loading fine epsilon_{ k + q } slice from binary file..."
        end if
    
        call read_epsilon_kplusq_data_binary_slice(binary_filename, i_k_start, i_k_end)
    
    end subroutine read_epsilon_kplusq_data_binary_distributed



    ! ====================================== !
    ! Subroutines for DOS / Number of States !
    ! ====================================== !

    function dN_States_per_UC_over_dE_i(E, T, mu, close2E_F) result(dn_de)
        logical , intent(in) :: close2E_F
        real(dp), intent(in) :: E, T, mu

        complex(dp) :: denominator, denominator2, G_ret_mk, Sigma
        integer     :: i_k_local, i_k_plus_q, i_q, m, n, nu
        real(dp)    :: bose_factor, coupling_prefactor, dn_de, epsilon_mk, epsilon_kplusq, fermi, local_sum_val, &
                       numerator1, numerator2, prefactor, sum_val
        

        ! The summation extends over two different k-point grids, depending on whether
        ! we are close to the Fermi level (E = 0) or not.

        ! Note: Prefactor of 2 because of spin-degeneracy
        prefactor = -2.0 / (pi * real(Nk_cube))
        local_sum_val = 0.0

        ! ------------------------------------- !
        ! Part 1: Non-Interacting Contributions !
        ! ------------------------------------- !

        call setup_mpi_kbounds("Non-Interacting")

        ! MPI parallelization over i_k_local, OpenMP over inner loops
        do i_k_local = i_k_start, i_k_end
            do m = 1, N_bnd
                epsilon_mk = epsilon_arr(m, i_k_local)
                Sigma = 0.0_dp

                ! Check if epsilon_mk is close to E_F, i.e. if | epsilon_mk - mu | < fsthick
                if (close2E_F .and. abs(epsilon_mk - mu) < fsthick) cycle
                
                G_ret_mk = 1.0 / (cmplx(E, eta, kind=dp) - (epsilon_mk + Sigma - mu))
                local_sum_val = local_sum_val + aimag(G_ret_mk)

            end do
        end do

        ! --------------------------------- !
        ! Part 2: Interacting Contributions !
        ! --------------------------------- !

        if (close2E_F) then 
            
            call setup_mpi_kbounds("Self-Energy")


            ! MPI parallelization over i_k_local, OpenMP over inner loops
            do i_k_local = i_k_start, i_k_end
                do m = 1, N_bnd

                    epsilon_mk = SE_epsilon_arr(m, i_k_local)

                    ! Compute the Fock-contribution Sigma
                    Sigma = 0.0
                    if (interacting .and. close2E_F) then 
                        !$OMP PARALLEL DO PRIVATE(i_q, n, nu, coupling_prefactor, i_k_plus_q, &
                        !$OMP                     bose_factor, fermi, numerator1, numerator2, &
                        !$OMP                     denominator, denominator2, epsilon_kplusq) &
                        !$OMP             REDUCTION(+:Sigma) SCHEDULE(STATIC)
                        do i_q = 1, Nq_tot_SE
                            do n = 1, N_bnd
                                do nu = 1, n_branches

                                    ! Skip iteration if coupling is zero
                                    if (SE_g_arr(n, m, nu, i_k_local - i_k_start + 1, i_q) < threshold_couplings .or. &
                                                                     SE_omega_arr(nu, i_q) < threshold_omegas) then
                                        cycle
                                    end if

                                    ! Read in the coupling g_{m,n,nu}(p,q) using local index  
                                    coupling_prefactor = (SE_g_arr(m, n, nu, i_k_local - i_k_start + 1, i_q) * 1.0e-3_dp)**2 &
                                                          / real(Nq_SE_cube)

                                    epsilon_kplusq = SE_epsilon_kplusq_arr(n, i_k_local - i_k_start + 1, i_q)

                                    ! Bose-Einstein distribution
                                    bose_factor = bose_function(SE_omega_arr(nu, i_q) * 1.0e-3_dp, T)

                                    ! First term
                                    fermi = fermi_function(epsilon_kplusq - mu, T)
                                    numerator1 = bose_factor + fermi
                                    denominator = cmplx(E, eta_q, kind=dp) - (epsilon_kplusq - mu) &
                                                                         + SE_omega_arr(nu, i_q) * 1.0e-3_dp

                                    ! Second term
                                    numerator2 = bose_factor + 1.0 - fermi
                                    denominator2 = cmplx(E, eta_q, kind=dp) - (epsilon_kplusq - mu) &
                                                                          - SE_omega_arr(nu, i_q) * 1.0e-3_dp
                                    Sigma = Sigma + coupling_prefactor * (numerator1 / denominator + numerator2 / denominator2)

                                end do
                            end do
                        end do
                        !$OMP END PARALLEL DO
                    end if 

                    G_ret_mk = 1.0 / (cmplx(E, eta, kind=dp) - (epsilon_mk + Sigma - mu))
                    local_sum_val = local_sum_val + aimag(G_ret_mk)

                end do
            end do

        end if 

        ! Sum contributions from all MPI processes
        call MPI_Allreduce(local_sum_val, sum_val, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_ierr)

        dn_de = prefactor * sum_val
    end function dN_States_per_UC_over_dE_i


    subroutine calculate_no_of_States_per_UC(T, mu)
        real(dp), intent(in)  :: T, mu

        real(dp), allocatable :: integrand(:)
        integer               :: file_unit, ios = -999
        character(len=50)     :: filename
        
        allocate(integrand(n_energy))

        ! Set up Energy-array for integration of DOS
        call setup_Energy_arr_DOS()

        ! Only master process handles file I/O in the case that we are debugging
        if (mpi_rank == 0 .and. debug_mode) then
            ! Create filename with temperature and chemical potential info for states
            write(filename, '(A,F6.1,A,F8.4,A)') trim(file_tree)//"GRIT_dN_states_per_UC_dE_T", T, "_mu", mu, ".txt" 

            ! Open file for writing
            open(newunit=file_unit, file=filename, status='replace', action='write', iostat=ios)
            if (ios /= 0) then
                write(*,*) "Warning: Could not open file for writing: ", trim(filename)
            end if

            ! Write header to file
            if (ios == 0) then
                write(file_unit, '(A)') "# Energy (eV)    dN_states_per_UC/dE"
            end if

        end if

        ! Sequential loop over energies (no OpenMP here)
        do i_E = 1, n_energy
            
            if (mpi_rank == 0) then
                ! Print only every 10% of iterations (at least every 10th iteration)
                if (mod(i_E-1, max(1, n_energy/10)) == 0 .or. i_E == n_energy) then
                    write(*, '(A,I5,A,I5,A,F10.3,A)') "   ===> Calculating integrand for i_E = ", i_E, " / ", n_energy, &
                                                                      " at energy E_i = ", E_arr(i_E), " eV"
                end if
            end if

            if (.not. interacting) then 

                integrand(i_E) = dN_States_per_UC_over_dE_i(E_arr(i_E), T, mu, .false.)

            else if (interacting .and. .not. Fermi_centered_k) then

                integrand(i_E) = dN_States_per_UC_over_dE_i(E_arr(i_E), T, mu, .true.)
            
            else if (interacting .and. Fermi_centered_k) then

                ! Case 1: Energy below the Fermi-window
                if (i_E .le. n_energy_1) then

                    integrand(i_E) = dN_States_per_UC_over_dE_i(E_arr(i_E), T, mu, .false.)

                ! Case 2: Energy within the Fermi-window
                else if (i_E .gt. n_energy_1 .and. i_E .le. n_energy_1 + n_energy_2) then

                    integrand(i_E) = dN_States_per_UC_over_dE_i(E_arr(i_E), T, mu, .true.)

                ! Case 3: Energy above the Fermi-window
                else if (i_E .gt. n_energy_1 + n_energy_2) then

                    integrand(i_E) = dN_States_per_UC_over_dE_i(E_arr(i_E), T, mu, .false.)

                end if
            
            end if
        end do

        ! Only master process writes results to file
        if (mpi_rank == 0 .and. debug_mode) then
            if (ios == 0) then
                do i_E = 1, n_energy
                    write(file_unit, '(F12.6,2X,ES14.6)') E_arr(i_E), integrand(i_E)
                end do
            end if

            ! Close file
            if (ios == 0) then
                close(file_unit)
                write(*,*)
                write(*,*) "- dN_{ (0) , states per UC } / dE data written to: ", trim(filename)
            end if
        end if

        N_states_per_UC = trapz(E_arr, integrand)
        deallocate(E_arr, integrand)
    end subroutine calculate_no_of_States_per_UC



    ! ================================ !
    ! Subroutines for Particle Numbers !
    ! ================================ !

    function dN_per_UC_el_over_dE(E, T, mu, close2E_F) result(dn_de)
        logical , intent(in) :: close2E_F
        real(dp), intent(in) :: E, T, mu

        complex(dp) :: denominator, denominator2, G_ret_mk, Sigma
        integer     :: i_k_local, i_k_plus_q, i_q, m, n, nu
        real(dp)    :: bose_factor, coupling_prefactor, dn_de, epsilon_mk, epsilon_kplusq, fermi, local_sum_val, &
                       numerator1, numerator2, prefactor, sum_val
        

        ! The summation extends over two different k-point grids, depending on whether
        ! we are close to the Fermi level (E = 0) or not.

        ! Note: Prefactor of 2 because of spin-degeneracy
        prefactor = -2.0 * fermi_function(E, T) / (pi * real(Nk_cube))
        local_sum_val = 0.0

        ! ------------------------------------- !
        ! Part 1: Non-Interacting Contributions !
        ! ------------------------------------- !

        call setup_mpi_kbounds("Non-Interacting")

        ! MPI parallelization over i_k_local, OpenMP over inner loops
        do i_k_local = i_k_start, i_k_end
            do m = 1, N_bnd
                epsilon_mk = epsilon_arr(m, i_k_local)
                Sigma = 0.0_dp

                ! Check if epsilon_mk is close to E_F, i.e. if | epsilon_mk - mu | < fsthick
                if (close2E_F .and. abs(epsilon_mk - mu) < fsthick ) cycle
                
                G_ret_mk = 1.0 / (cmplx(E, eta, kind=dp) - (epsilon_mk + Sigma - mu))
                local_sum_val = local_sum_val + aimag(G_ret_mk)

            end do
        end do

        ! --------------------------------- !
        ! Part 2: Interacting Contributions !
        ! --------------------------------- !

        if (close2E_F) then 
            
            call setup_mpi_kbounds("Self-Energy")


            ! MPI parallelization over i_k_local, OpenMP over inner loops
            do i_k_local = i_k_start, i_k_end
                do m = 1, N_bnd
                    epsilon_mk = SE_epsilon_arr(m, i_k_local)

                    ! Compute the Fock-contribution Sigma
                    Sigma = 0.0
                    if (interacting .and. close2E_F) then 
                        !$OMP PARALLEL DO PRIVATE(i_q, n, nu, coupling_prefactor, i_k_plus_q, &
                        !$OMP                     bose_factor, fermi, numerator1, numerator2, &
                        !$OMP                     denominator, denominator2, epsilon_kplusq) &
                        !$OMP             REDUCTION(+:Sigma) SCHEDULE(STATIC)
                        do i_q = 1, Nq_tot_SE
                            do n = 1, N_bnd
                                do nu = 1, n_branches

                                    if (SE_g_arr(n, m, nu, i_k_local - i_k_start + 1, i_q) < threshold_couplings) then
                                        cycle
                                    end if

                                    if (SE_omega_arr(nu, i_q) < threshold_omegas) then
                                        cycle
                                    end if

                                    ! Skip iteration if coupling is zero
                                    if (SE_g_arr(n, m, nu, i_k_local - i_k_start + 1, i_q) < threshold_couplings .or. &
                                                                     SE_omega_arr(nu, i_q) < threshold_omegas) then
                                        cycle
                                    end if

                                    ! Read in the coupling g_{m,n,nu}(p,q) using local index  
                                    coupling_prefactor = (SE_g_arr(m, n, nu, i_k_local - i_k_start + 1, i_q) * 1.0e-3_dp)**2 &
                                                          / real(Nq_SE_cube)

                                    epsilon_kplusq = SE_epsilon_kplusq_arr(n, i_k_local - i_k_start + 1, i_q)

                                    ! Bose-Einstein distribution
                                    bose_factor = bose_function(SE_omega_arr(nu, i_q) * 1.0e-3_dp, T)

                                    ! First term
                                    fermi = fermi_function(epsilon_kplusq - mu, T)
                                    numerator1 = bose_factor + fermi
                                    denominator = cmplx(E, eta_q, kind=dp) - (epsilon_kplusq - mu) &
                                                                         + SE_omega_arr(nu, i_q) * 1.0e-3_dp

                                    ! Second term
                                    numerator2 = bose_factor + 1.0 - fermi
                                    denominator2 = cmplx(E, eta_q, kind=dp) - (epsilon_kplusq - mu) &
                                                                          - SE_omega_arr(nu, i_q) * 1.0e-3_dp
                                    Sigma = Sigma + coupling_prefactor * (numerator1 / denominator + numerator2 / denominator2)

                                end do
                            end do
                        end do
                        !$OMP END PARALLEL DO
                    end if 

                    G_ret_mk = 1.0 / (cmplx(E, eta, kind=dp) - (epsilon_mk + Sigma - mu))
                    local_sum_val = local_sum_val + aimag(G_ret_mk)

                end do
            end do

        end if 

        ! Sum contributions from all MPI processes
        call MPI_Allreduce(local_sum_val, sum_val, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_ierr)

        dn_de = prefactor * sum_val
    end function dN_per_UC_el_over_dE


    function calculate_particle_number_per_UC(T, mu, showprogress) result(N_particles_per_UC)
        logical , intent(in)  :: showprogress
        real(dp), intent(in)  :: T, mu
        real(dp), allocatable :: integrand(:)
        
        real(dp)          :: N_particles_per_UC 
        integer           :: file_unit, ios
        character(len=50) :: filename

        ! Set up Energy-array for integration of DOOS
        call setup_Energy_arr_DOOS()

        allocate(integrand(n_energy))

        ! Only master process handles file I/O
        if (mpi_rank == 0 .and. debug_mode) then
            ! Create filename with temperature and chemical potential info
            write(filename, '(A,F6.1,A,F8.4,A)') trim(file_tree)//"GRIT_dN_per_UC_over_dE_T", T, "_mu", mu, ".txt"

            ! Open file for writing
            open(newunit=file_unit, file=filename, status='replace', action='write', iostat=ios)
            if (ios /= 0) then
                write(*,*) "Warning: Could not open file for writing: ", trim(filename)
            end if

            ! Write header to file
            if (ios == 0) then
                write(file_unit, '(A)') "# Energy (eV)    Integrand"
            end if

        end if


        ! Sequential loop over energies (no OpenMP here)
        do i_E = 1, n_energy

            if (mpi_rank == 0 .and. showprogress) then
                ! Print only every 10% of iterations (at least every 10th iteration)
                if (mod(i_E-1, max(1, n_energy/10)) == 0 .or. i_E == n_energy) then
                    write(*, '(A,I5,A,I5,A,F10.3,A)') "   ===> Calculating integrand for i_E = ", i_E, " / ", n_energy, &
                                                                        " at energy E_i = ", E_arr(i_E), " eV"
                end if
            end if

            if (.not. interacting) then 

                integrand(i_E) = dN_per_UC_el_over_dE(E_arr(i_E), T, mu, .false.)

            else if (interacting .and. .not. Fermi_centered_k) then

                integrand(i_E) = dN_per_UC_el_over_dE(E_arr(i_E), T, mu, .true.)
            
            else if (interacting .and. Fermi_centered_k) then

                ! Case 1: Energy below the Fermi-window
                if (i_E .le. n_energy_1) then

                    integrand(i_E) = dN_per_UC_el_over_dE(E_arr(i_E), T, mu, .false.)

                ! Case 2: Energy within the Fermi-window
                ! Note: Case 3 (Energy above the Fermi-window)
                !  ---> Doesn't matter cause integrand is flattened by Fermi-function.

                else if (i_E .gt. n_energy_1 .and. i_E .le. n_energy_1 + n_energy_2) then

                    integrand(i_E) = dN_per_UC_el_over_dE(E_arr(i_E), T, mu, .true.)
                
                ! Case 3: Energy above the Fermi-window
                else if (i_E .gt. n_energy_1 + n_energy_2) then

                    integrand(i_E) = dN_per_UC_el_over_dE(E_arr(i_E), T, mu, .false.)

                end if 

            end if
        end do

        ! Only master process writes results to file
        if (mpi_rank == 0 .and. debug_mode) then
            if (ios == 0) then
                do i_E = 1, n_energy
                    write(file_unit, '(F12.6,2X,ES14.6)') E_arr(i_E), integrand(i_E)
                end do
            end if

            ! Close file
            if (ios == 0) then
                close(file_unit)
                write(*,*)
                write(*,*) "- dN_{ (0) , el per UC } / dE data written to: ", trim(filename)
            end if
        end if

        N_particles_per_UC = trapz(E_arr, integrand)
        deallocate(E_arr, integrand)


    end function calculate_particle_number_per_UC



    ! ================================================ !
    ! Subroutines for Chemical Potential Determination !
    ! ================================================ !

    function zero_function(mu, T) result(zero_val)
        real(dp), intent(in) :: mu, T
        real(dp)             :: zero_val
        
        zero_val = calculate_particle_number_per_UC(T, mu, .false.) - n_UC

    end function zero_function


    function find_root_bisection(T, low, high) result(root)
        real(dp), intent(in) :: T
        real(dp), intent(in) :: low, high

        real(dp) :: root
        real(dp) :: mu_min, mu_max, mu_i
        real(dp) :: a, b, zero_a, zero_b, zero_mid
        integer            :: i_BIS
        integer, parameter :: max_iter = 100

        
        ! Set up Energy-array for integration of DOOS
        call setup_Energy_arr_DOOS()

        mu_min = low 
        mu_max = high
        a = low
        b = high

        ! Check if the chemical potential bounds are compatible with the energy array
        if (mpi_rank == 0) then
            write(*,'(A,F8.3,A,F8.3,A)') " - Chemical potential search bounds: [mu_min, mu_max] = [", mu_min, ", ", mu_max, "] eV"
            
            if (.true.) then
                write(*,*) "  Note 1: [epsilon_min - mu_min , epsilon_max - mu_min] = [", &
                            minval(epsilon_arr) - mu_min, ", ", maxval(epsilon_arr) - mu_min, "] eV"
                write(*,*) "  Note 2: [epsilon_min - mu_max , epsilon_max - mu_max] = [", &
                            minval(epsilon_arr) - mu_max, ", ", maxval(epsilon_arr) - mu_max, "] eV"
                write(*,'(A,F8.3,A,F8.3,A)') "   Note 3: E_arr = [", minval(E_arr), ", ", maxval(E_arr), "] eV"
            end if 
        end if

        zero_a = zero_function(a, T)
        zero_b = zero_function(b, T)

        ! Print zero_a and zero_b values
        if (mpi_rank == 0) then
            write(*,*)
            write(*,'(A,F8.3,A,F12.6)') "zero(mu_min = ", mu_min, " eV) = |<n_per_UC>(mu_min) - n_per_UC| = ", zero_a
            write(*,'(A,F8.3,A,F12.6)') "zero(mu_max = ", mu_max, " eV) = |<n_per_UC>(mu_max) - n_per_UC| = ", zero_b
            write(*,*)
        end if

        if (zero_a * zero_b > 0.0) then
            if (mpi_rank == 0) then
                write(*,*) "Error: Function values at bounds have same sign"
                write(*,'(A,F8.3,A,F8.6)') "zero(mu_min = ", mu_min, " eV) = ", zero_a
                write(*,'(A,F8.3,A,F8.6)') "zero(mu_max = ", mu_max, " eV) = ", zero_b
                write(*,*) "STOPPING EXECUTION - Bisection method cannot proceed"
            end if

            ! Finalize MPI before stopping
            call mpi_finalize(mpi_ierr)
            stop "Bisection method failed: function values at bounds have same sign"
        end if

        do i_BIS = 1, max_iter
            mu_i = 0.5 * (mu_min + mu_max)
            zero_mid = zero_function(mu_i, T)

            if (mpi_rank == 0) then
                write(*,'(A,I3,A,F10.6,A,F10.6)') "Iter ", i_BIS, ": mu_i = ", mu_i, &
                " eV, |zero(mu_i)| = |<n_per_UC>(mu_i) - n_per_UC| = ", abs(zero_mid)
            end if

            if (abs(zero_mid) < threshold_bisection) then
                if (mpi_rank == 0) then
                    write(*,*)
                    write(*,'(A,I3,A)') "Converged after ", i_BIS, " iterations"
                    write(*,*)
                end if
                root = mu_i
                return
            end if

            if (zero_mid < 0.0 .and. abs(mu_i - mu_min) > 1.0e-5) then
                mu_min = mu_i
            else if (zero_mid > 0.0 .and. abs(mu_i - mu_max) > 1.0e-5) then
                mu_max = mu_i
            end if
        end do

        if (mpi_rank == 0) then
            write(*,*) "Warning: Bisection method did not converge"
        end if
        root = mu_i

        deallocate(E_arr)
    end function find_root_bisection


    function calculate_mu_T(T, mu_T_guess) result(mu_T)
        real(dp), intent(in)  :: T, mu_T_guess
        real(dp) :: mu_T, N_per_UC_calc

        if (mpi_rank == 0) then
            write(*,'(A,F8.2,A)') " - Temperature: ", T, " K."
        end if

        mu_T = find_root_bisection(T, mu_T_guess - delta_mu, mu_T_guess + delta_mu)

        N_per_UC_calc = calculate_particle_number_per_UC(T, mu_T, .false.)

        if (mpi_rank == 0) then
            write(*,'(A,F18.12,A)') " - The chemical potential was computed by demanding of charge neutrality to mu_T = ", &
                              mu_T, " eV."
            write(*,'(A,F10.7,A)') "   (Compare this to the Fermi energy within the LDA-approx. at T = 0 K as E_F = ", E_F, " eV)."
            write(*,'(A,F8.3,A)') " - The Particle number per unit cell is given as N_per_UC = ", &
                   N_per_UC_calc, " (electrons)."
            write(*,*)
            write(*,*)
            write(*,*)
        end if

    end function calculate_mu_T



    ! ====================================== !
    ! Subroutines for Specific Heat Capacity !
    ! ====================================== !

    ! ============================= !
    ! Chemical Potential Derivative !
    ! ============================= !

    function partialT_F_over_dE_i(E, T, mu, close2E_F) result(partialT_F_over_dE)
        real(dp), intent(in) :: E, T, mu 
        logical , intent(in) :: close2E_F

        complex(dp) :: denominator1, denominator2, G_ret_mk, kBT2_partial_T_Sigma, Sigma
        integer     :: i_q, i_k_local, i_k_plus_q, m, n, nu
        real(dp)    :: beta, bose_factor, coupling_prefactor, epsilon_mk, epsilon_kplusq, factor, fermi, &
                       local_sum_val, numerator1, numerator2, partialT_F_over_dE, prefactor,             &
                       pTSE_numerator1, pTSE_numerator2, sum_val


        ! Note: Prefactor of 2 because of spin-degeneracy
        beta = 1.0 / (k_B * T)
        prefactor = -2.0 * beta / (pi * T * real(Nk_cube)) * fermi_function(E, T)
        factor = (1.0 - fermi_function(E, T)) * E

        sum_val = 0.0_dp
        local_sum_val = 0.0

        ! ------------------------------------- !
        ! Part 1: Non-Interacting Contributions !
        ! ------------------------------------- !

        call setup_mpi_kbounds("Non-Interacting")

        ! MPI parallelization over i_k (for each band), OpenMP over inner loops
        do i_k_local = i_k_start, i_k_end
            do m = 1, N_bnd
                epsilon_mk = epsilon_arr(m, i_k_local)

                Sigma = 0.0
                kBT2_partial_T_Sigma = 0.0
                
                ! Check if epsilon_mk is close to E_F, i.e. if | epsilon_mk - mu | < fsthick
                if (close2E_F .and. abs(epsilon_mk - mu) < fsthick ) cycle

                G_ret_mk = 1.0 / (cmplx(E, eta, kind=dp) - (epsilon_mk + Sigma - mu))
                local_sum_val = local_sum_val + (factor * aimag(G_ret_mk) - aimag(G_ret_mk**2 * kBT2_partial_T_Sigma))

            end do
        end do
        
        ! --------------------------------- !
        ! Part 2: Interacting Contributions !
        ! --------------------------------- !

        if (close2E_F) then 
            
            call setup_mpi_kbounds("Self-Energy")

            ! MPI parallelization over i_k (for each band), OpenMP over inner loops
            do i_k_local = i_k_start, i_k_end
                do m = 1, N_bnd
                    epsilon_mk = SE_epsilon_arr(m, i_k_local)

                    ! Compute the Fock-contribution Sigma
                    Sigma = 0.0
                    kBT2_partial_T_Sigma = 0.0
                    if (interacting .and. close2E_F) then 
                        !$OMP PARALLEL DO PRIVATE(n, nu, i_q, coupling_prefactor, i_k_plus_q, &
                        !$OMP                     bose_factor, fermi, numerator1, numerator2, &
                        !$OMP                     epsilon_kplusq, &
                        !$OMP                     denominator1, denominator2, &
                        !$OMP                     pTSE_numerator1, pTSE_numerator2) &
                        !$OMP             REDUCTION(+:Sigma, kBT2_partial_T_Sigma) SCHEDULE(STATIC)
                        do i_q = 1, Nq_tot_SE
                            do n = 1, N_bnd
                                do nu = 1, n_branches

                                    ! Skip iteration if coupling is zero
                                    if (SE_g_arr(n, m, nu, i_k_local - i_k_start + 1, i_q) < threshold_couplings .or. &
                                                                     SE_omega_arr(nu, i_q) < threshold_omegas) then
                                        cycle
                                    end if

                                    ! ================ !
                                    ! Self-Energy Part !
                                    ! ================ !

                                    ! Read in the coupling g_{n,m,nu}(k,q) using local index
                                    coupling_prefactor = (SE_g_arr(n, m, nu, i_k_local - i_k_start + 1, i_q) * 1.0e-3_dp)**2 &
                                                          / real(Nq_SE_cube)

                                    epsilon_kplusq = SE_epsilon_kplusq_arr(n, i_k_local - i_k_start + 1, i_q)

                                    ! Bose-Einstein distribution
                                    bose_factor = bose_function(SE_omega_arr(nu, i_q) * 1.0e-3_dp, T)

                                    ! First term
                                    fermi = fermi_function(epsilon_kplusq - mu, T)
                                    numerator1 = bose_factor + fermi
                                    denominator1 = cmplx(E, eta_q, kind=dp) - (epsilon_kplusq - mu) &
                                                                            + SE_omega_arr(nu, i_q) * 1.0e-3_dp

                                    ! Second term
                                    numerator2 = bose_factor + 1.0 - fermi
                                    denominator2 = cmplx(E, eta_q, kind=dp) - (epsilon_kplusq - mu) &
                                                                            - SE_omega_arr(nu, i_q) * 1.0e-3_dp

                                    Sigma = Sigma + coupling_prefactor * (numerator1 / denominator1 + numerator2 / denominator2)

                                    ! ============================ !
                                    ! k_B T^2 partial_T Sigma Part !
                                    ! ============================ !

                                    pTSE_numerator1 = SE_omega_arr(nu, i_q) * 1.0e-3_dp * &
                                                      csch( beta * SE_omega_arr(nu, i_q) * 1.0e-3_dp / 2.0 )**2
                                    pTSE_numerator2 = (epsilon_kplusq - mu) * &
                                                       sech( beta * (epsilon_kplusq - mu) / 2.0 )**2

                                    kBT2_partial_T_Sigma = kBT2_partial_T_Sigma &
                                                         + coupling_prefactor / 4.0 * (pTSE_numerator1 + pTSE_numerator2) &
                                                         / denominator1 &
                                                         + coupling_prefactor / 4.0 * (pTSE_numerator1 - pTSE_numerator2) &
                                                         / denominator2 

                                end do
                            end do
                        end do
                        !$OMP END PARALLEL DO
                    end if 

                    G_ret_mk = 1.0 / (cmplx(E, eta, kind=dp) - (epsilon_mk + Sigma - mu))
                    local_sum_val = local_sum_val + (factor * aimag(G_ret_mk) - aimag(G_ret_mk**2 * kBT2_partial_T_Sigma))

                end do
            end do
        
        end if 

        ! Sum contributions from all MPI processes
        call MPI_Allreduce(local_sum_val, sum_val, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_ierr)

        partialT_F_over_dE = prefactor * sum_val

    end function partialT_F_over_dE_i


    function partialmu_F_over_dE_i(E, T, mu, close2E_F) result(partialmu_F_over_dE)
        real(dp), intent(in) :: E, T, mu 
        logical , intent(in) :: close2E_F

        complex(dp) :: G_ret_mk, denominator1, denominator2, Sigma, sigma_mk
        integer     :: i_k_local, i_k_plus_q, i_q, m, n, nu
        real(dp)    :: beta, bose_factor, coupling_prefactor, epsilon_mk, epsilon_kplusq, factor1, factor2, &
                       fermi, local_sum_val, numerator1, numerator2, partialmu_F_over_dE, prefactor,        &
                       pTSE_numerator1, pTSE_numerator2, sum_val


        ! Note: Prefactor of 2 because of spin-degeneracy
        beta = 1.0 / (k_B * T)
        prefactor = -2.0 * beta / ( 4 * pi * real(Nk_cube) )
        factor1 = sech( beta * E / 2.0 )**2
        factor2 = fermi_function(E, T)

        sum_val = 0.0_dp
        local_sum_val = 0.0

        ! ------------------------------------- !
        ! Part 1: Non-Interacting Contributions !
        ! ------------------------------------- !

        call setup_mpi_kbounds("Non-Interacting")

        ! MPI parallelization over i_k (for each band), OpenMP over inner loops
        do i_k_local = i_k_start, i_k_end
            do m = 1, N_bnd
                epsilon_mk = epsilon_arr(m, i_k_local)

                if (close2E_F .and. abs(epsilon_mk - mu) < fsthick ) cycle

                ! Compute the Fock-contribution Sigma
                Sigma = 0.0
                sigma_mk = 0.0
                
                G_ret_mk = 1.0 / (cmplx(E, eta, kind=dp) - (epsilon_mk + Sigma - mu))

                local_sum_val = local_sum_val + (factor1 * aimag(G_ret_mk) - factor2 * aimag(G_ret_mk**2 * sigma_mk))
            end do
        end do


        ! --------------------------------- !
        ! Part 2: Interacting Contributions !
        ! --------------------------------- !

        if (close2E_F) then 
            
            call setup_mpi_kbounds("Self-Energy") 

            ! MPI parallelization over i_k (for each band), OpenMP over inner loops
            do i_k_local = i_k_start, i_k_end
                do m = 1, N_bnd
                    epsilon_mk = SE_epsilon_arr(m, i_k_local)

                    ! Compute the Fock-contribution Sigma
                    Sigma = 0.0
                    sigma_mk = 0.0
                    if (interacting .and. close2E_F) then 
                        !$OMP PARALLEL DO PRIVATE(n, nu, i_q, coupling_prefactor, i_k_plus_q, &
                        !$OMP                     bose_factor, fermi, numerator1, numerator2, &
                        !$OMP                     epsilon_kplusq, &
                        !$OMP                     denominator1, denominator2, &
                        !$OMP                     pTSE_numerator1, pTSE_numerator2) &
                        !$OMP             REDUCTION(+:Sigma, sigma_mk) SCHEDULE(STATIC)
                        do i_q = 1, Nq_tot_SE
                            do n = 1, N_bnd
                                do nu = 1, n_branches

                                    ! Skip iteration if coupling is zero
                                    if (SE_g_arr(n, m, nu, i_k_local - i_k_start + 1, i_q) < threshold_couplings .or. &
                                                                     SE_omega_arr(nu, i_q) < threshold_omegas) then
                                        cycle
                                    end if

                                    ! ================ !
                                    ! Self-Energy Part !
                                    ! ================ !

                                    ! Read in the coupling g_{n,m,nu}(k,q) using local index
                                    coupling_prefactor = (SE_g_arr(n, m, nu, i_k_local - i_k_start + 1, i_q) * 1.0e-3_dp)**2 &
                                                          / real(Nq_SE_cube)

                                    epsilon_kplusq = SE_epsilon_kplusq_arr(n, i_k_local - i_k_start + 1, i_q)

                                    ! Bose-Einstein distribution
                                    bose_factor = bose_function(SE_omega_arr(nu, i_q) * 1.0e-3_dp, T)

                                    ! First term
                                    fermi = fermi_function(epsilon_kplusq - mu, T)
                                    numerator1 = bose_factor + fermi
                                    denominator1 = cmplx(E, eta_q, kind=dp) - (epsilon_kplusq - mu) &
                                                                            + SE_omega_arr(nu, i_q) * 1.0e-3_dp

                                    ! Second term
                                    numerator2 = bose_factor + 1.0 - fermi
                                    denominator2 = cmplx(E, eta_q, kind=dp) - (epsilon_kplusq - mu) &
                                                                            - SE_omega_arr(nu, i_q) * 1.0e-3_dp
                                    Sigma = Sigma + coupling_prefactor * (numerator1 / denominator1 + numerator2 / denominator2)

                                    ! ============= !
                                    ! sigma_mk Part !
                                    ! ============= !

                                    sigma_mk = sigma_mk &
                                             + sech( beta * (epsilon_kplusq - mu) / 2.0 )**2 * coupling_prefactor * &
                                                   ( 1 / denominator1 &
                                                   - 1 / denominator2 ) 

                                end do
                            end do
                        end do
                        !$OMP END PARALLEL DO
                    end if 

                    G_ret_mk = 1.0 / (cmplx(E, eta, kind=dp) - (epsilon_mk + Sigma - mu))

                    local_sum_val = local_sum_val + (factor1 * aimag(G_ret_mk) - factor2 * aimag(G_ret_mk**2 * sigma_mk))
                end do
            end do
        
        end if

        ! Sum contributions from all MPI processes
        call MPI_Allreduce(local_sum_val, sum_val, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_ierr)

        partialmu_F_over_dE = prefactor * sum_val
    end function partialmu_F_over_dE_i


    function calc_dT_mu(T, mu) result(dT_mu)
        
        real(dp), intent(in)  :: T, mu
        real(dp), allocatable :: integrand1(:), integrand2(:)
        real(dp)              :: dT_mu, partialmu_F, partialT_F

        allocate(integrand1(n_energy), integrand2(n_energy))
        
        call setup_Energy_arr_SH() 

        ! Sequential loop over energies (no OpenMP here)
        do i_E = 1, n_energy

            if (.not. interacting) then

                integrand1(i_E) = partialT_F_over_dE_i(E_arr(i_E), T, mu, .false.)
                integrand2(i_E) = partialmu_F_over_dE_i(E_arr(i_E), T, mu, .false.)

            else if (interacting) then

                integrand1(i_E) = partialT_F_over_dE_i(E_arr(i_E), T, mu, .true.)
                integrand2(i_E) = partialmu_F_over_dE_i(E_arr(i_E), T, mu, .true.)

            end if
        end do

        partialT_F  = trapz(E_arr, integrand1)
        partialmu_F = trapz(E_arr, integrand2)

        deallocate(integrand1, integrand2)

        dT_mu = - partialT_F / partialmu_F
        

    end function calc_dT_mu



    ! ================================================= !
    ! Electronic Contributions (Temperature Derivative) !
    ! ================================================= !

    function partial_T_dH_el_per_UC_over_dE_i(E, T, mu, close2E_F) result(pT_dH_el_over_dE)
        real(dp), intent(in) :: E, T, mu
        logical , intent(in) :: close2E_F
        
        complex(dp) :: G_ret_mk, denominator1, denominator2, Sigma, kBT2_dT_Sigma
        integer     :: m, n, nu, i_k_plus_q, i_q, i_k_local
        real(dp)    :: beta, bose_factor, coupling_prefactor, epsilon_mk, epsilon_kplusq, factor1, fermi,  & 
                       local_sum_val, numerator1, numerator2, prefactor, pTSE_numerator1, pTSE_numerator2, &
                       pT_dH_el_over_dE, sum_val 


        ! Note: Prefactor of 2 because of spin-degeneracy
        beta = 1.0 / (k_B * T)
        prefactor = -2.0 / (k_B * pi * T ** 2 * real(Nk_cube)) * fermi_function(E, T)
        factor1 = (1.0 - fermi_function(E, T)) * E
        
        sum_val = 0.0_dp
        local_sum_val = 0.0

        ! ------------------------------------- !
        ! Part 1: Non-Interacting Contributions !
        ! ------------------------------------- !
                                    
        call setup_mpi_kbounds("Non-Interacting")

        ! MPI parallelization over i_k (for each band), OpenMP over inner loops
        do i_k_local = i_k_start, i_k_end
            do m = 1, N_bnd
                epsilon_mk = epsilon_arr(m, i_k_local)

                if (close2E_F .and. abs(epsilon_mk - mu) < fsthick ) cycle

                ! Compute the Fock-contribution Sigma
                Sigma = 0.0
                kBT2_dT_Sigma = 0.0
                
                G_ret_mk = 1.0 / (cmplx(E, eta, kind=dp) - (epsilon_mk + Sigma - mu))

                local_sum_val = local_sum_val &
                              + epsilon_mk * (factor1 * aimag(G_ret_mk) &
                                            - aimag(G_ret_mk**2 * kBT2_dT_Sigma))
            end do
        end do

        ! --------------------------------- !
        ! Part 2: Interacting Contributions !
        ! --------------------------------- !

        if (close2E_F) then

            call setup_mpi_kbounds("Self-Energy")

            ! MPI parallelization over i_k (for each band), OpenMP over inner loops
            do i_k_local = i_k_start, i_k_end
                do m = 1, N_bnd
                    epsilon_mk = SE_epsilon_arr(m, i_k_local)

                    ! Compute the Fock-contribution Sigma
                    Sigma = 0.0
                    kBT2_dT_Sigma = 0.0
                    if (interacting) then 
                        !$OMP PARALLEL DO PRIVATE(n, nu, i_q, coupling_prefactor, i_k_plus_q, &
                        !$OMP                     bose_factor, fermi, numerator1, numerator2, &
                        !$OMP                     denominator1, denominator2, &
                        !$OMP                     epsilon_kplusq, &
                        !$OMP                     pTSE_numerator1, pTSE_numerator2) &
                        !$OMP             REDUCTION(+:Sigma, kBT2_dT_Sigma) SCHEDULE(STATIC)
                        do i_q = 1, Nq_tot_SE
                            do n = 1, N_bnd
                                do nu = 1, n_branches

                                    ! Skip iteration if coupling is zero
                                    if (SE_g_arr(n, m, nu, i_k_local - i_k_start + 1, i_q) < threshold_couplings .or. &
                                                                     SE_omega_arr(nu, i_q) < threshold_omegas) then
                                        cycle
                                    end if

                                    ! ================ !
                                    ! Self-Energy Part !
                                    ! ================ !

                                    ! Read in the coupling g_{n,m,nu}(k,q) using local index
                                    coupling_prefactor = (SE_g_arr(n, m, nu, i_k_local - i_k_start + 1, i_q) * 1.0e-3_dp)**2 & 
                                                          / real(Nq_SE_cube)

                                    epsilon_kplusq = SE_epsilon_kplusq_arr(n, i_k_local - i_k_start + 1, i_q)

                                    ! Bose-Einstein distribution
                                    bose_factor = bose_function(SE_omega_arr(nu, i_q) * 1.0e-3_dp, T)

                                    ! First term
                                    fermi = fermi_function(epsilon_kplusq - mu, T)
                                    numerator1 = bose_factor + fermi
                                    denominator1 = cmplx(E, eta_q, kind=dp) - (epsilon_kplusq - mu) &
                                                                            + SE_omega_arr(nu, i_q) * 1.0e-3_dp

                                    ! Second term
                                    numerator2 = bose_factor + 1.0 - fermi
                                    denominator2 = cmplx(E, eta_q, kind=dp) - (epsilon_kplusq - mu) &
                                                                            - SE_omega_arr(nu, i_q) * 1.0e-3_dp
                                    Sigma = Sigma + coupling_prefactor * (numerator1 / denominator1 + numerator2 / denominator2)

                                    ! ============================ !
                                    ! k_B T^2 partial_T Sigma Part !
                                    ! ============================ !

                                    pTSE_numerator1 = SE_omega_arr(nu, i_q) * 1.0e-3_dp * &
                                                      csch( beta * SE_omega_arr(nu, i_q) * 1.0e-3_dp / 2.0 )**2
                                    pTSE_numerator2 = (epsilon_kplusq - mu) * &
                                                       sech( beta * (epsilon_kplusq - mu) / 2.0 )**2

                                    kBT2_dT_Sigma = kBT2_dT_Sigma &
                                                  + coupling_prefactor / 4.0 * (pTSE_numerator1 + pTSE_numerator2) / denominator1 &
                                                  + coupling_prefactor / 4.0 * (pTSE_numerator1 - pTSE_numerator2) / denominator2 

                                end do
                            end do
                        end do
                        !$OMP END PARALLEL DO
                    end if 

                    G_ret_mk = 1.0 / (cmplx(E, eta, kind=dp) - (epsilon_mk + Sigma - mu))

                    local_sum_val = local_sum_val &
                                  + epsilon_mk * (factor1 * aimag(G_ret_mk) &
                                                - aimag(G_ret_mk**2 * kBT2_dT_Sigma))
                end do
            end do

        end if 

        ! Sum contributions from all MPI processes
        call MPI_Allreduce(local_sum_val, sum_val, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_ierr)

        pT_dH_el_over_dE = prefactor * sum_val
    end function partial_T_dH_el_per_UC_over_dE_i


    function calc_partial_T_H_el_per_UC(T, mu) result(dT_H_el_per_UC)
        real(dp), intent(in) :: T, mu

        real(dp), allocatable :: integrand(:)
        real(dp)              :: dT_H_el_per_UC
        
        character(len=500)    :: filename

        integer               :: file_unit
        integer               :: ios = -1

        ! Set up Energy-array for integration of DOOS 
        call setup_Energy_arr_SH()

        allocate(integrand(n_energy))

        ! Only master process handles file I/O
        if (mpi_rank == 0) then
            ! Create filename with temperature and chemical potential info
            write(filename, '(A,ES12.5,A)') trim(file_tree)//"GRIT_partial_T_dH_el_per_UC_over_dE_i_T", T, ".txt"
            
            ! Open file for writing
            open(newunit=file_unit, file=filename, status='replace', action='write', iostat=ios)
            if (ios /= 0) then
                write(*,*) "Warning: Could not open file for writing: ", trim(filename)
            end if

            ! Write header to file
            if (ios == 0) then
                write(file_unit, '(A)') "# Energy (eV)    Integrand"
            end if

        end if


        ! Sequential loop over energies (no OpenMP here)
        do i_E = 1, n_energy

            if (.not. interacting) then
                integrand(i_E) = partial_T_dH_el_per_UC_over_dE_i(E_arr(i_E), T, mu, .false.)

            else if (interacting) then
                integrand(i_E) = partial_T_dH_el_per_UC_over_dE_i(E_arr(i_E), T, mu, .true.)
            end if 

        end do


        ! Only master process writes results to file
        if (mpi_rank == 0) then
            if (ios == 0) then
                
                do i_E = 1, n_energy
                    write(file_unit, '(F12.6,2X,ES12.5)') E_arr(i_E), integrand(i_E)
                end do

                close(file_unit)
                write(*,*)
                
                if (debug_mode) then 
                    write(*,*) "- partial_T dH_{ el per UC }/dE data written to: ", trim(filename)
                end if 
            end if

        end if

        dT_H_el_per_UC = trapz(E_arr, integrand)
        deallocate(integrand)
    end function calc_partial_T_H_el_per_UC



    ! ======================================================== !
    ! Electronic Contributions (Chemical Potential Derivative) !
    ! ======================================================== !

    function partial_mu_dH_el_per_UC_over_dE_i(E, T, mu, close2E_F) result(pmu_dH_el_over_dE)
        real(dp), intent(in) :: E, T, mu 
        logical, intent(in)  :: close2E_F

        complex(dp) :: G_ret_mk, denominator1, denominator2, Sigma, sigma_mk
        integer     :: i_k_local, i_k_plus_q, i_q, m, n, nu
        real(dp)    :: beta, bose_factor, coupling_prefactor, epsilon_mk, epsilon_kplusq, factor1, factor2, &
                       fermi, local_sum_val, numerator1, numerator2, pmu_dH_el_over_dE, prefactor,          &
                       pTSE_numerator1, pTSE_numerator2, sum_val


        ! Note: Prefactor of 2 because of spin-degeneracy
        beta = 1.0 / (k_B * T)
        prefactor = -2.0 * beta / ( 4 * pi * real(Nk_cube) )
        factor1 = sech( beta * E / 2.0 )**2
        factor2 = fermi_function(E, T)

        sum_val = 0.0_dp
        local_sum_val = 0.0_dp

        ! ------------------------------------- !
        ! Part 1: Non-Interacting Contributions !
        ! ------------------------------------- !

        call setup_mpi_kbounds("Non-Interacting")

        ! MPI parallelization over i_k (for each band), OpenMP over inner loops
        do i_k_local = i_k_start, i_k_end
            do m = 1, N_bnd
                epsilon_mk = epsilon_arr(m, i_k_local)

                if (close2E_F .and. abs(epsilon_mk - mu) < fsthick ) cycle

                ! Compute the Fock-contribution Sigma
                Sigma = 0.0
                sigma_mk = 0.0
                G_ret_mk = 1.0 / (cmplx(E, eta, kind=dp) - (epsilon_mk + Sigma - mu))

                local_sum_val = local_sum_val &
                              + epsilon_mk * (factor1 * aimag(G_ret_mk) - factor2 * aimag(G_ret_mk**2 * sigma_mk))
            end do
        end do

        ! --------------------------------- !
        ! Part 2: Interacting Contributions !
        ! --------------------------------- !

        if (close2E_F) then

            call setup_mpi_kbounds("Self-Energy")

            ! MPI parallelization over i_k (for each band), OpenMP over inner loops
            do i_k_local = i_k_start, i_k_end
                do m = 1, N_bnd
                    epsilon_mk = SE_epsilon_arr(m, i_k_local)

                    ! Compute the Fock-contribution Sigma
                    Sigma = 0.0
                    sigma_mk = 0.0
                    if (interacting) then 
                        !$OMP PARALLEL DO PRIVATE(n, nu, i_q, coupling_prefactor, i_k_plus_q, &
                        !$OMP                     bose_factor, fermi, numerator1, numerator2, &
                        !$OMP                     epsilon_kplusq, &
                        !$OMP                     denominator1, denominator2, &
                        !$OMP                     pTSE_numerator1, pTSE_numerator2) &
                        !$OMP             REDUCTION(+:Sigma, sigma_mk) SCHEDULE(STATIC)
                        do i_q = 1, Nq_tot_SE
                            do n = 1, N_bnd
                                do nu = 1, n_branches

                                    ! Skip iteration if coupling is zero
                                    if (SE_g_arr(n, m, nu, i_k_local - i_k_start + 1, i_q) < threshold_couplings .or. &
                                                                     SE_omega_arr(nu, i_q) < threshold_omegas) then
                                        cycle
                                    end if

                                    ! ================ !
                                    ! Self-Energy Part !
                                    ! ================ !

                                    ! Read in the coupling g_{n,m,nu}(k,q) using local index
                                    coupling_prefactor = (SE_g_arr(n, m, nu, i_k_local - i_k_start + 1, i_q) * 1.0e-3_dp)**2 &
                                                          / real(Nq_SE_cube)

                                    epsilon_kplusq = SE_epsilon_kplusq_arr(n, i_k_local - i_k_start + 1, i_q)

                                    ! Bose-Einstein distribution
                                    bose_factor = bose_function(SE_omega_arr(nu, i_q) * 1.0e-3_dp, T)

                                    ! First term
                                    fermi = fermi_function(epsilon_kplusq - mu, T)
                                    numerator1 = bose_factor + fermi
                                    denominator1 = cmplx(E, eta_q, kind=dp) - (epsilon_kplusq - mu) &
                                                                            + SE_omega_arr(nu, i_q) * 1.0e-3_dp

                                    ! Second term
                                    numerator2 = bose_factor + 1.0 - fermi
                                    denominator2 = cmplx(E, eta_q, kind=dp) - (epsilon_kplusq - mu) &
                                                                            - SE_omega_arr(nu, i_q) * 1.0e-3_dp
                                    Sigma = Sigma + coupling_prefactor * (numerator1 / denominator1 + numerator2 / denominator2)

                                    ! ============= !
                                    ! sigma_mk Part !
                                    ! ============= !

                                    sigma_mk = sigma_mk &
                                             + sech( beta * (epsilon_kplusq - mu) / 2.0 )**2 * coupling_prefactor * &
                                                   ( 1 / denominator1 &
                                                   - 1 / denominator2 ) 

                                end do
                            end do
                        end do
                        !$OMP END PARALLEL DO
                    end if 

                    G_ret_mk = 1.0 / (cmplx(E, eta, kind=dp) - (epsilon_mk + Sigma - mu))

                    local_sum_val = local_sum_val &
                                  + epsilon_mk * (factor1 * aimag(G_ret_mk) - factor2 * aimag(G_ret_mk**2 * sigma_mk))
                end do
            end do

        end if 

        ! Sum contributions from all MPI processes
        call MPI_Allreduce(local_sum_val, sum_val, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_ierr)

        pmu_dH_el_over_dE = prefactor * sum_val
    end function partial_mu_dH_el_per_UC_over_dE_i


    function calc_partial_mu_H_el_per_UC(T, mu) result(d_mu_H_el_per_UC)
        real(dp), intent(in) :: T, mu

        real(dp), allocatable :: integrand(:)
        real(dp)              :: d_mu_H_el_per_UC
        
        character(len=500)    :: filename

        integer               :: file_unit
        integer               :: ios = -1

        ! Set up Energy-array for integration  
        call setup_Energy_arr_SH()
        

        allocate(integrand(n_energy))

        ! Only master process handles file I/O
        if (mpi_rank == 0) then
            ! Create filename with temperature and chemical potential info
            write(filename, '(A,ES12.5,A)') trim(file_tree)//"GRIT_partial_mu_dH_el_per_UC_over_dE_i_T", T, ".txt"
            
            ! Open file for writing
            open(newunit=file_unit, file=filename, status='replace', action='write', iostat=ios)
            if (ios /= 0) then
                write(*,*) "Warning: Could not open file for writing: ", trim(filename)
            end if

            ! Write header to file
            if (ios == 0) then
                write(file_unit, '(A)') "# Energy (eV)    Integrand"
            end if

        end if

        ! Sequential loop over energies (no OpenMP here)
        do i_E = 1, n_energy

            if (.not. interacting) then 
                integrand(i_E) = partial_mu_dH_el_per_UC_over_dE_i(E_arr(i_E), T, mu, .false.)
            else if (interacting) then
                integrand(i_E) = partial_mu_dH_el_per_UC_over_dE_i(E_arr(i_E), T, mu, .true.)
            end if 

        end do

        ! Only master process writes results to file
        if (mpi_rank == 0) then
            if (ios == 0) then
                
                do i_E = 1, n_energy
                    write(file_unit, '(F12.6,2X,ES12.5)') E_arr(i_E), integrand(i_E)
                end do

                close(file_unit)
                write(*,*)
                
                if (debug_mode) then
                    write(*,*) "- partial_mu dH_{ el per UC }/dE data written to: ", trim(filename)
                end if 
            end if

        end if

        d_mu_H_el_per_UC = trapz(E_arr, integrand)
        deallocate(integrand)
    end function calc_partial_mu_H_el_per_UC



    ! ============================================================== !
    ! Non-Interacting Phononic Contribution with Gaussian Broadening !
    ! ============================================================== !

    ! function partial_T_dH_ph0_per_UC_over_dE_i(E, T) result(pT_dH_el_over_dE)
    !     real(dp), intent(in) :: E, T
    !     logical , intent(in) :: close2E_F
    !     
    !     complex(dp) :: G_ret_mk, denominator1, denominator2, Sigma, kBT2_dT_Sigma
    !     integer     :: m, n, nu, i_k_plus_q, i_q, i_k_local
    !     real(dp)    :: beta, bose_factor, coupling_prefactor, epsilon_mk, epsilon_kplusq, factor1, fermi,  & 
    !                    local_sum_val, numerator1, numerator2, prefactor, pTSE_numerator1, pTSE_numerator2, &
    !                    pT_dH_el_over_dE, sum_val 
! 
! 
    !     ! Note: Prefactor of 2 because of spin-degeneracy
    !     beta = 1.0 / (k_B * T)
    !     prefactor = -2.0 / (k_B * pi * T ** 2 * real(Nk_cube)) * fermi_function(E, T)
    !     factor1 = (1.0 - fermi_function(E, T)) * E
    !     
    !     sum_val = 0.0_dp
    !     local_sum_val = 0.0
! 
    !     do i_q = 1, Nq_tot_SE
    !         do nu = 1, n_branches
    !             bose_factor = bose_function(SE_omega_arr(nu, i_q) * 1.0e-3_dp, T)
! 
! 
    !             G_ret_mk = 1.0 / (cmplx(E, eta, kind=dp) - (epsilon_mk + Sigma - mu))
! 
    !             local_sum_val = local_sum_val &
    !                             + epsilon_mk * (factor1 * aimag(G_ret_mk) &
    !                                         - aimag(G_ret_mk**2 * kBT2_dT_Sigma))
    !         end do
    !     end do
! 
! 
    !     ! Sum contributions from all MPI processes
    !     call MPI_Allreduce(local_sum_val, sum_val, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_ierr)
! 
    !     pT_dH_el_over_dE = prefactor * sum_val
    ! end function partial_T_dH_ph0_per_UC_over_dE_i


    function calc_partial_T_H_ph0_per_UC(T) result(dT_H_ph0_per_UC)
        real(dp), intent(in) :: T

        real(dp), allocatable :: integrand(:)
        real(dp)              :: dT_H_ph0_per_UC
        
        character(len=500)    :: filename

        integer               :: file_unit
        integer               :: ios = -1

        ! Set up Energy-array for integration of DOOS 
        call setup_Energy_arr_SH()

        allocate(integrand(n_energy))

        ! Only master process handles file I/O
        if (mpi_rank == 0) then
            ! Create filename with temperature and chemical potential info
            write(filename, '(A,ES12.5,A)') trim(file_tree)//"GRIT_partial_T_dH_ph0_per_UC_over_dE_i_T", T, ".txt"
            
            ! Open file for writing
            open(newunit=file_unit, file=filename, status='replace', action='write', iostat=ios)
            if (ios /= 0) then
                write(*,*) "Warning: Could not open file for writing: ", trim(filename)
            end if

            ! Write header to file
            if (ios == 0) then
                write(file_unit, '(A)') "# Energy (eV)    Integrand"
            end if

        end if


        ! Sequential loop over energies (no OpenMP here)
        do i_E = 1, n_energy

            ! integrand(i_E) = partial_T_dH_ph0_per_UC_over_dE_i(E_arr(i_E), T)

        end do


        ! Only master process writes results to file
        if (mpi_rank == 0) then
            if (ios == 0) then
                
                do i_E = 1, n_energy
                    write(file_unit, '(F12.6,2X,ES12.5)') E_arr(i_E), integrand(i_E)
                end do

                close(file_unit)
                write(*,*)
                
                if (debug_mode) then 
                    write(*,*) "- partial_T dH_{ ph, 0 per UC }/dE data written to: ", trim(filename)
                end if 
            end if

        end if

        dT_H_ph0_per_UC = trapz(E_arr, integrand)
        deallocate(integrand)
    end function calc_partial_T_H_ph0_per_UC

 


    ! ====================== !
    ! Phononic Contributions !
    ! ====================== !

    function calc_partial_T_H_ph_0_per_UC(T) result(pT_H_ph_0_per_UC)
        real(dp), intent(in) :: T
        real(dp) :: pT_H_ph_0_per_UC, prefactor
        integer  :: nu, i_q

        pT_H_ph_0_per_UC = 0.0_dp
        prefactor = 1.0_dp / (4.0_dp * k_B * T**2 * real(Nq_cube, dp))

        ! TBD Parallelizing this. or only one MPI rank has to perform this calculation since no coupling data is needed

        ! No parallelization needed here since omega_arr is small
        do i_q = 1, Nq_tot_0
            do nu = 1, n_branches
                                
                if (omega_arr(nu, i_q) < threshold_omegas) then
                    cycle
                end if

                pT_H_ph_0_per_UC = pT_H_ph_0_per_UC + &
                            (omega_arr(nu, i_q) * 1.0e-3_dp * csch(0.5_dp * omega_arr(nu, i_q) * 1.0e-3_dp / (k_B * T)))**2

            end do
        end do

        pT_H_ph_0_per_UC = prefactor * pT_H_ph_0_per_UC
    end function calc_partial_T_H_ph_0_per_UC


    function calc_partial_T_H_ph_corr(T, mu) result(pT_H_ph_corr)
        real(dp), intent(in) :: T, mu

        integer  :: m, n, nu, i_k_local, i_k_plus_q, i_q
        real(dp) :: beta, denominator_factor, epsilon_kplusq, Fermi_nk, Fermi_mkpq, global_prefactor,     &
                    interior_loop_global_prefactor, local_pT_H_ph_corr, pT_H_ph_corr, pT_H_ph_corr_q,     &
                    term1, term1_prefactor, term1_numerator1_1, term1_numerator1_2, term1_denominator1,   &
                    term1_numerator2, term1_denominator2,                                                 &
                    term2, term2_prefactor, term2_numerator1, term2_numerator2, term2_denominator,        &
                    term3, term3_prefactor, term3_numerator, term3_denominator, term3_parentheses_factor, &
                    xi_diff, xi_nk, xi_mkpq


        call setup_mpi_kbounds("Uniform")
        
        beta = 1.0 / (k_B * T)
        global_prefactor = -2.0_dp / ( real(Ni_uniform_cube, dp) * 4.0_dp * k_B * T**2 ) ! Factor of 2 for spin-degeneracy
        pT_H_ph_corr = 0.0
        local_pT_H_ph_corr = 0.0 ! Initialize local accumulator

        ! MPI parallelization over i_k (for each band), OpenMP over inner loops
        do i_k_local = i_k_start, i_k_end
            do n = 1, N_bnd
                xi_nk = uniform_epsilon_arr(n, i_k_local) - mu

                Fermi_nk = fermi_function(xi_nk, T)
                
                do m = 1, N_bnd
                    ! Compute the q-contribution to the sum
                    pT_H_ph_corr_q = 0.0

                    ! Parallelize the sum over i_q (and nu) with OpenMP.
                    ! pT_H_ph_corr_q is a reduction variable; temporaries are PRIVATE.
                    !$OMP PARALLEL DO DEFAULT(shared) & 
                    !$OMP    PRIVATE(nu, i_q, interior_loop_global_prefactor, i_k_plus_q, xi_mkpq, Fermi_mkpq, &
                    !$OMP            xi_diff, denominator_factor, term1_prefactor, term1_numerator1_1, &
                    !$OMP            epsilon_kplusq, &
                    !$OMP            term1_numerator1_2, term1_denominator1, &
                    !$OMP            term1_numerator2, term1_denominator2, term1, &
                    !$OMP            term2_prefactor, term2_numerator1, term2_numerator2, term2_denominator, term2, &
                    !$OMP            term3_prefactor, term3_numerator, term3_denominator, term3_parentheses_factor, term3) &
                    !$OMP    REDUCTION(+:pT_H_ph_corr_q) SCHEDULE(STATIC)

                    do i_q = 1, Ni_uniform_cube
                        do nu = 1, n_branches
                            ! Skip iteration if coupling is zero
                            if (uniform_omega_arr(nu, i_q) < threshold_omegas .or. &
                                uniform_g_arr(m, n, nu, i_k_local - i_k_start + 1, i_q) < threshold_couplings) then
                                cycle
                            end if

                            ! Read in the coupling g_{n,m,nu}(k,q) using local index
                            interior_loop_global_prefactor &
                            = &
                            uniform_omega_arr(nu, i_q) * 1.0e-3_dp &
                            * (uniform_g_arr(m, n, nu, i_k_local - i_k_start + 1, i_q) * 1.0e-3_dp)**2 &
                            / real(Ni_uniform_cube)

                            i_k_plus_q = get_i_kplusq(i_k_local, i_q)
                            epsilon_kplusq = uniform_epsilon_arr(m, i_k_plus_q)

                            xi_mkpq = epsilon_kplusq - mu
                            Fermi_mkpq = fermi_function(xi_mkpq, T)

                            xi_diff = xi_nk - xi_mkpq

                            if (abs(xi_diff) < threshold_xis) then
                                cycle
                            end if

                            denominator_factor = uniform_omega_arr(nu, i_q) * 1.0e-3_dp + xi_diff

                            if (abs(denominator_factor) < threshold_denominator) then
                                cycle
                            end if


                            ! ========== !
                            ! First term !
                            ! ========== !

                            term1_prefactor = xi_nk * sech(0.5 * xi_nk / (k_B * T))**2 - &
                                              xi_mkpq * sech(0.5 * xi_mkpq / (k_B * T))**2

                            ! First term, first term inside the parentheses
                            term1_numerator1_1 = bose_function(uniform_omega_arr(nu, i_q) * 1.0e-3_dp, T) 
                            term1_numerator1_2 = bose_function(- xi_diff, T)
                            term1_denominator1 = (uniform_omega_arr(nu, i_q) * 1.0e-3_dp + xi_diff) ** 2

                            ! First term, second term inside the parentheses
                            term1_numerator2 = beta * (csch(0.5 * beta * uniform_omega_arr(nu, i_q) * 1.0e-3_dp))**2
                            term1_denominator2 = 4 * (uniform_omega_arr(nu, i_q) * 1.0e-3_dp + xi_diff)

                            term1 = term1_prefactor * (term1_numerator1_1 - term1_numerator1_2) / term1_denominator1 &
                                  + term1_prefactor * term1_numerator2 / term1_denominator2


                            
                            ! =========== !
                            ! Second term !
                            ! =========== !
                            
                            term2_prefactor = Fermi_nk - Fermi_mkpq
                            term2_numerator1 &
                            = &
                            uniform_omega_arr(nu, i_q) * 1.0e-3_dp * (csch(0.5 * uniform_omega_arr(nu, i_q) &
                                                       * 1.0e-3_dp / (k_B * T)))**2
                            
                            term2_numerator2 = (xi_diff) * (csch(- 0.5 * xi_diff / (k_B * T)))**2
                            term2_denominator = (uniform_omega_arr(nu, i_q) * 1.0e-3_dp + xi_diff) ** 2

                            term2 = term2_prefactor * (term2_numerator1 + term2_numerator2) / term2_denominator



                            ! ========== !
                            ! Third term !
                            ! ========== !

                            term3_prefactor = Fermi_nk - Fermi_mkpq
                            term3_numerator = (csch(0.5 * uniform_omega_arr(nu, i_q) * 1.0e-3_dp / (k_B * T)))**2
                            term3_denominator = (uniform_omega_arr(nu, i_q) * 1.0e-3_dp + xi_diff)
                            term3_parentheses_factor = beta * uniform_omega_arr(nu, i_q) * 1.0e-3_dp * &
                                                       coth(0.5 * beta * uniform_omega_arr(nu, i_q) * 1.0e-3_dp) - 1.0

                            term3 = term3_prefactor * term3_numerator / term3_denominator * term3_parentheses_factor

                            ! For debugging
                            if (.not. ieee_is_finite(term1)) then
                                write(*,*) "Error: term1 is not finite!"
                                term1 = 0.0
                            end if
                            if (.not. ieee_is_finite(term2)) then
                                write(*,*) "Error: term2 is not finite!"
                                term2 = 0.0
                            end if
                            if (.not. ieee_is_finite(term3)) then
                                write(*,*) "Error: term3 is not finite!"
                                term3 = 0.0
                            end if




                            ! ======================== !
                            ! Collecting Contributions !
                            ! ======================== !

                            pT_H_ph_corr_q = pT_H_ph_corr_q + interior_loop_global_prefactor * (term1 + term2 + term3)
                            
                        end do
                    end do
                    !$OMP END PARALLEL DO
                     
                    local_pT_H_ph_corr = local_pT_H_ph_corr + pT_H_ph_corr_q

                end do
            end do
        end do

        ! Sum contributions from all MPI processes
        call MPI_Allreduce(local_pT_H_ph_corr, pT_H_ph_corr, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_ierr)

        pT_H_ph_corr = global_prefactor * pT_H_ph_corr

    end function calc_partial_T_H_ph_corr


    function calc_partial_T_H_ph_corr2(T, mu) result(pT_H_ph_corr2)
        real(dp), intent(in) :: T, mu

        integer  :: n, nu, i_k_local
        real(dp) :: A_nu, B_nu, pT_B_nu, Fermi_nk, local_prefactor1, local_prefactor2, local_B_nu, local_pT_B_nu, &
                    xi_nk, sech2_nk, interior_loop_prefactor, pT_H_ph_corr2


        call setup_mpi_kbounds("Uniform")
        
        pT_H_ph_corr2 = 0.0

        local_prefactor1 = + 2.0_dp / ( real(Ni_uniform_cube, dp) ) ! Factor of 2 for spin-degeneracy
        local_prefactor2 = + 2.0_dp / ( real(Ni_uniform_cube, dp) * 4.0_dp * k_B * T**2 ) ! Factor of 2 for spin-degeneracy

        do nu = 1, n_branches

            if (uniform_omega_arr(nu, 1) < threshold_omegas) then
                cycle
            else
                A_nu = 1 / (uniform_omega_arr(nu, 1) * 1.0e-3_dp)
            end if

            B_nu = 0.0
            pT_B_nu = 0.0

            local_B_nu = 0.0
            local_pT_B_nu = 0.0
            
            ! MPI parallelization over i_k (for each band), OpenMP over inner loops
            do i_k_local = i_k_start, i_k_end
                do n = 1, N_bnd

                    xi_nk = uniform_epsilon_arr(n, i_k_local) - mu

                    Fermi_nk = fermi_function(xi_nk, T)

                    sech2_nk = sech(0.5 * xi_nk / (k_B * T))**2

                    ! Skip iteration if coupling is zero
                    if (uniform_omega_arr(nu, 1) < threshold_omegas .or. &
                        uniform_g_arr(n, n, nu, i_k_local - i_k_start + 1, 1) < threshold_couplings) then
                        cycle
                    end if

                    ! Read in the coupling g_{n,m,nu}(k,1) using local index
                    interior_loop_prefactor &
                    = &
                    uniform_g_arr(n, n, nu, i_k_local - i_k_start + 1, 1) * 1.0e-3_dp

                    local_B_nu &
                    = &
                    local_B_nu + local_prefactor1 * interior_loop_prefactor * Fermi_nk

                    local_pT_B_nu &
                    = &
                    local_pT_B_nu + local_prefactor2 * interior_loop_prefactor * xi_nk * sech2_nk

                end do
            end do

            call MPI_Allreduce(local_B_nu, B_nu, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_ierr)
            call MPI_Allreduce(local_pT_B_nu, pT_B_nu, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_ierr)

            pT_H_ph_corr2 = pT_H_ph_corr2 + A_nu * B_nu + pT_B_nu

        end do

    end function calc_partial_T_H_ph_corr2


    function calc_partial_mu_H_ph_corr(T, mu) result(pmu_H_ph_corr)
        real(dp), intent(in) :: T, mu

        integer  :: m, n, nu, i_k_local, i_k_plus_q, i_q
        real(dp) :: beta, denominator_factor, epsilon_kplusq, global_prefactor,                           &
                    interior_loop_global_prefactor, local_pmu_H_ph_corr, pmu_H_ph_corr, pmu_H_ph_corr_q,  &
                    term1, term1_prefactor, term1_numerator1_1, term1_numerator1_2, term1_denominator1,   &
                    term1_numerator2, term1_denominator2,                                                 &
                    term2, term2_prefactor, term2_numerator1, term2_numerator2, term2_denominator,        &
                    term3, term3_prefactor, term3_numerator, term3_denominator, term3_parentheses_factor, &
                    xi_diff, xi_nk, xi_mkpq, sech2_nk, sech2_mkpq

        call setup_mpi_kbounds("Uniform")
        
        beta = 1.0 / (k_B * T)
        global_prefactor = - 2.0 / (Ni_uniform_cube * 4.0 * k_B * T) ! Factor of 2 for spin-degeneracy
        pmu_H_ph_corr = 0.0
        local_pmu_H_ph_corr = 0.0 ! Initialize local accumulator

        ! MPI parallelization over i_k (for each band), OpenMP over inner loops
        do i_k_local = i_k_start, i_k_end
            do n = 1, N_bnd
                xi_nk = uniform_epsilon_arr(n, i_k_local) - mu
                sech2_nk = sech(0.5 * xi_nk / (k_B * T))**2
                
                do m = 1, N_bnd
                    ! Compute the q-contribution to the sum
                    pmu_H_ph_corr_q = 0.0

                    ! Parallelize the sum over i_q (and nu) with OpenMP.
                    ! pmu_H_ph_corr_q is a reduction variable; temporaries are PRIVATE.
                    !$OMP PARALLEL DO DEFAULT(shared) & 
                    !$OMP    PRIVATE(nu, i_q, interior_loop_global_prefactor, i_k_plus_q, xi_mkpq, sech2_mkpq, &
                    !$OMP            xi_diff, denominator_factor, term1_prefactor, term1_numerator1_1, &
                    !$OMP            epsilon_kplusq, &
                    !$OMP            term1_numerator1_2, term1_denominator1, &
                    !$OMP            term1_numerator2, term1_denominator2, term1, &
                    !$OMP            term2_prefactor, term2_numerator1, term2_numerator2, term2_denominator, term2, &
                    !$OMP            term3_prefactor, term3_numerator, term3_denominator, term3_parentheses_factor, term3) &
                    !$OMP    REDUCTION(+:pmu_H_ph_corr_q) SCHEDULE(STATIC)
                    do i_q = 1, Ni_uniform_cube
                        do nu = 1, n_branches
                            ! Skip iteration if coupling is zero
                            if (uniform_omega_arr(nu, i_q) < threshold_omegas .or. &
                                uniform_g_arr(m, n, nu, i_k_local - i_k_start + 1, i_q) < threshold_couplings) then
                                cycle
                            end if

                            ! Read in the coupling g_{n,m,nu}(k,q) using local index
                            interior_loop_global_prefactor &
                            = &
                            uniform_omega_arr(nu, i_q) * 1.0e-3_dp &
                            * (uniform_g_arr(m, n, nu, i_k_local - i_k_start + 1, i_q) * 1.0e-3_dp)**2 &
                            / real(Ni_uniform_cube)
                            
                            i_k_plus_q = get_i_kplusq(i_k_local, i_q)
                            epsilon_kplusq = uniform_epsilon_arr(m, i_k_plus_q)

                            xi_mkpq = epsilon_kplusq - mu
                            sech2_mkpq = sech(0.5 * xi_mkpq / (k_B * T))**2

                            xi_diff = xi_nk - xi_mkpq

                            if (abs(xi_diff) < threshold_xis) then
                                cycle
                            end if

                            denominator_factor = uniform_omega_arr(nu, i_q) * 1.0e-3_dp + xi_diff

                            if (abs(denominator_factor) < threshold_denominator) then
                                cycle
                            end if



                            ! ========== !
                            ! First term !
                            ! ========== !

                            term1_prefactor = sech2_nk - sech2_mkpq

                            ! First term, first term inside the parentheses
                            term1_numerator1_1 = bose_function(uniform_omega_arr(nu, i_q) * 1.0e-3_dp, T) 
                            term1_numerator1_2 = bose_function(- xi_diff, T)
                            term1_denominator1 = (uniform_omega_arr(nu, i_q) * 1.0e-3_dp + xi_diff) ** 2

                            ! First term, second term inside the parentheses
                            term1_numerator2   = beta * (csch(0.5 * beta * uniform_omega_arr(nu, i_q) * 1.0e-3_dp))**2
                            term1_denominator2 = 4 * (uniform_omega_arr(nu, i_q) * 1.0e-3_dp + xi_diff)
                            term1 = term1_prefactor * (term1_numerator1_1 - term1_numerator1_2) / term1_denominator1 &
                                  + term1_prefactor * term1_numerator2 / term1_denominator2
                            

                            ! For debugging
                            if (.not. ieee_is_finite(term1)) then
                                write(*,*) "Error: term1 is not finite!"
                                term1 = 0.0
                            end if
                            if (.not. ieee_is_finite(term2)) then
                                write(*,*) "Error: term2 is not finite!"
                                term2 = 0.0
                            end if
                            if (.not. ieee_is_finite(term3)) then
                                write(*,*) "Error: term3 is not finite!"
                                term3 = 0.0
                            end if



                            ! ======================== !
                            ! Collecting Contributions !
                            ! ======================== !

                            pmu_H_ph_corr_q = pmu_H_ph_corr_q + interior_loop_global_prefactor * term1
                            
                        end do
                    end do
                    !$OMP END PARALLEL DO
                     
                    local_pmu_H_ph_corr = local_pmu_H_ph_corr + pmu_H_ph_corr_q
                end do
            end do
        end do

        ! Sum contributions from all MPI processes
        call MPI_Allreduce(local_pmu_H_ph_corr, pmu_H_ph_corr, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_ierr)

        pmu_H_ph_corr = global_prefactor * pmu_H_ph_corr

    end function calc_partial_mu_H_ph_corr


    function calc_partial_mu_H_ph_corr2(T, mu) result(pmu_H_ph_corr2)
        real(dp), intent(in) :: T, mu

        integer  :: n, nu, i_k_local
        real(dp) :: A_nu, B_nu, pmu_B_nu, Fermi_nk, local_prefactor1, local_prefactor2, local_B_nu, local_pmu_B_nu, &
                    xi_nk, sech2_nk, interior_loop_prefactor, pmu_H_ph_corr2


        call setup_mpi_kbounds("Uniform")
        
        pmu_H_ph_corr2 = 0.0

        local_prefactor1 = + 2.0_dp / ( real(Ni_uniform_cube, dp) ) ! Factor of 2 for spin-degeneracy
        local_prefactor2 = + 2.0_dp / ( real(Ni_uniform_cube, dp) * 4.0_dp * k_B * T ) ! Factor of 2 for spin-degeneracy

        do nu = 1, n_branches

            if (uniform_omega_arr(nu, 1) < threshold_omegas) then
                cycle
            else
                A_nu = 1 / (uniform_omega_arr(nu, 1) * 1.0e-3_dp)
            end if

            B_nu = 0.0
            pmu_B_nu = 0.0

            local_B_nu = 0.0
            local_pmu_B_nu = 0.0
            
            ! MPI parallelization over i_k (for each band), OpenMP over inner loops
            do i_k_local = i_k_start, i_k_end
                do n = 1, N_bnd

                    xi_nk = uniform_epsilon_arr(n, i_k_local) - mu

                    Fermi_nk = fermi_function(xi_nk, T)

                    sech2_nk = sech(0.5 * xi_nk / (k_B * T))**2

                    ! Skip iteration if coupling is zero
                    if (uniform_omega_arr(nu, 1) < threshold_omegas .or. &
                        uniform_g_arr(n, n, nu, i_k_local - i_k_start + 1, 1) < threshold_couplings) then
                        cycle
                    end if

                    ! Read in the coupling g_{n,m,nu}(k,1) using local index
                    interior_loop_prefactor &
                    = &
                    uniform_g_arr(n, n, nu, i_k_local - i_k_start + 1, 1) * 1.0e-3_dp

                    local_B_nu &
                    = &
                    local_B_nu + local_prefactor1 * interior_loop_prefactor * Fermi_nk

                    local_pmu_B_nu &
                    = &
                    local_pmu_B_nu + local_prefactor2 * interior_loop_prefactor * sech2_nk

                end do
            end do

            call MPI_Allreduce(local_B_nu, B_nu, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_ierr)
            call MPI_Allreduce(local_pmu_B_nu, pmu_B_nu, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_ierr)

            pmu_H_ph_corr2 = pmu_H_ph_corr2 + A_nu * B_nu + pmu_B_nu

        end do

    end function calc_partial_mu_H_ph_corr2


    function calculate_c_V_per_uc(T, mu_T) result(c_V_molar_T)
        real(dp), intent(in) :: T, mu_T

        real(dp) :: dT_mu_T,                                            &
                    partial_mu_H_el_per_UC,                             &
                    partial_T_H_el_per_UC     ,  c_V_molar_el,          &
                    partial_T_H_ph_0_per_UC   ,  c_V_molar_ph_0,        &
                    partial_T_H_ph_per_UC_corr,  c_V_molar_ph_pT_corr,  &
                    partial_mu_H_ph_per_UC_corr, c_V_molar_ph_pmu_corr, &
                                                 c_V_molar_ph,          &
                                                 c_V_molar_int,         &
                                                 c_V_molar_T


        ! ========= !
        ! Electrons !
        ! ========= !

        dT_mu_T                = calc_dT_mu(T, mu_T) 
        partial_mu_H_el_per_UC = calc_partial_mu_H_el_per_UC(T, mu_T)
        partial_T_H_el_per_UC  = calc_partial_T_H_el_per_UC(T, mu_T)

        c_V_molar_el = ((partial_T_H_el_per_UC + partial_mu_H_el_per_UC * dT_mu_T) / f_factor) * eV_to_J



        ! ======= !
        ! Phonons !
        ! ======= !

        partial_T_H_ph_0_per_UC = calc_partial_T_H_ph_0_per_UC(T)

        if (.not. interacting) then 
        
            partial_T_H_ph_per_UC_corr  = 0.0
            partial_mu_H_ph_per_UC_corr = 0.0

        else if (interacting) then

            partial_T_H_ph_per_UC_corr  = calc_partial_T_H_ph_corr(T, mu_T)
            ! partial_T_H_ph_per_UC_corr = partial_T_H_ph_per_UC_corr + calc_partial_T_H_ph_corr2(T, mu_T) ! Adding the second correction term

            partial_mu_H_ph_per_UC_corr = calc_partial_mu_H_ph_corr(T, mu_T)
            ! partial_mu_H_ph_per_UC_corr = partial_mu_H_ph_per_UC_corr + calc_partial_mu_H_ph_corr2(T, mu_T) ! Adding the second correction term

        end if

        c_V_molar_ph_0        = (partial_T_H_ph_0_per_UC               / f_factor) * eV_to_J
        c_V_molar_ph_pT_corr  = (partial_T_H_ph_per_UC_corr            / f_factor) * eV_to_J
        c_V_molar_ph_pmu_corr = (partial_mu_H_ph_per_UC_corr * dT_mu_T / f_factor) * eV_to_J
        c_V_molar_ph          = c_V_molar_ph_0 + c_V_molar_ph_pT_corr + c_V_molar_ph_pmu_corr



        ! ============ !
        ! Interactions !
        ! ============ !

        c_V_molar_int = 0.0 ! Is only non-zero for systems with optical branches.



        ! ======================== !
        ! Summing up Contributions !
        ! ======================== !
        c_V_molar_T = c_V_molar_el + c_V_molar_ph + c_V_molar_int

        if (mpi_rank == 0) then
            write(*,*)
            write(*,'(A,F8.2,A)') " - Temperature: ", T, " K"
            write(*,'(A)') repeat('=', 115)
            write(*,'(A)') "Electrons:"
            write(*,'(A,F8.2,A)') " - The 1st electronic contribution  is                     \partial_T H_el      / n_moles = ", &
                     (partial_T_H_el_per_UC / f_factor * eV_to_J * 1.0e3), " mJ / K mol"
            write(*,'(A,F8.2,A)') " - The 2nd electronic contribution  is            dT mu * \partial_mu H_el      / n_moles = ", &
                     (partial_mu_H_el_per_UC * dT_mu_T / f_factor * eV_to_J * 1.0e3), " mJ / K mol"
            write(*,'(A,F8.2,A)') " - The tot electronic contribution  is c_V_molar_el     =        d/dT H_el      / n_moles = ", &
                     c_V_molar_el* 1.0e3, " mJ / K mol"
            write(*,*)
            write(*,'(A)') "Phonons:"
            write(*,'(A,F8.2,A)') " - The free phononic contribution   is c_V_molar_ph_0   = \partial_T  H_ph_0    / n_moles = ", &
                   c_V_molar_ph_0* 1.0e3, " mJ / K mol" 
            write(*,'(A,F8.2,A)') " - The 1st  phononic correction     is                    \partial_T  H_ph_corr / n_moles = ", &
                c_V_molar_ph_pT_corr* 1.0e3, " mJ / K mol"
            write(*,'(A,F8.2,A)') " - The 2nd  phononic correction     is            dT mu * \partial_mu H_ph_corr / n_moles = ", &
                c_V_molar_ph_pmu_corr* 1.0e3, " mJ / K mol"
            write(*,'(A,F8.2,A)') " - The tot  phononic contribution   is c_V_molar_ph     =        d/dT H_ph      / n_moles = ", &
                     c_V_molar_ph* 1.0e3, " mJ / K mol"
            write(*,*)
            write(*,'(A)') "Interactions:"
            write(*,'(A,F8.2,A)') " - The interaction    contribution  is c_V_molar_int     =       d/dT H_int     / n_moles = ", &
                    c_V_molar_int* 1.0e3, " mJ / K mol"
            write(*,*)
            write(*,'(A)') "Total Specific Heat Capacity:"
            write(*,'(A,F8.2,A)') " - The total specific heat capacity is c_V_molar_T       =       d/dT H_tot     / n_moles = ", &
                      c_V_molar_T* 1.0e3, " mJ / K mol"
            
            if (calc_gamma) then
                write(*,*)
                write(*,*) "Remarks: "
                write(*,'(A,F8.2,A)') " - One finds ... c_V_molar_el  / T   = (d/dT H_el) / (T   * n_moles) = ", &
                      (c_V_molar_el / T) * 1.0e3, " mJ / K^2 mol"
                write(*,'(A,F8.2,A)') " - One finds ... c_V_molar_ph0 / T^3 = (d/dT H_ph0) / (T^3 * n_moles) = ", &
                      (c_V_molar_ph_0 / T ** 3) * 1.0e3, " mJ / K^4 mol"
            end if 
            
            write(*,'(A)') repeat('=', 115)
            write(*,*)
        end if
    end function calculate_c_V_per_uc


    subroutine linear_fit(x, y, n, a, b)
        integer,  intent(in) :: n
        real(dp), intent(in) :: x(n), y(n)
        real(dp), intent(out):: a, b 

        real(dp) :: sx, sy, sxx, sxy
        
        integer :: i

        ! Find the linear fit y = a*x + b to the data (x,y) with n points

        if (mpi_rank /= 0) return

        sx  = 0.0
        sy  = 0.0
        sxx = 0.0
        sxy = 0.0   
        
        do i = 1, n
           sx  = sx  + x(i)
           sy  = sy  + y(i)
           sxx = sxx + x(i)*x(i)
           sxy = sxy + x(i)*y(i)
        end do  
        
        a = (n*sxy - sx*sy) / (n*sxx - sx*sx)
        b = (sy - a*sx) / n

    end subroutine linear_fit

end module utils
