module modmain
    use iso_fortran_env, only : real64

    implicit none

    ! ================ !
    ! Input-Parameters !
    ! ================ !
    
    ! --------- !
    ! Structure !
    ! --------- !
    
    ! character(len=20) :: structure_name = "Ag_Fm3m"
    ! character(len=20) :: structure_name = "Cs_Im3m"
    character(len=20) :: structure_name = "Pb_Fm3m"

    ! character(len=20) :: structure_name = "Cu_Fm3m"
    ! character(len=20) :: structure_name = "Au_Fm3m"
    ! character(len=20) :: structure_name = "K_Im3m"
    ! character(len=20) :: structure_name = "Rb_Im3m"

    ! character(len=20) :: structure_name = "Diamond_Fd3m"
    ! character(len=20) :: structure_name = "Si_Fd3m1"
    
    ! ----- !
    ! Modes !
    ! ----- !
    logical, parameter :: debug_mode        = .false.   ! Enable or disable debug mode
    logical, parameter :: interacting       = .false.   ! Turn electron-phonon interaction on / off - Requires Fermi- and Gamma-centered grids
    logical, parameter :: Fermi_centered_k  = .false.   ! Use fine k-grid around the Fermi-surface for better accuracy
    logical, parameter :: Gamma_centered_q  = .false.   ! Use a very fine grid around the Gamma-point.

    ! ----- !
    ! Tasks !
    ! ----- !
    logical, parameter :: create_binary                            = .true.   ! Create binary files from text files (only needed once).
    logical, parameter :: symmetrize                               = .true.   ! Symmetrize the uniform couplings "by hand" to fix small numerical asymmetries (only needed if create_binary = .true.)
    logical, parameter :: calc_N_states_per_UC_mu_eq_E_F_T_eq_zero = .true.   ! Calculate number of states from mu = E_F at T = 0 K
    logical, parameter :: calc_N_per_UC_mu_eq_E_F_T_eq_zero        = .true.   ! Calculate particle number per UC from mu = E_F at T = 0 K
    logical, parameter :: calc_mu_vs_T                             = .true.   ! Calculate chemical potential vs Temperature
    logical, parameter :: calc_c_V_per_UC                          = .true.   ! Calculate specific heat capacity vs Temperature
    logical, parameter :: calc_gamma                               = .false.  ! Calculate electronic coefficient of specific heat capacity
    
    ! -------------- !
    ! Precision kind !
    ! -------------- !
    integer, parameter  :: dp = real64
    
    ! ------------------------------------ !
    ! Temperature parameters and variables !
    ! ------------------------------------ !
    integer  :: i_T             ! Index for temperature
    integer  :: n_temp
    real(dp) :: T_min
    real(dp) :: T_max
    


    ! ================================ !
    ! Threshold parameters within GRIT !
    ! ================================ !

    real(dp), parameter :: threshold_bisection   = 1.0e-4_dp
    real(dp), parameter :: threshold_omegas      = 1.0e-10_dp ! (We want to avoid omega = 0 in our calculations.)
    real(dp), parameter :: threshold_couplings   = 1.0e-5_dp  ! (We want to avoid     g = 0 in our calculations.)
    real(dp), parameter :: threshold_denominator = 1.0e-10_dp
    real(dp), parameter :: threshold_xis         = 1.0e-10_dp



    ! ================================== !
    ! Convergence parameters within GRIT !
    ! ================================== !

    real(dp)            :: eta
    real(dp), parameter :: eta_q       = 1.0e-3_dp    ! Broadening for retarded Green function in eV
    real(dp), parameter :: E_Length_SH = 5.0e-3_dp    ! Energy Integration Range for [S]pecific [H]eat integration in eV
    real(dp), parameter :: E_tails     = 7.5_dp       ! Energy tails for integration in eV
    real(dp), parameter :: delta_mu    = 3.5_dp       ! Value around E_F to search for mu at finite T in eV
                                                      ! Notice that E_tails should be larger than delta_mu

    ! Energy integration parameters
    integer            :: i_E                             ! Index for energy integration
    integer, parameter :: n_energy    = 1000              ! Number of energy points for integration
                                                          ! Note: n_energy must be dividible by 4
    integer, parameter :: n_energy_1  = 7 * n_energy / 20 ! Number of energy points in lower tail
    integer, parameter :: n_energy_2  = 6 * n_energy / 20 ! Number of energy points in Fermi-window
    integer, parameter :: n_energy_3  = 7 * n_energy / 20 ! Number of energy points in upper tail
                                                          ! Note: n_energy_1 + n_energy_2 + n_energy_3 = n_energy
    


    ! =================== !
    ! Parameters from EPW !
    ! =================== !

    ! -------------------------------------- !
    ! Grids for the Non-Interacting Energies !
    ! -------------------------------------- !
    integer         :: Nk
    integer         :: Nk_cube
    integer         :: Nq
    integer         :: Nq_cube
    integer, target :: Nq_tot_0              ! Deviates from Nq_cube only if Gamma_centered_q = .true.

    ! ---------------------------------------------------- !
    ! Grids and Parameters for the Self-Energy-Calculation !
    ! ---------------------------------------------------- !
    real(dp)        :: fsthick                             ! Thickness of the fine region around E_F in eV
    integer         :: Nk_SE                               ! Number of fine energies within the Fermi-window; Deviates from Nk_cube only if Fermi_centered_k = .true.
    integer         :: Nq_SE
    integer         :: Nq_SE_cube
    integer, target :: Nq_tot_SE                           ! Deviates from Nq_SE_cube only if Gamma_centered_q = .true.

    ! -------------------------------- !
    ! Grids for Phonon-Correction-Term !
    ! -------------------------------- !
    integer :: Ni_uniform
    integer :: Ni_uniform_cube

    ! ------------------------------ !
    ! Structure-dependent parameters !
    ! ------------------------------ !
    integer  :: r                    ! = number of atoms in the UC
    integer  :: N_bnd                ! = number of electronic bands in the reduced Wannier-System
    integer  :: n_UC                 ! = number of electrons per UC
    integer  :: n_branches           ! = r * dx = Number of phonon branches
    real(dp) :: N_electrons          ! = n_UC * N_p Particle number for charge neutrality
    real(dp) :: E_F                  ! Fermi energy in eV

    ! ------------------ !
    ! Physical constants !
    ! ------------------ !
    integer , parameter :: dx       = 3                    ! Number of spatial dimensions
    real(dp), parameter :: pi       = 3.141592653589793_dp
    real(dp), parameter :: hbar     = 6.58211e-16_dp       ! in eV s
    real(dp), parameter :: k_B      = 8.61733e-5_dp        ! Boltzmann constant in eV / K
    real(dp), parameter :: f_factor = 1.66e-24_dp          ! Number of moles in 1 unit cell (N_p atoms)
    real(dp), parameter :: eV_to_J  = 1.60218e-19_dp       ! eV to Joules conversion factor

    ! ---------- !
    ! File paths !
    ! ---------- !
    character(len=10) :: Nk_str         ! Will hold string version of Nk
    character(len=10) :: Nq_str         ! Will hold string version of Nq
    character(len=10) :: Nq_SE_str      ! Will hold string version of Nq_SE
    character(len=10) :: fsthick_str    ! Will hold string version of fsthick
    character(len=10) :: Ni_uniform_str ! Will hold string version of Ni_uniform
    

    ! ----------------------------------------------- !
    ! File-Paths for the Non-Interacting Hamiltonians !
    ! ----------------------------------------------- !
    character(len=200) :: file_tree      ! Base filepath for all data files
    character(len=250) :: epsilons_file  ! = trim(file_tree)//"epsilon_data.txt"
    character(len=250) :: omegas_file    ! = trim(file_tree)//"omega_data.txt"

    ! ------------------------------------------- !
    ! File-Paths for the Self-Energy-Computations !
    ! ------------------------------------------- !
    character(len=250) :: SE_epsilons_file        ! = trim(file_tree)//"SE_epsilon_data.txt"
    character(len=250) :: SE_epsilons_kplusq_file ! = trim(file_tree)//"SE_epsilon_kplusq_data.txt"
    character(len=250) :: SE_couplings_file       ! = trim(file_tree)//"SE_g_data.txt"
    character(len=250) :: SE_omegas_file          ! = trim(file_tree)//"SE_omega_data.txt"

    ! -------------------------------------------------------- !
    ! File-Paths for the Specific-Heat-Correction-Computations !
    ! -------------------------------------------------------- !
    character(len=250) :: uniform_epsilons_file        ! = trim(file_tree)//"uniform_epsilon_data.txt"
    character(len=250) :: uniform_couplings_file       ! = trim(file_tree)//"uniform_g_data.txt"
    character(len=250) :: uniform_omegas_file          ! = trim(file_tree)//"uniform_omega_data.txt"
    character(len=250) :: uniform_k_points_file        ! = trim(file_tree)//"uniform_k_data.txt"
    
    ! ----------------------------- !
    ! File-Paths for Input / Output !
    ! ----------------------------- !
    character(len=200) :: mu_vs_T_filename
    character(len=200) :: T_filename
    character(len=200) :: c_V_per_uc_filename
    integer            :: unit, ios


    ! ------------------------------------- !
    ! Arrays for the Non-Interacting System !
    ! ------------------------------------- !
    real(dp), allocatable, target :: epsilon_arr(:,:) ! (N_bnd, N_p)                             in  eV
    real(dp), allocatable, target :: omega_arr(:,:)   ! (n_branches, Nq_fine_cube / Nq_tot_0)   in meV
    
    ! ---------------------------------------------------------------- !
    ! Arrays for the Interacting System in the Self-Energy-Calculation !
    ! ---------------------------------------------------------------- !
    real(dp), allocatable, target :: SE_epsilon_arr(:,:)          ! (N_bnd, Nk_SE)                                            in  eV
    real(dp), allocatable, target :: SE_epsilon_kplusq_arr(:,:,:) ! (N_bnd, Nk_SE, Nq_SE_cube / Nq_tot_SE)                    in  eV
    real(dp), allocatable, target :: SE_omega_arr(:,:)            ! (n_branches, Nq_SE_cube / Nq_tot_SE)                      in meV
    real(dp), allocatable, target :: SE_g_arr(:,:,:,:,:)          ! (N_bnd, N_bnd, n_branches, Nk_SE, Nq_SE_cube / Nq_tot_SE) in meV

    real(dp), allocatable, target :: uniform_epsilon_arr(:,:)          ! (N_bnd, Nk_uniform)                                  in  eV
    real(dp), allocatable, target :: uniform_omega_arr(:,:)            ! (n_branches, Nq_uniform)                             in meV
    real(dp), allocatable, target :: uniform_g_arr(:,:,:,:,:)          ! (N_bnd, N_bnd, n_branches, Nk_uniform, Nq_uniform)   in meV
    real(dp), allocatable         :: uniform_wavevector_arr(:,:)       ! (3, Nq_uniform)


    ! ----------------------- !
    ! Arrays for calculations !
    ! ----------------------- !
    real(dp), allocatable :: E_arr(:)      ! Energy array for integration              in eV
    real(dp), allocatable :: mu_arr(:)     ! Chemical potential array vs Temperature   in eV
    real(dp), allocatable :: c_V_per_uc(:) ! Specific heat capacity per unit cell vs T in eV/K
    real(dp), allocatable :: temp_arr(:)   ! Temperature array                         in    K

    ! Other variables  
    real(dp) :: E_min, E_max                 ! Min and max energies of the electronic structure 
    real(dp) :: E_min_rel, E_max_rel         ! E_min - E_F and E_max - E_F
    real(dp) :: N_states_per_UC              ! Computed Number of states per UC from the reduced Wannier-system
    real(dp) :: N_mu_eq_E_F_T_eq_zero        ! Particle number computed from mu = E_F at T = 0 K
    real(dp) :: N_per_UC_mu_eq_E_F_T_eq_zero ! Particle number per UC computed from mu = E_F at T = 0 K 

    real(dp) :: gamma, const ! For specific heat coefficient calculation
    real(dp) :: elapsed_time
    integer  :: start_time   ! For timing the first calculation
    integer  :: end_time
    integer  :: clock_rate


    ! MPI variables
    integer :: mpi_rank, mpi_size, mpi_ierr

    ! MPI distribution variables
    integer :: i_k_start, i_k_end, Nk_local_MPI

contains 

    subroutine set_structure_parameters()
        ! eta = Broadening for retarded Green function in eV; Note: 0.13 eV = 0.01 Ry
        select case(trim(structure_name))
        case("Au_Fm3m")
            r = 1; N_bnd = 1; n_UC = 1
            E_F = 19.3283
            eta        =   2.0 * 0.13_dp
            Nk         =  40
            Nq         =  60
            Nq_SE      =  6
            fsthick    =  50.0_dp
            Ni_uniform =  10
            T_min      =  10.0
            T_max      = 400.0
            n_temp     =  5
            !T_min      =  0.7
            !T_max      =  1.5
            !n_temp     =  5
        case("Cu_Fm3m")
            r = 1; N_bnd = 1; n_UC = 1
            E_F = 15.9764
            eta        =   2.0 * 0.13_dp
            Nk         =  40
            Nq         =  60
            Nq_SE      =  6
            fsthick    =  50.0_dp
            Ni_uniform =  10
            T_min      =  1.0
            T_max      = 1000.0
            n_temp     =  40
        case("Si_Fd3m1")
            r = 2; N_bnd = 8; n_UC = 8
            E_F = 5.7011
            eta        =   2.0 * 0.13_dp
            Nk         =  40
            Nq         =  60
            Nq_SE      =  6
            fsthick    =  50.0_dp
            Ni_uniform =  10
            T_min      =  20.0
            T_max      = 500.0
            n_temp     =  20
        case("Diamond_Fd3m")
            r = 2; N_bnd = 8; n_UC = 8
            E_F = 15.1478
            eta        =   2.0 * 0.13_dp
            Nk         =  40
            Nq         =  60
            Nq_SE      =  6
            fsthick    =  50.0_dp
            Ni_uniform =  10
            T_min      =  10.0
            T_max      = 300.0
            n_temp     =  20
        case("Rb_Im3m")
            r = 1; N_bnd = 1; n_UC = 1
            E_F = 2.0062
            eta        =   2.0 * 0.13_dp
            Nk         =  40
            Nq         =  60
            Nq_SE      =  6
            fsthick    =  50.0_dp
            Ni_uniform =  10
            T_min      =  20.0
            T_max      = 125.0
            n_temp     =  20
        case("Ag_Fm3m")
            r = 1; N_bnd = 1; n_UC = 1
            E_F = 15.6739
            eta        =   2.0 * 0.13_dp
            Nk         =  40
            Nq         =  60
            Nq_SE      =  6
            fsthick    =  50.0_dp
            Ni_uniform =  10
            T_min      =  10.0
            T_max      = 300.0
            n_temp     =  20
        case("K_Im3m")
            r = 1; N_bnd = 1; n_UC = 1
            E_F        =   1.5503
            ! eta        =  25.0e-3_dp
            eta        =   1.0 * 0.13_dp
            Nk         =  40
            Nq         =  60
            Nq_SE      =   6
            fsthick    =   50.0_dp
            Ni_uniform =  10
            ! T_min      =   0.5
            ! T_max      =   1.41
            T_min      =   10.0
            T_max      =   150.0
            n_temp     =   20
        case("Cs_Im3m")
            r = 1; N_bnd = 1; n_UC = 1
            E_F        =  1.98770
            eta        = 1.0 * 0.13_dp
            Nk         = 40
            Nq         = 60
            Nq_SE      = 6
            fsthick    = 50.0_dp
            Ni_uniform =  10
            T_min      =    2.0
            T_max      =  100.0
            n_temp     =  30
        case("Pb_Fm3m")
            r = 1; N_bnd = 3; n_UC = 2
            E_F        =  12.6068 
            eta        =   0.13_dp
            Nk         =  40
            Nq         =  60
            Nq_SE      =  6
            fsthick    =  50.0_dp
            Ni_uniform =  10
            T_min      =  10.0
            T_max      = 400.0
            n_temp     = 30
            ! n_temp     = 5
        case("Sn_Fd3m")
            r = 2; N_bnd = -1; n_UC = -1
            E_F = -1 
        case default
            write(*,*) "Unknown structure: ", trim(structure_name)
            stop
        end select
        Nk_cube         = Nk    ** 3
        Nq_cube         = Nq    ** 3
        Nq_SE_cube      = Nq_SE ** 3
        Ni_uniform_cube = Ni_uniform ** 3
        n_branches      =   r  * dx
        N_electrons     = real(n_UC * Nk_cube)
    end subroutine set_structure_parameters


    subroutine set_filepaths()
        write(Nk_str  , '(I0)') Nk             ! Convert Nk        (i.e. an integer number) to string
        write(Nq_str, '(I0)') Nq               ! Convert Nq        (i.e. an integer number) to string
        write(Nq_SE_str, '(I0)') Nq_SE ! Convert Nq_SE (i.e. an integer number) to string
        if (fsthick > 9.9) then 
            fsthick_str = "INF"
        else 
            ! Parse now fsthick * 1000 to an integer and subsequently write it to a string
            write(fsthick_str, '(I0)') int(fsthick * 1000.0_dp)
        end if 
        write(Ni_uniform_str, '(I0)') Ni_uniform ! Convert Ni_uniform (i.e. an integer number) to string

        file_tree = "../B_QE_n_EPW_Files/Step3_EPW/" // trim(structure_name) // &
                    "/Dispersions_" // trim(Nk_str) // "_" // trim(Nq_str) // "_" // trim(Nq_SE_str) // &
                    "_" // trim(fsthick_str) // "_" // trim(Ni_uniform_str) // "/"

        epsilons_file           = trim(file_tree)//"epsilon_data.txt"
        omegas_file             = trim(file_tree)//"omega_data.txt"

        SE_epsilons_file        = trim(file_tree)//"SE_epsilon_data.txt"
        SE_epsilons_kplusq_file = trim(file_tree)//"SE_epsilon_kplusq_data.txt"
        SE_couplings_file       = trim(file_tree)//"SE_g_data.txt"
        SE_omegas_file          = trim(file_tree)//"SE_omega_data.txt"

        uniform_epsilons_file  = trim(file_tree)//"uniform_epsilon_data.txt"
        uniform_couplings_file = trim(file_tree)//"uniform_g_data.txt"
        uniform_omegas_file    = trim(file_tree)//"uniform_omega_data.txt"
        uniform_k_points_file  = trim(file_tree)//"uniform_k_data.txt"

    end subroutine set_filepaths


    subroutine print_elapsed_time(task_name, elapsed_seconds)
        character(len=*), intent(in) :: task_name
        real(dp), intent(in) :: elapsed_seconds
        real(dp) :: days, hours, minutes, seconds
        
        ! Only print from master process
        if (mpi_rank /= 0) return
        
        ! Calculate days, hours, minutes, seconds (keeping fractional parts)
        days = elapsed_seconds / 86400.0    ! 86400 seconds in a day
        hours = elapsed_seconds / 3600.0    ! 3600 seconds in an hour
        minutes = elapsed_seconds / 60.0    ! 60 seconds in a minute
        seconds = elapsed_seconds
        
        ! Print formatted time - always show all units
        write(*,*)
        write(*,*)
        write(*,'(A,A,A,F8.4,A,F8.4,A,F8.4,A,F8.4,A)') &
            " - ", trim(task_name), " took ", &
            days, " days, ", &
            hours, " hours, ", &
            minutes, " minutes, ", &
            seconds, " seconds"
        write(*,'(A)') " - Approximated runtime for the remaining tasks: TBImplemented."
        write(*,*)
        write(*,*)
    end subroutine print_elapsed_time


    subroutine linspace(start, end, n, arr)
        real(dp), intent(in)    :: start, end
        real(dp), intent(out)   :: arr(:)
        real(dp)                :: step
        integer, intent(in) :: n
        integer             :: i

        if (n <= 1) then
            arr(1) = start
            return
        end if

        step = (end - start) / real(n - 1)
        do i = 1, n
            arr(i) = start + real(i - 1) * step
        end do

    end subroutine linspace


    subroutine print_with_border(message)
        character(len=*), intent(in) :: message
        integer :: len_msg

        ! Only master process prints
        if (mpi_rank /= 0) return

        len_msg = len_trim(message)
        write(*,*)
        write(*,'(A)') "+"//repeat("=", len_msg + 4)//"+"
        write(*,'(A)') "|| "//trim(message)//" ||"
        write(*,'(A)') "+"//repeat("=", len_msg + 4)//"+"
        write(*,*)
    end subroutine print_with_border


    subroutine print_with_border2(message)
        character(len=*), intent(in) :: message
        integer :: len_msg

        ! Only master process prints
        if (mpi_rank /= 0) return

        len_msg = len_trim(message)
        write(*,*)
        write(*,'(A)') "+"//repeat("-", len_msg + 4)//"+"
        write(*,'(A)') "|  "//trim(message)//"  |"
        write(*,'(A)') "+"//repeat("-", len_msg + 4)//"+"
        write(*,*)
    end subroutine print_with_border2


    subroutine print_intro()

        write(*,*)
        write(*,'(A)') "     ================================================================================ "
        write(*,'(A)') "    ||                                                                              || "
        write(*,'(A)') "    ||  GGGGGGGGGGGGGG    RRRRRRRRRRRRR       IIIIIIIIIIIIII    TTTTTTTTTTTTTTTTTT  || "
        write(*,'(A)') "    ||  GGGGGGGGGGGGGG    RRRRRRRRRRRRRR      IIIIIIIIIIIIII    TTTTTTTTTTTTTTTTTT  || "
        write(*,'(A)') "    ||  GGG               RRRR       RRRR          IIII                TTTT         || "
        write(*,'(A)') "    ||  GGG               RRRR       RRRR          IIII                TTTT         || "
        write(*,'(A)') "    ||  GGG      GGGG     RRRRRRRRRRRRR            IIII                TTTT         || "
        write(*,'(A)') "    ||  GGG      GGGGG    RRRRRRRRRRRRR            IIII                TTTT         || "
        write(*,'(A)') "    ||  GGG          GG   RRRR      RRRR           IIII                TTTT         || "
        write(*,'(A)') "    ||  GGG          GG   RRRR       RRRR          IIII                TTTT         || "
        write(*,'(A)') "    ||  GGGGGGGGGGGGGG    RRRR        RRRR    IIIIIIIIIIIIII           TTTT         || "
        write(*,'(A)') "    ||  GGGGGGGGGGGGGG    RRRR         RRRR   IIIIIIIIIIIIII           TTTT         || "
        write(*,'(A)') "    ||                                                                              || "
        write(*,'(A)') "     ================================================================================ "
        write(*,*)
        write(*,'(A)') "              [Gr]een-Functions for [I]nteracting [T]hermodynamics"
        write(*,'(A)') "              (Post-Processing Tool for the Electron-Phonon-Calculations)"
        write(*,*)
        write(*,*)

    end subroutine print_intro


    subroutine setup_temp_arr()

        if (.not. allocated(temp_arr)) then
            allocate(temp_arr(n_temp))
        end if
        
        call linspace(T_min, T_max, n_temp, temp_arr)
    end subroutine setup_temp_arr
    
    
    subroutine setup_Energy_arr_DOS()
        ! In the case of Fermi_centered_k == .true., we create a denser energy grid around E_F

        ! Variables for Fermi_centered_k = .true. 
        real(dp), allocatable :: E_arr_1(:), E_arr_2(:), E_arr_3(:)

        if (.not. Fermi_centered_k) then 
        
            allocate(E_arr(n_energy))
            call linspace(E_min_rel - E_tails, E_max_rel + E_tails, n_energy, E_arr)
        
        else if (Fermi_centered_k) then 

            ! Allocate partitions

            allocate(E_arr(n_energy))
            allocate(E_arr_1(n_energy_1))
            allocate(E_arr_2(n_energy_2))
            allocate(E_arr_3(n_energy_3))

            ! Lower-energy tail
            call linspace( E_min_rel - E_tails, &
                           0 - fsthick,         &
                           n_energy_1,          &
                           E_arr_1 )

            ! Fermi-region
            call linspace( 0 - fsthick, &
                           0 + fsthick, &
                           n_energy_2,  &
                           E_arr_2 )

            ! Upper-energy tail
            call linspace( 0 + fsthick,         &
                           E_max_rel + E_tails, &
                           n_energy_3,          &
                           E_arr_3 )

            ! Concatenate
            E_arr(1                           : n_energy_1)                           = E_arr_1
            E_arr(n_energy_1 + 1              : n_energy_1 + n_energy_2)              = E_arr_2
            E_arr(n_energy_1 + n_energy_2 + 1 : n_energy_1 + n_energy_2 + n_energy_3) = E_arr_3

            deallocate(E_arr_1, E_arr_2, E_arr_3)

        end if

    end subroutine setup_Energy_arr_DOS

    
    subroutine setup_Energy_arr_DOOS()

        ! Variables for Fermi_centered_k = .true. 
        real(dp), allocatable :: E_arr_1(:), E_arr_2(:), E_arr_3(:)

        if (.not. Fermi_centered_k) then 
        
            if (.not. allocated(E_arr)) allocate(E_arr(n_energy))
            call linspace(E_min_rel - E_tails, E_max_rel + E_tails, n_energy, E_arr)
        
        else if (Fermi_centered_k) then 
            ! In the case of Fermi_centered_k == .true., we create a denser energy grid around E_F

            ! Allocate partitions

            if (.not. allocated(E_arr)) allocate(E_arr(n_energy))
            allocate(E_arr_1(n_energy_1))
            allocate(E_arr_2(n_energy_2))
            allocate(E_arr_3(n_energy_3))

            ! Lower-energy tail
            call linspace( E_min_rel - E_tails, &
                           0 - fsthick,         &
                           n_energy_1,          &
                           E_arr_1 )

            ! Fermi-region
            call linspace( 0 - fsthick, &
                           0 + fsthick, &
                           n_energy_2,  &
                           E_arr_2 )

            ! Upper-energy tail
            call linspace( 0 + fsthick,         &
                           E_max_rel + E_tails, &
                           n_energy_3,          &
                           E_arr_3 )

            ! Concatenate
            E_arr(1                           : n_energy_1)                           = E_arr_1
            E_arr(n_energy_1 + 1              : n_energy_1 + n_energy_2)              = E_arr_2
            E_arr(n_energy_1 + n_energy_2 + 1 : n_energy_1 + n_energy_2 + n_energy_3) = E_arr_3

            deallocate(E_arr_1, E_arr_2, E_arr_3)

        end if

    end subroutine setup_Energy_arr_DOOS


    subroutine setup_Energy_arr_SH()
        call linspace( - E_Length_SH, + E_Length_SH, n_energy, E_arr)
    end subroutine setup_Energy_arr_SH


    function fermi_function(E, T) result(f)
        real(dp), intent(in) :: E, T
        real(dp)             :: f, beta, x

        if (T <= 1.0e-12) then
            if (E < 0.0) then
                f = 1.0
            else
                f = 0.0
            end if
            return
        end if

        beta = 1.0 / (k_B * T)
        x = beta * E

        f = 1.0 / (1.0 + exp(x))
        
    end function fermi_function


    function bose_function(E, T) result(b)
        real(dp), intent(in) :: E, T
        real(dp)             :: b, beta, x

        beta = 1.0 / (k_B * T)
        x = beta * E

        if (x > 40.0) then
            b = 0.0
        else if (abs(x) < 1.0e-6) then
            ! Correct series expansion: n_B(E) ≈ kT/E - 1/2 + E/(12kT) + ...
            b = 1.0 / x - 0.5 + x / 12.0
        else
            b = 1.0 / (exp(x) - 1.0)
        end if
    end function bose_function


    function trapz(x, y) result(integral)
        real(dp), intent(in) :: x(:), y(:)
        real(dp)             :: integral
        integer              :: i, n

        n = size(x)
        integral = 0.0

        do i = 1, n-1
            integral = integral + 0.5 * (y(i) + y(i+1)) * (x(i+1) - x(i))
        end do
    end function trapz


    function csch(x) result(csch_val)
        implicit none
        real(dp), intent(in) :: x
        real(dp)             :: csch_val
        csch_val = 1.0 / sinh(x)
    end function csch


    function sech(x) result(sech_val)
        implicit none
        real(dp), intent(in) :: x
        real(dp)             :: sech_val
        sech_val = 1.0 / cosh(x)
    end function sech


    function coth(x) result(coth_val)
        implicit none
        real(dp), intent(in) :: x
        real(dp)             :: coth_val
        coth_val = cosh(x) / sinh(x)
    end function coth

end module modmain
