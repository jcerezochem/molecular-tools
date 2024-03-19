program prepare_numder_coms


    !==============================================================
    ! This code uses of MOLECULAR_TOOLS (version 1.0/March 2012)
    !==============================================================
    !
    ! Description:
    ! -----------
    ! Program to analyse vibrations in term of internal coordinates.
    !
    ! Compilation instructions (for mymake script):
    !make$ echo "COMPILER: $FC"; sleep 1; $FC ../modules/alerts.f90 ../modules/structure_types_v3.f90 ../modules/line_preprocess.f90 ../modules/ff_build_module_v3.f90 ../modules/gro_manage_v2.f90 ../modules/pdb_manage_v2.f90 ../modules/constants_mod.f90 ../modules/atomic_geom_v2.f90 ../modules/gaussian_manage_v2.f90 ../modules/gaussian_fchk_manage_v2.f90 ../modules/symmetry_mod.f90 ../modules/MatrixMod.f90 internal_SR_v6.f90 internal_duschinski_v5.f90 -llapack -o internal_duschinski_v5.exe -cpp -DDOUBLE
    !
    ! Change log:
    !
    ! TODO:
    ! ------
    !
    ! History
    ! V3: Internal analysis is based on internal_duschinski_v5 (never finished...)
    ! V4: Internal analysis is based on internal_duschinski_v7
    !  V4b (not in the main streamline!): includes the generation of scans calc. for Gaussian.
    !  V4c: the same as 4b. Bug fixes on guess_connect (ff_build module_v3 was buggy)
    !
    !Addapted to v4 release (distribution upgrade). Feb '14
    !v4 releases:
    !v4.0.1:
    ! -use redundant coordinates
    !v4.0.1.1:
    ! - include UnSym (MOLCAS) file as input (freqs and nm)
    !v4.0.1.2:
    ! - if the calculation is only a internal scan, do not need the hessian (so do not try to read it)
    !
    !============================================================================    


    !*****************
    !   MODULE LOAD
    !*****************
    !============================================
    !   Generic
    !============================================
    use alerts
    use line_preprocess
    use constants 
    use verbosity
    use matrix
    use matrix_print
    !============================================
    !   Structure types module
    !============================================
    use structure_types
    !============================================
    !   File readers
    !============================================
    use generic_io
    use generic_io_molec
    use xyz_manage
    use gaussian_manage
    !============================================
    !  Structure-related modules
    !============================================
    use molecular_structure
    use atomic_geom
    use symmetry
    !============================================
    !  Vibrational info
    !============================================
    use vibrational_analysis, only: Hlt_to_Hess,Hess_to_Hlt


    implicit none

    integer,parameter :: NDIM = 600

    !====================== 
    !Options 
    logical :: overwrite=.false.
    !======================

    !====================== 
    !System variables
    type(str_resmol)  :: molecule, molec_aux
    character(len=50) :: calc_type,basis
    character(len=60) :: method
    integer           :: charge, mult
    integer,dimension(1:NDIM) :: iord
    integer :: Nat, Nvib
    character(len=5) :: PG
    ! Energies
    real(8) :: E_scf, E_td, E_tot
    ! FCHK data
    real(8),dimension(1:NDIM) :: Grad
    ! DELTA
    real(8) :: delta=1.d-3
    !======================

    !=====================
    ! Available data swithes
    logical :: have_gradient, have_hessian,&
               ! Energies
               have_SCF,have_TD,have_TOT
    !===================== 

    !====================== 
    !MATRICES
    !Other matrices
    real(8) :: R, E
    !====================== 

    !====================== 
    !Read fchk auxiliars
    real(8),dimension(:),allocatable :: A, GEOMS, GRADS
    integer,dimension(:),allocatable :: IA
    character(len=1) :: dtype
    integer :: error, N, lenght
    !====================== 

    ! IRC vars
    integer :: ninfo, nsteps

    !====================== 
    !Auxiliars for LAPACK matrix nversion
    integer :: info
    integer,dimension(NDIM) :: ipiv
    real(8),dimension(NDIM,NDIM) :: work
    !====================== 

    !====================== 
    !Auxiliar variables
    character(1) :: null
    character(len=50) :: dummy_char
    character(len=240):: line
    integer :: iat_new, iat_orig, j_orig, j_new, nswap
    !====================== 

    !=============
    !Counters
    integer :: i,j,k,l, ii,jj,kk, iat, nn, imin, imax, iii, i_scan
    !=============

    !================
    !I/O stuff 
    !units
    integer :: I_INP =10,  &
               O_STR =20,  &
               O_GAU =21,  &
               O_FCHK=22
    !files
    character(len=10)  :: filetype="guess", filetype_out="guess"
    character(len=200) :: inpfile ="input.fchk", &
                          outfile="default"
    character(len=200) :: basefile, title, path
    character(len=10)  :: extension, label
    !status
    integer :: IOstatus
    !===================

    !===================
    !CPU time 
    real(8) :: ti, tf
    !===================

! (End of variables declaration) 
!==================================================================================
    call cpu_time(ti)

    !===========================
    ! Allocate atoms (default)
    call allocate_atoms(molecule)
    call allocate_atoms(molec_aux)
    !===========================

    ! 0. GET COMMAND LINE ARGUMENTS
    call parse_input(inpfile,filetype,outfile,filetype_out,overwrite,delta)

    ! Output names are labeled with steps
    call split_line_back(outfile,".",basefile,extension)
    ! Files are generated on the folder where they are called, so remove relative path if any
    if (index(basefile,'/') /= 0) then
        call split_line_back(basefile,"/",path,basefile)
    else
        path='.'
    endif
 
    !================
    ! READ DATA
    !================
    if (adjustl(filetype) == "guess") &
    call split_line_back(inpfile,".",null,filetype)
    if (filetype /= "fchk") call alert_msg("fatal","Only FCHK are currently supported")
    if (adjustl(filetype_out) == "guess") &
    call split_line_back(outfile,".",null,filetype_out)

    ! STRUCTURE FILE
    print'(X,A)', "READING STRUCTURE..."
    open(I_INP,file=inpfile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(inpfile)) )
    call generic_strmol_reader(I_INP,filetype,molecule)
    !Shortcuts
    Nat = molecule%natoms
    Nvib = 3*Nat-6
    print'(X,A,/)', "Done"

    ! JOB INFO
    print'(X,A)', "READING JOB INFO..."
    rewind(I_INP)
    call read_gauss_job(I_INP,filetype,calc_type,method,basis)
    print'(X,A)', " Job type: "//trim(adjustl(calc_type))
    print'(X,A)', " Method  : "//trim(adjustl(method))
    print'(X,A)', " Basis   : "//trim(adjustl(basis))
    rewind(I_INP)
    call read_gauss_chargemult(I_INP,filetype,charge,mult)
    print'(X,A,I0)', " Charge  : ", charge
    print'(X,A,I0)', " Mult.   : ", mult
    print'(X,A,/)', "Done"

    
    ! Write steps
    do i=1,molecule%natoms
        ! ** X **
        ! Bw
        molecule%atom(i)%x = molecule%atom(i)%x - delta
        write(outfile,'(A,I0,A)') trim(adjustl(basefile))//"_at", i, "_xyz1_bw.com"
        write(title,'(A,I0,A)') "Backward step for atom=", i, " and xyz=1"
        open(O_GAU,file=outfile)
        call write_gcom(O_GAU,molecule,&
                              !Optional args
                              chkname=outfile, &! 
                              calc='SP NoSymm',&!calc_type,  &! e.g. SP, Freq, Opt...
                              method=method,   &!
                              basis=basis,     &!
                              title=title,     &!
                              free_format=.false.)
        close(O_GAU)
        molecule%atom(i)%x = molecule%atom(i)%x + delta
        ! Fw
        molecule%atom(i)%x = molecule%atom(i)%x + delta
        write(outfile,'(A,I0,A)') trim(adjustl(basefile))//"_at", i, "_xyz1_fw.com"
        write(title,'(A,I0,A)') "Forward step for atom=", i, " and xyz=1"
        open(O_GAU,file=outfile)
        call write_gcom(O_GAU,molecule,&
                              !Optional args
                              chkname=outfile, &! 
                              calc='SP NoSymm',&!calc_type,  &! e.g. SP, Freq, Opt...
                              method=method,   &!
                              basis=basis,     &!
                              title=title,     &!
                              free_format=.false.)
        close(O_GAU)
        molecule%atom(i)%x = molecule%atom(i)%x - delta
        ! ** Y **
        ! Bw
        molecule%atom(i)%y = molecule%atom(i)%y - delta
        write(outfile,'(A,I0,A)') trim(adjustl(basefile))//"_at", i, "_xyz2_bw.com"
        write(title,'(A,I0,A)') "Backward step for atom=", i, " and xyz=2"
        open(O_GAU,file=outfile)
        call write_gcom(O_GAU,molecule,&
                              !Optional args
                              chkname=outfile, &! 
                              calc='SP NoSymm',&!calc_type,  &! e.g. SP, Freq, Opt...
                              method=method,   &!
                              basis=basis,     &!
                              title=title,     &!
                              free_format=.false.)
        close(O_GAU)
        molecule%atom(i)%y = molecule%atom(i)%y + delta
        ! Fw
        molecule%atom(i)%y = molecule%atom(i)%y + delta
        write(outfile,'(A,I0,A)') trim(adjustl(basefile))//"_at", i, "_xyz2_fw.com"
        write(title,'(A,I0,A)') "Forward step for atom=", i, " and xyz=2"
        open(O_GAU,file=outfile)
        call write_gcom(O_GAU,molecule,&
                              !Optional args
                              chkname=outfile, &! 
                              calc='SP Nosymm',&!calc_type,  &! e.g. SP, Freq, Opt...
                              method=method,   &!
                              basis=basis,     &!
                              title=title,     &!
                              free_format=.false.)
        close(O_GAU)
        molecule%atom(i)%y = molecule%atom(i)%y - delta
        ! ** Z **
        ! Bw
        molecule%atom(i)%z = molecule%atom(i)%z - delta
        write(outfile,'(A,I0,A)') trim(adjustl(basefile))//"_at", i, "_xyz3_bw.com"
        write(title,'(A,I0,A)') "Backward step for atom=", i, " and xyz=3"
        open(O_GAU,file=outfile)
        call write_gcom(O_GAU,molecule,&
                              !Optional args
                              chkname=outfile, &! 
                              calc='SP NoSymm',&!calc_type,  &! e.g. SP, Freq, Opt...
                              method=method,   &!
                              basis=basis,     &!
                              title=title,     &!
                              free_format=.false.)
        close(O_GAU)
        molecule%atom(i)%z = molecule%atom(i)%z + delta
        ! Fw
        molecule%atom(i)%z = molecule%atom(i)%z + delta
        write(outfile,'(A,I0,A)') trim(adjustl(basefile))//"_at", i, "_xyz3_fw.com"
        write(title,'(A,I0,A)') "Forward step for atom=", i, " and xyz=3"
        open(O_GAU,file=outfile)
        call write_gcom(O_GAU,molecule,&
                              !Optional args
                              chkname=outfile, &! 
                              calc='SP NoSymm',&!calc_type,  &! e.g. SP, Freq, Opt...
                              method=method,   &!
                              basis=basis,     &!
                              title=title,     &!
                              free_format=.false.)
        close(O_GAU)
        molecule%atom(i)%z = molecule%atom(i)%z - delta
    
    enddo    

    call cpu_time(tf)
    write(0,'(/,A,X,F12.3,/)') "CPU time (s)", tf-ti

    stop


    !==============================================
    contains
    !=============================================

    subroutine parse_input(inpfile,filetype,outfile,filetype_out,overwrite,delta)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: inpfile,filetype,outfile,filetype_out
        logical,intent(inout)          :: overwrite
        real(8),intent(inout)          :: delta
        ! Local
        logical :: argument_retrieved,  &
                   need_help = .false.
        integer:: i
        character(len=200) :: arg
        ! iargc type must be specified with implicit none (strict compilation)
        integer :: iargc

        argument_retrieved=.false.
        do i=1,iargc()
            if (argument_retrieved) then
                argument_retrieved=.false.
                cycle
            endif
            call getarg(i, arg) 
            select case (adjustl(arg))
                case ("-f") 
                    call getarg(i+1, inpfile)
                    argument_retrieved=.true.
                case ("-ft") 
                    call getarg(i+1, filetype)
                    argument_retrieved=.true.  
                case("-o")
                    call getarg(i+1, outfile)
                    argument_retrieved=.true.
                case ("-fto") 
                    call getarg(i+1, filetype_out)
                    argument_retrieved=.true.  
                case("-ow")
                    overwrite=.true.
                case("-delta")
                    call getarg(i+1, arg)
                    read(arg,*) delta
                    argument_retrieved=.true.
                case ("-h")
                    need_help=.true.

                case default
                    call alert_msg("fatal","Unkown command line argument: "//adjustl(arg))
            end select
        enddo 

        !The default output is now xyz
        if (adjustl(outfile) == "default") then
            call split_line_back(inpfile,".",outfile,null)
            outfile=trim(adjustl(outfile))//".xyz"
        endif

       !Print options (to stdx)
        write(0,'(/,A)') '========================================================'
        write(0,'(/,A)') '             P R E P A R E   N U M D E R  C O M S '    
        write(0,'(/,A)') '       Prepare input files (.com) for numerical ' 
        write(0,'(A)')   '       derivatives displacing the initila structure'
        write(0,'(A)')   '               +/- delta Angs for each atom'
        call print_version()
        write(0,'(/,A)') '========================================================'
        write(0,'(/,A)') '-------------------------------------------------------------------'
        write(0,'(A)')   ' Flag         Description                      Value'
        write(0,'(A)')   '-------------------------------------------------------------------'
        write(0,*)       '-f           Input file                       ', trim(adjustl(inpfile))
        write(0,*)       '-ft          \_ FileTyep                      ', trim(adjustl(filetype))
        write(0,*)       '-o           Output file                      ', trim(adjustl(outfile))
        write(0,*)       '-fto         \_ FileTyep                      ', trim(adjustl(filetype_out))
        write(0,'(X,A,ES9.3)')    '-delta       Delta (angstrom)                 ', delta
        write(0,*)       '-ow          Force overwrite output          ',  overwrite
        write(0,*)       '-h           This help                       ',  need_help
        write(0,*)       '-------------------------------------------------------------------'
        if (need_help) call alert_msg("fatal", 'There is no manual (for the moment)' )

        return
    end subroutine parse_input
       

end program prepare_numder_coms

