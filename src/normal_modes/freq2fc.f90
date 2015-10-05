program freq2fc


    !==============================================================
    ! This code uses of MOLECULAR_TOOLS (version 1.0/March 2012)
    !==============================================================
    !
    ! Description:
    ! -----------
    ! Program to convert between Force Constant (atomic units) and Frequencies (cm^-1) (both directions)
    !
    ! Compilation instructions (for mymake script): now using make
    !
    ! Change log:
    !
    ! TODO:
    ! ------
    !
    ! History
    !============================================================================    

!*****************
!   MODULE LOAD
!*****************
!============================================
!   Generic (structure_types independent)
!============================================
    use alerts
    use constants

    implicit none

    real(8) :: freq, fc

    integer :: i,j


! (End of variables declaration) 
!==================================================================================

    ! 0. GET COMMAND LINE ARGUMENTS
    freq=-99999.d0
    fc  =-99999.d0
    call parse_input(freq,fc)

    if (freq == -99999.d0) then
        !convert fc to freq
        freq= sign(dsqrt(abs(fc)*HARTtoJ/BOHRtoM**2/AUtoKG)/2.d0/pi/clight/1.d2,&
                   fc)
        print'(A,G15.7)', "FREQ(cm-1) =", freq
    else if (fc == -99999.d0) then
        !convert freq to fc
        fc = (2.d0*pi*clight*1.d2*freq)**2 * BOHRtoM**2*AUtoKG/HARTtoJ
        print'(A,G15.7)', "FC(a.u.) =", fc
        fc = fc * HtoKCALM * CALtoJ
        print'(A,G15.7)', "FC(kJ/mol bohr^-1) =", fc
        fc = plank * clight * freq*1.d2 /HARTtoJ * HtoeV
        print'(A,G15.7)', "h\nu(eV) = ", fc
    else
        call alert_msg("fatal","You need to input either a -fq or a -fc")
    endif


    print*, ""

    stop


    !==============================================
    contains
    !=============================================

    subroutine parse_input(freq,fc)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        real(8),intent(inout) :: freq,fc
        ! Local
        logical :: argument_retrieved,  &
                   need_help = .false.
        integer:: i
        character(len=200) :: arg
        !Dummy arguments, for printing
        character(len=200) :: cfreq="", cfc=""

        argument_retrieved=.false.
        do i=1,iargc()
            if (argument_retrieved) then
                argument_retrieved=.false.
                cycle
            endif
            call getarg(i, arg) 
            select case (adjustl(arg))
                case ("-fq") 
                    call getarg(i+1, cfreq)
                    read(cfreq,*) freq
                    argument_retrieved=.true.

                case ("-fc") 
                    call getarg(i+1, cfc)
                    read(cfc,*) fc
                    argument_retrieved=.true.
        
                case ("-h")
                    need_help=.true.

                case default
                    print*, "Unkown command line argument: "//adjustl(arg)
                    stop
                    call alert_msg("fatal","Unkown command line argument: "//adjustl(arg))
            end select
        enddo 

       !Print options (to stderr)
        write(6,'(/,A)') '--------------------------------------------------'
        write(6,'(/,A)') '           FREQ to FC converter '          
        write(6,'(/,A)') '--------------------------------------------------'
        write(6,*) '-fq            ', trim(adjustl(cfreq))
        write(6,*) '-fc            ', trim(adjustl(cfc))
        write(6,*) '-h            ',  need_help
        write(6,*) '--------------------------------------------------'
        if (need_help) then
            write(6,*) ""
            write(6,*) "SHORT DESCRIPTION"
            write(6,*) " This program converts force constants (fc, in Atomic Units) into"
            write(6,*) " frequencies (freq, in cm-1) and viceversa."
            write(6,*) "USAGE"
            write(6,*) " freq2fc -fq 3000"
            write(6,*) "  (returns the fc in A.U)"
            write(6,*) " freq2fc -fc 30"
            write(6,*) "  (returns the fq in cm-1)"
            write(6,*) ""
            stop
        endif
!call alert_msg("fatal", 'There is no manual (for the moment)' )

        return
    end subroutine parse_input
       

end program freq2fc

