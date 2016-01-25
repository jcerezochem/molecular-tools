module alerts

    !==============================================================
    ! This code is part of MOLECULAR_TOOLS 
    !==============================================================
    ! Description
    !  Subroutine to print out alert messages
    !
    ! Dependences
    !
    ! Notes
    !  The module also includes some variables (acting as
    !  a common block for them)
    !-----------------------------------------------------

    integer, save :: n_notes=0 , &
                     n_errors=0, &
                     alert_unt=0
    integer, save :: nmax_notes=20
    ! Silent will mute the NOTES
    logical, save :: silent_notes = .true.

    contains

    subroutine alert_msg(attype,SENTENCE)

        implicit none

        character(len=*),intent(in):: attype, SENTENCE

        if (.not.silent_notes .and. n_notes == nmax_notes) then
            write(alert_unt,'(/,A,/)') "NOTE: Too many notes. Desactivating notes on output"
            n_notes=n_notes+1
            silent_notes=.true.
        endif

        select case (adjustl(attype))
            case ("note")
                if (.not.silent_notes) &
                write(alert_unt,'(/,A,A,/)') "NOTE: ",SENTENCE
                n_notes = n_notes + 1

            !The following is deprecated (maintained for backward compatibility)
            case ("error")
                write(alert_unt,'(/,A,A,/)') "WARNING: ",SENTENCE
                n_errors = n_errors + 1

            case ("warning")
                write(alert_unt,'(/,A,A,/)') "WARNING: ",SENTENCE
                n_errors = n_errors + 1

            case ("fatal")
                write(alert_unt,'(/,A,/)') "============================================"
                write(alert_unt,'(A,A,/)') "FATAL ERROR: ",trim(SENTENCE)
                write(alert_unt,'(A)')     "============================================"
                write(alert_unt,'(A,/)') "Exiting..."
                stop 1

            case default
                write(alert_unt,'(/,A,A,A,/)') attype," ", SENTENCE
        end select

        return

    end subroutine alert_msg


    subroutine summary_alerts

        if (n_notes>0) &
         print'(/,A,I0,A)', "There were ",n_notes," notes in this run"
        if (n_errors>0) &
         print'(/,A,I0,A)', "There were ",n_errors," warnings in this run"

    end subroutine summary_alerts

end module alerts
