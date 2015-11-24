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

    integer, save :: n_notes=0,&
                     n_errors=0
    ! Silent will mute the NOTES
    logical, save :: silent_notes = .true.

    contains

    subroutine alert_msg(attype,SENTENCE)

        implicit none

        character(len=*),intent(in):: attype, SENTENCE

        select case (adjustl(attype))
            case ("note")
                if (.not.silent_notes) &
                write(0,'(/,A,A,/)') "NOTE: ",SENTENCE
                n_notes = n_notes + 1

            !The following is deprecated (maintained for backward compatibility)
            case ("error")
                write(0,'(/,A,A,/)') "WARNING: ",SENTENCE
                n_errors = n_errors + 1

            case ("warning")
                write(0,'(/,A,A,/)') "WARNING: ",SENTENCE
                n_errors = n_errors + 1

            case ("fatal")
                write(0,'(/,A,/)') "============================================"
                write(0,'(A,A,/)') "FATAL ERROR: ",trim(SENTENCE)
                write(0,'(A)')     "============================================"
                write(0,'(A,/)') "Exiting..."
                stop

            case default
                write(0,'(/,A,A,A,/)') attype," ", SENTENCE
        end select

        return

    end subroutine alert_msg

end module alerts
