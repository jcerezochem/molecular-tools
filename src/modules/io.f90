module io

    !==============================================================
    ! This code is part of MOLECULAR_TOOLS 
    !==============================================================
    ! Description
    !  Subroutine customizing io tasks, as open/read/write files (TODO)
    !  read command line options
    !  This is also a good place to store general output things
    !  such as the type of output (stdout or stderr)
    !
    ! Dependences
    !
    ! Notes
    !  The module also includes some variables (acting as
    !  a common block for them)
    !-----------------------------------------------------

    contains

    subroutine get_input(i,arg,input_command)

        integer,intent(in)                      :: i
        character(len=*),intent(out)            :: arg 
        character(len=*),intent(inout),optional :: input_command 

        call getarg(i, arg)

        if (present(input_command)) &
         input_command=trim(adjustl(input_command))//" "//trim(adjustl(arg))

        return

    end subroutine get_input


end module io
