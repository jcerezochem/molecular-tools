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

    use line_preprocess
    use verbosity

    integer,save :: uout=6

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

    subroutine heading(unt,title,print_always)

        integer,intent(in)          :: unt
        character(len=*),intent(in) :: title
        logical,intent(in),optional :: print_always
        ! local
        character(len=len_trim(title)) :: line
        logical :: print_now

        if (present(print_always) .and. print_always) then
            print_now = .true.
        else
            print_now = .false.
        endif

        line=adjustl(title)
        call set_word_upper_case(line)

        if (verbose == 0 .and. .not. print_now) return

        write(unt,'(2/,X,A)') "==============================================="
        write(unt,'(2X,A)' )  line
        write(unt,'(X,A)')   "==============================================="

        return

    end subroutine heading

    subroutine subheading(unt,title,keep_case,upper_case)

        integer,intent(in)          :: unt
        character(len=*),intent(in) :: title
        logical,intent(in),optional :: keep_case
        logical,intent(in),optional :: upper_case
        ! local
        character(len=len_trim(title)) :: line
        integer :: i
        logical :: first_letter

        if (verbose == 0) return

        if (present(keep_case) .and. keep_case) then
            line=adjustl(title)
        elseif (present(upper_case) .and. upper_case) then
            line=adjustl(title)
            call set_word_upper_case(line)
        else
            line=adjustl(title)
            first_letter = .true.
            do i=1,len_trim(line)
                if (len_trim(line(i:i)) == 0) then
                    first_letter=.true.
                else if (first_letter) then
                    call set_upper_case(line(i:i))
                    first_letter=.false.
                else
                    call set_lower_case(line(i:i))
                endif
            enddo
        endif
            
        if (verbose == 0) return

        write(unt,'(/,X,A)') "------------------------------------------"
        write(unt,'(2X,A)' )  line
        write(unt,'(X,A)')   "------------------------------------------"

        return

    end subroutine subheading

    subroutine subsubheading(unt,title,keep_case,upper_case)

        integer,intent(in)          :: unt
        character(len=*),intent(in) :: title
        logical,intent(in),optional :: keep_case
        logical,intent(in),optional :: upper_case
        ! local
        character(len=len_trim(title)) :: line
        integer :: i
        logical :: first_letter

        if (verbose == 0) return

        if (present(keep_case) .and. keep_case) then
            line=adjustl(title)
        elseif (present(upper_case) .and. upper_case) then
            line=adjustl(title)
            call set_word_upper_case(line)
        else
            line=adjustl(title)
            first_letter = .true.
            do i=1,len_trim(line)
                if (len_trim(line(i:i)) == 0) then
                    first_letter=.true.
                else if (first_letter) then
                    call set_upper_case(line(i:i))
                    first_letter=.false.
                else
                    call set_lower_case(line(i:i))
                endif
            enddo
        endif
            
        if (verbose == 0) return

        write(unt,'(/,2X,A)' )  line
        write(unt,'(X,A)')     "------------------------------------------"

        return

    end subroutine subsubheading

    subroutine statement(unt,title,keep_case,upper_case)

        integer,intent(in)          :: unt
        character(len=*),intent(in) :: title
        logical,intent(in),optional :: keep_case
        logical,intent(in),optional :: upper_case
        ! local
        character(len=len_trim(title)) :: line
        integer :: i
        logical :: first_letter
        
        if (verbose == 0) return

        if (present(keep_case) .and. keep_case) then
            line=adjustl(title)
        elseif (present(upper_case) .and. upper_case) then
            line=adjustl(title)
            call set_word_upper_case(line)
        else
            line=adjustl(title)
            first_letter = .true.
            do i=1,len_trim(line)
                if (len_trim(line(i:i)) == 0) then
                    first_letter=.true.
                else if (first_letter) then
                    call set_upper_case(line(i:i))
                    first_letter=.false.
                else
                    call set_lower_case(line(i:i))
                endif
            enddo
        endif

        write(unt,'(3X,A,/)' )  line

        return

    end subroutine statement

end module io
