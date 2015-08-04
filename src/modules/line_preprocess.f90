module line_preprocess

    !==============================================================
    ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
    !==============================================================
    !
    ! Despription
    ! -----------
    ! Subroutines to parse and manage lines
    !
    ! History
    ! -------
    ! 29/11/13: split_line_back: changed to store the whole thing
    !           in line_b in case there is no match (as opposed to
    !           to normal split_line SR
    ! 13/02/14: added change case subroutines
    ! 17/02/14  added "string2vector" subroutine
    ! 27/02/14  added "string2vector_char" (there also is "string2vector_int")
    !==============================================================


    implicit none

    contains

   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine split_line(line,splitter,line_a,line_b)

        !Split a line from a given marker. If it is not present, it does not
        !split the line (the whole is preserved in line_a

        character(len=*),intent(in):: line,splitter
        character(len=*),intent(out):: line_a,line_b

        !local
        integer :: i,j
        !Auxiliar helps when line(input) is also one 
        !of the outputs, line_a or line_b
        character(len=(len(line_a))) :: aux_line_a

        i=INDEX(line,splitter)
        if ( i == 0 ) then
            line_a=line
            line_b=""
            return
        endif
        j=len_trim(splitter)
        
        aux_line_a=line(1:i-1)
        line_b=line(i+j:)
        line_a=aux_line_a

        return

    end subroutine split_line

   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine split_line_back(line,splitter,line_a,line_b)

        !Split a line from a given marker. If it is not present, it does not
        !split the line (the whole is preserved in line_a >> change: for back
        ! this is stored in line_b (so that if no match, the whole thing is in line_b
        ! the important thing for this SR

        character(len=*),intent(in):: line,splitter
        character(len=*),intent(out):: line_a,line_b

        !local
        integer :: i,j
        !Auxiliar helps when line(input) is also one 
        !of the outputs, line_a or line_b
        character(len=(len(line_a))) :: aux_line_a

        !INDEX with BACK=.true., search match from the end of the string (useful to get file extensions)
        i=INDEX(line,splitter,.true.)
        if ( i == 0 ) then
            line_a=""
            line_b=line
            return
        endif
        j=len_trim(splitter)
        
        aux_line_a=line(1:i-1)
        line_b=line(i+j:)
        line_a=aux_line_a

        return

    end subroutine split_line_back

   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine parse_line(line,narg,arg)

        !======================================================  
        ! Description
        !  Divides a line in the words that are separated by blank spaces
        ! Arguments
        !  line   INPUT   String to be evaluated
        !  narg   OUTPUT  Number of words
        !  arg(:) OUTPUT  String vector with all the words
        !======================================================  
    
        character(len=*),intent(in) :: line
        integer,intent(out) :: narg
        character(len=*),dimension(:),intent(out) :: arg
        !Local
        character(len=len_trim(line)) :: line_pp, args
        integer :: nchar, i, id

        line_pp = adjustl(line)
        nchar=len(line_pp)

        id=1
        args=""
        do i=1, nchar
            if ( line(i:i) == " " .and. args /= "" ) then
                arg(id)=adjustl(args)
                args=""
                id=id+1
                cycle
            elseif ( line(i:i) == " " .and. args == "" ) then
                cycle
            endif
            
            args=trim(args)//line(i:i)

        enddo
        arg(id)=adjustl(args)
        narg=id
            
        return

    end subroutine parse_line


    subroutine read_list_int(record,nelements,int_array)

        !Read an array of integers of unknownn size stored in from a character string

        character(len=*),intent(inout) :: record
        integer,intent(out) :: nelements
        integer,dimension(:),intent(out) :: int_array

        !local
        integer :: imax, i

        nelements=0
        record=adjustl(record)
        imax=len_trim(record)
        !Eliminate extra blanks
        do i=1,imax
            if (len_trim(record(i:i)) == 0 .and. len_trim(record(i-1:i-1)) == 0) cycle
            if (len_trim(record(i:i)) == 0) nelements = nelements+1
        enddo
        nelements = nelements+1
        read(record,*) int_array(1:nelements)

        return

    end subroutine read_list_int

    subroutine set_upper_case(letter)

        character(len=1),intent(inout) :: letter

        if ( (ichar(letter) >= 97) .and. &
             (ichar(letter) <= 122) ) then
            letter = char(ichar(letter)-32)
        endif
        
        return

    end subroutine set_upper_case

    subroutine set_lower_case(letter)

        character(len=1),intent(inout) :: letter

        if ( (ichar(letter) >= 65) .and. &
             (ichar(letter) <= 90) ) then
            letter = char(ichar(letter)+32)
        endif
        
        return

    end subroutine set_lower_case


    subroutine string2vector(raw_vector,array_vector,n_elem)

        !Description
        ! Tranforms a string of comma sepparated values into an
        ! array of such real vaues 

        character(len=*),intent(inout) :: raw_vector
#ifdef DOUBLE
        double precision,dimension(:),intent(out) :: array_vector
#else
        real,dimension(:),intent(out) :: array_vector
#endif
        integer,intent(out) :: n_elem

        !Local
        character(len=240) :: auxchar
        integer :: i
    
        
        !Read unknown length vector (comma sepparated)
        i=0
        do 
            i=i+1
            if ( INDEX(raw_vector,',') /= 0 ) then
                call split_line(raw_vector,',',auxchar,raw_vector)
                read(auxchar,*) array_vector(i)
            else 
                read(raw_vector,*) array_vector(i)
                exit
            endif
        enddo  
        n_elem=i

        return

    end subroutine string2vector

    subroutine string2vector_int(raw_vector,array_vector,n_elem)

        !Description
        ! Tranforms a string of comma sepparated values into an
        ! array of such integer vaues (integer version) 

        character(len=*),intent(inout) :: raw_vector
        integer,dimension(:),intent(out) :: array_vector
        integer,intent(out) :: n_elem

        !Local
        character(len=240) :: auxchar
        integer :: i
    
        
        !Read unknown length vector (comma sepparated)
        i=0
        do 
            i=i+1
            if ( INDEX(raw_vector,',') /= 0 ) then
                call split_line(raw_vector,',',auxchar,raw_vector)
                read(auxchar,*) array_vector(i)
            else 
                read(raw_vector,*) array_vector(i)
                exit
            endif
        enddo  
        n_elem=i

        return

    end subroutine string2vector_int

    subroutine string2vector_char(raw_vector,array_vector,n_elem,sep)

        !Description
        ! Tranforms a string of comma sepparated values into an
        ! array of such integer vaues (integer version) 

        character(len=*),intent(inout) :: raw_vector
        character(len=*),dimension(:),intent(out) :: array_vector
        integer,intent(out) :: n_elem
        character(len=*) :: sep !separador

        !Local
        character(len=240) :: auxchar
        integer :: i
    
        
        !Read unknown length vector (comma sepparated)
        i=0
        do 
            i=i+1
            if ( INDEX(raw_vector,sep) /= 0 ) then
                call split_line(raw_vector,sep,auxchar,raw_vector)
                read(auxchar,*) array_vector(i)
            else 
                read(raw_vector,*) array_vector(i)
                exit
            endif
        enddo  
        n_elem=i

        return

    end subroutine string2vector_char


    ! Functions that get character from numbers

    function int20char(i,length) result(c)

        integer,intent(in)    :: i
        integer,intent(in)    :: length
        character(len=length) :: c
        ! Local
        character(len=10) :: fmt
        integer           :: ilength, j
        character(len=10) :: dummy_char

        !If length<number of digits, rise an error
        ilength = int(log10(float(i)))+1
        if (ilength>length) then
            write(0,*) "Error in int20char: more digits in number than character size"
            stop
        endif
        ! otherwise fill head with zeroes
        dummy_char = ""
        do j=1,length-ilength
            dummy_char = trim(dummy_char)//"0"
        enddo

        ! Write format
        write(fmt,'(a,i0,a)') '(A,I',ilength,')'
        write(c,fmt) trim(dummy_char), i

        return

    end function int20char

    function int2char(i,length) result(c)

        integer,intent(in)    :: i
        integer,intent(in)    :: length
        character(len=length) :: c
        ! Local
        character(len=10) :: fmt
        integer           :: ilength

        !If length<number of digits, rise an error
        ilength = int(log10(float(i)))+1
        if (ilength>length) then
            write(0,*) "Error in int2char: more digits in number than character size"
            stop
        endif

        ! Write format 
        write(fmt,'(a,i0,a)') '(I',length,')'
        write(c,fmt) i

        return

    end function int2char

    function real2char(r,length,decimals) result(c)

        real(8),intent(in)    :: r
        integer,intent(in)    :: length, decimals
        character(len=length) :: c
        ! Local
        character(len=10) :: fmt

        ! Write format 
        write(fmt,'(a,i0,a,i0,a)') '(F',length,'.',decimals,')'
        write(c,fmt) r

        return

    end function real2char


end module line_preprocess
