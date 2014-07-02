program test_read_fchk

    !Compilation
    ! $FC ../modules/alerts.f90 ../modules/line_preprocess.f90 ../modules/gaussian_fchk_manage_v2.f90 ../modules/gaussian_manage_vNoTypes.f90 read_state_fchk.f90 -o read_state_fchk.exe -cpp -DDOUBLE 
    ! THIS IS A SELF CONTAINED VERSION (no needs molecular tools)

    !NOTES
    ! Normal modes matrix indeces: T(1:Nvib,1:3Nat)
    ! This might be the contrary as the usual convention
    ! (anywaym who cares, since we use a Tvector)

!     use gaussian_fchk_manage
!     use alerts
!     use line_preprocess

    implicit none

    real(8),parameter :: BOHRtoAMS = 5.2917720859D-1
    !CONVERSIONS
    double precision, parameter :: uma_to_Kg=1.66053873d-27,  &
                                   H_to_J=4.3597482d-18,      &
                                   bohr_to_m=5.291772083d-11, &
                                   A_to_m=1d-10

    !CONSTANTS
    double precision, parameter :: pi=4.d0*datan(1.d0),       &
                                   c_vac=2.99792458d8,        &
                                   plank=6.62606957d-34   

    !Interesting info..
    integer :: Nat, Nvib
    real(8),dimension(1:1000) :: GEOM, Gvector
    real(8),dimension(:),allocatable :: Tvector
    real(8),dimension(1:1000,1:1000) :: T
    

    !Auxiliars
    character(len=100) :: section
    integer :: N, N_T
    character :: dtype, cnull
    real(8),dimension(:),allocatable :: A
    integer,dimension(:),allocatable :: I
    integer :: error
    real(8) :: Norm
    !Counters
    integer :: j,k,l

    !I/O
    character(len=100) :: fchkfile, outfile
    integer :: ios
    integer :: I_FCHK=10,&
               O_STA =20

    ! Get filenames from commandline or read standard input/output
    if ( iargc() == 0 ) then
        I_FCHK = 5
        O_STA  = 6
    else
        call getarg(1, fchkfile)
        open(I_FCHK,file=fchkfile,status="old",iostat=ios)
        if ( ios /= 0 ) call alert_msg("fatal","could not open the file "//trim(adjustl(fchkfile)))
    endif

    !INFORMATION IN THE FCHK FILE
    ! Natoms
    call read_fchk(I_FCHK,"Number of atoms",dtype,N,A,I,error)
    if (error == 0) then
        Nat = I(1)
        deallocate(I)
    endif

    ! Nvib
    call read_fchk(I_FCHK,"Number of Normal Modes",dtype,N,A,I,error)
    if (error == 0) then
        Nvib = I(1)
        deallocate(I)
    else
        print*, "You likely missed 'SaveNM' keyword"
        print*, "Estimating Nvib assuming non-linear molecule"
        Nvib = 3*Nat - 6
    endif

    ! Geom
    call read_fchk(I_FCHK,"Current cartesian coordinates",dtype,N,A,I,error)
    if (error == 0) then
        GEOM(1:N) = A*BOHRtoAMS
        deallocate(A)
    endif

    ! L_cart matrix
    !from the fchk a vector is read where T(Nvib,3N) is stored as
    ! T(1,1:3N), T(2,1:3N), ..., T(Nvib,1:3N)
    ! but FCclasses needs a vector where the fast index is the first:
    ! T(1:Nvib,1), ..., T(1:Nvib,3N)
    call read_fchk(I_FCHK,"Vib-Modes",dtype,N,A,I,error)
    if (error == 0) then
        allocate( Tvector(1:N) )
        N_T = N
        Tvector = A
        deallocate(A)
    else
        print*, "You likely missed 'SaveNM' keyword"
        print*, "Try to get it from the log file..."
        print*, "Setting them to 1.000. The file is not valid for AH!"
        N=Nvib*3*Nat
        N_T = N
        allocate( Tvector(1:N) )
        Tvector(1:N) = 1.d0
    endif
    ! Reconstruct Lcart (non-symmetric)
    l=0
    do j=1,Nvib
        do k=1,3*Nat
            l=l+1
            T(j,k) = Tvector(l)
         enddo
    enddo
    ! Rebuild Tvector in FCclasses compliant way
    l=0
    do j=1,3*Nat
        do k=1,Nvib
            l=l+1
            Tvector(l) = T(k,j)
         enddo
    enddo

    call read_fchk(I_FCHK,"Vib-E2",dtype,N,A,I,error)
    if (error == 0) then
        Gvector(1:Nvib) = A(1:Nvib)
        deallocate(A)
    else
        print*, "You likely missed 'SaveNM' keyword"
        print*, "Try to get it from the log file..."
        print*, "Setting them to 1.000. The file is not valid for AH!"
        Gvector(1:Nvib) = 1.d0
    endif
    close(I_FCHK)


    !WRITE STATE FILE
    if (O_STA /= 6) then
        call split_line(fchkfile,".",outfile,cnull)
        outfile = "state_"//trim(adjustl(outfile))//"_fchk"
        open(O_STA,file=outfile,status="replace")
    endif
!     write(O_STA,*) "NEW"
    do j=1,3*Nat
!         write(O_STA,'(E20.10)') GEOM(j)
        write(O_STA,*) GEOM(j)
    enddo
    do j=1,N_T
!         write(O_STA,'(E20.10)') Tvector(j)
        write(O_STA,*) Tvector(j)
    enddo
    do j=1,Nvib
        write(O_STA,'(F10.4)') Gvector(j)
    enddo
    close(O_STA)


    stop

    contains

    subroutine alert_msg(attype,SENTENCE)

        implicit none

        common /ALERT_COUNT/ n_notes, n_errors

        character(len=*),intent(in):: attype, SENTENCE
        integer :: n_notes, n_errors

        select case (adjustl(attype))
            case ("note")
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

    subroutine read_fchk(unt,section,data_type,N,A,I,error_flag)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 1.0/March 2012)
        !==============================================================
        !Description
        ! Generic SR to read any section of the checkpoint
        ! Enter deallocated arrays
        !Arguments
        ! unt (int;in): unit number of the log file
        ! section(char,in): name of the section to be read
        ! data_type(char,out): Integer (I) or Real (R) data read
        ! A(real,dimension(:)): Real array to store real data
        ! I(integer,dimension(:)): Int array to store int data
        ! error_flag(integer,out): 0: success
        !                          1: section not found
        !==============================================================

        integer,intent(in) :: unt
        character(len=*),intent(in) :: section
        character(len=1),intent(out) :: data_type
        integer,intent(out) :: N
! #ifdef DOUBLE
        double precision, dimension(:), allocatable, intent(out) :: A
! #else
!         real, dimension(:), allocatable, intent(out) :: A
! #endif
        integer,dimension(:), allocatable, intent(out) :: I
        integer,intent(out) :: error_flag

        !Local stuff
        !=============
        character(len=240) :: line=""
        character(len=42) :: section_full
        character(len=1) :: is_array
        character(len=40) :: cdata
        !I/O
        integer :: IOstatus
        
        
        ! Search section
        error_flag = 0
        do
                read(unt,'(A)',IOSTAT=IOstatus) line
                ! Two possible scenarios while reading:
                ! 1) End of file
                if ( IOstatus < 0 ) then
                    call alert_msg("note","Section '"//trim(adjustl(section))//"' not present in the FCHK file.")
                    error_flag=1
                    rewind(unt)
                    return
                endif
                ! 2) Found what looked for!      
                if ( INDEX(line,trim(adjustl(section))) /= 0 ) then
                    read(line,'(A42)') section_full
                    if (adjustl(section_full) == adjustl(section)) exit
                endif
        enddo

        !Get info about section from line just read
        read(line,'(A42,X,A1,3X,A1,X,A)') section_full, data_type, is_array, cdata
        if (is_array /= "N") then
            !Is not an array
            N=1
            if ( data_type == "R" ) then 
                allocate( A(1:1) )
                read(cdata,*) A(1)
            elseif ( data_type == "I" ) then
                allocate( I(1:1) )
                read(cdata,*) I(1) 
            endif
        else
            read(cdata,*) N
            if ( data_type == "R" ) then
                allocate( A(1:N) )
                read(unt,*) A(1:N)
            elseif ( data_type == "I" ) then
                allocate( I(1:N) )
                read(unt,*) I(1:N)
            endif
        endif 

        rewind(unt)
        return

    end subroutine read_fchk

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

end program test_read_fchk

