program gen_oniom

    !==============================================================
    ! This code uses of MOLECULAR_TOOLS (version 1.0/March 2012)
    !==============================================================
    !
    ! Description:
    ! -----------
    !
    ! Compilation instructions (for mymake script):
    !make$ gfortran ../modules/alerts.f90 ../modules/structure_types_v2.f90 ../modules/line_preprocess.f90 ../modules/gro_manage_v3.f90 ../modules/gaussian_manage_v2.f90 ../modules/pdb_manage_v2.f90 ../modules/molpro_manage_v2.f90 gen_oniom_v3.f90 -o gen_oniom_v3.exe -cpp
    !
    !============================================================================    

    implicit none

    !constants
    real(8),parameter :: BOHRtoANGS= 5.2917720859D-1

    !================
    !Structural data
    character(len=5),dimension(1:400) :: AtName
    real(8),dimension(1:400)          :: Xat, Yat, Zat
    integer,dimension(1:400)          :: AtNum
    real(8),dimension(1:100000)       :: Xq, Yq, Zq, q    
    integer                           :: Nat, Nq
    !================

    !=============
    !Counters
    integer :: i,j,k, ii
    !=============

    !================
    !I/O stuff 
    integer :: I_INP=10
    !files
    character(len=200):: inpfile="input.gro"
    !reading 
    character(len=260)  :: line
    !status
    integer :: IOstatus

    !Elec pot things
    real(8) :: r, pot
    !===================

    ! 0. GET COMMAND LINE ARGUMENTS AND OPEN FILES
    call parse_input(inpfile)

    ! 1. READ DATA
    open(I_INP,file=inpfile,iostat=IOstatus,status="old")
    if (IOstatus /= 0) then
        write(0,*) "Error reading input file"
        stop
    endif
    do
        read(I_INP,'(A)') line
        if (index(line,"#") == 1) exit
    enddo
    if (index(line,"Charge")/=0.or.index(line,"charge")/=0) then
        continue
    else
        write(0,*) "Not a charge calculation"
        stop
    endif
    !Read blank
    read(I_INP,'(A)') line
    do
        read(I_INP,'(A)') line
        if (len_trim(line) == 0) exit
    enddo
    read(I_INP,'(A)') line
    !Now read the structure
    i=0
    do 
        read(I_INP,'(A)') line
        if (len_trim(line) == 0) exit
        i=i+1
        read(line,*) AtName(i), Xat(i), Yat(i), Zat(i)
    enddo
    Nat = i
    read(I_INP,'(A)') line
    !And now the charges
    i=0
    do 
        read(I_INP,'(A)') line
        if (len_trim(line) == 0) exit
        i=i+1
        read(line,*) Xq(i), Yq(i), Zq(i), q(i)
    enddo
    Nq = i
    
    ! 3. Take only relevant molecules (those in L layer)
    j=0
    k=0
    do i=1,Nat
        !Change units
        Xat(i) = Xat(i) / BOHRtoANGS
        Yat(i) = Yat(i) / BOHRtoANGS
        Zat(i) = Zat(i) / BOHRtoANGS
    enddo
    do i=1,Nq
        !Change units
        Xq(i) = Xq(i) / BOHRtoANGS
        Yq(i) = Yq(i) / BOHRtoANGS
        Zq(i) = Zq(i) / BOHRtoANGS
    enddo


    !Print file
    do i=1,Nat
        pot=0.d0
        do j=1,Nq
            r = dsqrt((Xat(i)-Xq(j))**2+&
                      (Yat(i)-Yq(j))**2+&
                      (Zat(i)-Zq(j))**2)
            pot = pot+q(j)/r
        enddo
        print*, pot
    enddo     

    stop

    contains

    subroutine parse_input(inpfile)

        implicit none

        character(len=*),intent(inout) :: inpfile
        ! Local
        logical :: need_help = .false., &
                   argument_retrieved
        integer :: i
        character(len=200) :: arg

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

                case ("-h")
                    need_help=.true.

                case default
                    write(0,*) "Error reading input file"//adjustl(arg)
                    stop
            end select
        enddo 
          

       !Print options (to stderr)
        write(0,'(/,A)') '--------------------------------------------------'
        write(0,'(/,A)') '              Charges_Potential '               
        write(0,'(/,A)') '--------------------------------------------------'
        write(0,*) '-f              ', trim(adjustl(inpfile))
        write(0,*) '-h             ',  need_help
        write(0,*) '--------------------------------------------------'
        if (need_help) stop

        return
    end subroutine parse_input

end program gen_oniom

