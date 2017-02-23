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
    real(8),dimension(1:400)          :: Xat, Yat, Zat, qlist
    integer,dimension(1:400)          :: AtNum
    real(8),dimension(1:100000)       :: Xq, Yq, Zq, q   
    integer                           :: Nat, Nq
    !================

    !=============
    !Counters
    integer :: i,j,k, ii
    integer :: ires, iold
    !=============

    !================
    !I/O stuff 
    integer :: I_INP=10
    !files
    character(len=200):: inpfile="input.gro"
    !reading 
    character(len=260)  :: line
    character :: null
    !status
    integer :: IOstatus

    !Elec pot things
    real(8) :: r, dq19, dq20, q1, q18, q19, q20
    real(8),dimension(1:40) :: pot
    !===================

    ! 0. GET COMMAND LINE ARGUMENTS AND OPEN FILES
    call parse_input(inpfile)

    ! 1. READ DATA
    open(I_INP,file=inpfile,iostat=IOstatus,status="old")
    if (IOstatus /= 0) then
        write(0,*) "Error reading input file"
        stop
    endif

    qlist(1:9) = (/-0.23710,  0.05950,  0.05950,  0.05950,  0.51200, -0.08060, -0.08060, -0.69740,  0.40520/) 


    do

        read(I_INP,*,iostat=IOstatus) line
        if (IOstatus /= 0) exit

        Nat=36
        if (adjustl(line) == "POSITION" ) then
            i=0
            do i=1,Nat
                read(I_INP,'(A)') line
                if (adjustl(line) == "END" ) exit
                read(line,*)  ii  ,                  & !resseq
                              null,                  & !resname
                              null,                  & 
                              ii,                    & !serial
                              Xat(i),                  &
                              Yat(i),                  &
                              Zat(i)
            enddo
            !Transform to BOHR
            Xat(1:Nat)=Xat(1:Nat)*10.d0/ BOHRtoANGS
            Yat(1:Nat)=Yat(1:Nat)*10.d0/ BOHRtoANGS
            Zat(1:Nat)=Zat(1:Nat)*10.d0/ BOHRtoANGS
            !Read charges
            i = 0 
            iold = 1
            do
                read(I_INP,'(A)') line
                if (adjustl(line) == "END" ) exit
                i=i+1
                read(line,*)  ires,                  & !resseq
                              null,                  & !resname
                              null,                  & 
                              ii,                    & !serial
                              Xq(i),                  &
                              Yq(i),                  &
                              Zq(i)
                if (ires /= iold) k=0 
                k=k+1
                q(i)=qlist(k)
                iold=ires
            enddo
            Nq = i
            !Transform to BOHR
            Xq(1:Nq)=Xq(1:Nq)*10.d0/ BOHRtoANGS
            Yq(1:Nq)=Yq(1:Nq)*10.d0/ BOHRtoANGS
            Zq(1:Nq)=Zq(1:Nq)*10.d0/ BOHRtoANGS
        endif

    enddo

    !Print file
    do i=1,Nat
        pot(i)=0.d0
        do j=1,Nq
            r = dsqrt((Xat(i)-Xq(j))**2+&
                      (Yat(i)-Yq(j))**2+&
                      (Zat(i)-Zq(j))**2)
            pot(i) = pot(i)+q(j)/r
        enddo
    enddo     

    dq19 = -1.1*pot(19) + 0.04
    dq20 = -1.1*pot(20) + 0.04

    q19 = -0.569172 + dq19
    q20 = -0.611033 + dq20
    q18 =  0.755627 - dq19
    q1  =  0.769597 - dq20

    print'(4F15.6)', q1, q18, q19, q20

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
          

!        !Print options (to stderr)
!         write(0,'(/,A)') '--------------------------------------------------'
!         write(0,'(/,A)') '              Charges_Potential '               
!         write(0,'(/,A)') '--------------------------------------------------'
!         write(0,*) '-f              ', trim(adjustl(inpfile))
!         write(0,*) '-h             ',  need_help
!         write(0,*) '--------------------------------------------------'
!         if (need_help) stop

        return
    end subroutine parse_input

end program gen_oniom

