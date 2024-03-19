program rotate_cube


    !==============================================================
    ! This code uses of MOLECULAR_TOOLS (version 1.0/March 2012)
    !==============================================================
    !
    ! Description:
    ! -----------
    ! Forms a gro structure with the one present in the fchk
    !
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

    implicit none

    !====================== 
    !System variables
    type(str_resmol) :: molec
    !====================== 

    real(8),dimension(:),allocatable :: X,Y,Z, Chrg
    integer,dimension(:),allocatable :: AtNum
    
    integer :: Nat, N1,N2,N3
    real(8) :: x0,  y0,  z0,   &
               xr,  yr,  zr,   &
               u1x,u1y,u1z, &
               u2x,u2y,u2z, &
               u3x,u3y,u3z

    !=============
    !Counters and dummies
    integer :: i,j,k,l, jj,kk, iat
    character(len=1) :: null
    logical :: overwrite    = .false. ,&
               make_connect = .false. ,&
               use_elements = .false. ,&
               remove_com   = .false.
    !Swap related counters
    integer :: iat_new, iat_orig, nswap, nbonds
    !Rotation/translation
    real(8),dimension(3)   :: Tr
    real(8),dimension(3,3) :: R
    character(len=100) :: rotcenter_str='COM'
    !=============

    !================
    !I/O stuff 
    !units
    integer :: I_INP=10, &
               I_SWP=11, &
               I_ROT=12, &
               I_TRA=13, &
               O_OUT=20, &
               O_XYZ=21
    !files
    character(len=200) :: title1,title2
    character(len=10) :: filetype_inp="guess",&
                         filetype_out="guess"
    character(len=200):: inpfile="input.cube",&
                         outfile="default"   ,&
                         rotfile="rot.dat"
    !status
    integer :: IOstatus
    character(len=7) :: stat="new" !do not overwrite when writting
    !===================
    
    ! 0. GET COMMAND LINE ARGUMENTS
    call parse_input(inpfile,filetype_inp,outfile,filetype_out,overwrite,&
                     rotfile,rotcenter_str)
 
    ! 1. READ INPUT (cube file)
    ! ---------------------------------
    open(I_INP,file=inpfile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(inpfile)) )

    if (adjustl(filetype_inp) == "guess") call split_line_back(inpfile,".",null,filetype_inp)
    ! Read cube
    if (filetype_inp /= 'cube') then
        call alert_msg('fatal','Only cube files supported')
    endif
    read(I_INP,'(A)') title1
    read(I_INP,'(A)') title2
    read(I_INP,*) Nat,  x0,    y0,   z0
    call allocate_atoms(molec,Nat)
    allocate(X(Nat),Y(Nat),Z(Nat),Chrg(Nat),AtNum(Nat))
    read(I_INP,*) N1, u1x,u1y,u1z
    read(I_INP,*) N2, u2x,u2y,u2z
    read(I_INP,*) N3, u3x,u3y,u3z
    molec%natoms = Nat
    do i=1,Nat
        read(I_INP,*) AtNum(i), Chrg(i), X(i), Y(i), Z(i)
        molec%atom(i)%AtNum = AtNum(i)
        molec%atom(i)%Mass  = atmass_from_atnum(AtNum(i))
        molec%atom(i)%Name  = atname_from_atnum(AtNum(i))
        molec%atom(i)%x = X(i)*BOHRtoANGS
        molec%atom(i)%y = Y(i)*BOHRtoANGS
        molec%atom(i)%z = Z(i)*BOHRtoANGS
    enddo

    ! 2. ROTATE
    ! -------------------------------
    ! Mulecule
    if (adjustl(rotfile)/="none") then
        print'(X,A)', "ROTATING:"
        open(I_ROT,file=rotfile,status="old")
        ! Read rotation matrix.
        do i=1,3
            read(I_ROT,*) R(i,1:3)
        enddo
        close(I_ROT)
        print'(X,A)', "  STRUCTURE..."
        if (rotcenter_str == 'COM') then
            print'(X,A)', "  (from the Center of Mass)"
            call get_com(molec)
            Tr = (/molec%comX,molec%comY,molec%comZ/)
            print'(X,3F8.2)',  Tr(1:3)
        elseif (rotcenter_str == 'COG') then
            print'(X,A)', "  (from the Center of Geometry)"
            call get_cog(molec)
            Tr = (/molec%comX,molec%comY,molec%comZ/)
            print'(X,3F8.2)',  Tr(1:3)
        else
            print'(X,A)', "  (from reference point: "//trim(adjustl(rotcenter_str))//")"
            call string2vector(rotcenter_str,Tr,i,',')
            if ( i/=3 ) call alert_msg('fatal','Reference rotation point is not well defined')
        endif
        print'(3F8.3)', molec%atom(1)%x, molec%atom(1)%y, molec%atom(1)%z
        call translate_molec(molec,-Tr)
        print'(3F8.3)', molec%atom(1)%x, molec%atom(1)%y, molec%atom(1)%z
        call rotate_molec(molec,R)
        print'(3F8.3)', molec%atom(1)%x, molec%atom(1)%y, molec%atom(1)%z
        call translate_molec(molec,Tr)
        print'(3F8.3)', molec%atom(1)%x, molec%atom(1)%y, molec%atom(1)%z
        print'(X,A,/)', "Done"
    endif
    ! Cube grid
    !------------
    ! Initial point of the grid
    x0 = x0 - Tr(1)
    y0 = y0 - Tr(2)
    z0 = z0 - Tr(3)
    ! Xr = R X0
    xr = R(1,1)*x0 + R(1,2)*y0 + R(1,3)*z0
    yr = R(2,1)*x0 + R(2,2)*y0 + R(2,3)*z0
    zr = R(3,1)*x0 + R(3,2)*y0 + R(3,3)*z0
    ! Reset center of rotation
    x0 = xr + Tr(1)
    y0 = yr + Tr(2)
    z0 = zr + Tr(3)
    !
    ! Grid vectors
    ! u1
    xr = R(1,1)*u1x + R(1,2)*u1y + R(1,3)*u1z
    yr = R(2,1)*u1x + R(2,2)*u1y + R(2,3)*u1z
    zr = R(3,1)*u1x + R(3,2)*u1y + R(3,3)*u1z
    u1x = xr 
    u1y = yr 
    u1z = zr 
    ! u2
    xr = R(1,1)*u2x + R(1,2)*u2y + R(1,3)*u2z
    yr = R(2,1)*u2x + R(2,2)*u2y + R(2,3)*u2z
    zr = R(3,1)*u2x + R(3,2)*u2y + R(3,3)*u2z
    u2x = xr 
    u2y = yr 
    u2z = zr 
    ! u3
    xr = R(1,1)*u3x + R(1,2)*u3y + R(1,3)*u3z
    yr = R(2,1)*u3x + R(2,2)*u3y + R(2,3)*u3z
    zr = R(3,1)*u3x + R(3,2)*u3y + R(3,3)*u3z
    u3x = xr 
    u3y = yr 
    u3z = zr 

    ! 3. WRITE OUTPUT
    ! ---------------------------------
    if (overwrite) stat="unknown"
    open(O_OUT,file=outfile,status=stat,iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Cannot write in "//trim(adjustl(outfile))//&
                                                            ". Already exists? Use -ow to overwrite" )
    open(O_XYZ,file='str.xyz')
    write(O_XYZ,*) Nat
    write(O_XYZ,*) 'Structure'

    write(O_OUT,'(A)') trim(title1)
    write(O_OUT,'(A)') trim(title2)
    write(O_OUT,'(I5,3F12.6,I5)') Nat, x0,y0,z0, 1
    write(O_OUT,'(I5,3F12.6)')    N1,  u1x,u1y,u1z
    write(O_OUT,'(I5,3F12.6)')    N2,  u2x,u2y,u2z
    write(O_OUT,'(I5,3F12.6)')    N3,  u3x,u3y,u3z
    do i=1,Nat
        X(i) = molec%atom(i)%x/BOHRtoANGS
        Y(i) = molec%atom(i)%y/BOHRtoANGS
        Z(i) = molec%atom(i)%z/BOHRtoANGS
        write(O_OUT,'(I5,4F12.6)') AtNum(i), Chrg(i), X(i), Y(i), Z(i)
        write(O_XYZ,'(I5,3F12.6)') AtNum(i), X(i)*BOHRtoANGS, Y(i)*BOHRtoANGS, Z(i)*BOHRtoANGS
    enddo
    !
    close(O_OUT)
    close(O_XYZ)

    stop


    !==============================================
    contains
    !=============================================

    subroutine parse_input(inpfile,filetype_inp,outfile,filetype_out,overwrite,&
                           rotfile,rotcenter_str)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: inpfile,outfile,&
                                          filetype_inp,filetype_out, &
                                          rotfile,rotcenter_str
        logical,intent(inout) :: overwrite
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
                case ("-fti") 
                    call getarg(i+1, filetype_inp)
                    argument_retrieved=.true.

                case ("-o") 
                    call getarg(i+1, outfile)
                    argument_retrieved=.true.
                case ("-fto") 
                    call getarg(i+1, filetype_out)
                    argument_retrieved=.true.
                    
                case ("-rot") 
                    call getarg(i+1, rotfile)
                    argument_retrieved=.true.
                    
                case ("-refrot")
                    call getarg(i+1, rotcenter_str)
                    argument_retrieved=.true.

                case ("-r")
                    overwrite=.true.

                case ("-ow")
                    overwrite=.true.
        
                case ("-h")
                    need_help=.true.

                case default
                    call alert_msg("fatal","Unkown command line argument: "//adjustl(arg))
            end select
        enddo 

        !The default output is now xyz
        if (adjustl(outfile) == "default") then
            call split_line_back(inpfile,".",outfile,null)
            outfile=trim(adjustl(outfile))//"_rot.cube"
        endif

        ! Some checks on the input
        !----------------------------
        !Print options (to stderr)
        write(0,'(/,A)') '========================================================'
        write(0,'(/,A)') '             R O T A T E   C U B E '    
        write(0,'(/,A)') '           Rotate cube files (header) '      
        call print_version()
        write(0,'(/,A)') '========================================================'
        write(0,'(/,A)') '-------------------------------------------------------------------'
        write(0,'(A)')   ' Flag         Description                      Value'
        write(0,'(A)')   '-------------------------------------------------------------------'
        write(0,*)       '-f           Input structure file             ', trim(adjustl(inpfile))
        write(0,*)       '-fti         \_ FileTyep                      ', trim(adjustl(filetype_inp))
        write(0,*)       '-o           Output structure file            ', trim(adjustl(outfile))
        write(0,*)       '-fto         \_ FileTyep                      ', trim(adjustl(filetype_out))
        write(0,*)       '-ow          Overwrite output if exists       ',  overwrite
        write(0,*)       '-rot         File with rotation matrix        ',  trim(adjustl(rotfile))
        write(0,*)       '-refrot      Location of the rotation point   ',  trim(adjustl(rotcenter_str))
        write(0,*)       '-h           This help                       ',  need_help
        write(0,*)       '-------------------------------------------------------------------'
        if (need_help) call alert_msg("fatal", 'There is no manual (for the moment)' )

        return
    end subroutine parse_input


end program rotate_cube
