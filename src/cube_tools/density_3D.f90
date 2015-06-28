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
    real(8),dimension(1:100000)       :: Xs, Ys, Zs !, q    
    integer                           :: Nat, Ns
    !================

    !=============
    !Counters
    integer :: i,j,k, ii
    !=============

    !================
    !I/O stuff 
    integer :: I_INP=10,  &
               O_CUB=20
    !files
    character(len=200):: inpfile , &
                         inpfile_bn="input.com", & 
                         outcube="out.cube" 
    integer :: ind0, indmax, ifile
                         
    !reading 
    character(len=260)  :: line
    !status
    integer :: IOstatus

    !Elec pot things
    real(8) :: xmax=-9999.d0, ymax=-9999.d0, zmax=-9999.d0, &
               x0= 9999.d0,   y0= 9999.d0,   z0= 9999.d0,   &
               d1x,d1y,d1z, &
               d2x,d2y,d2z, &
               d3x,d3y,d3z, &
               x,y,z, qq, dd, qsel
    integer :: n1,n2,n3, i1, i2, i3
    real(8),dimension(:,:,:),allocatable :: Dens
    !Binning
    integer :: bin1, bin2, bin3
    !===================

    !Defaults
    dd = 1.0
    ind0 = 0
    indmax = 1

    ! 0. GET COMMAND LINE ARGUMENTS AND OPEN FILES
    call parse_input(inpfile_bn,ind0,indmax,outcube,qsel,dd)

    !======= START WITH FILES ==================================
    do ifile=ind0,indmax
        ! 1. READ DATA
        write(line,'(I0)') ifile
        inpfile = trim(adjustl(inpfile_bn))//trim(adjustl(line))//".com"
        print*, "Reading "//trim(adjustl(inpfile))//"..."
        open(I_INP,file=inpfile,iostat=IOstatus,status="old")
        if (IOstatus /= 0) then
            write(0,*) "Error reading input file"
            stop
        endif
        do
            read(I_INP,'(A)') line
            if (index(line,"#") == 1) exit
        enddo
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
        !And now the get the atoms (selected based on their charge)
        i=0
        do 
            read(I_INP,'(A)') line
            if (len_trim(line) == 0) exit
            read(line,*) x, y, z, qq
            if ( qq == qsel ) then
                i=i+1
                Xs(i) = x
                Ys(i) = y
                Zs(i) = z
            endif
        enddo
        Ns = i

        close(I_INP)

        !Get AtNums
        call AtName2AtNum(Nat,AtName,AtNum)

        ! 3. Take only relevant molecules (those in L layer)
        j=0
        k=0
        do i=1,Nat
            !Change units
            Xat(i) = Xat(i) / BOHRtoANGS
            Yat(i) = Yat(i) / BOHRtoANGS
            Zat(i) = Zat(i) / BOHRtoANGS
            xmax = max(xmax,Xat(i))
            ymax = max(ymax,Yat(i))
            zmax = max(zmax,Zat(i))
            x0   = min(x0,Xat(i))
            y0   = min(y0,Yat(i))
            z0   = min(z0,Zat(i))
        enddo
        do i=1,Ns
            !Change units
            Xs(i) = Xs(i) / BOHRtoANGS
            Ys(i) = Ys(i) / BOHRtoANGS
            Zs(i) = Zs(i) / BOHRtoANGS
        enddo

        !For the first frame, get limits and allocate Dens
        if (ifile == ind0) then
            !Set the limits of the box using a buffer of 3AA
            x0   = x0 - 6.d0/BOHRtoANGS
            xmax = xmax + 6.d0/BOHRtoANGS
            y0   = y0 - 1.5d0/BOHRtoANGS
            ymax = ymax + 1.5d0/BOHRtoANGS
            z0   = z0 - 6.d0/BOHRtoANGS
            zmax = zmax + 6.d0/BOHRtoANGS

            !For the moment we use canonical vectors
            d1x = dd
            d1y = 0.d0
            d1z = 0.d0
            d2x = 0.d0
            d2y = dd
            d2z = 0.d0
            d3x = 0.d0
            d3y = 0.d0
            d3z = dd

            n1 = int((xmax-x0)/d1x)+1
            n2 = int((ymax-y0)/d2y)+1
            n3 = int((zmax-z0)/d3z)+1

            print*, "Will make a box with limits"
            print'(A,3(I0,X))', "  Size: ", n1,n2,n3
            print'(A,3F8.2)', "  X:", x0, xmax, x0+dd*(n1-1)
            print'(A,3F8.2)', "  Y:", y0, ymax, y0+dd*(n2-1)
            print'(A,3F8.2)', "  Z:", z0, zmax, z0+dd*(n3-1)
            !Redefine xmax
            xmax = x0+(n1-1)*dd
            ymax = y0+(n2-1)*dd
            zmax = z0+(n3-1)*dd

            allocate(Dens(1:n1,1:n2,1:n3),stat=IOstatus)
            if (IOstatus /= 0) then
                print*, "ERROR: Cannot allocate an array of size:", n1, n2, n3
                stop
            endif
            Dens(:,:,:) = 0.d0
        endif

        !Compute density
        do i=1,Ns
            !If molecule out of the grid, discard it
            if (Xs(i)-dd/2.d0 <= x0 .or. Xs(i) >= xmax+dd/2.d0 .or. &
                Ys(i)-dd/2.d0 <= y0 .or. Ys(i) >= ymax+dd/2.d0 .or. &
                Zs(i)-dd/2.d0 <= z0 .or. Zs(i) >= zmax+dd/2.d0) cycle

            ! Fill 3D bins
            bin1 = int(anint(((Xs(i)-x0-dd/2.d0)/(xmax-x0+dd)) * float(n1))) + 1
            bin2 = int(anint(((Ys(i)-y0-dd/2.d0)/(ymax-y0+dd)) * float(n2))) + 1
            bin3 = int(anint(((Zs(i)-z0-dd/2.d0)/(zmax-z0+dd)) * float(n3))) + 1

            if (bin1 < 1 .or. bin1 > n1) print*, "Incorrect bin1", bin1
            if (bin2 < 1 .or. bin2 > n2) print*, "Incorrect bin2", bin2
            if (bin3 < 1 .or. bin3 > n3) print*, "Incorrect bin3", bin3

            Dens(bin1,bin2,bin3) = Dens(bin1,bin2,bin3) + 1.d0

        enddo
    enddo
    Dens = Dens/float(indmax-ind0+1)/dd**3*1.d3
    !======= END WITH FILES ==================================

    !Print file
    open(O_CUB,file=outcube)
    write(O_CUB,*) "File generated from "//trim(adjustl(inpfile))
    write(O_CUB,*) "3D density map"
    write(O_CUB,'(I5,4F12.6)') Nat,x0,y0,z0
    write(O_CUB,'(I5,4F12.6)') n1, d1x,d1y,d1z
    write(O_CUB,'(I5,4F12.6)') n2, d2x,d2y,d2z
    write(O_CUB,'(I5,4F12.6)') n3, d3x,d3y,d3z
    do i=1,Nat
        write(O_CUB,'(I5,4F12.6)') AtNum(i), float(AtNum(i)),Xat(i),Yat(i),Zat(i)
    enddo     
    do i1=1,n1
    do i2=1,n2
        write(O_CUB,'(6F12.6)') Dens(i1,i2,1:n3)
    enddo
    enddo
    close(O_CUB)

    stop

    contains

    subroutine AtName2AtNum(Nat,AtName,AtNum)

        integer,intent(in)                       :: Nat
        character(len=*),dimension(:),intent(in) :: AtName
        integer,dimension(:),intent(out)         :: AtNum
        !Local
        integer :: iel
        character(len=5),dimension(103) :: atom_names_from_atnum
 
        !This should be elsewhere (constants_mod?)
        data atom_names_from_atnum(1:103) &
         /'H' ,                                                                                'He',&
          'Li','Be',                                                  'B' ,'C' ,'N' ,'O' ,'F' ,'Ne',&
          'Na','Mg',                                                  'Al','Si','P' ,'S' ,'Cl','Ar',&
          'K' ,'Ca','Sc','Ti','V' ,'Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',&
          'Rb','Sr','Y' ,'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I' ,'Xe',&
          'Cs','Ba','La',& !Lantanides:  
!                  ---------------------------------------------------
                    'Ce','Pr','Nd','Pm','Sm','Eu','Gd',&
                    'Tb','Dy','Ho','Er','Tm','Yb','Lu',&
!                  ---------------------------------------------------
                         'Hf','Ta','W' ,'Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn',&
         'Fr','Ra','Ac',& !Actinides:
!                  ---------------------------------------------------
                   'Th','Pa','U' ,'Np','Pu','Am','Cm',&
                   'Bk','Cf','Es','Fm','Md','No','Lr'&
!                  ---------------------------------------------------
         /

        AtNum(1:Nat) = 0
        do i=1,Nat
            do iel=1,103
                if (adjustl(AtName(i)) == &
                    adjustl(atom_names_from_atnum(iel))) then
                    AtNum(i) = iel
                    exit
                endif
            enddo
        enddo

        return

    end subroutine

    subroutine parse_input(inpfile_bn,ind0,indmax,outfile,qsel,delta)

        implicit none

        character(len=*),intent(inout) :: inpfile_bn,outfile
        real(8),intent(inout) :: delta, qsel
        integer,intent(inout) :: ind0,indmax
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
                case ("-fbn") 
                    call getarg(i+1, inpfile_bn)
                    argument_retrieved=.true.

                case ("-f0")
                    call getarg(i+1, arg)
                    argument_retrieved=.true.
                    read(arg,*) ind0

                case ("-fm")
                    call getarg(i+1, arg)
                    argument_retrieved=.true.
                    read(arg,*) indmax

                case ("-o") 
                    call getarg(i+1, outfile)
                    argument_retrieved=.true.

                case ("-delta")
                    call getarg(i+1, arg)
                    argument_retrieved=.true.
                    read(arg,*) delta

                case ("-qsel")
                    call getarg(i+1, arg)
                    argument_retrieved=.true.
                    read(arg,*) qsel

                case ("-h")
                    need_help=.true.

                case default
                    write(0,*) "Error reading input file"//adjustl(arg)
                    stop
            end select
        enddo 
          

       !Print options (to stderr)
        write(0,'(/,A)') '--------------------------------------------------'
        write(0,'(/,A)') '               3D-Density '               
        write(0,'(/,A)') '--------------------------------------------------'
        write(0,*) '-fbn            ', trim(adjustl(inpfile_bn))
        write(0,'(X,A,I0)') '-f0             ', ind0
        write(0,'(X,A,I0)') '-fm             ', indmax
        write(0,*) '-o              ', trim(adjustl(outfile))
        write(0,'(X,A,F8.3)') '-delta          ', delta
        write(0,'(X,A,F8.3)') '-qsel           ', qsel
        write(0,*) '-h             ',  need_help
        write(0,*) '--------------------------------------------------'
        if (need_help) stop

        return
    end subroutine parse_input

end program gen_oniom

