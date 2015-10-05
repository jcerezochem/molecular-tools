program normal_modes_Cartesian

    !Compilation
    ! $FC ../modules/constants_mod.f90 ../modules/alerts.f90 ../modules/line_preprocess.f90 ../modules/structure_types_v2.f90 ../modules/gro_manage_v2.f90 ../modules/gaussian_fchk_manage_v2.f90 ../modules/geom_meter_v2.f90 ../modules/xyz_manage.f90 normal_modes_Cartesian_v1b.f90 -o normal_modes_Cartesian_v1b.exe -cpp -DDOUBLE 

    !NOTES
    ! Normal modes matrix indeces: T(1:Nvib,1:3Nat)
    ! This might be the contrary as the usual convention
    ! (anywaym who cares, since we use a Tvector)

    !V1b: from normal_mode_animation_v1b. Naming systematic in compliance with internal version
    !
    !Addapted to v4 release (distribution upgrade). 25/02/14
    !Additional Changes
    ! -Added filetype support with generic_strfile_read 


!*****************
!   MODULE LOAD
!*****************
!============================================
!   Generic (structure_types independent)
!============================================
    use alerts
    use line_preprocess
    use constants
!   Matrix manipulation (i.e. rotation matrices)
    use MatrixMod
!============================================
!   Structure types module
!============================================
    use structure_types
!============================================
!   Structure dependent modules
!============================================
    use gro_manage
    use pdb_manage
    use gaussian_manage
    use gaussian_manage_notypes
    use gaussian_fchk_manage
    use xyz_manage
    use fcc_manage
!   Structural parameters
    use molecular_structure
    use ff_build
!   Bond/angle/dihed meassurement
    use atomic_geom

    type(str_resmol) :: molec

    !Interesting info..
    integer :: Nat, Nvib
    real(8),dimension(1:1000) :: GEOM, RedMass, Freq, Factor
    real(8),dimension(1:1000,1:1000) :: T
    character(len=2) :: atname 
    integer          :: Zat

    real(8) :: Amplitude, qcoord
    

    !Auxiliars
    character(1) :: null
    character(len=50) :: dummy_char
    real(8) :: Tmwc_ij
    ! Reading FCHK part
    character(len=100) :: section
    integer :: N, N_T
    character :: dtype, cnull
    real(8),dimension(:),allocatable :: A
    integer,dimension(:),allocatable :: IA
    integer :: error


    !Moving normal modes
    integer,dimension(1:1000) :: nm=0
    real(8) :: Qstep
    logical :: call_vmd = .false.
    character(len=10000) :: vmdcall

    !Counters
    integer :: i,j,k, jdh,idh,istep, iat, kk

    !I/O
    integer :: I_INP=10,  &
               O_GRO=20,  &
               O_G09=21,  & 
               O_Q  =22,  &
               O_NUM=23,  &
               O_LIS=24,  &
               S_VMD=30
    !files
    character(len=20) :: filetype="guess"
    character(len=20) :: addfile=""
    character(len=200):: inpfile ="input.fchk"
    character(len=100),dimension(1:1000) :: grofile
    character(len=100) :: nmfile='none', g09file,qfile,numfile,tmpfile
    !Control of stdout
    logical :: verbose=.false.

    !MOVIE things
    integer :: movie_cycles, movie_steps

    !===============================================

    !Defaults
    Amplitude = 1.d0
    movie_cycles = 0 !No movie

    !Get options from command line
    nm(1) = 0
    call parse_input(inpfile,addfile,nmfile,nm,Nsel,Amplitude,filetype,verbose,call_vmd,movie_cycles)


    ! 1. CARTESIAN VIBRATIONAL ANALYSIS -- done by G09 (just read it) 
 
    ! 1. READ DATA
    ! ---------------------------------
    open(I_INP,file=inpfile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(inpfile)) )

    !Read structure
    ! If fcc state file, we also need the input to get natoms, nvib and atom info
    if (filetype == "fcc") then
        !Need to read also the input file
        if (len_trim(addfile) == 0) then
            print*, "Need to set the fcc input file with the hidded -add flag"
            stop
        endif
        !-------- read input ------------
        open(I_ADD,file=addfile)
        read(I_ADD,*) molec%natoms
        read(I_ADD,*) Nvib
        do i=1,molec%natoms
            read(I_ADD,*) molec%atom(i)%mass
        enddo
        close(I_ADD)
        !----------------------------------
    endif
    call generic_strfile_read(I_INP,filetype,molec)
    !Shortcuts
    Nat = molec%natoms
    if (filetype /= "fcc") then
        !Caution if linear...
        Nvib = 3*Nat-6
    endif

    if (adjustl(filetype) == "log") then
        !Gaussian logfile (directly read vibronic analysis by Gaussian)
        N = 3*Nat * Nvib
        allocate( A(1:N) )
        call read_freq_NT(I_INP,Nvib,Nat,Freq(1:Nvib),RedMass(1:Nvib),A(1:N),error)
        if (error == 1) then
            call alert_msg("warning", "Normal modes (T) matrix in low precision. This will lead to poor results")
        elseif (error == -1) then
            call alert_msg("fatal","No frequency information in the log file")
        endif
        ! Reconstruct Lcart (non-symmetric) [in T]
        l=0
        do j=1,Nvib
            do k=1,3*molec%natoms
                l=l+1
                T(j,k) = A(l)
            enddo
        enddo
        deallocate(A)
       !Reduced masses (to account for the normalization factor in Tcart)
        RedMass(1:Nvib) = RedMass(1:Nvib)*UMAtoAU
!         !In the future, that would be good to perform the vib analysis here, instead
!         allocate(props)
!         call parse_summary(I_INP,molec,props,"read_hess")
!         ! RECONSTRUCT THE FULL HESSIAN
!         k=0
!         do i=1,3*Nat
!             do j=1,i
!                 k=k+1
!                 Hess(i,j) = props%H(k) 
!                 Hess(j,i) = Hess(i,j)
!             enddo
!         enddo
!         deallocate(props)zz
    !Read the Hessian: only two possibilities supported
    else if (adjustl(filetype) == "fcc") then
        N = 3*Nat * Nvib
        allocate( A(1:N) )
        call read_fccstate_freq(I_INP,Nvib,Nat,Freq(1:Nvib),A(1:N))
        ! Reconstruct Lcart (non-symmetric) [in T]
        l=0
        do k=1,3*molec%natoms
            do j=1,Nvib
                l=l+1
                T(j,k) = A(l)
            enddo
        enddo
        deallocate(A)
        !Compute reduced masses
        do k=1,Nvib
            kk=0
            RedMass(k)=0.d0
            do i=1,Nat
            do j=1,3
                kk=kk+1
                Tmwc_ij=T(k,kk)*dsqrt(molec%atom(i)%mass)
                RedMass(k)=RedMass(k)+Tmwc_ij**2
            enddo
            enddo
        enddo
        !Reduced masses (to account for the normalization factor in Tcart)
        RedMass(1:Nvib) = RedMass(1:Nvib)*UMAtoAU
        !Get atom names from masses
        do i=1,molec%natoms
            call atominfo_from_atmass(molec%atom(i)%mass,Zat,atname)
            molec%atom(i)%name    = atname
            molec%atom(i)%element = atname
            molec%atom(i)%AtNum   = Zat
        enddo

    else if (adjustl(filetype) == "fchk") then
        !FCHK file: READ VIB ANALYSIS
        ! L_cart matrix
        !from the fchk a vector is read where T(Nvib,3N) is stored as
        ! T(1,1:3N), T(2,1:3N), ..., T(Nvib,1:3N)
        ! but FCclasses needs a vector where the fast index is the first:
        ! T(1:Nvib,1), ..., T(1:Nvib,3N)
        call read_fchk(I_INP,"Vib-Modes",dtype,N,A,IA,error)
        if (error == 0) then
            ! Reconstruct Lcart (non-symmetric) [in T]
            l=0
            do j=1,Nvib
                do k=1,3*molec%natoms
                    l=l+1
                    T(j,k) = A(l)
                enddo
            enddo
            deallocate(A)
        else
            print*, "You likely missed 'SaveNM' keyword"
            print*, "Run 'freqchk -save' utility and come back"
            stop
        endif

        call read_fchk(I_INP,"Vib-E2",dtype,N,A,IA,error)
        if (error == 0) then
           !Reduced masses (to account for the normalization factor in Tcart)
            RedMass(1:Nvib) = A(Nvib+1:2*Nvib)*UMAtoAU

            !freq based factor (to make displacements equivalent in dimensionless units
            Freq(1:Nvib) = A(1:Nvib)
            deallocate(A)
        else
            print*, "You likely missed 'SaveNM' keyword"
            print*, "Run 'freqchk -save' utility and come back"
            stop
        endif
        close(I_INP)
    endif

    !Transoform "Gaussian normalized" L'cart to "real" Lcart
    !   Lcart = L'cart/sqrt(mu)
    do j=1,Nvib
        !Note is stored as the transpose, i.e. T(Nvib x 3Nat)
        !     RedMass in atomic units
        T(j,:) = T(j,:)/dsqrt(RedMass(j))/1.895d0
        !NOTE: the factor 1.89 has been empirically identified as necessary, but its origin is
        !      unclear. This should be carrefully revised!
    enddo

    !Define the Factor to convert shift into addimensional displacements
    ! from the shift in SI units:
    Factor(1:Nvib) = dsqrt(dabs(Freq(1:Nvib))*1.d2*clight*2.d0*PI/plankbar)
    ! but we need it from au not SI
    Factor(1:Nvib)=Factor(1:Nvib)*BOHRtoM*dsqrt(AUtoKG)
 

    !NORMAL MODES SELECTION SWITCH
    if (Nsel == 0) then
        !Then select them all
        Nsel = Nvib
        do i=1,Nsel
            nm(i) = i
        enddo
    endif
    if (Nsel > 1000) call alert_msg("fatal", "Too many normal modes. Dying")


    !=========================================================
    !MOVE NORMAL MODES (deltaQ = 0.01)
    !=========================================================
    ! This stuff should be introduced from command line
    ! Displacement need to be scaled for each mode (how??)
!     Qstep = 5.d-3
!     Nsteps = int(Amplitude*1.d-1/Qstep)
    !Do the opposite: fix Nsteps and modiffy the displacment
    Nsteps = 100
!     Qstep = Amplitude*20.d0/float(Nsteps)*249.8859d0
    Qstep = Amplitude/float(Nsteps)
    if ( mod(Nsteps,2) /= 0 ) Nsteps = Nsteps + 1
    do jj=1,Nsel
        kk=0
        qcoord=0.d0
        molec%atom(1:molec%natoms)%resname = "RES"
        j = nm(jj)
        write(grofile(jj),*) j
        molec%title = "Animation of normal mode "//trim(adjustl(grofile(jj)))
        g09file="Mode"//trim(adjustl(grofile(jj)))//"_cart.com" 
        numfile="Mode"//trim(adjustl(grofile(jj)))//"_cart_num.com"
        qfile="Mode"//trim(adjustl(grofile(jj)))//"_cart_steps.dat"
        grofile(jj) = "Mode"//trim(adjustl(grofile(jj)))//"_cart.gro"
        print*, "Writting results to: "//trim(adjustl(grofile(jj)))
        open(O_GRO,file=grofile(jj))
        open(O_Q  ,file=qfile)
        open(O_G09,file=g09file)

        !===========================
        !Start from equilibrium. 
        !===========================
        print*, "STEP:", kk
        call write_gro(O_GRO,molec)
        !===========================
        !Half Forward oscillation
        !===========================
        do istep = 1,nsteps/2
            kk=kk+1
            qcoord = qcoord + Qstep/Factor(j)
            print*, "STEP:", kk
            write(dummy_char,*) kk
            molec%title = trim(adjustl(grofile(jj)))//".step "//trim(adjustl(dummy_char))
            do i = 1, molec%natoms
                k=3*(i-1) + 1
                molec%atom(i)%x = molec%atom(i)%x + T(j,k)   * Qstep /Factor(j)
                molec%atom(i)%y = molec%atom(i)%y + T(j,k+1) * Qstep /Factor(j)
                molec%atom(i)%z = molec%atom(i)%z + T(j,k+2) * Qstep /Factor(j)
                !Note: itcould be reduced-mass weighted by mult * dsqrt(Freq(j)
            enddo
            call write_gro(O_GRO,molec)
            !Write the max amplitude step to G09 scan
            if (kk==nsteps/2) then
                call write_gcom(O_G09,molec)
                write(O_Q,*) qcoord, qcoord*Factor(j)
            endif
        enddo
        !=======================================
        ! Reached amplitude. Back oscillation
        !=======================================
        do istep = 1,nsteps/2-3
            kk=kk+1
            qcoord = qcoord - Qstep/Factor(j)
            print*, "STEP:", kk
            write(dummy_char,*) kk
            molec%title = trim(adjustl(grofile(jj)))//".step "//trim(adjustl(dummy_char))
            do i = 1, molec%natoms
                k=3*(i-1) + 1
                molec%atom(i)%x = molec%atom(i)%x - T(j,k) * Qstep   /Factor(j)
                molec%atom(i)%y = molec%atom(i)%y - T(j,k+1) * Qstep /Factor(j)
                molec%atom(i)%z = molec%atom(i)%z - T(j,k+2) * Qstep /Factor(j)
                !Note: itcould be reduced-mass weighted by mult * dsqrt(Freq(j)
            enddo
            call write_gro(O_GRO,molec)
            if (mod(kk,10) == 0) then
                molec%job%title=trim(adjustl(grofile(jj)))//".step "//trim(adjustl(dummy_char))
                molec%title=trim(adjustl(g09file))
                call write_gcom(O_G09,molec)
                write(O_Q,*) qcoord, qcoord*Factor(j)
            endif
        enddo
        !=======================================
        ! Reached equilibrium again (5 points for numerical second derivatives)
        !=======================================
        open(O_NUM,file=numfile,status="replace")
        do istep = 1,5
            kk=kk+1
            qcoord = qcoord - Qstep/Factor(j)
            print*, "STEP:", kk
            write(dummy_char,*) kk
            molec%title = trim(adjustl(grofile(jj)))//".step "//trim(adjustl(dummy_char))
            do i = 1, molec%natoms
                k=3*(i-1) + 1
                molec%atom(i)%x = molec%atom(i)%x - T(j,k) * Qstep   /Factor(j)
                molec%atom(i)%y = molec%atom(i)%y - T(j,k+1) * Qstep /Factor(j)
                molec%atom(i)%z = molec%atom(i)%z - T(j,k+2) * Qstep /Factor(j)
                !Note: itcould be reduced-mass weighted by mult * dsqrt(Freq(j)
            enddo
            call write_gro(O_GRO,molec)
            !This time write all five numbers
            molec%job%title=trim(adjustl(grofile(jj)))//".step "//trim(adjustl(dummy_char))
            molec%title=trim(adjustl(g09file))
            call write_gcom(O_G09,molec)
            write(dummy_char,*) qcoord
            molec%job%title = "Displacement = "//trim(adjustl(dummy_char))
            molec%title=trim(adjustl(numfile))
            call write_gcom(O_NUM,molec)
            write(O_Q,*) qcoord, qcoord*Factor(j)
        enddo
        !=======================================
        ! Continue Back oscillation
        !=======================================
        do istep = 1,nsteps/2-2
            kk=kk+1
            qcoord = qcoord - Qstep/Factor(j)
            print*, "STEP:", kk
            write(dummy_char,*) kk
            molec%title = trim(adjustl(grofile(jj)))//".step "//trim(adjustl(dummy_char))
            do i = 1, molec%natoms
                k=3*(i-1) + 1
                molec%atom(i)%x = molec%atom(i)%x - T(j,k) * Qstep   /Factor(j)
                molec%atom(i)%y = molec%atom(i)%y - T(j,k+1) * Qstep /Factor(j)
                molec%atom(i)%z = molec%atom(i)%z - T(j,k+2) * Qstep /Factor(j)
                !Note: itcould be reduced-mass weighted by mult * dsqrt(Freq(j)
            enddo
            call write_gro(O_GRO,molec)
            if (mod(kk,10) == 0) then
                call write_gcom(O_G09,molec)
                write(O_Q,*) qcoord, qcoord*Factor(j)
            endif
        enddo
        !=======================================
        ! Reached amplitude. Half Forward oscillation (till equilibrium)
        !=======================================
        do istep = 1,nsteps/2-1
            kk=kk+1
            print*, "STEP:", kk
            do i = 1, molec%natoms
                k=3*(i-1) + 1
                molec%atom(i)%x = molec%atom(i)%x + T(j,k) * Qstep   /Factor(j)
                molec%atom(i)%y = molec%atom(i)%y + T(j,k+1) * Qstep /Factor(j)
                molec%atom(i)%z = molec%atom(i)%z + T(j,k+2) * Qstep /Factor(j)
                !Note: itcould be reduced-mass weighted by mult * dsqrt(Freq(j)
            enddo
            call write_gro(O_GRO,molec)
        enddo
        close(O_GRO)
        close(O_G09)
        close(O_Q)
        close(O_NUM)
    enddo


    if (call_vmd) then
        open(S_VMD,file="vmd_conf.dat",status="replace")
        !Set general display settings (mimic gv)
        write(S_VMD,*) "color Display Background iceblue"
        write(S_VMD,*) "color Name {C} silver"
        write(S_VMD,*) "axes location off"
        !Set molecule representation
        do i=0,Nsel-1
            j = nm(i+1)
            write(S_VMD,*) "mol representation CPK"
            write(S_VMD,*) "mol addrep 0"
            if (i==0) then
                write(S_VMD,*) "molinfo ", i, " set drawn 1"
            else
                write(S_VMD,*) "molinfo ", i, " set drawn 0"
            endif
            write(S_VMD,*) "mol addrep ", i
            write(dummy_char,'(A,I4,X,F8.2,A)') "{Mode",j, Freq(j),"cm-1}"
            dummy_char=trim(adjustl(dummy_char))
            write(S_VMD,*) "mol rename ", i, trim(dummy_char)
            vmdcall = trim(adjustl(vmdcall))//" "//trim(adjustl(grofile(i+1)))
        enddo
        write(S_VMD,*) "display projection Orthographic"
        close(S_VMD)
        !Call vmd
        vmdcall = 'vmd -m '
        do i=1,Nsel
        vmdcall = trim(adjustl(vmdcall))//" "//trim(adjustl(grofile(i)))
        enddo
        vmdcall = trim(adjustl(vmdcall))//" -e vmd_conf.dat"
        call system(vmdcall)
    endif

    if (movie_cycles > 0) then
        open(S_VMD,file="vmd_movie.dat",status="replace")
        !Set general display settings (mimic gv)
        write(S_VMD,*) "color Display Background white"
        write(S_VMD,*) "color Name {C} silver"
        write(S_VMD,*) "axes location off"
        write(S_VMD,*) "display projection Orthographic"
        !Set molecule representation
        do i=0,Nsel-1
            j = nm(i+1)
            write(S_VMD,*) "mol representation CPK"
            write(S_VMD,*) "molinfo ", i, " set drawn 0"
            write(S_VMD,*) "mol addrep ", i
            write(dummy_char,'(A,I4,X,F8.2,A)') "{Mode",j, Freq(j),"cm-1}"
            dummy_char=trim(adjustl(dummy_char))
            write(S_VMD,*) "mol rename ", i, trim(dummy_char)
            vmdcall = trim(adjustl(vmdcall))//" "//trim(adjustl(grofile(i+1)))
        enddo
        write(S_VMD,'(A)') "#====================="
        write(S_VMD,'(A)') "# Start movies"
        write(S_VMD,'(A)') "#====================="
        !Set length of the movie
        movie_steps = movie_cycles*20
        do i=0,Nsel-1
            j = nm(i+1)
            write(S_VMD,'(A,I4)') "# Mode", j
            write(tmpfile,*) j
            tmpfile="Mode"//trim(adjustl(tmpfile))
            write(S_VMD,*) "molinfo ", i, " set drawn 1"
            write(S_VMD,*) "set figfile "//trim(adjustl(tmpfile))
            write(S_VMD,'(A,I3,A)') "for {set xx 0} {$xx <=", movie_steps,&
                                    "} {incr xx} {"
            write(S_VMD,*) "set x [expr {($xx-($xx/20)*20)*10}]"
            write(S_VMD,*) 'echo "step $x"'
            write(S_VMD,*) "animate goto $x"
            write(S_VMD,*) "render Tachyon $figfile-$xx.dat"
            write(S_VMD,'(A)') '"/usr/local/lib/vmd/tachyon_LINUX" -aasamples 12 '//& 
                           '$figfile-$xx.dat -format TARGA -o $figfile-$xx.tga'
            write(S_VMD,'(A)') 'convert -font URW-Palladio-Roman -pointsize 30 -draw '//&
                           '"text 30,70 '//"'"//trim(adjustl(tmpfile))//"'"//&
                           '" $figfile-$xx.tga $figfile-$xx.jpg'
            write(S_VMD,'(A)') "}"
            !Updated ffmpeg call. The output is now loadable from ipynb
            write(S_VMD,'(A)') 'ffmpeg -i $figfile-%d.jpg -vcodec libx264 -s 640x360 $figfile.mp4'
            write(S_VMD,*) "molinfo ", i, " set drawn 0"
        enddo
        write(S_VMD,*) "exit"
        close(S_VMD)
        !Call vmd
        vmdcall = 'vmd -m '
        do i=1,Nsel
        vmdcall = trim(adjustl(vmdcall))//" "//trim(adjustl(grofile(i)))
        enddo
        vmdcall = trim(adjustl(vmdcall))//" -e vmd_movie.dat -size 500 500"
        open(O_LIS,file="movie.cmd",status="replace")
        write(O_LIS,'(A)') trim(adjustl(vmdcall))
        write(O_LIS,'(A)') "rm Mode*jpg Mode*dat Mode*tga"
        close(O_LIS)
        print*, ""
        print*, "============================================================"
        print*, "TO GENERATE THE MOVIES (AVI) EXECUTE COMMANDS IN 'movie.cmd'"
        print*, "(you may want to edit 'vmd_movie.dat'  first)"
        print*, "============================================================"
        print*, ""
    endif


    stop

    !==============================================
    contains
    !=============================================

    subroutine parse_input(inpfile,addfile,nmfile,nm,Nsel,Amplitude,filetype,verbose,call_vmd,movie_cycles)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: inpfile,addfile,nmfile,filetype
        logical,intent(inout) ::  verbose, call_vmd
        integer,dimension(:),intent(inout) :: nm
        integer,intent(inout) :: Nsel, movie_cycles
        real(8),intent(out) :: Amplitude
        ! Local
        logical :: argument_retrieved,  &
                   need_help = .false.
        integer:: i
        character(len=200) :: arg
        real(8) :: maxd

        !Prelimirary defaults
        Nsel = 0

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
                case ("-ft") 
                    call getarg(i+1, filetype)
                    argument_retrieved=.true.
                case ("-add") 
                    call getarg(i+1, addfile)
                    argument_retrieved=.true.

                case ("-nm") 
                    call getarg(i+1, arg)
                    argument_retrieved=.true.
                    call string2vector_int(arg,nm,Nsel)

                case ("-nmf") 
                    call getarg(i+1, nmfile)
                    argument_retrieved=.true.

                case ("-maxd") 
                    call getarg(i+1, arg)
                    read(arg,*) maxd
                    !The whole Amplitude is twide the max displacement
                    Amplitude = maxd*2
                    argument_retrieved=.true.

                case ("-vmd")
                    call_vmd=.true.

                case ("-movie")
                    call getarg(i+1, arg)
                    read(arg,*) movie_cycles
                    argument_retrieved=.true.

                case ("-nosym")
                    write(0,*) trim(adjustl(arg))//" is not a valid option for Cartesian"
                case ("-sym")
                    write(0,*) trim(adjustl(arg))//" is not a valid option for Cartesian"

                case ("-sa")
                    write(0,*) trim(adjustl(arg))//" is not a valid option for Cartesian"
                case ("-nosa")
                    write(0,*) trim(adjustl(arg))//" is not a valid option for Cartesian"

                case ("-readz") 
                    write(0,*) trim(adjustl(arg))//" is not a valid option for Cartesian"

                case ("-zmat")
                    write(0,*) trim(adjustl(arg))//" is not a valid option for Cartesian"
                case ("-nozmat")
                    write(0,*) trim(adjustl(arg))//" is not a valid option for Cartesian"

                case ("-v")
                    verbose=.true.
        
                case ("-h")
                    need_help=.true.

                case default
                    call alert_msg("fatal","Unkown command line argument: "//adjustl(arg))
            end select
        enddo 

       !Print options (to stderr)
        write(0,'(/,A)') '--------------------------------------------------'
        write(0,'(/,A)') '            NORMAL MODES in CARTESIAN '    
        write(0,'(/,A)') '       Generator of normal mode animations' 
        write(0,'(/,A)') '         Revision: nm_cartesian-140320-1'
        write(0,'(/,A)') '--------------------------------------------------'
        write(0,*) '-f              ', trim(adjustl(inpfile))
        write(0,*) '-ft             ', trim(adjustl(filetype))
        write(0,*) '-nm            ', nm(1),"-",nm(Nsel)
!         write(0,*) '-nmf           ', nm(1:Nsel)
        write(0,*) '-vmd           ',  call_vmd
        write(0,*) '-movie (cycles)',  movie_cycles
        write(0,*) '-maxd          ',  Amplitude
        write(0,*) '-v             ', verbose
        write(0,*) '-h             ',  need_help
        write(0,*) '--------------------------------------------------'
        if (need_help) call alert_msg("fatal", 'There is no manual (for the moment)' )

        return
    end subroutine parse_input


    subroutine generic_strfile_read(unt,filetype,molec)

        integer, intent(in) :: unt
        character(len=*),intent(inout) :: filetype
        type(str_resmol),intent(inout) :: molec

        !local
        type(str_molprops) :: props

        if (adjustl(filetype) == "guess") then
        ! Guess file type
        call split_line_back(inpfile,".",null,filetype)
        select case (adjustl(filetype))
            case("gro")
             call read_gro(I_INP,molec)
             call atname2element(molec)
             call assign_masses(molec)
            case("pdb")
             call read_pdb_new(I_INP,molec)
             call atname2element(molec)
             call assign_masses(molec)
            case("log")
             call parse_summary(I_INP,molec,props,"struct_only")
             call atname2element(molec)
             call assign_masses(molec)
            case("fchk")
             call read_fchk_geom(I_INP,molec)
             call atname2element(molec)
!              call assign_masses(molec) !read_fchk_geom includes the fchk masses
            case default
             call alert_msg("fatal","Trying to guess, but file type but not known: "//adjustl(trim(filetype))&
                        //". Try forcing the filetype with -ft flag (available: log, fchk)")
        end select

        else
        ! Predefined filetypes
        select case (adjustl(filetype))
            case("gro")
             call read_gro(I_INP,molec)
             call atname2element(molec)
             call assign_masses(molec)
            case("fcc")
             call read_fccstate_atoms(I_INP,molec)
             molec%atom(1:molec%natoms)%name="X"
            case("pdb")
             call read_pdb_new(I_INP,molec)
             call atname2element(molec)
             call assign_masses(molec)
            case("log")
             call parse_summary(I_INP,molec,props,"struct_only")
             call atname2element(molec)
             call assign_masses(molec)
            case("fchk")
             call read_fchk_geom(I_INP,molec)
             call atname2element(molec)
!              call assign_masses(molec) !read_fchk_geom includes the fchk masses
            case default
             call alert_msg("fatal","File type not supported: "//filetype)
        end select
        endif


        return


    end subroutine generic_strfile_read


end program normal_modes_Cartesian

