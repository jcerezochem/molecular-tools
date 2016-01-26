module molecular_structure

    !==============================================================
    ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
    !==============================================================
    ! Description
    !  This MODULE contains subroutines to manage molecular_strucuture
    !  atributes, such as conectivity, atom/element naming, atomic 
    !  masses assignement...
    !
    ! NOTE: Now, most routines work with atom%element instead of
    !       of atom%name, although reading subroutines get atom%name
    !       Therefore, the use of atname2element is a must!
    !
    ! SUBROUTINES
    !    subroutine renum(system)
    !    subroutine sist2res(sistema,residuo)
    !    subroutine res2sist(residuo,nres,boxX,boxY,boxZ,sistema)
    !==============================================================

    !Common declarations:
    !===================
    use structure_types
    use alerts
    use line_preprocess
    use constants
    use verbosity
    implicit none

    contains

    subroutine set_geom_units(molec,units)

        type(str_resmol),intent(inout) :: molec
        character(len=*),intent(in)    :: units

        !Local
        integer :: Nat
        real(8) :: factor
        character(len=len(units)) :: units_local

        units_local = units
        call set_word_upper_case(units_local)
        call set_word_upper_case(molec%units)

        if (molec%units == units_local) return

        Nat = molec%natoms

        ! First set from input to ANGS
        select case(adjustl(molec%units))
            case("BOHR") !AU (bohr)
             factor=BOHRtoANGS
            case("NM") !
             factor=1.d-1
            case("ANGS") !
             factor=1.d0
            case default
             call alert_msg("fatal","Unknow units (input): "//molec%units)
        end select

       ! Now convert to the required units
        select case(adjustl(units_local))
            case("BOHR") 
             factor=factor/BOHRtoANGS
            case("NM") 
             factor=factor/1.d-1
            case("ANGS") 
             factor=factor/1.d0
            case default
             call alert_msg("fatal","Unknow units (output): "//units_local)
        end select

        molec%atom(1:Nat)%x = molec%atom(1:Nat)%x*factor
        molec%atom(1:Nat)%y = molec%atom(1:Nat)%y*factor
        molec%atom(1:Nat)%z = molec%atom(1:Nat)%z*factor

        molec%units = adjustl(units)

        return

    end subroutine set_geom_units


    subroutine assign_masses_molec(molec)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
        !==============================================================
        !Description
        ! Assign masses based on atom names
        !Arguments
        ! molec (str_system,inout): molecule
        !
        !Notes
        !(20/06/2012) Improved to get one or two letter symbols 
        !
        ! History
        ! 12/02/2014  Now atom%element is used instead of atom%name.
        !             With two-letter elements through "atname2element" subroutine
        !==============================================================

        use constants

        type(str_resmol),intent(inout) :: molec

        integer :: i
        character(len=2) :: atname

        do i=1,molec%natoms
            atname = adjustl(molec%atom(i)%element)
            molec%atom(i)%mass = atmass_from_atname(atname)
        enddo

        return

    end subroutine assign_masses_molec

    subroutine assign_atnum_molec(molec)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
        !==============================================================
        !Description
        ! Assign masses based on atom names
        !Arguments
        ! molec (str_system,inout): molecule
        !
        !Notes
        !Based on assign_masses_molec
        !==============================================================

        use constants

        type(str_resmol),intent(inout) :: molec

        integer :: i
        character(len=2) :: atname

        do i=1,molec%natoms
            atname = adjustl(molec%atom(i)%element)
            molec%atom(i)%Atnum = atnum_from_atname(atname)
        enddo

        return

    end subroutine assign_atnum_molec

   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine write_connect(unt,molec)

        implicit none
        integer,intent(in)             :: unt
        type(str_resmol),intent(inout) :: molec

        !Local
        integer :: i, nbonds

        write(unt,*) molec%natoms
        do i=1,molec%natoms
            nbonds=molec%atom(i)%nbonds
            write(unt,*) i, nbonds, molec%atom(i)%connect(1:nbonds)
        enddo

        return

    end subroutine write_connect

    subroutine read_connect(unt,molec)

        implicit none
        integer,intent(in)             :: unt
        type(str_resmol),intent(inout) :: molec

        !Local
        integer :: i,j, N, nbonds

        read(unt,*) N
        if (N /= molec%natoms) call alert_msg("note","Number of connections "//&
                                                     "does not match Natoms")
        molec%atom(1:molec%natoms)%nbonds = 0
        do i=1,N
            read(unt,*) j, nbonds, molec%atom(j)%connect(1:nbonds)
            molec%atom(j)%nbonds = nbonds
        enddo

        return

    end subroutine read_connect

    subroutine guess_connect(molec,inc_hbond)

        use constants

        implicit none
        type(str_resmol),intent(inout) :: molec
        logical,optional               :: inc_hbond
        !local
        real :: av_len, dist
        integer :: i,j, i_cnx
        logical :: include_hbond=.false.
        character(len=4) ::  current_units

        ! Save current_units and ensure Angs
        current_units=molec%units
        call set_geom_units(molec,"Angs")

        !InOut, better reassigned
        molec=molec

        !Set include_hbond
        if (present(inc_hbond)) include_hbond=inc_hbond

        ! If not done, get AtNum (if not yet assign)
        do i=1,molec%natoms
            if (molec%atom(i)%AtNum == 0) &
             molec%atom(i)%AtNum = atnum_from_atname(molec%atom(i)%name)
        enddo

        !Loop over all atoms. This is not optimal, distance matrix is symmetric
        !(it is not costly, anyway)
        !TODO: Change it to retrieve symmetric from already calc.
        i_cnx=0
        do i=1,molec%natoms
            molec%atom(i)%connect(:)=0
            i_cnx=0
            do j=1,molec%natoms
                if (i == j) cycle

                !Calculate the distance:
                dist = (molec%atom(i)%x - molec%atom(j)%x)**2 &
                     + (molec%atom(i)%y - molec%atom(j)%y)**2 &
                     + (molec%atom(i)%z - molec%atom(j)%z)**2 
                dist = sqrt(dist)

                !Get default bond length from database
                av_len=bond_length_db(molec%atom(i)%AtNum,molec%atom(j)%AtNum,include_hbond)
                ! Criterium: dist < av. length +10% --increased from 5% (17/12/12) due to H2O2
                ! a more sophiticated data base might include hibridization
                av_len=av_len*1.1

                if (dist < av_len) then
                    i_cnx=i_cnx+1
                    molec%atom(i)%connect(i_cnx)=j
                endif

            enddo
            molec%atom(i)%nbonds=i_cnx
        enddo

        ! Reset input units before leaving
        call set_geom_units(molec,adjustl(current_units))

        return

    end subroutine guess_connect


    function bond_length_db(iat1,iat2,inc_hbond) result(av_len)

        !=======================================================================
        !Description
        !-------------
        ! Seting up a database is not an obvious task. Here you have one try:
        ! (ok, you'll find a nice description here, soon) 
        !
        ! Pros and cons
        !---------------
        ! Pros: easy extension to external DB
        !       clarity of the implmented DB
        ! Cons: Memory waste storing DB (not that much, though)
        !
        ! History
        ! 12/02/2014  Now atom%element is used instead of atom%name.
        !             With two-letter elements through "atname2element" subroutine
        !=======================================================================

        integer,intent(in) :: iat1,iat2
        logical,optional   :: inc_hbond
        real(8)            :: av_len

        !Local
        real,dimension(1:118,1:3) :: CovRad
        integer :: i,j,k
        logical :: include_hbond

        ! DATABASE GENERATION
        CovRad(1  ,1:3) = (/0.32, -1., -1./)
        CovRad(2  ,1:3) = (/0.46, -1., -1./)
        CovRad(3  ,1:3) = (/1.33,1.24, -1./)
        CovRad(4  ,1:3) = (/1.02,0.90,0.85/)
        CovRad(5  ,1:3) = (/0.85,0.78,0.73/)
        CovRad(6  ,1:3) = (/0.75,0.67,0.60/)
        CovRad(7  ,1:3) = (/0.71,0.60,0.54/)
        CovRad(8  ,1:3) = (/0.63,0.57,0.53/)
        CovRad(9  ,1:3) = (/0.64,0.59,0.53/)
        CovRad(10 ,1:3) = (/0.67,0.96, -1./)
        CovRad(11 ,1:3) = (/1.55,1.60, -1./)
        CovRad(12 ,1:3) = (/1.39,1.32,1.27/)
        CovRad(13 ,1:3) = (/1.26,1.13,1.11/)
        CovRad(14 ,1:3) = (/1.16,1.07,1.02/)
        CovRad(15 ,1:3) = (/1.11,1.02,0.94/)
        CovRad(16 ,1:3) = (/1.03,0.94,0.95/)
        CovRad(17 ,1:3) = (/0.99,0.95,0.93/)
        CovRad(18 ,1:3) = (/0.96,1.07,0.96/)
        CovRad(19 ,1:3) = (/1.96,1.93, -1./)
        CovRad(20 ,1:3) = (/1.71,1.47,1.33/)
        CovRad(21 ,1:3) = (/1.48,1.16,1.14/)
        CovRad(22 ,1:3) = (/1.36,1.17,1.08/)
        CovRad(23 ,1:3) = (/1.34,1.12,1.06/)
        CovRad(24 ,1:3) = (/1.22,1.11,1.03/)
        CovRad(25 ,1:3) = (/1.19,1.05,1.03/)
        CovRad(26 ,1:3) = (/1.16,1.09,1.02/)
        CovRad(27 ,1:3) = (/1.11,1.03,0.96/)
        CovRad(28 ,1:3) = (/1.10,1.01,1.01/)
        CovRad(29 ,1:3) = (/1.12,1.15,1.20/)
        CovRad(30 ,1:3) = (/1.18,1.20, -1./)
        CovRad(31 ,1:3) = (/1.24,1.17,1.21/)
        CovRad(32 ,1:3) = (/1.21,1.11,1.14/)
        CovRad(33 ,1:3) = (/1.21,1.14,1.06/)
        CovRad(34 ,1:3) = (/1.16,1.07,1.07/)
        CovRad(35 ,1:3) = (/1.14,1.09,1.10/)
        CovRad(36 ,1:3) = (/1.17,1.21,1.08/)
        CovRad(37 ,1:3) = (/2.10,2.02, -1./)
        CovRad(38 ,1:3) = (/1.85,1.57,1.39/)
        CovRad(39 ,1:3) = (/1.63,1.30,1.24/)
        CovRad(40 ,1:3) = (/1.54,1.27,1.21/)
        CovRad(41 ,1:3) = (/1.47,1.25,1.16/)
        CovRad(42 ,1:3) = (/1.38,1.21,1.13/)
        CovRad(43 ,1:3) = (/1.28,1.20,1.10/)
        CovRad(44 ,1:3) = (/1.25,1.14,1.03/)
        CovRad(45 ,1:3) = (/1.25,1.10,1.06/)
        CovRad(46 ,1:3) = (/1.20,1.17,1.12/)
        CovRad(47 ,1:3) = (/1.28,1.39,1.37/)
        CovRad(48 ,1:3) = (/1.36,1.44, -1./)
        CovRad(49 ,1:3) = (/1.42,1.36,1.46/)
        CovRad(50 ,1:3) = (/1.40,1.30,1.32/)
        CovRad(51 ,1:3) = (/1.40,1.33,1.27/)
        CovRad(52 ,1:3) = (/1.36,1.28,1.21/)
        CovRad(53 ,1:3) = (/1.33,1.29,1.25/)
        CovRad(54 ,1:3) = (/1.31,1.35,1.22/)
        CovRad(55 ,1:3) = (/2.32,2.09, -1./)
        CovRad(56 ,1:3) = (/1.96,1.61,1.49/)
        CovRad(57 ,1:3) = (/1.80,1.39,1.39/)
        CovRad(58 ,1:3) = (/1.63,1.37,1.31/)
        CovRad(59 ,1:3) = (/1.76,1.38,1.28/)
        CovRad(60 ,1:3) = (/1.74,1.37, -1./)
        CovRad(61 ,1:3) = (/1.73,1.35, -1./)
        CovRad(62 ,1:3) = (/1.72,1.34, -1./)
        CovRad(63 ,1:3) = (/1.68,1.34, -1./)
        CovRad(64 ,1:3) = (/1.69,1.35,1.32/)
        CovRad(65 ,1:3) = (/1.68,1.35, -1./)
        CovRad(66 ,1:3) = (/1.67,1.33, -1./)
        CovRad(67 ,1:3) = (/1.66,1.33, -1./)
        CovRad(68 ,1:3) = (/1.65,1.33, -1./)
        CovRad(69 ,1:3) = (/1.64,1.31, -1./)
        CovRad(70 ,1:3) = (/1.70,1.29, -1./)
        CovRad(71 ,1:3) = (/1.62,1.31,1.31/)
        CovRad(72 ,1:3) = (/1.52,1.28,1.22/)
        CovRad(73 ,1:3) = (/1.46,1.26,1.19/)
        CovRad(74 ,1:3) = (/1.37,1.20,1.15/)
        CovRad(75 ,1:3) = (/1.31,1.19,1.10/)
        CovRad(76 ,1:3) = (/1.29,1.16,1.09/)
        CovRad(77 ,1:3) = (/1.22,1.15,1.07/)
        CovRad(78 ,1:3) = (/1.23,1.12,1.10/)
        CovRad(79 ,1:3) = (/1.24,1.21,1.23/)
        CovRad(80 ,1:3) = (/1.33,1.42, -1./)
        CovRad(81 ,1:3) = (/1.44,1.42,1.50/)
        CovRad(82 ,1:3) = (/1.44,1.35,1.37/)
        CovRad(83 ,1:3) = (/1.51,1.41,1.35/)
        CovRad(84 ,1:3) = (/1.45,1.35,1.29/)
        CovRad(85 ,1:3) = (/1.47,1.38,1.38/)
        CovRad(86 ,1:3) = (/1.42,1.45,1.33/)
        CovRad(87 ,1:3) = (/2.23,2.18, -1./)
        CovRad(88 ,1:3) = (/2.01,1.73,1.59/)
        CovRad(89 ,1:3) = (/1.86,1.53,1.40/)
        CovRad(90 ,1:3) = (/1.75,1.43,1.36/)
        CovRad(91 ,1:3) = (/1.69,1.38,1.29/)
        CovRad(92 ,1:3) = (/1.70,1.34,1.18/)
        CovRad(93 ,1:3) = (/1.71,1.36,1.16/)
        CovRad(94 ,1:3) = (/1.72,1.35, -1./)
        CovRad(95 ,1:3) = (/1.66,1.35, -1./)
        CovRad(96 ,1:3) = (/1.66,1.36, -1./)
        CovRad(97 ,1:3) = (/1.68,1.39, -1./)
        CovRad(98 ,1:3) = (/1.68,1.40, -1./)
        CovRad(99 ,1:3) = (/1.65,1.40, -1./)
        CovRad(100,1:3) = (/1.67, -1., -1./)
        CovRad(101,1:3) = (/1.73,1.39, -1./)
        CovRad(102,1:3) = (/1.76, -1., -1./)
        CovRad(103,1:3) = (/1.61,1.41, -1./)
        CovRad(104,1:3) = (/1.57,1.40,1.31/)
        CovRad(105,1:3) = (/1.49,1.36,1.26/)
        CovRad(106,1:3) = (/1.43,1.28,1.21/)
        CovRad(107,1:3) = (/1.41,1.28,1.19/)
        CovRad(108,1:3) = (/1.34,1.25,1.18/)
        CovRad(109,1:3) = (/1.29,1.25,1.13/)
        CovRad(110,1:3) = (/1.28,1.16,1.12/)
        CovRad(111,1:3) = (/1.21,1.16,1.18/)
        CovRad(112,1:3) = (/1.22,1.37,1.30/)
        CovRad(113,1:3) = (/1.36, -1., -1./)
        CovRad(114,1:3) = (/1.43, -1., -1./)
        CovRad(115,1:3) = (/1.62, -1., -1./)
        CovRad(116,1:3) = (/1.75, -1., -1./)
        CovRad(117,1:3) = (/1.65, -1., -1./)
        CovRad(118,1:3) = (/1.57, -1., -1./)

        ! Compute the lenght
        av_len = CovRad(iat1,1) + CovRad(iat2,1)

        ! Take O-O bond into account
        if (iat1 == 8 .and. iat2 == 8) av_len=av_len+0.8

        if (present(inc_hbond)) include_hbond=inc_hbond
        if (include_hbond) then
            ! Only between H and N/O
            if (i==1.or.j==1) then
                k=i+j
                if (k==8.or.k==9) then
                    av_len = av_len + 0.8
                endif
            endif
        endif

        return
    end function bond_length_db

    subroutine atname2element(molec)

    ! Subroutine to generate element names from atom names

        type(str_resmol),intent(inout) :: molec

        !local
        integer :: i
        character(len=3) :: dummy_char

        do i=1,molec%natoms
            call get_element_from_name(molec%atom(i))

            !Conflicting cases (for the moment not using additional info)
            if (molec%atom(i)%element == "Nx") then
                molec%atom(i)%element = "N"
                write(dummy_char,'(I3)') i 
                call alert_msg("note","NA atom name treated as nitrogen for atom "//dummy_char)
            elseif (molec%atom(i)%element == "Cx") then
                molec%atom(i)%element = "C"
                write(dummy_char,'(I3)') i 
                call alert_msg("note","CA atom name treated as carbon for atom "//dummy_char)
            endif
        enddo


        return

        contains

        subroutine get_element_from_name(atom)

            type(str_atom),intent(inout) :: atom

            !local
            character(len=5) :: atname

            !1. Process atom name 
            atname = atom%name
            atname = adjustl(atname)

            ! Sometimes, H atoms start with a number in PDB, GRO..
            if (atname(1:1) == "1" .or. &
                atname(1:1) == "2" .or. &
                atname(1:1) == "3" .or. &
                atname(1:1) == "4" .or. &
                atname(1:1) == "5" .or. &
                atname(1:1) == "6" .or. &
                atname(1:1) == "7" .or. &
                atname(1:1) == "8" .or. &
                atname(1:1) == "9") then
                if (atname(2:2) == "H" .or. &
                    atname(2:2) == "h") then
                    atom%element="H"
                    return
                else
                    !Then we simply remove the number and go on
                    atname(1:4) = atname(2:5)
                    atname(5:5) = ""
                endif
            endif

            !Set first letter to upper case
            call set_upper_case(atname(1:1))

            !First solve conflicts with one-letter elements
            select case (atname(1:1))
               !==========
                case ("H")
               !==========
                !It can be H, He, Hf (not considered the lanthanide: Ho)
                    ! We consider that:
                    !  HE is hidrogen labeled as "E" (strange, though)
                    !  He is helium
                    !  HF is hidrogen labeled as "F" (strange, though)
                    !  Hf is hafnium
                    select case (atname(2:2))
                        case ("e")
                            atom%element = "He"
                            call alert_msg("warning","He taken as helium")
                        case ("f")
                            atom%element = "Hf"
                            call alert_msg("warning","He taken as hafnium")
                        case default
                            atom%element = "H"
                            if (adjustl(atom%element) /= adjustl(atom%name) ) &
                             call alert_msg("note",trim(adjustl(atom%name))//" taken as hydrogen")
                    end select
                    return
               !==========
                case ("B")
               !==========
                !It can be B, Be, Br, Ba
                    select case (atname(2:2))
                        case ("a")
                            atom%element = "Ba"
                        case ("A")
                            atom%element = "Ba"
                            call alert_msg("note","BA taken as barium")
                        case ("e")
                            atom%element = "Be"
                        case ("E")
                            atom%element = "Be"
                            call alert_msg("note","BE taken as berium")
                        case ("r")
                            atom%element = "Br"
                        case ("R")
                            atom%element = "Br"
                            call alert_msg("note","BR taken as bromine")
                        case default
                            atom%element = "B"
                            if (adjustl(atom%element) /= adjustl(atom%name) ) &
                             call alert_msg("warning",trim(adjustl(atom%name))//" taken as borium")
                    end select
                    return
               !==========
                case ("C")
               !==========
                    !C is a nightmare... It can be C Cl Cd Ca Cr Cs (not considered the lanthanide/actinides: Ce, Cm, Cf)
                    ! We consider that:
                    !  CD is carbon labeled as "D"
                    !  Cd is Cadmium
                    !  CL and Cl are chlorine (there is not usually an "L" label)
                    !  CR and Cr are chromium (WARNING: chirality label?)
                    !  CS and Cs are cesium (WARNING: chirality label?)
                    !  CA is carbon labeled as "A" or or calcium: use more info later
                    !  Ca is calcium
                    select case (atname(2:2))
                        case ("d")
                            atom%element = "Cd"
                            call alert_msg("warning","Cd taken as cadmium")
                        case ("r")
                            atom%element = "Cr"
                            call alert_msg("warning","Cd taken as chromium")
                        case ("R")
                            atom%element = "Cr"
                            call alert_msg("warning","CR taken as chromium")
                        case ("s")
                            atom%element = "Cs"
                            call alert_msg("warning","Cs taken as cesium")
                        case ("S")
                            atom%element = "Cs"
                            call alert_msg("warning","CS taken as cesium")
                        case ("l")
                            atom%element = "Cl"
                            call alert_msg("warning","Cl taken as chlorine")
                        case ("L")
                            atom%element = "Cl"
                            call alert_msg("warning","CL taken as chlorine")
                        case ("a")
                            ! If it has additional labels (e.g. Ca1), 
                            ! this is probably not Ca but Carbon
                            if (len_trim(atname) > 2) then
                                atom%element = "C"
                                call alert_msg("warning",trim(atname)//" taken as carbon")
                            else
                                atom%element = "Ca"
                                call alert_msg("warning",trim(atname)//" taken as calcium")
                            endif
                        case ("A")
                            !This case can be either C"A" or Ca. Mark with x to check later
                            atom%element = "C"
                            call alert_msg("note","CA taken as carbone")
                        case default
                            atom%element = "C"
                            if (adjustl(atom%element) /= adjustl(atom%name) ) &
                             call alert_msg("note",trim(adjustl(atom%name))//" taken as carbone")
                    end select
                    return
               !==========
                case ("N")
               !==========
                !It can be N, Na, Ni, Nb (not considered the lanthanide/actinides: Nd, Np, No)
                    ! We consider that:
                    !  NB is carbon labeled as "B"
                    !  Nb is niobium
                    !  Ni and NI are nickel (there is not usually an "I" label)
                    !  NA is nitrogen labeled as "A" or or sodium: use more info later
                    !  Na is sodium
                    select case (atname(2:2))
                        case ("b")
                            atom%element = "Nb"
                        case ("i")
                            atom%element = "Ni"
                        case ("I")
                            atom%element = "Ni"
                        case ("a")
                            atom%element = "Na"
                        case ("A")
                            !This case can be either C"A" or Ca. Mark with x to check later
                            atom%element = "Nx"
                        case default
                            atom%element = "N"
                    end select
                    return
               !==========
                case ("O")
               !==========
                !It can be O, Os
                    ! We consider that:
                    !  OS is carbon labeled as "S" (strange, although Os is more strange)
                    !  Os is osmium
                    select case (atname(2:2))
                        case ("s")
                            atom%element = "Os"
                        case default
                            atom%element = "O"
                    end select
                    return
               !==========
                case ("F")
               !==========
                !It can be F, Fe
                    ! We consider that:
                    !  Fe and FE are iron
                    select case (atname(2:2))
                        case ("e")
                            atom%element = "Fe"
                            call alert_msg("warning","Fe taken as iron")
                        case ("E")
                            atom%element = "Fe"
                            call alert_msg("warning","FE taken as iron")
                        case default
                            atom%element = "F"
                            if (adjustl(atom%element) /= adjustl(atom%name) ) &
                             call alert_msg("note",trim(adjustl(atom%name))//" taken as fluorine")
                    end select
                    return
               !==========
                case ("P")
               !==========
                !It can be P, Pb, Po
                    ! We consider that:
                    !  Pb and PB are lead
                    !  Po is polonium
                    !  PO is P labeled "O"
                    select case (atname(2:2))
                        case ("o")
                            atom%element = "Po"
                        case ("O")
                            atom%element = "Po"
                        case ("t")
                            atom%element = "Pt"
                        case ("T")
                            atom%element = "Pt"
                        case default
                            atom%element = "P"
                    end select
                    return
               !==========
                case ("S")
               !==========
                !It can be S, Sr, Se, Sn, Si
                    ! We consider that:
                    !  Sb is antimonium 
                    !  SB sulfur labeled as "B"
                    select case (atname(2:2))
                        case ("i")
                            atom%element = "Si"
                            call alert_msg("warning","Si taken as silicon")
                        case ("I")
                            atom%element = "Si"
                            call alert_msg("warning","SI taken as silicon")
                        case ("r")
                            atom%element = "Sr"
                            call alert_msg("warning","Sr taken as strontium")
                        case ("R")
                            atom%element = "Sr"
                            call alert_msg("warning","SR taken as strontium")
                        case ("n")
                            atom%element = "Sn"
                            call alert_msg("warning","Sn taken as tin (Sn)")
                        case ("N")
                            atom%element = "Sn"
                            call alert_msg("warning","SN taken as tin (Sn)")
                        case ("b")
                            atom%element = "Sb"
                            call alert_msg("warning","Sb taken as antimony")
                        case default
                            atom%element = "S"
                            if (adjustl(atom%element) /= adjustl(atom%name) ) &
                             call alert_msg("note",trim(adjustl(atom%name))//" taken as sulfur")
                    end select
                    return
               !==========
                case ("K")
               !==========
                !It can be K, Kr
                    select case (atname(2:2))
                        case ("r")
                            atom%element = "Kr"
                        case ("R")
                            atom%element = "Kr"
                        case default
                            atom%element = "K"
                    end select
                    return
               !==========
                case ("V")
               !==========
                !It can only be V
                    atom%element = "V"
                    return
               !==========
                case ("W")
               !==========
                !It can only be W
                    atom%element = "W"
                    return
               !==========
                case ("Y")
               !==========
                !It can only be Y
                    atom%element = "Y"
                    return
               !==========
                case ("U")
               !==========
                !It can only be U
                    atom%element = "U"
                    return
               !==========
                case ("I")
               !==========
                !It can be I, Ir
                    ! We consider that:
                    !  Ir and IR are iridium
                    select case (atname(2:2))
                        case ("r")
                            atom%element = "Ir"
                        case ("R")
                            atom%element = "Ir"
                        case default
                            atom%element = "I"
                    end select
                    return
            end select

            !Once one-letter conflicts are solved, the rest are trivial
            call set_lower_case(atname(2:2))
            atom%element = atname(1:2)

            return

        end subroutine get_element_from_name

    end subroutine atname2element

!======================================

    subroutine element2AtNum(molec)

    ! Subroutine to generate Atomic Numbers from element name

        type(str_resmol),intent(inout) :: molec

        !local
        integer :: i, iel

        do i=1,molec%natoms
            do iel=1,103
                if (adjustl(molec%atom(i)%element) == &
                    adjustl(atname_from_atnum(iel))) then
                    molec%atom(i)%AtNum = iel
                    exit
                endif
            enddo
        enddo

        return

    end subroutine element2AtNum


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine get_cog(molec)

        !--------------------------------------------------------------------------
        ! Compute the center of geometry and store it in molec%cog
        !--------------------------------------------------------------------------
        
        type(str_resmol),intent(inout)::molec

        !Local
        integer :: i

        molec%cogX = 0.d0
        molec%cogY = 0.d0
        molec%cogZ = 0.d0
        do i=1,molec%natoms
            molec%cogX = molec%cogX + molec%atom(i)%x
            molec%cogY = molec%cogY + molec%atom(i)%y
            molec%cogZ = molec%cogZ + molec%atom(i)%z
        enddo
        molec%cogX = molec%cogX/molec%natoms
        molec%cogY = molec%cogY/molec%natoms
        molec%cogZ = molec%cogZ/molec%natoms

        return

    end subroutine get_cog

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine get_com(molec)

        !--------------------------------------------------------------------------
        ! Compute the center of mass and store it in molec%com
        !--------------------------------------------------------------------------
        
        type(str_resmol),intent(inout)::molec

        !Local
        integer :: i
#ifdef DOUBLE
        real(8) :: Mass
#else
        real :: Mass
#endif
        

        molec%comX = 0.d0
        molec%comY = 0.d0
        molec%comZ = 0.d0
        Mass       = 0.d0
        do i=1,molec%natoms
            molec%comX = molec%comX + molec%atom(i)%x*molec%atom(i)%mass
            molec%comY = molec%comY + molec%atom(i)%y*molec%atom(i)%mass
            molec%comZ = molec%comZ + molec%atom(i)%z*molec%atom(i)%mass
            Mass = Mass + molec%atom(i)%mass
        enddo
        molec%comX = molec%comX/Mass
        molec%comY = molec%comY/Mass
        molec%comZ = molec%comZ/Mass

        return

    end subroutine get_com

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine inertia(molec,IM)

        type(str_resmol),intent(in) :: molec
        double precision,dimension(1:3,1:3),intent(out) :: IM 

        !Co6ers
        integer :: i, j


        !Moment of inertia
        IM=0.d0
        do i=1,molec%natoms
            IM(1,1)=IM(1,1)+molec%atom(i)%mass*(molec%atom(i)%y**2+molec%atom(i)%z**2)
            IM(2,2)=IM(2,2)+molec%atom(i)%mass*(molec%atom(i)%x**2+molec%atom(i)%z**2)
            IM(3,3)=IM(3,3)+molec%atom(i)%mass*(molec%atom(i)%x**2+molec%atom(i)%y**2)

            IM(2,1)=IM(2,1)-molec%atom(i)%mass*(molec%atom(i)%y*molec%atom(i)%x)
            IM(3,1)=IM(3,1)-molec%atom(i)%mass*(molec%atom(i)%z*molec%atom(i)%x)
            IM(3,2)=IM(3,2)-molec%atom(i)%mass*(molec%atom(i)%z*molec%atom(i)%y)
       enddo
       do i=1,3
          do j=1,i-1
              IM(j,i) = IM(i,j)
          enddo
        enddo

        return

    end subroutine inertia


    SUBROUTINE ROTATA1(molec,molecRef,Rot)

        !============================================================
        ! Description
        ! ------------
        !  routine that rotates one molecRot to minimize the
        !  RMSD with molecRef based on quaternion formalism
        !
        ! This version is addapted from FCclasses2, interfaced 
        ! with structure_types (v4)
        !============================================================

        use matrix
        use structure_types
        
        type(str_resmol),intent(in)            :: molec
        type(str_resmol),intent(in)            :: molecRef
        real(8),dimension(1:3,1:3),intent(out) :: rot
        !Local
        ! scalar
        integer :: i, j, k, ii, jj, kk, l, m, kin, imax, jmax, kmax, ivm
        integer :: n, n3
        real(8) :: aa, dist, dist0, pos, diff, diffmax, distmin, distmin0, distmin1, ppp
        real(8) :: aaa, aaam
        ! static arrays
        real(8),dimension(1:4,1:4) :: ar1, al2, aaq, aaqs, cvec
        real(8),dimension(1:4) :: e
        ! allocatable
        real(8),dimension(:),   allocatable :: disvet
        real(8),dimension(:,:), allocatable :: dist1, dist2
        ! interface specific
        real(8),dimension(:), allocatable :: geo1,geo2


        !interface with structure_types
        n= molecRef%natoms

        !allocation
        allocate(disvet(1:N))
        allocate(dist1(1:N,1:N), dist2(1:N,1:N))
        allocate(geo1(1:3*N),geo2(1:3*N))
       
        ! geo1 -> molec
        ! geo2 -> molecRef
        do i=1,3*N,3
            j = (i-1)/3 + 1
            geo1(i  ) = molec%atom(j)%x
            geo1(i+1) = molec%atom(j)%y
            geo1(i+2) = molec%atom(j)%z
            geo2(i  ) = molecRef%atom(j)%x
            geo2(i+1) = molecRef%atom(j)%y
            geo2(i+2) = molecRef%atom(j)%z
        enddo


        do k=1,n
            disvet(k)=0.d0
            kk=3*(k-1)
            do i=1,n
                ii=3*(i-1)
                aa=0.d0
                do j=1,3
                    aa=aa+(geo1(kk+j)-geo1(ii+j))**2
                enddo
                dist1(i,k)=dsqrt(aa)
            enddo
        enddo
        do k=1,n
            kk=3*(k-1)
            do i=1,n
                ii=3*(i-1)
                aa=0.d0
                do j=1,3
                    aa=aa+(geo2(kk+j)-geo2(ii+j))**2
                enddo
                dist2(i,k)=dsqrt(aa)
            enddo
        enddo

        k=0
        dist=0.d0
        do i=1,n
        do j=1,3
            k=k+1
            dist=dist+(geo1(k)-geo2(k))**2
        enddo
        enddo

            dist0=dist
            do k=1,4
            do kk=1,k
            aaq(kk,k)=0.d0
            enddo
            enddo
            do i=1,n
                k=3*(i-1)
                ar1(1,1)=0.d0
                ar1(2,1)=geo1(k+1)
                ar1(3,1)=geo1(k+2)
                ar1(4,1)=geo1(k+3)
                ar1(1,2)=-geo1(k+1)
                ar1(2,2)=0.d0
                ar1(3,2)=-geo1(k+3)
                ar1(4,2)=geo1(k+2)
                ar1(1,3)=-geo1(k+2)
                ar1(2,3)=geo1(k+3)
                ar1(3,3)=0.d0
                ar1(4,3)=-geo1(k+1)
                ar1(1,4)=-geo1(k+3)
                ar1(2,4)=-geo1(k+2)
                ar1(3,4)=geo1(k+1)
                ar1(4,4)=0.d0
! c     
                al2(1,1)=0.d0
                al2(2,1)=geo2(k+1)
                al2(3,1)=geo2(k+2)
                al2(4,1)=geo2(k+3)
                al2(1,2)=-geo2(k+1)
                al2(2,2)=0.d0
                al2(3,2)=geo2(k+3)
                al2(4,2)=-geo2(k+2)
                al2(1,3)=-geo2(k+2)
                al2(2,3)=-geo2(k+3)
                al2(3,3)=0.d0
                al2(4,3)=geo2(k+1)
                al2(1,4)=-geo2(k+3)
                al2(2,4)=geo2(k+2)
                al2(3,4)=-geo2(k+1)
                al2(4,4)=0.d0
! c     
                do kk=1,4
                do l=1,4
                do m=1,4
                aaq(kk,l)=aaq(kk,l)+al2(kk,m)*ar1(m,l)
                enddo
                enddo
                enddo 
! c     
            enddo
            do kk=1,4
            do l=1,4
            aaq(kk,l)=-aaq(kk,l)
            aaqs(kk,l)=aaq(kk,l)*1.d-6
            enddo
            enddo

        ! Diagonalize (original eigen call chaged to:)
        call diagonalize_full(aaq,4,cvec,e,"lapack")

        do i=1,4
        do j=1,4
        ar1(i,j)=0.d0
        do k=1,4
        ar1(i,j)=ar1(i,j)+aaq(i,k)*cvec(k,j)
        enddo
        enddo
        enddo
        do i=1,4
        do j=1,4
        aaq(i,j)=0.d0
        do k=1,4
        aaq(i,j)=aaq(i,j)+cvec(k,i)*ar1(k,j)
        enddo
        aaqs(i,j)=aaq(i,j)*1.d-6
        enddo
        enddo

        distmin0=0.d0
        do k=1,3*n
        distmin0=distmin0+geo1(k)**2+geo2(k)**2
        enddo
        
        distmin=distmin0-2.d0*e(4)  
        if (verbose>2) &  
         write(0,*) 'minimal distance =',distmin        
        if (-e(1)-e(4).gt.1.d-6) then
            distmin1=distmin0+2.d0*e(1)
            if (verbose>1) then
                write(0,*) 'should consider reflection'
                write(0,*) 'the minimal distance would be', distmin1
            endif
        endif
        rot(1,1)=cvec(1,4)**2+cvec(2,4)**2-cvec(3,4)**2-cvec(4,4)**2
        rot(2,1)=2.d0*(cvec(2,4)*cvec(3,4)+cvec(1,4)*cvec(4,4))
        rot(3,1)=2.d0*(cvec(2,4)*cvec(4,4)-cvec(1,4)*cvec(3,4))
        rot(1,2)=2.d0*(cvec(2,4)*cvec(3,4)-cvec(1,4)*cvec(4,4))
        rot(2,2)=cvec(1,4)**2-cvec(2,4)**2+cvec(3,4)**2-cvec(4,4)**2
        rot(3,2)=2.d0*(cvec(3,4)*cvec(4,4)+cvec(1,4)*cvec(2,4))
        rot(1,3)=2.d0*(cvec(2,4)*cvec(4,4)+cvec(1,4)*cvec(3,4))
        rot(2,3)=2.d0*(cvec(3,4)*cvec(4,4)-cvec(1,4)*cvec(2,4))
        rot(3,3)=cvec(1,4)**2-cvec(2,4)**2-cvec(3,4)**2+cvec(4,4)**2
             
        return
      end subroutine ROTATA1


      subroutine rotate_molec(molec,Rot)

        ! X' = R * X
        ! (from FClasses code)

        type(str_resmol),intent(inout)         :: molec
        real(8),dimension(1:3,1:3),intent(in)  :: rot
        !Local
        ! scalar
        integer :: i, j, k, kin
        integer :: n
        ! interface specific
        real(8),dimension(:), allocatable :: geo1,geo1r

        !interface with structure_types
        n= molec%natoms

        !allocation
        allocate(geo1(1:3*N),geo1r(1:3*N))

        ! molec -> geo1
        do i=1,3*N,3
            j = (i-1)/3 + 1
            geo1(i  ) = molec%atom(j)%x
            geo1(i+1) = molec%atom(j)%y
            geo1(i+2) = molec%atom(j)%z
        enddo
          
          DO  I=1,N
           kin=(i-1)*3
            DO  J=1,3
             geo1r(kin+j)=0.d0
              do k=1,3
               geo1r(kin+j)=geo1r(kin+j)+ROT(j,k)*geo1(kin+k)
              enddo
            enddo
           enddo  

        ! geo1r -> molec
        do i=1,3*N,3
            j = (i-1)/3 + 1
            molec%atom(j)%x = geo1r(i  )
            molec%atom(j)%y = geo1r(i+1)
            molec%atom(j)%z = geo1r(i+2)
        enddo

          return

      end subroutine rotate_molec

end module molecular_structure
