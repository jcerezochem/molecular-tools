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
    implicit none

    contains

    subroutine assign_masses(molec)

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

        type(str_resmol),intent(inout) :: molec

        integer :: i
        character(len=2) :: atname

        do i=1,molec%natoms
            atname = adjustl(molec%atom(i)%element)
            select case (atname)
               case ("H")
                molec%atom(i)%mass=1.0078250
               case ("He")
                molec%atom(i)%mass=4.0026033
               case ("Li")
                molec%atom(i)%mass= 7.0160045
               case ("Be")
                molec%atom(i)%mass= 9.0121825
               case ("B")
                molec%atom(i)%mass=11.0093053
               case ("C")
                molec%atom(i)%mass=12.000000
               case ("N")
                molec%atom(i)%mass=14.0030740
               case ("O")
                molec%atom(i)%mass=15.9949146
               case ("F")
                molec%atom(i)%mass=18.9984033
               case ("Ne")
                molec%atom(i)%mass=19.9924391
               case ("Na")
                molec%atom(i)%mass=22.9897697
               case ("Mg")
                molec%atom(i)%mass=23.9850450
               case ("Al")
                molec%atom(i)%mass=26.9815413
               case ("Si")
                molec%atom(i)%mass=27.9769284
               case ("P")
                molec%atom(i)%mass=30.9737634
               case ("S")
                molec%atom(i)%mass=31.9720718
               case ("Cl")
                molec%atom(i)%mass=34.9688527
               case ("Ar")
                molec%atom(i)%mass=39.9623831
               case ("K")
                molec%atom(i)%mass=38.9637079
               case ("Ca")
                molec%atom(i)%mass=39.9625907
               case ("Sc")
                molec%atom(i)%mass=44.9559136
               case ("Ti")
                molec%atom(i)%mass=47.9479467
               case ("V")
                molec%atom(i)%mass=50.9439625
               case ("Cr")
                molec%atom(i)%mass=51.9405097
               case ("Mn")
                molec%atom(i)%mass=54.9380463
               case ("Fe")
                molec%atom(i)%mass=55.9349393
               case ("Co")
                molec%atom(i)%mass=58.9331978
               case ("Ni")
                molec%atom(i)%mass=57.9353471
               case ("Cu")
                molec%atom(i)%mass=62.9295992
               case ("Zn")
                molec%atom(i)%mass=63.9291454
               case ("Ga")
                molec%atom(i)%mass=68.9255809
               case ("Ge")
                molec%atom(i)%mass=73.9211788
               case ("As")
                molec%atom(i)%mass=74.9215955
               case ("Se")
                molec%atom(i)%mass=79.9165205
               case ("Br")
                molec%atom(i)%mass=78.9183361
               case ("Kr")
                molec%atom(i)%mass=83.9115064
               !Desordenados
               case ("Pd")
                molec%atom(i)%mass=105.9032000
               case ("Pt")
                molec%atom(i)%mass=194.9648000
               case ("I")
                molec%atom(i)%mass=126.9004000
               !Default
               case default
                call alert_msg("warning","Don't know how to assign mass to "//atname//" Set to zero.")
                molec%atom(i)%mass=0.00
            end select
        enddo

        return

    end subroutine assign_masses


   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine guess_connect(molec)

        implicit none
        type(str_resmol),intent(inout) :: molec
        !local
        real :: av_len, dist
        integer :: i,j, i_cnx

        !InOut, better reassigned
        molec=molec

        !Loop over all atoms. This is not optimal, distance matrix is symmetric
        !(nor it is not costly, anyway)
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
                
                av_len=bond_length_db(molec%atom(i),molec%atom(j))
                ! Criterium: dist < av. length +10% --increased from 5% (17/12/12) due to H2O2
                ! a more sophiticated data base might include hibridization
                av_len=av_len*1.1

                if (dist < av_len) then
! print*, "Connection"
! print*, molec%atom(i)%name,molec%atom(j)%name
! print*, i,j,av_len
! print*, ""
                    i_cnx=i_cnx+1
                    molec%atom(i)%connect(i_cnx)=j
                endif

            enddo
            molec%atom(i)%nbonds=i_cnx
        enddo
   
        return

    end subroutine guess_connect


    function bond_length_db(atom1,atom2) result(av_len)

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

        !DB entries
        character(len=2) :: elementA, elementB
        real :: db_length

        type(str_atom),intent(in) :: atom1, atom2
        !Local atoms
        type(str_atom) :: Latom1, Latom2
        real :: av_len
        !local
        character(len=14),dimension(100) :: database
        integer :: n_entries, i
        logical :: external_DB=.false.

        !Set local atoms
        Latom1%element = adjustl(atom1%element)
        Latom2%element = adjustl(atom2%element)


        ! DATABASE GENERATION
        if (external_DB) then
            print*, "External DB not yet supported"
        else
            n_entries=28
            ! Default implementation
            ! Taken from CRC Handbook of Chemistry and Physics
            !(that's a good plan, for the moment, they're gv defaults)
            database(1:n_entries)=            &
                    (/                        &
                    ! A    B    bond length(\AA)
                    "C    C    1.54",      & !1
                    "C    O    1.43",      & !2
                    "C    S    1.78",      & !3
                    "C    H    1.07",      & !4
                    "O    O    1.34",      & !5 was 1.32, but was not enough for H2O2 (1.47 needed 1.337)
                    "O    S    1.67",      & !6
                    "O    H    1.80",      & !7 !to include Hbonds use 1.80, else use 1.07
                    "S    S    2.02",      & !8
                    "S    H    1.31",      & !9
                    "H    H    0.60",      & !10
                    "N    N    1.40",      & !11
                    "N    C    1.47",      & !12
                    "N    O    1.36",      & !13
                    "N    S    1.71",      & !14
                    "N    H    1.00",      & !15
                    "F    H    0.88",      & !16
                    "F    C    1.35",      & !17
                    "F    N    1.28",      & !18
                    "F    O    1.24",      & !19
                    "F    S    1.59",      & !20
                    "F    F    1.16",      & !21
                    "Mg   N    2.15",      & !22  Larger than gv standard (2.06)
                    "Ni   O    2.10",      & !23  To be revised
                    "Ni   N    2.10",      & !24  To be revised
                    "Cl   C    1.80",      & !25  Larger than gv standard (1.76)
                    "P    O    1.71",      & !26
                    "P    H    1.35",      & !27
                    "Pt   N    2.20"       & !28
                    /)                       
        endif

        ! SELECTION LOOPS 
        av_len=0.
        do i=1,n_entries 
            !Get database entry
            read(database(i),*) elementA, elementB, db_length
            !Compare with input
            if ( ( Latom1%element(1:2) == adjustl(elementA) .and. &
                   Latom2%element(1:2) == adjustl(elementB) )     &
                .or.                                            &
                 ( Latom1%element(1:2) == adjustl(elementB) .and. &
                   Latom2%element(1:2) == adjustl(elementA) )     &
                ) then
                   av_len=db_length
!                    print*, "Found ", db_entry%A," ", db_entry%B," ", Latom1%name(1:1)," ", Latom2%name(1:1)
            endif
       enddo     

! if ( av_len == 0.) then
!     call alert_msg("fatal",Latom1%name(1:1)//" and "//Latom2%name(1:1)//" missing in the DB")
! endif

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
                        case ("f")
                            atom%element = "Hf"
                        case default
                            atom%element = "H"
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
                        case ("e")
                            atom%element = "Be"
                        case ("E")
                            atom%element = "Be"
                        case ("r")
                            atom%element = "Br"
                        case ("R")
                            atom%element = "Br"
                        case default
                            atom%element = "B"
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
                        case ("r")
                            atom%element = "Cr"
                        case ("R")
                            atom%element = "Cr"
                        case ("s")
                            atom%element = "Cs"
                        case ("S")
                            atom%element = "Cs"
                        case ("l")
                            atom%element = "Cl"
                        case ("L")
                            atom%element = "Cl"
                        case ("a")
                            atom%element = "Ca"
                        case ("A")
                            !This case can be either C"A" or Ca. Mark with x to check later
                            atom%element = "Cx"
                        case default
                            atom%element = "C"
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
                        case ("E")
                            atom%element = "Fe"
                        case default
                            atom%element = "F"
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
                        case ("e")
                            atom%element = "Fe"
                        case ("E")
                            atom%element = "Fe"
                        case ("o")
                            atom%element = "Po"
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
                        case ("I")
                            atom%element = "Si"
                        case ("r")
                            atom%element = "Sr"
                        case ("R")
                            atom%element = "Sr"
                        case ("n")
                            atom%element = "Sn"
                        case ("N")
                            atom%element = "Sn"
                        case ("b")
                            atom%element = "Sb"
                        case default
                            atom%element = "S"
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

        do i=1,molec%natoms
            do iel=1,103
                if (adjustl(molec%atom(i)%name) == &
                    adjustl(atom_names_from_atnum(iel))) then
                    molec%atom(i)%AtNum = iel
                    exit
                endif
            enddo
        enddo

        return

    end subroutine element2AtNum

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!   RUTINES FOR MOLECULE <-> RESIDUE MANAGEMENT (previously was in gro_manage


    subroutine renum(system)

        !====================================
        ! Renumber residues in a gro system
        !====================================

        type(str_resmol),intent(inout)::system

        !local
        integer::i, i_res, prev_res


        i_res=1
        prev_res=system%atom(1)%resseq

        do i=1,system%natoms

            if (system%atom(i)%resseq /= prev_res) i_res=i_res+1
            prev_res=system%atom(i)%resseq
            system%atom(i)%resseq=i_res

        enddo

        return

    end subroutine renum

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine sist_nres(sistema)

        !--------------------------------------------------------------------------

        type(str_resmol),intent(inout)::sistema

        !Contadores
         integer::i,ii,iat,ires

        !criterios para busqueda de residuos
        integer::resI_test
        character::resA_test*5


        !Como es intetn(inout), damos valor a todas las direcciones
        !de la estructura dentro de la subrutina.
        sistema = sistema


        !Encontramos los residuos
        !Para ello aplicamos conjuntamente los dos siguientes criterios:
        resI_test=-1     !criterio 1: cambio en res.seq
        resA_test='???'  !criterio 2: cambio en res.name
        ires=0
        do i=1,sistema%natoms
!            if (indice /= sistema%atom(i)%res_index) ires = ires + 1
!       La siguiente condición es más "robusta" para localizar residuos
            if (resI_test /= sistema%atom(i)%resseq .or. resA_test /= sistema%atom(i)%resname) ires = ires + 1
            resI_test=sistema%atom(i)%resseq
            resA_test=sistema%atom(i)%resname
        enddo

        sistema%nres=ires
 
        return

    end subroutine sist_nres

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine sist2res(sistema,residuo)

        ! TO BE DEPRECATED

        !--------------------------------------------------------------------------
        ! Subrutina que localiza y guarda en una estructura (tipo derivado)
        ! residuos en un sistema leido de PDB o GRO. No es trivial, ya que los
        ! ficheros pueden ser algo desordenados. El uso de dos criterios asegura
        ! encontrar los residuos siempre que las estructuras no esten partidas 
        ! y dispersas entre el PDB o GRO. En ese caso la única posibilidad sería
        ! usar la topología del sistema leída de ficheros separados
        !--------------------------------------------------------------------------

        type(str_resmol),intent(inout)::sistema
        type(str_resmol),intent(inout)::residuo(:)

        !Contadores
         integer::i,ii,iat,ires

        !criterios para busqueda de residuos
        integer::resI_test
        character::resA_test*5


        !Como es intetn(inout), damos valor a todas las direcciones
        !de la estructura dentro de la subrutina.
        sistema = sistema
        residuo(:)%frst_atom=-1


        !Encontramos los residuos
        !Para ello aplicamos conjuntamente los dos siguientes criterios:
        resI_test=-1     !criterio 1: cambio en res.seq
        resA_test='???'  !criterio 2: cambio en res.name
        ires=0
        do i=1,sistema%natoms
!            if (indice /= sistema%atom(i)%res_index) ires = ires + 1
!       La siguiente condición es más "robusta" para localizar residuos
            if (resI_test /= sistema%atom(i)%resseq .or. resA_test /= sistema%atom(i)%resname) ires = ires + 1
            resI_test=sistema%atom(i)%resseq
            resA_test=sistema%atom(i)%resname
            residuo(ires)%name = sistema%atom(i)%resname
            if ( residuo(ires)%frst_atom == -1 ) residuo(ires)%frst_atom = i
            residuo(ires)%lst_atom = i
        enddo

        sistema%nres=ires

        do i=1,sistema%nres
            residuo(i)%natoms = residuo(i)%lst_atom - residuo(i)%frst_atom + 1
        enddo


        !Y les asignamos sus átomos
        iat = 0
        do i=1,sistema%nres
            do ii=1,residuo(i)%natoms
                iat = iat + 1
                residuo(i)%atom(ii) = sistema%atom(iat)
            enddo
        enddo

print*, sistema%nres

       return
    end subroutine sist2res


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine sist2res_new(sistema,residuo)

        !EN PRUEBAS!!!

        !--------------------------------------------------------------------------
        ! La subrutina antigua es un poco paranoica, no creo que sea muy normal que
        ! en un pdb/gro no cambie el numero de residuo cuando debería porque cambie
        ! el nombre... (pero...   no ocurria para BCR???)
        !--------------------------------------------------------------------------

        type(str_resmol),intent(inout)::sistema
        type(str_resmol),intent(inout)::residuo(:)

        !Contadores
         integer::i,ii,iat,ires

        !auxiliares para busqueda de residuos
        integer::nres


        !Como es intetn(inout), damos valor a todas las direcciones
        !de la estructura dentro de la subrutina.
        sistema = sistema

        !Llenamos la array "residuo()"
        residuo(:)%natoms=0
        nres=0
        do i=1,sistema%natoms
            ires=sistema%atom(i)%resseq
            residuo(ires)%name = sistema%atom(i)%resname
            residuo(ires)%natoms=residuo(ires)%natoms+1
            residuo(ires)%atom(residuo(ires)%natoms)=sistema%atom(i)
            !Lo siguiente nos ayudará a tener una secuencia continua:
            nres=max(ires,nres)
        enddo

        !Ahora comprobamos que la numeración de los residuos es continua
        ires=0
        do i=1,nres
            ires=ires+1
            if ( residuo(ires)%natoms == 0 ) then
                ires = ires - 1
                cycle
            endif
            residuo(ires)=residuo(i)
        enddo

       return
    end subroutine sist2res_new

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine sist2res1_new(sistema,residuo)

        !EN PRUEBAS!!!

        !--------------------------------------------------------------------------
        ! La subrutina antigua es un poco paranoica, no creo que sea muy normal que
        ! en un pdb/gro no cambie el numero de residuo cuando debería porque cambie
        ! el nombre... (no ocurria para BCR???)
        ! VERSION FOR ONE-RESIDUE SYSTEMS -- should be stable for this case!
        !
        ! This version manages only one residue (not an array)
        !--------------------------------------------------------------------------

        type(str_resmol),intent(inout)::sistema
        type(str_resmol),intent(inout)::residuo

        !Contadores
         integer::i,ii,iat,ires

        !auxiliares para busqueda de residuos
        integer::nres


        !Como es intetn(inout), damos valor a todas las direcciones
        !de la estructura dentro de la subrutina.
        sistema = sistema

        !Consideramos que solo hay un residuo (aunque haya más... WARNING?)
        residuo%natoms=0
        do i=1,sistema%natoms
            ires=sistema%atom(i)%resseq
            residuo%name = sistema%atom(i)%resname
            residuo%natoms=residuo%natoms+1
            residuo%atom(residuo%natoms)=sistema%atom(i)
        enddo

       return
    end subroutine sist2res1_new


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine res2sist(residuo,nres,boxX,boxY,boxZ,sistema)

        !--------------------------------------------------------------------------
        ! The inverse as befeore (this is trivial)
        !--------------------------------------------------------------------------

        type(str_resmol),intent(in)::residuo(:)
        integer, intent(in)::nres
        type(str_resmol),intent(inout)::sistema
#ifdef DOUBLE
        real(8),intent(in) :: boxX,boxY,boxZ
#else
        real,intent(in) :: boxX,boxY,boxZ
#endif

        !Contadores
        integer::i,ii,iat


        iat=0
        do i=1,nres
            do ii=1,residuo(i)%natoms
                iat=iat+1
                sistema%atom(iat)=residuo(i)%atom(ii)
            enddo
        enddo

        sistema%natoms=iat
        sistema%nres=nres

        sistema%boxX=boxX
        sistema%boxY=boxY
        sistema%boxZ=boxZ
            
        return

    end subroutine res2sist

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine res1_2_sist(residuo,nres,boxX,boxY,boxZ,sistema)

        !--------------------------------------------------------------------------
        ! The inverse as befeore (this is trivial)
        !--------------------------------------------------------------------------

        type(str_resmol),intent(in)::residuo
        integer, intent(in)::nres
        type(str_resmol),intent(inout)::sistema
#ifdef DOUBLE
        real(8),intent(in) :: boxX,boxY,boxZ
#else
        real,intent(in) :: boxX,boxY,boxZ
#endif

        !Contadores
        integer::i,ii,iat


        iat=0
        do i=1,1 !nres
            do ii=1,residuo%natoms
                iat=iat+1
                sistema%atom(iat)=residuo%atom(ii)
            enddo
        enddo

        sistema%natoms=iat
        sistema%nres=nres

        sistema%boxX=boxX
        sistema%boxY=boxY
        sistema%boxZ=boxZ
            
        return

    end subroutine res1_2_sist

end module molecular_structure
