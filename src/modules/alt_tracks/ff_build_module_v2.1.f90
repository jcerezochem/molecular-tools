module ff_build

    !==============================================================
    ! This code is part of MOLECULAR_TOOLS (version 1.0/March 2012)
    !==============================================================
    ! MODULE with subroutines to associate atomtypes to atom names based on 
    ! the element name and the connectivity
    ! Atom names are the name of the element
    ! We give names based on that, adding a numbered dummy_char2 at the end, i.e.
    ! C becomes C1 the first time it appears. Hidrogens are named after the
    ! the element they are linked to: HC11, HC12... in the previus case
    !
    ! Subroutines and funcitons in this module (all are independent):
    !    subroutine fftype_db(atom,numH, fftype,H_fftype, I_DB, serial)
    !    subroutine refine_charges(residue) 
    !    subroutine build_ff(residue)
    !    subroutine gen_bonded(residue)
    !    subroutine assign_masses(molec)
    !    subroutine guess_connect(molec)
    !    subroutine guess_connect_res(residue)
    !    function bond_length_db(atom1,atom2) result(av_len)
    !    subroutine aa2ua(molecAA,molecUA)
    !    subroutine aa2ua_res(residueAA,residueUA)
    !
    ! V2.1: contains the symmetrize_charges subroutine (out of track)
    !
    !===========================================================================
    use alerts
    use structure_types
    implicit none

    contains

   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine fftype_db(atom,numH, fftype,H_fftype, I_DB, serial)

    !=======================================================================
    !Description
    !-------------
    ! Database assignement of the atom_type. It's based on the atom name
    ! and environment.
    ! External database support: must be opened in main, the unit nr is 
    ! passed as I_DB. If I_DB=0, hard-coded defaults are used.
    !
    ! About this database approach
    !-----------------------------
    ! Seting up a database is not an obvious task. Here you have one try:
    ! (ok, you'll find a nice description here, soon) 
    !
    ! Pros and cons
    !---------------
    ! Pros: easy extension to external DB
    !       clarity of the implmented DB
    ! Cons: Memory waste storing DB (not that much, though)
    !
    !Changelog 
    !---------
    ! Added external database feature (6/06/12): Involved changes in the 
    ! number of arguments passed to the SR!!
    ! Added other database native support using specific I_DB (22/06/2012):
    !   I_DB > 0     : Read external DB
    !   I_DB = -1    : Hybridization database
    !   I_DB = other : CHARMM database (e.g. if not set in the main program, I_DB=0)
    !
    !Whislist
    !--------
    ! serial should no be here. The error call is more suitable from main,
    ! where this SR is called. So, it'd better to add a ierr flag. (whised
    ! 6/06/12)
    ! Wildcards "-1" does not work properly. To be extended for nb adn nH
    !
    !=======================================================================

        use line_preprocess

        !-----------------------------------------
        type DB_atom_types
            character(len=2) :: A
            integer :: nb, nH
            character(len=4) :: fftype, H_fftype         
        end type DB_atom_types
        !-----------------------------------------

        type(str_atom),intent(in) :: atom
        integer,intent(in) :: numH, serial
        character(len=*),intent(out) :: fftype, H_fftype 
        !local
        character(len=5) :: Latom_name
        integer :: i
        character(len=21),dimension(100) :: database
        !new structure not really needed (not hurting anyway)
        type(DB_atom_types) :: db_entry
        integer :: n_entries
        logical :: found
        !aux
        character(len=100) :: dummy_char, dummy_char2, dummy_char3
        !I/O
        integer,intent(in) :: I_DB
        integer :: iostatus

        !Get local atom name
        Latom_name = adjustl(atom%name)
        if (Latom_name(1:1) == "M") then
            Latom_name = adjustl(Latom_name(1:2))
        else
            Latom_name = adjustl(Latom_name(1:1))
        endif

        ! DATABASE GENERATION (units > 0)
        if (I_DB > 0) then
            !TODO: Store the database to avoid reading everytime!!
            rewind(I_DB)
            read(I_DB,*,IOStat=iostatus) n_entries
            if (iostatus /= 0) call alert_msg("fatal","Cannot read the external data base")
            i=1
            do
                read(I_DB,'(A)',IOSTAT=iostatus) dummy_char
                if (iostatus /= 0) call alert_msg("fatal","Unexpected error while reading the external data base")
                call split_line(dummy_char,"!",dummy_char2,dummy_char3)
                if (len_trim(dummy_char2) /= 0) then
                    database(i)=adjustl(dummy_char2)
                    i=i+1
                endif
                if ( i > n_entries) exit
            enddo
        else if (I_DB == -1) then
           !Hybrid types database
            n_entries=11
            database(1:n_entries)=                   &
                    (/                               &
                     !Elem  nb nH fftype H_fftype
                       "C   4 -1   SP3   XX  ",      & !1
                       "C   3 -1   SP2   XX  ",      & !2
                       "C   2 -1   SP    XX  ",      & !3
                       "O   1 -1   SP2   XX  ",      & !4
                       "O   2 -1   SP3   XX  ",      & !5 
                       "S   2 -1   SP3   XX  ",      & !6 
                       "N   3  0   SP2   XX  ",      & !7 ;planar
                       "N   3  1   SP2   XX  ",      & !8 ;planar
                       "N   3  2   SP2   XX  ",      & !9 ;planar
                       "N   3  3   SP3   XX  ",      & !10
                       "MG -1 -1   SP3   XX  "       & !11
                     /)                       
        else
            ! Default implementation
            ! CHARMM27
            n_entries=19
            database(1:n_entries)=                   &
                    (/                               &
                     !Elem  nb nH fftype H_fftype
                       "C   4  0   CT    XX  ",      & !1
                       "C   4  1   CT1   HA  ",      & !2
                       "C   4  2   CT2   HA  ",      & !3
                       "C   4  3   CT3   HA  ",      & !4
!                        "C   4  1   CT1x  HA1 ",      & !2
!                        "C   4  2   CT2x  HA2 ",      & !3
!                        "C   4  3   CT3x  HA3 ",      & !4
                       "C   3  0   CA    XX  ",      & !5 ;Aromatic
!                        "C   3  1   CE1   HE1 ",      & !6
!                        "C   3  2   CE2   HE2 ",      & !7
                       "C   3  1   CA    HP  ",      & !6 ;Aromatic
                       "C   3  2   CA    HP  ",      & !7 ;Aromatic
                       "O   1  0   O     XX  ",      & !8
                       "O   2  0   OS    XX  ",      & !9
                       "O   2  1   OH1   H   ",      & !10
                       "O   2  2   OT    HT  ",      & !11
                       "S   2  0   S     XX  ",      & !12
                       "S   2  1   SS    HT  ",      & !13
                       "N   3  0   NPH   XX  ",      & !14 (HEME pirrole ring)
                       "N   3  1   NH1   H   ",      & !15 (peptide nitrogen)
                       "N   3  2   NH2   H   ",      & !16
                       "N   4  3   NH3   H   ",      & !17
                       "F   1  0   F1    XX  ",      & !18
                       "MG -1  0   MG    XX  "       & !19
                     /)                       
        endif

        ! SELECTION LOOPS
        found=.false.
        do i=1,n_entries
            !Get database entry
            read(database(i),*) db_entry%A, db_entry%nb, db_entry%nH, db_entry%fftype, db_entry%H_fftype
            !Compare with input
            if ( Latom_name /= adjustl(db_entry%A) ) cycle
            if ( atom%nbonds /= db_entry%nb .and. db_entry%nb /= -1) cycle
            if ( numH /= db_entry%nH .and. db_entry%nH /= -1 ) cycle
            !The following complex "if" structure also works (amazing)
!            if ( Latom_name == adjustl(db_entry%A) .and. &
!                 ( atom%nbonds == db_entry%nb            &
!                   .or. db_entry%nb == -1 )       .and.  &
!                 numH == db_entry%nH                     & 
!               ) then
                fftype=db_entry%fftype
                H_fftype=db_entry%H_fftype
                found=.true.
                exit
!            endif
        enddo 


        write(dummy_char,*)  numH
        write(dummy_char2,*) atom%nbonds
        write(dummy_char3,*) serial
        !Change it to note? Avoid stopping if the gro can still be generated properly. Now it the error makes clear, however.
        if (.not.found) then
            call alert_msg("warning","Atom "//trim(adjustl(dummy_char3))//" not found in data base: "//&
                                       trim(adjustl(atom%name))//" with "//trim(adjustl(dummy_char))//"H/"//    &
                                       trim(adjustl(dummy_char2))//"bonds" )    
            H_fftype = "HX"
            fftype = "XX"
        endif

        return
    end subroutine fftype_db


    subroutine hybrid_db(atom,numH, fftype,H_fftype, I_DB, serial)

    !=======================================================================
    !Description
    !-------------
    ! Database assignement of the atom_type. It's based on the atom name
    ! and environment.
    ! External database support: must be opened in main, the unit nr is 
    ! passed as I_DB. If I_DB=0, hard-coded defaults are used.
    !
    ! About this database approach
    !-----------------------------
    ! Seting up a database is not an obvious task. Here you have one try:
    ! (ok, you'll find a nice description here, soon) 
    !
    ! Pros and cons
    !---------------
    ! Pros: easy extension to external DB
    !       clarity of the implmented DB
    ! Cons: Memory waste storing DB (not that much, though)
    !
    !Changelog 
    !---------
    ! Added external database feature (6/06/12): Involved changes in the 
    ! number of arguments passed to the SR!!
    !
    !Whislist
    !--------
    ! serial should no be here. The error call is more suitable from main,
    ! where this SR is called. So, it'd better to add a ierr flag. (whised
    ! 6/06/12)
    ! Wildcards "-1" does not work properly. To be extended for nb adn nH
    !
    !=======================================================================

        use line_preprocess

        !-----------------------------------------
        type DB_atom_types
            character(len=2) :: A
            integer :: nb, nH
            character(len=4) :: fftype, H_fftype         
        end type DB_atom_types
        !-----------------------------------------

        type(str_atom),intent(in) :: atom
        integer,intent(in) :: numH, serial
        character(len=*),intent(out) :: fftype, H_fftype 
        !local
        character(len=5) :: Latom_name
        integer :: i
        character(len=21),dimension(100) :: database
        !new structure not really needed (not hurting anyway)
        type(DB_atom_types) :: db_entry
        integer :: n_entries
        logical :: found
        !aux
        character(len=100) :: dummy_char, dummy_char2, dummy_char3
        !I/O
        integer,intent(in) :: I_DB
        integer :: iostatus

        !Get local atom name
        Latom_name = adjustl(atom%name)
        if (Latom_name(1:1) == "M") then
            Latom_name = adjustl(Latom_name(1:2))
        else
            Latom_name = adjustl(Latom_name(1:1))
        endif

        ! DATABASE GENERATION
        if (I_DB /= 0) then
            !TODO: Store the database to avoid reading everytime!!
            rewind(I_DB)
            read(I_DB,*,IOStat=iostatus) n_entries
            if (iostatus /= 0) call alert_msg("fatal","Cannot read the external data base")
            i=1
            do
                read(I_DB,'(A)',IOSTAT=iostatus) dummy_char
                if (iostatus /= 0) call alert_msg("fatal","Unexpected error while reading the external data base")
                call split_line(dummy_char,"!",dummy_char2,dummy_char3)
                if (len_trim(dummy_char2) /= 0) then
                    database(i)=adjustl(dummy_char2)
                    i=i+1
                endif
                if ( i > n_entries) exit
            enddo
        else
           n_entries=11
            ! Default implementation
            database(1:n_entries)=                   &
                    (/                               &
                     !Elem  nb nH fftype H_fftype
                       "C   4 -1   SP3   XX  ",      & !1
                       "C   3 -1   SP2   XX  ",      & !2 
                       "C   2 -1   SP    XX  ",      & !3 
                       "O   1 -1   SP2   XX  ",      & !4
                       "O   2 -1   SP3   XX  ",      & !5
                       "S   2 -1   SP3   XX  ",      & !6
                       "N   3  0   SP2   XX  ",      & !7 (plane, but could not be..)
                       "N   3  1   SP2   XX  ",      & !8 (plane, but could not be..)
                       "N   3  2   SP2   XX  ",      & !9 (plane, but could not be..)
                       "N   3  3   SP3   XX  ",      & !10
                       "MG -1 -1   SP3   XX  "       & !11
                     /)                       
        endif

        ! SELECTION LOOPS
        found=.false.
        do i=1,n_entries
            !Get database entry
            read(database(i),*) db_entry%A, db_entry%nb, db_entry%nH, db_entry%fftype, db_entry%H_fftype
            !Compare with input
            if ( Latom_name /= adjustl(db_entry%A) ) cycle
            if ( atom%nbonds /= db_entry%nb .and. db_entry%nb /= -1) cycle
            if ( numH /= db_entry%nH .and. db_entry%nH /= -1 ) cycle
            !The following complex "if" structure also works (amazing)
!            if ( Latom_name == adjustl(db_entry%A) .and. &
!                 ( atom%nbonds == db_entry%nb            &
!                   .or. db_entry%nb == -1 )       .and.  &
!                 numH == db_entry%nH                     & 
!               ) then
                fftype=db_entry%fftype
                H_fftype=db_entry%H_fftype
                found=.true.
                exit
!            endif
        enddo 


        write(dummy_char,*)  numH
        write(dummy_char2,*) atom%nbonds
        write(dummy_char3,*) serial
        !Change it to note? Avoid stopping if the gro can still be generated properly. Now it the error makes clear, however.
        if (.not.found) then
            call alert_msg("warning","Atom "//trim(adjustl(dummy_char3))//" not found in data base: "//&
                                       trim(adjustl(atom%name))//" with "//trim(adjustl(dummy_char))//"H/"//    &
                                       trim(adjustl(dummy_char2))//"bonds" )    
            H_fftype = "HX"
            fftype = "XX"
        endif

        return
    end subroutine hybrid_db

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine refine_charges(residue) 

        !New type
        type list
            character(len=5) :: name
            integer :: items
            integer,dimension(100) :: atom_list 
            character(len=6),dimension(7) :: env
            real :: charge
         end type list

        type(str_residue),intent(inout) :: residue

        !Fixed charges
        real :: q_Hal=0.09,  &  ! aliphatic hydrogen
                q_Car=-0.115, &  ! aromatic carbon attached to H
                q_Har=0.115     ! aromatic hydrogen (NOTE: proposed value
                                ! in GenFF pare is +0.115/-0.115 for H/C
                                ! however, ParamChem gives +0.150/-0.150
                                ! Maybe paper value update

        type(list),dimension(50) :: fftype
        logical :: existing_fftype
        real :: q_tot_prev, q_tot, delta_q

        !Counters
        integer :: i,j, jj, n_fftypes, n_extra_charges, n_compens

        character(len=50) ::  dummy_char,dummy_char2,dummy_char3


        ! =========================================
        ! MAKE SAME ATOMTYPES WITH THE SAME CHARGE  
        ! =========================================

!         ! 1. Identify atom types with atom numbers and calculata average charges
!         n_fftypes=0
!         q_tot_prev=0.
!         do i=1,system%natoms
! 
!             ! Look for the atom type in the list
!             existing_fftype=.false.
!             do j=1,n_fftypes
!                 if ( adjustl(system%atom(i)%fftype) == adjustl(fftype(j)%name) ) then
!                     existing_fftype=.true.
!                     fftype(j)%charge=fftype(j)%charge + system%atom(i)%q
!                     fftype(j)%items=fftype(j)%items + 1
!                     fftype(j)%atom_list(fftype(j)%items) = i
!                     exit
!                 endif
!             enddo
! 
!             !If not in the list, update it
!             if (.not.existing_fftype) then
!                 n_fftypes=n_fftypes+1
!                 fftype(n_fftypes)%name = adjustl(system%atom(i)%fftype)
!                 !Initialize this type
!                 fftype(n_fftypes)%charge = system%atom(i)%q
!                 fftype(n_fftypes)%items = 1
!                 fftype(j)%atom_list(fftype(j)%items) = i
!             endif
! 
!             q_tot_prev = q_tot_prev + system%atom(i)%q
! 
!         enddo
        ! 1. Identify atom types with atom numbers and calculate average charges
        !    Atoms are groupped by type&environment
        n_fftypes=0
        q_tot_prev=0.
        fftype(n_fftypes)%env=""
        fftype(n_fftypes)%name=""
        do i=1,residue%natoms

            ! Look for the atom type in the list
!             existing_fftype=.false.
            do j=1,n_fftypes
                existing_fftype=.true.
                if ( residue%atom(i)%fftype == fftype(j)%name ) then
                    do jj=1,residue%atom(i)%nbonds
                        if ( residue%atom(i)%env(jj) /= fftype(j)%env(jj) ) &
                          existing_fftype=.false.
                    enddo 
                else 
                    existing_fftype=.false.
                endif
                if (existing_fftype) then
                    fftype(j)%charge=fftype(j)%charge + residue%atom(i)%q
                    fftype(j)%items=fftype(j)%items + 1
                    fftype(j)%atom_list(fftype(j)%items) = i
                    exit
                endif
            enddo

            !If not in the list, update it
            if (.not.existing_fftype) then
                n_fftypes=n_fftypes+1
                !Initialize this type
                fftype(n_fftypes)%env = residue%atom(i)%env(1:residue%atom(i)%nbonds)
                fftype(n_fftypes)%name = adjustl(residue%atom(i)%fftype)
                fftype(n_fftypes)%charge = residue%atom(i)%q
                fftype(n_fftypes)%items = 1
                fftype(n_fftypes)%atom_list(fftype(n_fftypes)%items) = i
            endif

            q_tot_prev = q_tot_prev + residue%atom(i)%q

        enddo
        !Average charges
        do j=1,n_fftypes
            fftype(j)%charge = fftype(j)%charge/fftype(j)%items
        enddo


        !2. Reassing charges to have the same for the same atomtypes. If the actual
        !   charge is very different from the old one, it is not reassigned (warning)
        !   The loop is based in atomtypes groups
        q_tot=0.
        do j=1,n_fftypes

            do jj=1,fftype(j)%items
                i=fftype(j)%atom_list(jj)
                if ( abs(residue%atom(i)%q - fftype(j)%charge) > 0.07 ) then
                    write(dummy_char,*) i
                    write(dummy_char2,*) residue%atom(i)%q
                    write(dummy_char3,*) fftype(j)%charge
                    call alert_msg("error", "Charge on atom "//trim(adjustl(dummy_char))// &
                                           " (" //trim(adjustl(dummy_char2))//             &
                                           ") deviates more than 0.07 a.u. from the "//       &
                                           "average charge for this atom type (" //        &
                                           trim(adjustl(dummy_char3))//"). Hence, it is"// &
                                           " not updated"                                  &
                                    )
                    !Just round to the 3rd decimal
                    residue%atom(i)%q = FLOAT(INT(residue%atom(i)%q * 1000.0 + 0.5)) / 1000.0
                elseif ( abs(residue%atom(i)%q - fftype(j)%charge) > 0.05 ) then
                    write(dummy_char,*) i
                    write(dummy_char2,'(F7.3)') residue%atom(i)%q
                    write(dummy_char3,'(F7.3)') fftype(j)%charge
                    call alert_msg("note", "Charge on atom "//trim(adjustl(dummy_char))// &
                                           " (" //trim(adjustl(dummy_char2))//            &
                                           ") deviated more than 0.05 a.u. from the "//   &
                                           "average charge for this atom type (" //       &
                                           trim(adjustl(dummy_char3))//")."  &
                                   )
                     !Update and round to the 3rd decimal
                     residue%atom(i)%q = FLOAT(INT(fftype(j)%charge * 1000.0 + 0.5)) / 1000.0

                else 
                    !Update and round to the 3rd decimal
                    residue%atom(i)%q = FLOAT(INT(fftype(j)%charge * 1000.0 + 0.5)) / 1000.0
                endif
                q_tot = q_tot + residue%atom(i)%q
            enddo

        enddo


        ! 3. Compensate rounding errors and differences after averaging charges
        !    The difference is allocated over atoms with larger charges (relative
        !    error is lower)
        ! Reorder atomtype list according to absolute charges (take adv. of unused array space fro auxiliars)
        do i=1,n_fftypes-1
            do j=i+1,n_fftypes
                if ( (abs(fftype(i)%charge) - abs(fftype(j)%charge)) < 0. ) then
                    fftype(n_fftypes+1)=fftype(i)
                    fftype(i)=fftype(j)
                    fftype(j)=fftype(n_fftypes+1)
                 endif
            enddo
        enddo

        ! Number of extra charges to compensate
        delta_q=q_tot - q_tot_prev
        n_extra_charges = int( abs(delta_q*1000) + 0.5 )
        if (n_extra_charges > residue%natoms ) &
          call alert_msg("error", "Rounding error larger that 0.001a.u. per atom when reassigning charges")

        ! Assing extra charges over largelier charged atoms
        q_tot=0.
        n_compens=0
        do while ( n_compens < n_extra_charges )

            do j=1,n_fftypes
                do jj=1,fftype(j)%items
                    i=fftype(j)%atom_list(jj)
                    if ( n_compens < n_extra_charges ) then
                        residue%atom(i)%q = residue%atom(i)%q - 0.001 * delta_q/abs(delta_q)
                        n_compens=n_compens+1
                    endif
                    q_tot = q_tot + residue%atom(i)%q            
                enddo
            enddo

        enddo

        ! Fix aliphatic H charges to a given value (CHARMM correction)
        do i=1,residue%natoms
            residue%atom(i)%fftype=adjustl(residue%atom(i)%fftype)
            !Aliphatic hydrogen
            if ( residue%atom(i)%fftype(1:2) == "HA" ) then
                delta_q = residue%atom(i)%q - q_Hal
                residue%atom(i)%q = q_Hal  
                !Transfer extra charge to bonded carbon atom      
                j= residue%atom(i)%connect(1)
                residue%atom(j)%q = residue%atom(j)%q + delta_q
            !Aromatic hydrogen (also include aquenil H)
            else if ( (residue%atom(i)%fftype(1:2) == "HP") .or. &
                      (residue%atom(i)%fftype(1:2) == "HE") ) then
                delta_q = residue%atom(i)%q - q_Har
                residue%atom(i)%q = q_Har  
                !Transfer extra charge to bonded carbon atom      
                j= residue%atom(i)%connect(1)
                residue%atom(j)%q = residue%atom(j)%q + delta_q
            endif
        enddo
                    
        return               

    end subroutine refine_charges

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine symmetrize_charges(molec)

    !Done with molec not residue. Este jaleo tiene que acabar!! 

    use symmetry_mod

    implicit none

    common /GLOBAL_OPTS/ do_refine_charges

    type(str_system) :: molec
    integer :: ifreq=0, i, j
    integer,dimension(1:500) :: isym
    double precision :: av_charge
    logical :: do_refine_charges=.false.

    open(10,file="assignement.dat")

    call read_log(5,"elec_pot",ifreq,molec)
    call symm_atoms(molec,isym)

    !Read list of assignements (manual)
    do i=1,molec%natoms
        read(10,*) j, isym(j)
    enddo

    do i=1,molec%natoms
        print'(I5,2X,F8.3)', i, molec%atom(i)%q
    enddo
print*, "---------------------"

    do i=1,molec%natoms
        print*, i, isym(i)
        av_charge=(molec%atom(i)%q+molec%atom(isym(i))%q)/2.d0
        molec%atom(i)%q=av_charge
        molec%atom(isym(i))%q=av_charge
    enddo

    do i=1,molec%natoms
        print'(2I5,2X,F8.3)', i, isym(i), molec%atom(i)%q
    enddo

    return

    end subroutine symmetrize_charges

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine build_ff(residue)

        common /GLOBAL_OPTS/ do_refine_charges, I_DB, rename_atoms

        type(str_residue),intent(inout) :: residue

        !Mantained for easy life
        integer,dimension(1000,6) :: connect
        integer,dimension(1000) :: nbonds

        !locales
        integer,dimension(1000) :: non_H, list_H
        integer,dimension(6) :: hbond
        integer :: n_nonH, n_H, imp_bonds, numH
        integer :: I_DB
        !Contadores
        integer :: i,j,k, jj
        !Aux
        integer :: id
        character(len=150) ::  dummy_char, dummy_char2
        !We use atom name to save the adjustl character
        character(len=6) :: H_fftype, atom_name

        !Switching features on/off
        logical ::  do_refine_charges, &
                    rename_atoms 


        !Form the connexion matrix and nbonds vector
        do i=1,residue%natoms
            nbonds(i) = residue%atom(i)%nbonds
            connect(i,:) = residue%atom(i)%connect(:)
        enddo

        ! Identify non-H atoms first
        id=0
        do i=1,residue%natoms
            !Avoid making adjustl directly on the variable!
            !residue%atom(i)%name = adjustl(residue%atom(i)%name)
            atom_name = adjustl(residue%atom(i)%name)
            if ( atom_name(1:1) /= "H") then
!            if (trim(adjustl(residue%atom(i)%name)) /= "H") then
                id=id+1
                non_H(id)=i
            endif
        enddo
        n_nonH=id

        !
        ! The job is done focussed on non-H atoms. For each of them, the name is changed
        ! to Element_Name//id, where is is the number of non-H atom. Then linked H are 
        ! treated at the same time
        !
        do i=1,n_nonH
            id=non_H(i)

            !COUNT BONDED H
            numH=0               
            do j=1,nbonds(id)
                !Avoid making adjustl directly on the variable!
                !residue%atom(connect(id,j))%name=adjustl(residue%atom(connect(id,j))%name)
                atom_name = adjustl(residue%atom(connect(id,j))%name)
                if ( atom_name(1:1) == "H" ) then
                    numH=numH+1
                    hbond(numH)=j
                endif
            enddo

            !ASSIGN ATOMTYPES (before the atomname has been changed!
            call fftype_db(residue%atom(id),numH,              &
                           residue%atom(id)%fftype,H_fftype,   &
                           I_DB, id )
            if (numH /= 0) &
              residue%atom(connect(id,hbond(1:numH)))%fftype=H_fftype          

            !RENAME ATOMS 
            if (rename_atoms) then
                write(dummy_char,*) id !dummy_char holds "id" as character
                residue%atom(id)%name=trim(adjustl(residue%atom(id)%name))//trim(adjustl(dummy_char))
                !Also rename bound H            
                do j=1,numH
                    !dummy_char2 holds "h_index" as character
                    write(dummy_char2,*) j
                    residue%atom(connect(id,hbond(j)))%name="H"//trim(adjustl(dummy_char))//trim(adjustl(dummy_char2))
                enddo
            endif

        enddo

        !RETRIEVE ENVIRONMENT (according to neighbouring atomtypes)
        do i=1,residue%natoms
            do j=1,residue%atom(i)%nbonds
                jj=residue%atom(i)%connect(j)
                residue%atom(i)%env(j) = adjustl(residue%atom(jj)%fftype)
            enddo
        enddo

        !Order the env matrix (to allow trustful comparisons)
        do i=1,residue%natoms
           do j=1,residue%atom(i)%nbonds-1
                do k=j+1,residue%atom(i)%nbonds
                    ! Lexical Greater Than
                    if ( lgt(residue%atom(i)%env(k),residue%atom(i)%env(j)) ) then
                        residue%atom(i)%env(residue%atom(i)%nbonds+1) = residue%atom(i)%env(j)
                        residue%atom(i)%env(j) = residue%atom(i)%env(k)
                        residue%atom(i)%env(k) = residue%atom(i)%env(residue%atom(i)%nbonds+1)
                        residue%atom(i)%env(residue%atom(i)%nbonds+1) = ""
                    endif
                enddo
            enddo
        enddo

        if (do_refine_charges) call refine_charges(residue)

    end subroutine build_ff

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine gen_bonded(residue)

        type(str_residue),intent(inout) ::  residue
        
        !Counters
        integer :: i_1, i_2, i_3, i_4, &
                   j, k, l
        integer :: n_angle, n_dihed

        !Ordering stuff
        integer :: i, min_dihed_index, imin
        integer,dimension(1:4) :: aux

        ! Stepwise walk along the connection network to retrieve
        ! bonding relations (angle, dihedral and 1-4 interaction)
        ! 13/02/12: angles, dihedrals and pairs do not belong to str_atom anymore
        n_angle=0
        n_dihed=0
        do i_1=1,residue%natoms 
            ! i_1 is the epicenter-atom

            !Step 1 (bond)
            do j=1,residue%atom(i_1)%nbonds
                i_2=residue%atom(i_1)%connect(j)

                !Step 2 (angle)
                do k=1,residue%atom(i_2)%nbonds
                    i_3=residue%atom(i_2)%connect(k)
                    if ( i_3 == i_1 ) cycle !do not go back!
                    n_angle=n_angle + 1
                    !Get the angle
                    residue%geom%angle(n_angle,1)=i_1
                    residue%geom%angle(n_angle,2)=i_2
                    residue%geom%angle(n_angle,3)=i_3
                    !Remove repeated angles (take just those with incremental index)
                    if ( residue%geom%angle(n_angle,1) > residue%geom%angle(n_angle,3) ) &
                       n_angle=n_angle - 1

                    !Step 3 (for dihedral and 1-4 int)
                    do l=1,residue%atom(i_3)%nbonds
                        i_4=residue%atom(i_3)%connect(l)
                        if ( i_4 == i_2 ) cycle !do not go back!
                        n_dihed=n_dihed + 1 ! = n_pair
                        !Get the dihedral&pair
                        residue%geom%dihed(n_dihed,1)=i_1
                        residue%geom%dihed(n_dihed,2)=i_2
                        residue%geom%dihed(n_dihed,3)=i_3
                        residue%geom%dihed(n_dihed,4)=i_4
                        residue%geom%pair(n_dihed,1)=i_1 !1-4 int
                        residue%geom%pair(n_dihed,2)=i_4 !1-4 int
                    if ( residue%geom%dihed(n_dihed,2) > residue%geom%dihed(n_dihed,3) ) &
                       n_dihed=n_dihed - 1

                    enddo !... leaving step 3

                enddo !... leaving step 2


            enddo !... leaving step 1

        enddo !... leaving step 0 (I'm free!)
        residue%geom%nangles=n_angle
        residue%geom%ndihed=n_dihed
        residue%geom%npairs=n_dihed

        !==================================
        !Order dihedrals by rotated bond:
        !==================================
        do i=1,n_dihed-1
            min_dihed_index=residue%geom%dihed(i,2)
            imin=i
            do j=i+1,n_dihed
                if (residue%geom%dihed(j,2) < min_dihed_index) then
                    min_dihed_index=residue%geom%dihed(j,2)
                    imin=j
                endif
            enddo
        
            aux=residue%geom%dihed(i,1:4)
            residue%geom%dihed(i,1:4)=residue%geom%dihed(imin,1:4)
            residue%geom%dihed(imin,1:4)=aux
        enddo   

        do i=1,n_dihed-1
            min_dihed_index=residue%geom%dihed(i,3)
            imin=i
            do j=i+1,n_dihed
                if (residue%geom%dihed(j,2) /= residue%geom%dihed(i,2)) exit
                if (residue%geom%dihed(j,3) < min_dihed_index) then
                    min_dihed_index=residue%geom%dihed(j,3)
                    imin=j
                endif
             enddo
        
             aux=residue%geom%dihed(i,1:4)
             residue%geom%dihed(i,1:4)=residue%geom%dihed(imin,1:4)
             residue%geom%dihed(imin,1:4)=aux

        enddo   
                    
        return

    end subroutine gen_bonded

    subroutine assign_masses(molec)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 1.0/March 2012)
        !==============================================================
        !Description
        ! Assign masses based on atom names
        !Arguments
        ! molec (str_system,inout): molecule
        !
        !Notes
        !(20/06/2012) Improved to get one or two letter symbols 
        !==============================================================

        type(str_system),intent(inout) :: molec

        integer :: i
        character(len=2) :: atname

        do i=1,molec%natoms
                atname = adjustl(molec%atom(i)%name)
                !TODO: call change_case("upper",atname(1:2),atname(1:2))
                !First two-letter symbols
                if (atname(1:2) == "HE") then
                    molec%atom(i)%mass=4.003
                elseif (atname(1:2) == "LI") then
                    molec%atom(i)%mass=6.941
                elseif (atname(1:2) == "BE") then
                    molec%atom(i)%mass=9.012
                elseif (atname(1:2) == "MG") then
                    molec%atom(i)%mass=24.30500
                ! To be completed...
                !elseif (atname(1:2) == "") then
                !    molec%atom(i)%mass=
                !Then one-letter symbols
                elseif (atname(1:1) == "H") then
                    molec%atom(i)%mass=1.008
                elseif (atname(1:1) == "B") then
                    molec%atom(i)%mass=10.611
                elseif (atname(1:1) == "C") then
                    molec%atom(i)%mass=12.011
                elseif (atname(1:1) == "N") then
                    molec%atom(i)%mass=14.00700
                elseif (atname(1:1) == "O") then
                    molec%atom(i)%mass=15.9994
                elseif (atname(1:1) == "S") then
                    molec%atom(i)%mass=32.060
                elseif (atname(1:1) == "F") then
                    molec%atom(i)%mass=18.99840
                ! To be completed...
                !elseif (atname(1:1) == "") then
                !    molec%atom(i)%mass=
                else
                    call alert_msg("warning","Don't know how to assign mass to "//atname//" Set to zero.")
                    molec%atom(i)%mass=0.00
                endif
        enddo

        return

    end subroutine assign_masses


   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine guess_connect(molec)

        implicit none
        type(str_system),intent(inout) :: molec
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
                    i_cnx=i_cnx+1
                    molec%atom(i)%connect(i_cnx)=j
                endif

            enddo
            molec%atom(i)%nbonds=i_cnx
        enddo
   
        return

    end subroutine guess_connect

    subroutine guess_connect_res(residue)

    !The same, but worknig on residues (so, now the function "bond_length_db" is independent

        implicit none
        type(str_residue),intent(inout) :: residue
        !local
        real :: av_len, dist
        integer :: i,j, i_cnx

        !InOut, better reassigned
        residue=residue

        !Loop over all atoms. This is not optimal, distance matrix is symmetric
        !(nor it is not costly, anyway)
        !TODO: Change it to retrieve symmetric from already calc.
        i_cnx=0
        do i=1,residue%natoms
            residue%atom(i)%connect(:)=0
            i_cnx=0

            do j=1,residue%natoms
                if (i == j) cycle

                !Calculate the distance:
                dist = (residue%atom(i)%x - residue%atom(j)%x)**2 &
                     + (residue%atom(i)%y - residue%atom(j)%y)**2 &
                     + (residue%atom(i)%z - residue%atom(j)%z)**2
                dist = sqrt(dist)

                !Get default bond length from database

                av_len=bond_length_db(residue%atom(i),residue%atom(j))
                ! Criterium: dist < av. length +10% --increased from 5% (17/12/12) due to H2O2
                ! a more sophiticated data base might include hibridization
                av_len=av_len*1.1
                if (dist < av_len) then
                    i_cnx=i_cnx+1
                    residue%atom(i)%connect(i_cnx)=j
                endif

            enddo
            residue%atom(i)%nbonds=i_cnx
        enddo

        return

    end subroutine guess_connect_res


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
        !=======================================================================

        implicit none
        type DB_bond_length
            character(len=1) :: A, B
            real :: length 
        end type DB_bond_length

        type(str_atom),intent(in) :: atom1, atom2
        !Local atoms
        type(str_atom) :: Latom1, Latom2
        real :: av_len
        !local
        character(len=14),dimension(100) :: database
        !new structure not really needed (not hurting anyway)
        type(DB_bond_length) :: db_entry
        integer :: n_entries, i
        logical :: external_DB=.false.

        !Set local atoms
        Latom1%name = adjustl(atom1%name)
        Latom2%name = adjustl(atom2%name)


        ! DATABASE GENERATION
        if (external_DB) then
            print*, "External DB not yet supported"
        else
            n_entries=30
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
                    "O    H    1.07",      & !7
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
                    "MG   N    2.15",      & !22  Larger than gv standard (2.06)
                    "Ni   O    2.10",      & !23  To be revised
                    "Ni   N    2.10",      & !24  To be revised
                    "Cl   C    1.80",      & !25  Larger than gv standard (1.76)
                    "P    O    2.71",      & !26
                    "P    H    1.35",      & !27
                    "P    N    2.71",      & !28 ! Same as P-O
                    "Pt   O    2.50",      & !29 !A ojo
                    "Pt   N    2.50"       & !30
                    /)                       
        endif

        ! SELECTION LOOPS
        av_len=0.
        do i=1,n_entries 
            !Get database entry
            read(database(i),*) db_entry%A, db_entry%B, db_entry%length
            !Compare with input
            if ( ( Latom1%name(1:1) == db_entry%A .and. &
                   Latom2%name(1:1) == db_entry%B )     &
                .or.                                            &
                 ( Latom1%name(1:1) == db_entry%B .and. &
                   Latom2%name(1:1) == db_entry%A )     &
                ) av_len=db_entry%length
       enddo     

! if ( av_len == 0.) then
!     call alert_msg("fatal",Latom1%name(1:1)//" and "//Latom2%name(1:1)//" missing in the DB")
! endif

        return
    end function bond_length_db

    subroutine aa2ua(molecAA,molecUA)

        type(str_system),intent(in)  :: molecAA
        type(str_system),intent(out) :: molecUA

        !Local
        integer :: i, j, k, atom2
        character(len=6) :: atname, atname2

        !This first assignement prevents a segmentation fault arrising from molecUA
        !allocation. Maybe an explicit "interface" is strongly adviced here
        molecUA = molecAA

        k=0
        do i=1,molecAA%natoms

            atname = adjustl(molecAA%atom(i)%name)
            !Non polar hydrogens are not included in the UA model
            if (atname(1:1) == "H") then
                atom2 = molecAA%atom(i)%connect(1)
                atname2 = adjustl(molecAA%atom(atom2)%name)
                if (atname2(1:1) == "C" ) cycle
            endif

            k = k + 1
            molecUA%atom(k) = molecAA%atom(i)

            !Charge and mass from non polar hydrogens go to the corresponding carbon
            if (atname(1:1) == "C") then
                do j=1,molecAA%atom(i)%nbonds
                    atom2 =  molecAA%atom(i)%connect(i)
                    atname2 = adjustl(molecAA%atom(atom2)%name)
                    if (atname2(1:1) == "H" ) &
                     molecUA%atom(k)%q    = molecUA%atom(k)%q    + molecAA%atom(atom2)%q
                     molecUA%atom(k)%mass = molecUA%atom(k)%mass + molecAA%atom(atom2)%mass
                enddo
            endif

        enddo
        molecUA%natoms = k
            

        return

    end subroutine aa2ua

    subroutine aa2ua_res(residueAA,residueUA)

        type(str_residue),intent(in)  :: residueAA
        type(str_residue),intent(out) :: residueUA

        !Local
        integer :: i, k, atom2
        character(len=6) :: atname, atname2

        !This first assignement prevents a segmentation fault arrising from molecUA
        !allocation. Maybe an explicit "interface" is strongly adviced here
        residueUA = residueAA

        k=0
        do i=1,residueAA%natoms

            atname = adjustl(residueAA%atom(i)%name)
            if (atname(1:1) == "H") then
                atom2 = residueAA%atom(i)%connect(1)
                atname2 = adjustl(residueAA%atom(atom2)%name)
                if (atname2(1:1) == "C" ) cycle
            endif

            k = k + 1
            residueUA%atom(k) = residueAA%atom(i)

        enddo
        residueUA%natoms = k


        return

    end subroutine aa2ua_res

end module ff_build

