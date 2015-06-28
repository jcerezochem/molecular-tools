module psi4_manage

    !==============================================================
    ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
    !==============================================================
    ! Description
    !  This MODULE contains subroutines to get molecular information 
    ! Note: all subroutines rewind the file after using it
    ! Notes
    !  Since variables molec, residue... are allocatable, they should
    !  be passed to the sr even if they are output to get the allocation.
    !==============================================================

    !Common modules for all subroutines
    use structure_types
    use constants
    use alerts

    contains

    subroutine read_psi_geom(unt,molec)

        use line_preprocess

        integer,intent(in)::unt
        type(str_resmol),intent(inout)::molec

        !local
        integer::i, natoms
        character :: cnull
        character(len=10) :: DIST_UNITS
        character(len=280) :: line

        !Rewind the file in case it was not in the gegining
        rewind(unt)

        ! Search the Geometry in the output 
        do 
            read(unt,'(X,A)',IOSTAT=IOstatus) line
            ! Two possible scenarios while reading:
            ! 1) End of file
            if ( IOstatus < 0 ) call alert_msg("fatal","Unexpected end of file while searching Geometry")
            ! 2) Geom section reached
            if ( INDEX(line,"==> Geometry <==") /= 0 ) then
                exit
            endif
        enddo

        ! Overpass lines until reaching the target table
        do while (INDEX(line,"Center") == 0)
            read(unt,'(X,A)',IOSTAT=IOstatus) line
            if ( INDEX(line,"Full point group:") /= 0 ) then
                call split_line(line,":",cnull,line)
                read(line,*) molec%PG
            else if ( INDEX(line,"Geometry (in Angstrom)") /= 0 ) then
                DIST_UNITS="Angstrom"
            else if ( INDEX(line,"Geometry (in Bohr)") /= 0 ) then
                DIST_UNITS="Bohr"
            endif
        enddo
        !Read Table lines
        read(unt,'(X,A)',IOSTAT=IOstatus) line
        !Start reading geometry
        i=0
        read(unt,'(X,A)',IOSTAT=IOstatus) line
        do while (len_trim(line) /= 0)
            i=i+1
            read(line,*) molec%atom(i)%name, &
                         molec%atom(i)%x,    &
                         molec%atom(i)%y,    &
                         molec%atom(i)%z
            read(unt,'(X,A)',IOSTAT=IOstatus) line
            if ( IOstatus < 0 ) call alert_msg("fatal","Unexpected end of file while reading Geometry")
        enddo
        molec%natoms = i
        
        !Manage lenght units
        if (adjustl(DIST_UNITS) == "Bohr") then
            molec%atom(1:i)%x = molec%atom(1:i)%x * BOHRtoANGS
            molec%atom(1:i)%y = molec%atom(1:i)%y * BOHRtoANGS
            molec%atom(1:i)%z = molec%atom(1:i)%z * BOHRtoANGS
        endif

        rewind(unt)
        return

    end subroutine read_psi_geom

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine read_psi_hess(unt,N,Hess,error)

        use line_preprocess

        integer,intent(in)::unt
        integer,intent(in)::N
#ifdef DOUBLE
        real(8),dimension(:,:),intent(out)::Hess
#else
        real(4),dimension(:,:),intent(out)::Hess
#endif 
        integer,intent(out) :: error

        !local
        integer:: i, j, k,             &
                  iblock, jini, jfin,  &
                  Nblocks
        character :: cnull
        character(len=280) :: line

        error=1
        !Rewind the file in case it was not in the gegining
        rewind(unt)

        ! Search the Geometry in the output 
        do 
            read(unt,'(X,A)',IOSTAT=IOstatus) line
            ! Two possible scenarios while reading:
            ! 1) End of file
            if ( IOstatus < 0 ) call alert_msg("fatal","Unexpected end of file while searching Hessian")
            ! 2) Geom section reached
            if ( INDEX(line,"Force Constants in cartesian coordinates.") /= 0 ) then
                exit
            endif
        enddo


        !Hessian elements arranged in blocks of 5 columns each
        Nblocks = N/5
        if (N /= Nblocks*5) Nblocks=Nblocks+1
        do iblock=1,Nblocks
            !Three first lines useless
            read(unt,'(A)') line
            read(unt,'(A)') line
            read(unt,'(A)') line
            jini=(iblock-1)*5+1
            jfin=min(N,iblock*5)
            do i=1,N
                read(unt,*) k, Hess(i,jini:jfin)
            enddo
        enddo

        error=0
        return

    end subroutine read_psi_hess

end module psi4_manage
