module fcc_manage

    !==============================================================
    ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
    !==============================================================
    ! Description
    !-----------------------------------------------------

    !Common modules for all subroutines
    use constants
    implicit none

    contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine read_fccstate_atoms(unt,system)

        !Input is in \AA.

        use structure_types

        integer,intent(in)::unt
        type(str_resmol),intent(inout)::system

        ! local
        integer :: i, natoms
        real(8) :: natoms_real 
        ! read
        character(len=250) :: line
        integer :: IOstatus


        if (system%natoms == -1) then
            !Guess the number of atoms here. Assume Nvib=3N-6
            i=0
            do 
                read(unt,'(A)',iostat=IOstatus) line
                if (len_trim(line) == 0 .or. IOstatus /= 0) exit
                i=i+1
            enddo

            natoms_real = (12.d0+dsqrt(12.0**2+36.d0*(6.d0+float(i))))/18.d0
            if (natoms_real - int(natoms_real) /= 0) then
                print*, "ERROR: cannot guess the number of atoms from the state", natoms, natoms_real
                stop
            endif
            natoms = int(natoms_real)
            system%natoms=natoms
            rewind(unt)
        endif

        do i=1,system%natoms
            read(unt,*) system%atom(i)%x         
            read(unt,*) system%atom(i)%y  
            read(unt,*) system%atom(i)%z
        enddo

        return

    end subroutine read_fccstate_atoms

    subroutine read_fccstate_geom(unt,Nat,X,Y,Z)

        integer,intent(in) :: unt
        integer,intent(in) :: Nat
        real(8),dimension(:),intent(out)  :: X,Y,Z
        ! Local
        integer :: i, j

        rewind(unt)

        do i=1,Nat 
            read(unt,*) X(i)
            read(unt,*) Y(i)
            read(unt,*) Z(i)
        enddo

        return

    end subroutine read_fccstate_geom

    subroutine read_fccstate_nm(unt,Nvib,Nat,Freq,L)

        integer,intent(in) :: unt
        integer,intent(in) :: Nvib,Nat
        real(8),dimension(:),intent(out)  :: Freq
        real(8),dimension(:,:),intent(out):: L
        ! Local
        integer :: i, j

        rewind(unt)

        do i=1,3*Nat 
            read(unt,*) 
        enddo

        do i=1,3*Nat
            do j=1,Nvib
                read(unt,*) L(i,j)
            enddo
        enddo

        do i=1,Nvib
            read(unt,*) Freq(i)
        enddo

        return

    end subroutine read_fccstate_nm

end module fcc_manage
