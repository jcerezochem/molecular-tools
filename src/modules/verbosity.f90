module verbosity

    !==============================================================
    ! This code is part of MOLECULAR_TOOLS 
    !==============================================================
    ! Description
    !  This module acts as a COMMON block with a common structure
    !  to set verbosity level in all the programs of the distribution
    ! 
    ! Level       Behaviour
    ! 0           Do not write anything
    ! 1           Basic verbose level (default)
    ! 2           Debug verbose level
    ! 3+          Pedantic verbose level
    !-----------------------------------------------------

    integer, save :: verbose = 1
    integer, save :: verbose_current
    integer, save :: verbose_mute_depth = 0

    contains

    subroutine verbose_mute()

        if (verbose>1) return
        if (verbose_mute_depth==0) verbose_current = verbose
        verbose=0
        verbose_mute_depth = verbose_mute_depth + 1

         return

    end subroutine verbose_mute

    subroutine verbose_continue()

        if (verbose>1) return
        verbose_mute_depth = verbose_mute_depth - 1
        if (verbose_mute_depth==0) verbose = verbose_current

        return

    end subroutine verbose_continue

end module verbosity
