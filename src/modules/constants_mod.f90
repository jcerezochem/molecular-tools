module constants

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
        !==============================================================
        !Description
        ! This module contains mathematical and physical constants given
        ! in SI units. Taken from http://physics.nist.gov/constants
        !==============================================================

        !CONSTANTS
#ifdef DOUBLE
        double precision, parameter :: &
#else
        real, parameter :: &
#endif
                           PI      = 4.0d0*datan(1.0d0),   &
                           clight  = 2.99792458D8,    &
                           plank   = 6.62606957D-34,  &
                           plankbar= 1.054571726D-34, &
                           boltz   = 1.3806488D-23,   &
                           NAv     = 6.02214129D23,   &
                           atmass  = 1.660538921D-27  


        !CONVERSION FACTORS
#ifdef DOUBLE
        double precision, parameter :: &
#else
        real, parameter :: &
#endif
                           BOHRtoAMS = 5.2917720859D-1, &
                           BOHRtoANGS= 5.2917720859D-1, &
                           UMAtoKG   = 1.66053873d-27,  &
                           UMAtoAU   = 1.82288839d3,    &
                           AUtoKG    = 9.10938291d-31,  &
                           BOHRtoM   = 5.291772083d-11, &
                           AMStoM    = 1.d-10,          &
                           ANGStoM   = 1.d-10,          &
                           HARTtoJ   = 4.3597482d-18,   &
                           HtoKCALM  = 627.5095d0,      &
                           CALtoJ    = 4.184,           &
                           HtoeV     = 27.2114




end module constants
