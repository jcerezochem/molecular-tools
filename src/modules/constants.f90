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
                           PI      = 3.14159265358979323846d0, & !
                           NAv     = 6.02214129D23,      & ! Avogadro number
                           clight  = 2.99792458D8,       & ! Speed of light
                           SL      = 2.99792458D8,       & ! (synonym to above)
                           plank   = 6.62606957D-34,     & ! Planck constant
                           plankbar= 1.054571726D-34,    & ! Planck constant over 2PI
                           kboltz  = 1.3806488D-23,      & ! Boltzman constant
                           boltz   = 1.3806488D-23,      & ! (synonym to above)
                           atmass  = 1.660538921D-27,    & ! Atomic mass
                           cvelau  = 137.0369


        !CONVERSION FACTORS
#ifdef DOUBLE
        double precision, parameter :: &
#else
        real, parameter :: &
#endif
                           BOHRtoAMS = 5.2917720859D-1, &
                           BOHRtoANGS= 5.2917720859D-1, &
                           UMAtoKG   = 1.66053873d-27,  &
                           UMAtoAU   = 1.82288839d3,    & ! deprecated (keep for comopatibility)
                           AMUtoAU   = 1.82288839d3,    &
                           AMUtoKG   = 1.66053873d-27,  &
                           AUtoKG    = 9.10938291d-31,  &
                           BOHRtoM   = 5.291772083d-11, &
                           BOHRtoNM  = 5.291772083d-2,  &
                           AMStoM    = 1.d-10,          &
                           ANGStoM   = 1.d-10,          &
                           HARTtoJ   = 4.3597482d-18,   &
                           HtoKCALM  = 627.5095d0,      &
                           CALtoJ    = 4.184,           &
                           HtoeV     = 27.2114,         &
                           autown    = 2.1947463068d5    !From FCclasses Freq from AU to cm-1



        !AtNum to element name conversion
        character(len=5),dimension(103) :: atname_from_atnum
        !This should be elsewhere (constants_mod?)
        data atname_from_atnum(1:103) &
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

        !AtNum to atom mass conversion
        ! Masses from IUPAC site
        real(8),dimension(118),parameter :: atmass_from_atnum = &
           (/1.008       ,  & !H   -- 1  
             4.002602    ,  & !He  -- 2  
             6.94        ,  & !Li  -- 3  
             9.0121831   ,  & !Be  -- 4  
             10.81       ,  & !B   -- 5  
             12.011      ,  & !C   -- 6  
             14.007      ,  & !N   -- 7  
             15.999      ,  & !O   -- 8  
             18.998403163,  & !F   -- 9  
             20.1797     ,  & !Ne  -- 10 
             22.98976928 ,  & !Na  -- 11 
             24.305      ,  & !Mg  -- 12 
             26.9815385  ,  & !Al  -- 13 
             28.085      ,  & !Si  -- 14 
             30.973761998,  & !P   -- 15 
             32.06       ,  & !S   -- 16 
             35.45       ,  & !Cl  -- 17 
             39.948      ,  & !Ar  -- 18 
             39.0983     ,  & !K   -- 19 
             40.078      ,  & !Ca  -- 20 
             44.955908   ,  & !Sc  -- 21 
             47.867      ,  & !Ti  -- 22 
             50.9415     ,  & !V   -- 23 
             51.9961     ,  & !Cr  -- 24 
             54.938044   ,  & !Mn  -- 25 
             55.845      ,  & !Fe  -- 26 
             58.933194   ,  & !Co  -- 27 
             58.6934     ,  & !Ni  -- 28 
             63.546      ,  & !Cu  -- 29 
             65.38       ,  & !Zn  -- 30 
             69.723      ,  & !Ga  -- 31 
             72.63       ,  & !Ge  -- 32 
             74.921595   ,  & !As  -- 33 
             78.971      ,  & !Se  -- 34 
             79.904      ,  & !Br  -- 35 
             83.798      ,  & !Kr  -- 36 
             85.4678     ,  & !Rb  -- 37 
             87.62       ,  & !Sr  -- 38 
             88.90584    ,  & !Y   -- 39 
             91.224      ,  & !Zr  -- 40 
             92.90637    ,  & !Nb  -- 41 
             95.95       ,  & !Mo  -- 42 
             97.0        ,  & !Tc  -- 43 
             101.07      ,  & !Ru  -- 44 
             102.9055    ,  & !Rh  -- 45 
             106.42      ,  & !Pd  -- 46 
             107.8682    ,  & !Ag  -- 47 
             112.414     ,  & !Cd  -- 48 
             114.818     ,  & !In  -- 49 
             118.71      ,  & !Sn  -- 50 
             121.76      ,  & !Sb  -- 51 
             127.6       ,  & !Te  -- 52 
             126.90447   ,  & !I   -- 53 
             131.293     ,  & !Xe  -- 54 
             132.90545196,  & !Cs  -- 55 
             137.327     ,  & !Ba  -- 56 
             138.90547   ,  & !La  -- 57 
             140.116     ,  & !Ce  -- 58 
             140.90766   ,  & !Pr  -- 59 
             144.242     ,  & !Nd  -- 60 
             145.0       ,  & !Pm  -- 61 
             150.36      ,  & !Sm  -- 62 
             151.964     ,  & !Eu  -- 63 
             157.25      ,  & !Gd  -- 64 
             158.92535   ,  & !Tb  -- 65 
             162.5       ,  & !Dy  -- 66 
             164.93033   ,  & !Ho  -- 67 
             167.259     ,  & !Er  -- 68 
             168.93422   ,  & !Tm  -- 69 
             173.045     ,  & !Yb  -- 70 
             174.9668    ,  & !Lu  -- 71 
             178.49      ,  & !Hf  -- 72 
             180.94788   ,  & !Ta  -- 73 
             183.84      ,  & !W   -- 74 
             186.207     ,  & !Re  -- 75 
             190.23      ,  & !Os  -- 76 
             192.217     ,  & !Ir  -- 77 
             195.084     ,  & !Pt  -- 78 
             196.966569  ,  & !Au  -- 79 
             200.592     ,  & !Hg  -- 80 
             204.38      ,  & !Tl  -- 81 
             207.2       ,  & !Pb  -- 82 
             208.9804    ,  & !Bi  -- 83 
             209.0       ,  & !Po  -- 84 
             210.0       ,  & !At  -- 85 
             222.0       ,  & !Rn  -- 86 
             223.0       ,  & !Fr  -- 87 
             226.0       ,  & !Ra  -- 88 
             227.0       ,  & !Ac  -- 89 
             232.0377    ,  & !Th  -- 90 
             231.03588   ,  & !Pa  -- 91 
             238.02891   ,  & !U   -- 92 
             237.0       ,  & !Np  -- 93 
             244.0       ,  & !Pu  -- 94 
             243.0       ,  & !Am  -- 95 
             247.0       ,  & !Cm  -- 96 
             247.0       ,  & !Bk  -- 97 
             251.0       ,  & !Cf  -- 98 
             252.0       ,  & !Es  -- 99 
             257.0       ,  & !Fm  -- 100
             258.0       ,  & !Md  -- 101
             259.0       ,  & !No  -- 102
             262.0       ,  & !Lr  -- 103
             267.0       ,  & !Rf  -- 104
             270.0       ,  & !Db  -- 105
             269.0       ,  & !Sg  -- 106
             270.0       ,  & !Bh  -- 107
             270.0       ,  & !Hs  -- 108
             278.0       ,  & !Mt  -- 109
             281.0       ,  & !Ds  -- 110
             281.0       ,  & !Rg  -- 111
             285.0       ,  & !Cn  -- 112
             286.0       ,  & !Nh  -- 113
             289.0       ,  & !Fl  -- 114
             289.0       ,  & !Mc  -- 115
             293.0       ,  & !Lv  -- 116
             293.0       ,  & !Ts  -- 117
             294.0          & !Og  -- 118
            /)

    contains

    !Helper functions to access mass, AtNum and names

    function atnum_from_atname(name) result(Zat)

        character(len=*),intent(in) :: name
        integer                     :: Zat 
        !local  
        integer                     :: i, n

        n = size(atname_from_atnum)
        do i=1,n
            if (adjustl(name)==atname_from_atnum(i)) then
                Zat = i
                exit
            endif
        enddo

        return

    end function atnum_from_atname


    function atmass_from_atname(name) result(mass)

        character(len=*),intent(in) :: name
        real(8)                     :: mass 
        !local  
        integer                     :: i, n

        mass=0.d0
        n = size(atmass_from_atnum)
        do i=1,n
            if (adjustl(name)==atname_from_atnum(i)) then
                mass = atmass_from_atnum(i)
                exit
            endif
        enddo
        if (mass == 0.d0) then
            print*, "ERROR: mass cannot be determined for: "//name
            stop
        endif

        return

    end function atmass_from_atname

    subroutine atominfo_from_atmass(mass,Zat,name) 
        !================================================
        ! Description
        ! From mass value, get atom%name and atom%Zat
        ! So we can still use old FCclasses input format
        ! We need atomnames to use guess_connect
        !================================================

        real(8),intent(in)           :: mass
        integer,intent(out)          :: Zat
        character(len=*),intent(out) :: name
        !Local
        integer :: i, n

        name = "XX"
        n = size(atmass_from_atnum)
        do i=1,n
            if (abs(mass-atmass_from_atnum(i)) < 1.d-1) then
                name = atname_from_atnum(i)
                Zat  = i
            endif
        enddo
        if (adjustl(name) == "XX") then
            write(0,'(A,F12.6)') "ERROR: atom info could not be set from mass: ",mass
            stop
        endif

        return

    end subroutine atominfo_from_atmass


end module constants
