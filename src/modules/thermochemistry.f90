module thermochemistry
!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012!

    !==============================================================
    ! This code is part of MOLECULAR_TOOLS 
    !==============================================================
    ! Description
    !  This MODULE contains subroutines to perform the vibrational
    !  analysis, by properly Diagonalizing the Hessian
    ! 
    !==============================================================

    !Common declarations:
    !===================
    use matrix
    use matrix_print
    use verbosity
    use constants
    use alerts

    implicit none

    contains

    subroutine thermo(Nat,Nvib,X,Y,Z,Mass,Freq,Temp)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS
        !==============================================================
        !Description
        ! Diagonalizes a Hessian after mass-weighing and translation to 
        ! the internal frame defined by satifying the Eckart-Saytvez conditions.
        !
        !Arguments
        ! Nat     (inp) int /scalar   Number of atoms
        ! Nvib    (inp) int /scalar   Number of vibrational degrees of freedom
        ! X,Y,Z   (inp) real/vectors  Coordinate vectors (ANGSTRONG)
        ! Mass    (inp) real/vector   Atomic masses (AMU)
        ! Freq    (inp) real/vector   Frequencies (cm-1)
        !
        !
        !==============================================================

        !Approximate zero
        real(kind=8),parameter :: ZERO=1.d-10

        integer,intent(in)                   :: Nat
        integer,intent(in)                   :: Nvib
        real(kind=8),dimension(:),intent(in) :: X,Y,Z
        real(kind=8),dimension(:),intent(in) :: Mass
        real(kind=8),dimension(:),intent(in) :: Freq
        real(kind=8),intent(in)              :: Temp

        !Local
        ! Counters
        integer :: i, j, k, ii, jj
        ! Auxiliar scalars and arrays
        real(kind=8)                :: MassTot
        real(kind=8),dimension(3)   :: R, RCOM, MIp, Arot
        real(kind=8),dimension(3,3) :: MI, Xrot

        ! Thermodynamic properties
        real(8) :: Pressure, kbwn, kbau
        real(8) :: ThetaRx, ThetaRy, ThetaRz, ThetaV
        real(8) :: qt, qr, qv
        real(8) :: St, Sr, Sv, S
        real(8) :: Eint, G, H, zpe

        print*, "====================="
        print*, " THERMOCHEMISTRY     "
        print*, "====================="

        !-----------------------------------
        ! Initial manipulations
        !-----------------------------------
        ! Boltzmann constant in cm-1
        kbwn=kboltz/clight/plank * 1.d-2
        kbau=kboltz/HARTtoJ
        Pressure = 1.01325d5 

        print'(2X,A,F8.2)', "Temperature (K):", Temp
        print'(2X,A,F8.2)', "Pressure (atm) :", Pressure/1.01325d5
        print*, ""

        !Get COM 
        RCOM(1:3) = 0.d0
        MassTot   = 0.d0
        do i=1,Nat
            RCOM(1) = RCOM(1) + X(i)*Mass(i)
            RCOM(2) = RCOM(2) + Y(i)*Mass(i)
            RCOM(3) = RCOM(3) + Z(i)*Mass(i)
            MassTot = MassTot + Mass(i)
        enddo
        RCOM(1:3) = RCOM(1:3)/MassTot
        
        !Get moment of intertia
        MI=0.d0
        do i=1,Nat
            R=(/X(i)-RCOM(1),Y(i)-RCOM(2),Z(i)-RCOM(3)/)
            !diag
            MI(1,1)=MI(1,1)+Mass(i)*(R(2)**2+R(3)**2)
            MI(2,2)=MI(2,2)+Mass(i)*(R(1)**2+R(3)**2)
            MI(3,3)=MI(3,3)+Mass(i)*(R(1)**2+R(2)**2)
            !off-diag
            MI(2,1)=MI(2,1)-Mass(i)*(R(2)*R(1))
            MI(3,1)=MI(3,1)-Mass(i)*(R(3)*R(1))
            MI(3,2)=MI(3,2)-Mass(i)*(R(3)*R(2))
        enddo
        do i=1,3
          do j=1,i-1
              MI(j,i) = MI(i,j)
          enddo
        enddo
        ! Diagonalize to get the rotation to the principal axes
        call diagonalize_full(MI(1:3,1:3),3,Xrot(1:3,1:3),MIp(1:3),"lapack")
        ! Principal moment of inertia
        MIp(1) = MIp(1) * ANGStoM**2 * AMUtoKG
        MIp(2) = MIp(2) * ANGStoM**2 * AMUtoKG
        MIp(3) = MIp(3) * ANGStoM**2 * AMUtoKG
        ! Rotational constants
        Arot(1) = plank/8.d0/pi**2/MIp(1)
        Arot(2) = plank/8.d0/pi**2/MIp(2)
        Arot(3) = plank/8.d0/pi**2/MIp(3)
        ! Rotational temperatures
        ThetaRx = plank**2/8.d0/pi**2/MIp(1)/kboltz
        ThetaRy = plank**2/8.d0/pi**2/MIp(2)/kboltz
        ThetaRz = plank**2/8.d0/pi**2/MIp(3)/kboltz

        print'(2X,A)', "Rotational constants (GHz)"
        print'(3F12.5)', Arot(1)*1.d-9, Arot(2)*1.d-9, Arot(3)*1.d-9
        print'(2X,A)', "Rotational temperatures (K)"
        print'(3F12.5)', ThetaRx, ThetaRy,ThetaRz
        print*, ""

        ! Tr and Rot partition function
        qt = (dsqrt(2.d0*pi*MassTot*AMUtoKG*kboltz*Temp)**3/plank**3) * kboltz*Temp/Pressure
        qr = dsqrt(pi) * dsqrt(Temp**3) / dsqrt(ThetaRx*ThetaRy*ThetaRz)
        ! Tr and Rot entropy
        St = (dlog(qt) + 2.5d0)
        Sr = (dlog(qr) + 1.5d0)

        qv=1.d0
        zpe  = 0.d0
        Eint = 0.d0
        Sv   = 0.d0
        do i=1,Nvib
            if (Freq(i) < 0) then
                print*, i, Freq(i)
                call alert_msg("warning","Imaginary frequency ignored")
                cycle
            endif
            ThetaV = Freq(i)/kbwn
            qv   = qv * dexp(-ThetaV/2.d0/Temp)/(1.d0-dexp(-ThetaV/Temp))
            Eint = Eint + ThetaV*(0.5d0 + 1.d0/(dexp(ThetaV/Temp) - 1.d0))
            Sv   = Sv + ThetaV/Temp/(dexp(ThetaV/Temp) - 1.d0) - dlog(1.d0-dexp(-ThetaV/Temp))
            zpe  = zpe + 0.5d0 * Freq(i)
        enddo
        ! Sum tr+rot+vib (all quantitues need to be multiplied by kB)
        S    = St+Sr+Sv 
        H    = (Eint + 3.d0*Temp) + Temp
        G    = H - S*Temp

        print'(2X,A)', "Partition Functions"
        print'(2X,A,G15.4)', "qt =", qt
        print'(2X,A,G15.4)', "qr =", qr
        print'(2X,A,G15.4)', "qv =", qv
        print*, ""
        print'(2X,A)', "Contributions to entropy"
        print'(2X,A,F10.4)', "St (cm-1/K) =", St*kbwn
        print'(2X,A,F10.4)', "Sr (cm-1/K) =", Sr*kbwn
        print'(2X,A,F10.4)', "Sv (cm-1/K) =", Sv*kbwn
        print*, ""
        print'(2X,A)', "Thermal corrections (AU)"
        print'(2X,A,G15.6)', "Zero Point Energy =", zpe/autown
        print'(2X,A,G15.6)', "Internal Energy   =", (Eint+3.d0*Temp)*kbau 
        print'(2X,A,G15.6)', "Enthalpy          =", H*kbau         
        print'(2X,A,G15.6)', "Gibbs Free Energy =", G*kbau    
        print*, ""       
        print'(2X,A)', "Thermal corrections (AU) - vibrational only"
        print'(2X,A,I8)', "Number of vibrations:", Nvib
        H    = Eint + Temp
        G    = H - Sv*Temp
        print'(2X,A,G15.6)', "V-internal Energy   =", Eint*kbau
        print'(2X,A,G15.6)', "V-enthalpy          =", H*kbau         
        print'(2X,A,G15.6)', "V-gibbs Free Energy =", G*kbau   

        return

    end subroutine thermo


end module thermochemistry