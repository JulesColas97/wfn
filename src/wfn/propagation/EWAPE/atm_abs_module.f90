module atm_abs_mod
    use kind_mod, only: wp

    public :: atm_abs_air, atm_abs_water

contains 
    function atm_abs_air(T0, p0, rh, f) result(alpha)
        ! function return the atmospheric attenuation of sound
        ! due to the thermo-viscous effects and relaxation of oxygen and nitrogen.
        ! References:   Salomons p.109-111

        real(wp), intent(in) :: T0 ! temperature K
        real(wp), intent(in) :: p0 ! static pressure Pascal 
        real(wp), intent(in) :: rh ! relative humidy in %
        real(wp), intent(in) :: f !  frequency

        real(wp) :: alpha ! atm_abs coef Delta L = -alpha R


        real(wp) ::  p0_ref=1.01325e05, T0_tpw=273.15, T0_ref=293.15
        real(wp) :: rho, tau, Csat, p0_sat, h, frN, fr0, b1, b2


        rho = p0/p0_ref
        tau = T0/T0_ref

        ! calculate saturation pressure
        Csat = -6.8346*(T0_tpw/T0)**1.261 + 4.6151
        p0_sat = p0_ref*10**Csat
        h = rh*p0_sat/p0  ! absolute humidity

        ! Scaled relaxation frequency for Nitrogen
        frN = rho*tau**(-1/2)*(9 + 280*h*exp(-4.17*(tau**(-1/3)-1)))

        ! scaled relaxation frequency for Oxygen
        frO = rho*(24 + 40400*h*(0.02+h)/(0.391+h))

        ! attenuation coefficient in dB/m
        b1 = 0.1068*exp(-3352/T0)/(frN + f**2/frN)
        b2 = 0.01275*exp(-2239.1/T0)/(frO + f**2/frO)
        alpha = 8.686*f**2*tau**(1/2)*(1.84e-11/rho + tau**(-3)*(b1 + b2))
    end function


    function atm_abs_water(f) result(alpha)
        ! from A simplified formula for viscous and chemical absorption in sea water
        ! Michael A. Ainslie; James G. McColm
        ! J. Acoust. Soc. Am. 103, 1671–1672 (1998)
        ! https://doi.org/10.1121/1.421258

        ! simplify , thorp's formula: https://dl.acm.org/doi/epdf/10.1145/1347364.1347373
        real(wp), intent(in) :: f
        real(wp) :: alpha, alpha_dB_km

        f_2= (f/1000.0_wp)**2
        if (f > 300) then 
            alpha_dB_km =  0.11 *  f_2/(1+f_2) + 44.0_wp * f_2/(4100 +f_2) + &
                2.75e-4_wp * f_2 + 0.003
        else 
            alpha_dB_km = 0.002_wp + 0.11*f_2/(1+f_2) +0.011 * f_2
        endif

        alpha = alpha_dB_km /1000
    end function

end module 
