module impedance_model

    public :: compute_beta 
    public :: compute_beta_porous_layer
    public :: compute_beta_variable_porosity

contains 

    subroutine compute_beta(mode, f, imp, var0)
        use types
        integer, intent(in)            :: mode
        real, intent(in)               :: f
        type(impedance), intent(inout) :: imp
        type(variables_0), intent(in)  :: var0

        if (mode == 0) then 
            imp%beta = 0
        else if (mode == 1) then 
            call compute_beta_variable_porosity(f, imp%sigmae, imp%alphae, var0%c, &
            var0%gamma, var0%rho, imp%beta)
        else if (mode == 2) then 
            call compute_beta_porous_layer(f, imp%Rs, imp%Om, imp%q2, imp%e, &
                var0%T, var0%P, var0%gamma, imp%beta)
        end if 
    end subroutine

    subroutine compute_beta_variable_porosity(f, sigmae, alphae, c0,  gamma, rho0, beta)
        real, intent(in) :: f, sigmae, alphae, c0, gamma, rho0
        complex, intent(out) :: beta
        real :: omega, pi=3.141592653589793
        complex :: Zs, eye=(0.0, 1.0) 
        
        omega = 2*pi*f
        Zs = sqrt(4*sigmae/(-eye*omega*gamma*rho0)) + c0*alphae/(-eye*omega*4*gamma)
        beta = 1/Zs
    end subroutine


    subroutine compute_beta_porous_layer(f,Rs,Om,q2,e, T0, P0, gamma, beta)
        real, intent(in) :: f, Rs, Om, q2, e, gamma, T0, P0
        complex, intent(out) :: beta
        real :: omega,k0, pi=3.141592653589793
        real :: Np, fu, ft, rGP, c0, rho0
        complex :: Zc, pg, kg, kc, rhoc, n1, m1,   eye=(0.0, 1.0) 
    
    ! onstantes
    Np = 0.7 ! nombre de Prandtl
    rGP = 287
    
    c0 = sqrt(gamma*rGP*T0)
    rho0 = P0/(rGP*T0)

    k0 = 2*pi*f/c0
    omega = 2*pi*f
    
    !! Modele Hamet & Berengier
    fu = Om*Rs/(2*pi*rho0*q2) ! effet de la viscosite
    ft = Rs/(2*pi*rho0*Np) ! effet thermique
    pg = rho0*q2*(1 + eye*fu/f) ! densite dynamique complexe
    kg = (gamma*P0)/(1+(gamma-1)/(1-eye*f/ft)) ! module de compressibilitÃ©
    
    Zc = sqrt(pg*kg)/Om ! impedance caracteristique
    kc = omega*sqrt(pg/kg) ! nombre d'onde dans la couche
    
     ! Reflexion etendue avec effet d epaisseur
    rhoc = kc*Zc/omega
    n1 = kc/k0
    m1 = rho0/rhoc
    
    ! Admittance
    beta = -eye*m1*sqrt(n1**2-1)*tan(k0*e*sqrt(n1**2-1))!with hard-backed layer
    end subroutine

end module
