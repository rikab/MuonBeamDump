# Code to calculate cross sections for the process N+mu -> N' + mu' + X
# Using the Weizacker-Williams approximation and the "Improved" Weizacker-Williams approximation


from tqdm import tqdm
import numpy as np
import scipy.integrate as integrate
from functools import lru_cache
import warnings
import os.path


warnings.filterwarnings('ignore')


# # Save config file
# os.makedirs(os.path.dirname(config["Config Directory"]), exist_ok=True)
# np.save(config["Config Directory"], config)


def calculate_cross_sections(config, force_rerun=False):

    # #########################################
    # ########## STEP 0: Config File ##########
    # #########################################

    E_0 = config["E_0"]

    # # Allow for optional command line input for the energy
    # args = sys.argv[1:]
    # if len(args) > 0:
    #     E_0 = int(args[0])
    #     config["E_0"] = E_0
    #     if E_0 > 10:
    #         config["WW"] = False
    #     reset_config_directories()

    m_lepton = config["Lepton Mass"]
    theta_max = config["theta_max"]
    x_max = 1 - m_lepton / E_0

    xbins = config["xbins"]
    xs = np.linspace(0, 1, xbins)

    m_Xs = np.array(config["m_X"])

    Z = config["Target Z"]
    A = config["Target A"]

    # Options
    run_IWW = config["IWW"]
    run_WW = config["WW"]
    cross_section_directory = config["Cross Sections File"]

    cases = config["cases"]

    # Constants
    m_mu = 0.105
    m_e = 0.000511
    m_tau = 1.776
    m_p = 0.937
    alpha_EM = 1.0/137

    # #################################################
    # ######### STEP 1: Nuclear Form Factors ##########
    # #################################################

    a = 111.0 * np.power(Z, -1/3) / m_e
    d = 0.164 * np.power(A, -2/3)

    a_prime = 773.0 * np.power(Z, -2/3) / m_e
    mu_p = 2.79

    def g2_elastic(t):

        # Electron screening term
        term1 = np.power(a**2 * t / (1 + a**2 * t), 2)

        # Nuclear size term
        term2 = np.power(1.0 / (1 + t/d), 2)

        return term1 * term2 * Z**2

    def g2_inelastic(t):

        # Inelastic atomic form factor
        term1 = np.power(a_prime**2 * t / (1 + a_prime**2 * t), 2)

        # Inelastic nuclear form factor
        term2 = np.power((1 + t/(4*m_p**2) * (mu_p**2 - 1)) / (1 + t / 0.71)**4, 2)

        return term1 * term2 * Z

    # Total Form Factor

    def g2(t):
        return g2_elastic(t) + g2_inelastic(t)

    # Total photon flux

    def chi(t_min, t_max):

        def integrand(t):
            return (t - t_min) * g2(t) / t**2

        return integrate.quad(integrand, t_min, t_max, epsabs=1.5e-6, epsrel=1.5e-6, limit=600)[0]

    # ###############################################
    # ######### STEP 2: Squared Amplitudes ##########
    # ###############################################

    def beta(E, m):
        return np.nan_to_num(np.sqrt(1 - m**2 / E**2))

    # From [1609.06781]

    def A22_tmin_scalar(u, x, m_X):

        # Not-u dependent
        term1 = x**2 / (1-x)

        # u-dependent
        term2 = 2*(m_X**2 - 4*m_lepton**2) / u**2
        term3 = u*x + m_X**2 * (1-x) + m_lepton**2 * x**2

        return term1 + term2*term3

    # From [1705.01633]

    def A22_tmin_psuedoscalar(u, x, m_X):

        # Not-u dependent
        term1 = x**2 / (1-x)

        # u-dependent
        term2 = 2*(m_X**2) / u**2
        term3 = u*x + m_X**2 * (1-x) + m_lepton**2 * x**2

        return term1 + term2*term3

    def A22_tmin_vector(u, x, m_X):

        # Not-u dependent
        term1 = 2 * (2 - 2*x + x**2) / (1-x)

        # u-dependent
        term2 = 4*(m_X**2 + 2 * m_lepton**2) / u**2
        term3 = u*x + m_X**2 * (1-x) + m_lepton**2 * x**2

        return term1 + term2*term3

    def A22_tmin_Axialvector(u, x, m_X):

        # Not-u dependent
        term0 = 4 * m_lepton**2 * x ** 2 / (1-x) / m_X**2  # Possible typo in paper? m without a subscript.
        term1 = 2 * (2 - 2*x + x**2) / (1-x)

        # u-dependent
        term2 = 4*(m_X**2 - 4 * m_lepton**2) / u**2
        term3 = u*x + m_X**2 * (1-x) + m_lepton**2 * x**2

        return term0 + term1 + term2*term3

    def A22_tmin_primakoff(u, x, m_X):

        # Not-u dependent
        term0 = 16 * u**2 * ((x-2)*x + 2)
        term1 = (x-1) * (u * x - m_X**2 * (x-1))

        return term0 / term1

    A22_tmins = [A22_tmin_scalar, A22_tmin_psuedoscalar, A22_tmin_vector, A22_tmin_Axialvector, A22_tmin_primakoff]

    # ##################################################
    # ######### STEP 3: Phase Space Integrals ##########
    # ##################################################

    # IWW: Integrate squared amplitude over u

    @lru_cache(maxsize=1000)
    def integrated_A22_IWW(A, x, m_X):

        u_max = -m_X**2 * (1-x) / x - m_lepton**2 * x

        def integrand(u):
            return A(u, x, m_X) / u**2

        return integrate.quad(integrand, -np.infty, u_max, epsabs=1.5e-9, epsrel=1.5e-9, limit=500)[0]

    # IWW Combine integrated squared ampitude with prefactors and photon flux

    @lru_cache(maxsize=1000)
    def dsigma_dx_IWW(A, x, m_X, E_0):

        if x > x_max or x < m_X / E_0:
            return 0

        # Pretend epsilon is 1. Multiply by epsilon**2 if necessary
        t_min = (m_X**2 / (2 * E_0))**2
        t_max = m_X**2 + m_lepton**2
        chi_factor = (chi(t_min, t_max))

        # multiplying beta by x to get |k|/E_0 factor
        flux = 1**2 * alpha_EM**3 * chi_factor * beta(x * E_0, m_X) * x

        u_max = -m_X**2 * (1-x) / x - m_lepton**2 * x
        return flux * (1-x) / x * (integrated_A22_IWW(A, x, m_X))

    # Full WW integral (squared amplitude and photon flux convoluted with theta)

    @lru_cache(maxsize=1000)
    def dsigma_dx_WW(A, x, m_X, E_0):

        if x > x_max or x < m_X / E_0:
            return 0

        # Pretend epsilon is 1. Multiply by epsilon**2 if necessary
        def integrand(theta):

            u = -m_X**2 * (1-x) / x - (m_lepton**2 + theta**2 * E_0**2) * x
            t_min = (u / (2 * E_0) / (1-x))**2
            t_max = m_X**2 + m_lepton**2
            chi_factor = (chi(t_min, t_max))

            # sin(theta) due to converting dcostheta to dtheta sintheta
            return np.sin(theta) * A(u, x, m_X) / u**2 * chi_factor

        # multiplying beta by x to get |k|/E_0 factor
        flux = 2 * 1**2 * alpha_EM**3 * beta(x * E_0, m_X) * x * E_0**2 * (1-x)
        theta_0 = min(theta_max, 10 * np.sqrt(m_lepton**2 + (1-x)*m_X**2 / x**2) / E_0)
        return flux * integrate.quad(integrand, 0, theta_0, epsabs=1.5e-9, epsrel=1.5e-9, limit=500)[0]

    # ############################################
    # ########## STEP 4: Run everything ##########
    # ############################################

    print(f"Checking for existing IWW cross section data at {config['Cross Sections File']}.npy ...")

    # IWW Version
    if run_IWW:

        # Open previous dictionary file, if it exists.
        if os.path.isfile(f"{config['Cross Sections File']}.npy") and not force_rerun:
            cross_section_dict = np.load(f"{config['Cross Sections File']}.npy", allow_pickle=True)[()]
            print(f"Found existing IWW cross section data for {config['Target Name']} at {E_0} GeV, no need to recalculate!")
        else:
            cross_section_dict = {}
            print(f"Could not find existing IWW cross section data for {config['Target Name']} at {E_0} GeV, calculating now ...")
            for (i, case) in enumerate(cases):

                cross_sections = []

                for m_X in m_Xs[m_Xs < E_0]:

                    print(f"IWW, {E_0} GeV, {case}, m_X = {m_X} GeV", )
                    cross_section = []
                    for x in tqdm(xs):
                        cross_section.append(dsigma_dx_IWW(A22_tmins[i], x, m_X, E_0))
                    cross_section = np.array(cross_section)
                    cross_sections.append(cross_section)

                cross_sections = np.array(cross_sections)
                cross_section_dict[case] = cross_sections

            np.save(cross_section_directory, cross_section_dict, )

    # WW Version
    if run_WW:

        print("Running with WW ...")
        print(f"Checking for existing WW cross section data at {config['Cross Sections File']}_WW.npy ...")

        # Open previous dictionary file, if it exists.
        if os.path.isfile(f"{config['Cross Sections File']}_WW.npy") and not force_rerun:
            cross_section_dict = np.load(f"{config['Cross Sections File']}_WW.npy", allow_pickle=True)[()]
            print(f"Found existing WW cross section data for {config['Target Name']} at {E_0} GeV, no need to recalculate!")

        else:

            print(f"Could not find existing WW cross section data for {config['Target Name']} at {E_0} GeV, calculating now ...")
            cross_section_dict = {}

            def cross_section_x(A, m_X, E_0):
                cross_section = np.zeros_like(xs)
                for i in tqdm(range(len(xs))):
                    cross_section[i] = (dsigma_dx_WW(A, xs[i], m_X, E_0))
                return cross_section

            for (i, case) in enumerate(cases):

                cross_sections = []
                for m_X in m_Xs[m_Xs < E_0]:
                    print(f"WW, {E_0} GeV, {case}, m_X = {m_X} GeV", )
                    cross_sections.append(cross_section_x(A22_tmins[i], m_X, E_0))

                cross_sections = np.array(cross_sections)
                cross_section_dict[case] = cross_sections

            np.save(cross_section_directory + "_WW", cross_section_dict, )

    print(f"Cross sections for {config['Target Name']} at {E_0} GeV are available!")
    print()
