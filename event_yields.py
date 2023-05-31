# Code to calculate expected event yields for the process N+mu -> N' + mu' + X, given cross sections


from tqdm import tqdm
import numpy as np
import scipy.integrate as integrate
import warnings
import os.path


warnings.filterwarnings('ignore')


def calculate_event_yields(config):

    # #########################################
    # ########## STEP 0: Config File ##########
    # #########################################

    # ##### PARAMETERS #####
    E_0 = config["E_0"]
    N_mu = config["N_mu"]
    m_lepton = config["Lepton Mass"]
    theta_max = config["theta_max"]
    x_max = 1 - m_lepton / E_0

    xbins = config["xbins"]
    dx = 1.0/xbins
    xs = np.linspace(0, 1, xbins)

    m_Xs = np.array(config["m_X"])
    m_Xs = m_Xs[m_Xs < E_0]

    # Target Stuff
    target_name = config["Target Name"]
    target_A = config["Target A"]
    target_X0 = 1/(15e-6)  # g/cm2, https://pdg.lbl.gov/2015/AtomicNuclearProperties/MUB/muB_tungsten_W.pdf @ 1-2 TeV, b_tot inverse
    target_density = config["Target Density"]  # g/cm3

    # Experiment Parameters
    l_target = config["Target Length"]  # Target Length (m)
    l_detector = config["Detector Length"]  # Detector Length (m)
    l_shield = config["Shield Length"]  # Shield Length (m)

    # Format is [weight, mass]
    decay_name_string = config["Model String"]
    decay_dictionary = config["Decays"]

    # Options
    run_IWW = config["IWW"]
    run_WW = config["WW"]
    cross_section_file = config["Cross Sections File"]
    experiment_dir = config["Experiment Directory"]

    if run_IWW:
        cross_section_dict = np.load(f"{cross_section_file}.npy", allow_pickle=True)[()]
    if run_WW:
        cross_section_dict = np.load(f"{cross_section_file}_WW.npy", allow_pickle=True)[()]

    # Constants
    m_mu = 0.105
    m_e = 0.000511
    m_tau = 1.776
    m_p = 0.937
    alpha_EM = 1.0/137

    hbarc = .197  # GeV fm
    hbarc_cm_invGeV2 = hbarc * 1e-13
    N_A = 6e23  # Avagadros number (mol^-1)
    maximum_exponent = 250

    l_target = l_target * 100  # Target Length (cm)
    l_detector = l_detector * 100  # Detector Length (cm)
    l_shield = l_shield * 100  # Shield Length (cm)

    cases = config["cases"]
    colors = ["red", "yellow", "green", "blue"]

# ######################################################
# ######### STEP 1: Lifetime and Decay Length ##########
# ######################################################

    # 4 Different decay widths

    def scalar_decay_width(m_X, epsilon, m_l):

        prefactor = epsilon**2 * alpha_EM / 2
        kinematic_factor = m_X * np.nan_to_num(np.sqrt(1 - 4 * m_l**2 / m_X**2)**3)
        return prefactor * kinematic_factor

    def pseudoscalar_decay_width(m_X, epsilon, m_l):

        prefactor = epsilon**2 * alpha_EM / 2
        kinematic_factor = m_X * np.nan_to_num(np.sqrt(1 - 4 * m_l**2 / m_X**2))
        return prefactor * kinematic_factor

    def vector_decay_width(m_X, epsilon, m_l):

        prefactor = epsilon**2 * alpha_EM / 3
        kinematic_factor = m_X * (1 + 2*m_l**2 / m_X**2) * np.nan_to_num(np.sqrt(1 - 4 * m_l**2 / m_X**2))
        return prefactor * kinematic_factor

    def axial_vector_decay_width(m_X, epsilon, m_l):

        prefactor = epsilon**2 * alpha_EM / 3
        kinematic_factor = m_X * np.nan_to_num(np.sqrt(1 - 4 * m_l**2 / m_X**2)**3)
        return prefactor * kinematic_factor

    decay_widths = [scalar_decay_width, pseudoscalar_decay_width, vector_decay_width, axial_vector_decay_width]

    # Helper function to read PDG Data files

    def is_float(string):
        """ True if given string is float else False"""
        try:
            return float(string)
        except ValueError:
            return False

    # R-ratio for hadronic decays
    pdg_r = []
    with open('rpp_cleaned.dat', 'r') as f:
        d = f.readlines()
        for i in d:
            k = i.rstrip().split(" ")
            if "*" not in k:
                row = []
                for x in k:
                    if x != "" and is_float(x):
                        row.append(float(x))
                if len(row) > 6 and is_float(row[0]):
                    pdg_r.append(np.array(row[:7]))

    pdg_r = np.array(pdg_r)

    def Rpp(m_X_linspace):
        return np.interp(m_X_linspace, pdg_r[:, 0], pdg_r[:, 3])

    # Return the decay length of a particle with energy fraction x, mass m_X, coupling epsilon, in meters

    def decay_length(decay_width, x, m_X, epsilon):

        gamma = (x * E_0) / m_X

        total_decay_width = 0
        for decay in decay_dictionary:
            weight, mass = decay_dictionary[decay]
            if decay == "hadrons":
                weight = Rpp(m_X) * weight
            total_decay_width += weight * decay_width(m_X, epsilon, mass)

        return 1e-15 * hbarc * gamma / total_decay_width


# ####################################################
# ########## STEP 2: Calculate Event Yields ##########
# ####################################################

    # Function to interpolate (in log scale) the cross section between two pre-calculated mass points

    def interpolate_cross_section(x, m_X, cross_sections):

        if x > x_max:
            return 0

        # logarithmic interpolation
        return (x > m_X / E_0) * dx * np.interp(m_X, m_Xs, np.nan_to_num([cross_sections[i][int(x * xbins)-1] for i in range(len(cross_sections))]))

    # numpy-vectorized function to approximate log[exp(x)-1] when x is either very large or very small

    def approximate_log_exp_minus_1(x, maximum_exponent, sign=1):

        answer = np.nan_to_num(np.log(sign * (np.exp(x)-1)))

        # if x is large and positive
        large_positive_x_bools = (x > maximum_exponent) * (x > 0)
        answer = x*large_positive_x_bools + answer*(1-large_positive_x_bools)

        # if x is small, but positive
        small_positive_x_bools = (x < 1/maximum_exponent) * (x > 0)
        answer = np.nan_to_num(np.log(x)*small_positive_x_bools) + answer*(1-small_positive_x_bools)

        # if x is large and negative
        large_negative_x_bools = (-1*x > maximum_exponent) * (x < 0)
        answer = -np.exp(x)*large_negative_x_bools + answer*(1-large_negative_x_bools)

        # if x is small and negative
        small_negative_x_bools = (-1*x < 1/maximum_exponent) * (x < 0)
        answer = np.nan_to_num(np.log(-1 * x)*small_negative_x_bools) + answer*(1-small_negative_x_bools)

        return answer

    # Calculate (the log of) dN/dx

    def log_normalized_dNdx(N_mu, x, m_X, epsilon, Gamma, cross_sections):

        if x > x_max:
            return -99999

        length = decay_length(Gamma, x, m_X, epsilon) * 100  # cm

        term1 = np.nan_to_num(np.log(target_density * length / target_A * hbarc_cm_invGeV2**2))
        term2 = np.nan_to_num(np.log(interpolate_cross_section(x, m_X, cross_sections)))

        term3 = np.nan_to_num(approximate_log_exp_minus_1(l_target / length, maximum_exponent, 1))
        term4 = np.nan_to_num(-1 * (l_target + l_shield) / length)

        term5 = np.nan_to_num(approximate_log_exp_minus_1(-1 * l_detector / length, maximum_exponent, -1))

        numerical_check_3 = (length > 0)

        numerical_check = numerical_check_3
        return (numerical_check * (term1 + term2 + term3 + term4 + term5))

    # Calculate the total number of events by integrating (the exponential of) log(dN/dX)

    def compute_events(N_mu, m_X, epsilon, width, case):

        events = np.zeros_like(m_X)
        for (i, x) in tqdm(enumerate(xs)):

            length = decay_length(width, x, m_X, epsilon) * 100  # cm
            experiment_lengthscale = l_detector + l_shield + l_target

            # Check to make sure that the decay length is both positive and not too big
            numerical_check_3 = (length > 0)  # * (length < maximum_exponent * experiment_lengthscale)

            numerical_check = numerical_check_3
            dNdX = log_normalized_dNdx(N_mu, x, m_X, epsilon, width, cross_section_dict[case])
            events += numerical_check * np.exp((np.nan_to_num(dNdX)))
        return N_A * N_mu * events

    # ############################################
    # ########## STEP 3: Run Everything ##########
    # ############################################

    m_X_logspace = np.logspace(-2, 1.2, 250)
    epsilon_logspace = np.logspace(-10, -2, 250)

    m_X_linspace = np.linspace(-2, 1.2, 250)
    epsilon_linspace = np.linspace(-10, -2, 250)

    M, E = np.meshgrid(m_X_logspace, epsilon_logspace)
    Mlin, Elin = np.meshgrid(m_X_linspace, epsilon_linspace)

    print(f"Checking for existing event yield data at {config['Event Yields File']}.npy ...")

    if not os.path.isfile(f"{config['Event Yields File']}.npy"):

        print(f"Could not find existing event yield data for the specified experiment, calculating now ...")

        # Calculate everything
        scalar_events = E**2 * compute_events(N_mu, M, E, scalar_decay_width, cases[0])
        psuedoscalar_events = E**2 * compute_events(N_mu, M, E, pseudoscalar_decay_width, cases[1])
        vector_events = E**2 * compute_events(N_mu, M, E, vector_decay_width, cases[2])
        axial_vector_events = E**2 * compute_events(N_mu, M, E, axial_vector_decay_width, cases[3])

        events_dict = {
            "Scalar": scalar_events,
            "Pseudoscalar": psuedoscalar_events,
            "Vector": vector_events,
            "Axial Vector": axial_vector_events,
        }

        np.save(config["Event Yields File"], events_dict)

    else:
        print(f"Found existing event yield data for the specified experiment, no need to recalculate!")

    print(f"Event yields for {config['Target Name']} with target, shield, detector lengths of {int(l_target/100)}m, {int(l_shield/100)}m, {int(l_detector/100)}m, at {E_0} GeV are available!")
    print()
