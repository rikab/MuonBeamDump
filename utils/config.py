import numpy as np
import os.path

# Experiment config file


def build_full_config(run_WW, xbins, cases, m_Xs, cross_section_dir, experiment_dir, E_0, m_lepton, theta_max, target_name, target_Z, target_A, target_density, target_length, shield_length, detector_length, model_name_string, decay_dictionary, N_mu, mrange, erange,):

    config = {

        #   Parameters needed for cross section calculations

        "E_0": E_0,                      # Lepton Beam Energy in GeV
        "Lepton Mass": m_lepton,             # Lepton Mass in GeV
        "theta_max": theta_max,                # Maximum emission angle
        "xbins": xbins,                    # Number of x bins (0 <  x < 1)
        "Target Name": target_name,           # Name of Target
        "Target Z": target_Z,                   # Target Charge / Proton number
        "Target A": target_A,                # Target Mass (in g/mol OR GeV/nuclei)
        "IWW": True,                      # Compute cross sections using IWW
        "WW": run_WW,                       # Compute cross sections using WW. Automatically set to false if E_0 > 10
        "m_X": m_Xs,    # Range of X masses to compute
        "cases": cases,  # DO NOT CHANGE THIS

        # Parameters needed only for computing exclusion limits

        "Target Density": target_density,           # Target density in g/cm3
        "Target Length": target_length,          # Length of target in meters
        "Shield Length": shield_length,          # Length of shield in meters
        "Detector Length": detector_length,       # Length of detector in meters

        "Model String": model_name_string,
        "Decays": decay_dictionary,

        "N_mu": N_mu,
        "m_range": mrange,
        "e_range": erange,
    }

    # Directories
    os.makedirs(cross_section_dir, exist_ok=True)
    config["Cross Sections Directory"] = cross_section_dir
    config["Cross Sections File"] = cross_section_dir + f"cross_sections_{E_0}"
    config["Cross Sections Config"] = cross_section_dir + f"config_{E_0}"

    os.makedirs(experiment_dir, exist_ok=True)
    config["Experiment Directory"] = experiment_dir
    config["Event Yields File"] = experiment_dir + "event_counts"
    config["Event Yields Config"] = experiment_dir + "config"

    np.save(config["Event Yields Config"], config)
    return config


def build_cross_section_config(run_WW, xbins, cases, m_Xs, cross_section_dir, E_0, m_lepton, theta_max, target_name, target_Z, target_A):

    config = {

        #   Parameters needed for cross section calculations

        "E_0": E_0,                      # Lepton Beam Energy in GeV
        "Lepton Mass": m_lepton,             # Lepton Mass in GeV
        "theta_max": theta_max,                # Maximum emission angle
        "xbins": xbins,                    # Number of x bins (0 <  x < 1)
        "Target Name": target_name,           # Name of Target
        "Target Z": target_Z,                   # Target Charge / Proton number
        "Target A": target_A,                # Target Mass (in g/mol OR GeV/nuclei)
        "IWW": True,                      # Compute cross sections using IWW
        "WW": run_WW,                       # Compute cross sections using WW. Automatically set to false if E_0 > 10
        "m_X": m_Xs,    # Range of X masses to compute
        "cases": cases,  # DO NOT CHANGE THIS

    }

    # Directories
    os.makedirs(cross_section_dir, exist_ok=True)
    config["Cross Sections Directory"] = cross_section_dir
    config["Cross Sections File"] = cross_section_dir + f"cross_sections_{E_0}"
    config["Cross Sections Config"] = cross_section_dir + f"config_{E_0}"

    np.save(config["Cross Sections Config"], config)
    return config
