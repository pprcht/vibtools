import numpy as np


# Implement metrics SID and SIS from https://doi.org/10.1021/acs.jcim.1c00055
def SID(y_pred, y_target):
    # Ensure inputs are numpy arrays
    y_pred = np.asarray(y_pred)
    y_target = np.asarray(y_target)

    # Avoid division by zero and log of zero by adding a small epsilon
    epsilon = 1e-12
    y_pred = np.clip(y_pred, epsilon, None)
    y_target = np.clip(y_target, epsilon, None)

    # Calculate SID
    term1 = y_pred * np.log(y_pred / y_target)
    term2 = y_target * np.log(y_target / y_pred)
    
    return np.sum(term1 + term2)


def SIS(y_pred, y_target):
    # Ensure inputs are numpy arrays
    y_pred = np.asarray(y_pred)
    y_target = np.asarray(y_target)

    sid_calc = SID(y_pred, y_target)

    return 1.0/(1.0 + sid_calc)
    
