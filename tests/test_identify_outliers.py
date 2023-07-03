import numpy as np

def identify_outliers(array, threshold=0.8):
    """
    Identifies outliers in an array of floats using a majority consensus style veto.

    Args:
        array (numpy.ndarray or list): The input array of floats.
        threshold (float, optional): The threshold for outlier identification. Defaults to 0.8.

    Returns:
        numpy.ndarray: A boolean mask indicating the outliers in the input array.
    """
    array = np.array(array)
    mean = np.mean(array)
    std = np.std(array)
    z_scores = (array - mean) / std

    consensus_threshold = int(threshold * len(array))
    outlier_mask = np.abs(z_scores) > 1

    while np.sum(outlier_mask) > consensus_threshold:
        array = array[~outlier_mask]
        mean = np.mean(array)
        std = np.std(array)
        z_scores = (array - mean) / std
        outlier_mask = np.abs(z_scores) > 1

    final_mask = np.zeros_like(z_scores, dtype=bool)
    final_mask[np.abs(z_scores) > 1] = True

    return final_mask


print(identify_outliers(np.array([1, 1, 2])))

print(identify_outliers(np.array([1, 1])))

print(identify_outliers(np.array([1, 1, 1, 1, 1])))

print(identify_outliers(np.array([1, 1, 1, 1, 1.01])))
