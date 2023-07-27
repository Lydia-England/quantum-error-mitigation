######################################
###  functions_binning.py          ###
###  Author:        Lydia England  ###
###  Creation Date: July 27, 2023  ###
###  Last Modified: July 27, 2023  ###
######################################


#==== Define a function to create bins =======================#
def create_bins(lower_bound, width, quantity):
    """ create_bins returns an equal-width (distance) partitioning. 
        It returns an ascending list of tuples, representing the intervals.
        A tuple bins[i], i.e. (bins[i][0], bins[i][1])  with i > 0 
        and i < quantity, satisfies the following conditions:
            (1) bins[i][0] + width == bins[i][1]
            (2) bins[i-1][0] + width == bins[i][0] and
                bins[i-1][1] + width == bins[i][1]
    """
    bins    = []                                              # initialize array to store bins in 
    low_arr = []                                              # initialize array of lowest values in each interval
    count   = lower_bound                                     # initialize counter at given lower bound of bins
    for i in range(quantity):                                 # for loop through number of bins
       low_arr.append(count)                                  # add counter value to array of lowest values per interval
       bins.append((low_arr[i], low_arr[i]+width))            # add new bin to array of bins
       count = count + width                                  # increment counter so the lowest value of the next bin is count+width
    return bins                                               # return array of bins


#==== Define a function to put values in bins ================#
def find_bin(value, bins):
    """ bins is a list of tuples, like [(0,20), (20, 40), (40, 60)],
        binning returns the smallest index i of bins so that
        bin[i][0] <= value < bin[i][1]
    """
    for i in range(0, len(bins)):                             # for bin in range (number of bins)
        if bins[i][0] <= value < bins[i][1]:                  # if the value is within the bounds of this bin
            return i                                          # return bin index
    return -1                                                 # otherwise, return false



