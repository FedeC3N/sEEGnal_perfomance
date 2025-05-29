#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Use this script to load all the data and plot different information

Created on Thu 21/05/2025

@author: Fede
"""

# Imports
import numpy
import matplotlib.pyplot as plt

from scripts.miscellaneous.private import channel_information



#### PARAMETERS
# Select the database
database = 'AI_Mind_database'



# Plot deviation of the channels
"""deviation_all = channel_information.estimate_deviation(database)
x = range(deviation_all.shape[1])
plt.plot(x,deviation_all.transpose(),'o')
plt.show(block=True)"""

# Plot the type of badchannels
"""channel_information.badchannel_types(database)"""

# Plot correlation values
"""channel_information.plot_correlation_values(database)"""

# Compare correlation with the whole channels and with epochs
"""channel_information.compare_correlation_values(database)"""

# Print channel distance
"""channel_information.print_channels_distance()"""

# Print clean subjects
channel_information.plot_clean_subject(database)
