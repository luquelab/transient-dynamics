import numpy             as np
import pandas            as pd
import matplotlib.pyplot as plt

#Declare the array containing the series you want to plot. 
#For example:
time_series_array = np.sin(np.linspace(-np.pi, np.pi, 400)) + np.random.rand((400))
n_steps           = 15 #number of rolling steps for the mean/std.

#Compute curves of interest:
time_series_df = pd.DataFrame(time_series_array)
smooth_path    = time_series_df.rolling(n_steps).mean()
path_deviation = 2 * time_series_df.rolling(n_steps).std()
print(time_series_df)
print(smooth_path)
under_line     = (smooth_path-path_deviation)[0]
over_line      = (smooth_path+path_deviation)[0]

#Plotting:
plt.plot(smooth_path, linewidth=2) #mean curve.
plt.fill_between(path_deviation.index, under_line, over_line, color='b', alpha=.1)
plt.show()
