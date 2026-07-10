#!/usr/bin/python                                                                                                                                                                                                                                                                               
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
 
font = {'family' : 'serif',
        'size'   : 13}
 
plt.rc('font', **font)
 
DATA_LU0 = pd.read_csv("box_reactor_results_lu0.csv")
DATA_LU1 = pd.read_csv("box_reactor_results_lu1.csv")

fig = plt.figure(figsize=(10,6))
plt.plot(DATA_LU0["Time"], DATA_LU0["T_max"], 'o-', label="T_max (LU=0)")
plt.plot(DATA_LU0["Time"], DATA_LU0["T_min"], 'o-', label="T_min (LU=0)")
plt.plot(DATA_LU1["Time"], DATA_LU1["T_max"], 's-', label="T_max (LU=1)")
plt.plot(DATA_LU1["Time"], DATA_LU1["T_min"], 's-', label="T_min (LU=1)")
plt.xlabel('Time (s)')
plt.ylabel('Temperature (K)')
plt.legend()
plt.grid()
plt.savefig("zrk_h2_9sp_box_temp_trend.png")
 
plt.show()
