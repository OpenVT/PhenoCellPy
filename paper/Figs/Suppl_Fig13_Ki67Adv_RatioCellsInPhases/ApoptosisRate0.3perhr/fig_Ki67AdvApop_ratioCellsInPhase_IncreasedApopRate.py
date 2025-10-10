# plot results

import pandas as pd
from matplotlib import pyplot as plt

plt.rcParams["figure.figsize"] = [7.00, 3.50]
plt.rcParams["figure.autolayout"] = True
columns = ["hr", "Ki67neg", "Ki67pos_pre", "Ki67pos_post", "Apoptosis" ]

df = pd.read_csv("ratio_cells_in_phase_005.csv", usecols=columns)

plt.plot(df.hr, df.Ki67neg, label = "Ki67-")
plt.plot(df.hr, df.Ki67pos_pre, label = "Ki67+ pre")
plt.plot(df.hr, df.Ki67pos_post, label = "Ki67+ post")
plt.plot(df.hr, df.Apoptosis, linewidth = "2", ls = '--', label = "Apoptosis")
plt.xlabel('hours')
plt.xlim(0, 150)
plt.ylim(0,1.0)
plt.ylabel("Ratio of cells in each phase to total cells")
plt.title("Ki67 cell cycle with Apoptosis of 0.3 1/hr")
plt.grid()
plt.legend()
plt.savefig("figure_ratio_cells_Apop_0.3.pdf")
plt.show()
