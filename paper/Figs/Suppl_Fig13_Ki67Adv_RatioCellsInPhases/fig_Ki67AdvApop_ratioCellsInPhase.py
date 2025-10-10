# plot results

import pandas as pd
from matplotlib import pyplot as plt

plt.rcParams["figure.figsize"] = [7.00, 3.50]
plt.rcParams["figure.autolayout"] = True
columns = ["hr", "Ki67neg", "Ki67pos_pre", "Ki67pos_post", "Apoptosis" ]
df = pd.read_csv("ratio_cells_in_phase.csv", usecols=columns)
columns_physicell = ["x", "y"]
df_physicell_pos_pre = pd.read_csv("Fig11_a_ODE.csv", usecols=columns_physicell)
df_physicell_pos_post = pd.read_csv("Fig11_b_ODE.csv", usecols=columns_physicell)
df_physicell_apop = pd.read_csv("Fig11_c_ODE.csv", usecols=columns_physicell)

plt.plot(df.hr, df.Ki67neg, label = "Ki67-")
plt.plot(df.hr, df.Ki67pos_pre, label = "Ki67+ pre")
plt.plot(df.hr, df.Ki67pos_post, label = "Ki67+ post")
plt.plot(df.hr, df.Apoptosis, ls = "--", linewidth = "3", label = "Apoptosis")
plt.plot(df_physicell_pos_pre.x, df_physicell_pos_pre.y, label = "Analytical soln: Ki67+ pre")
plt.plot(df_physicell_pos_post.x, df_physicell_pos_post.y, label = "Analytical soln: Ki67+ post")
plt.plot(df_physicell_apop.x, df_physicell_apop.y, color = "grey", ls = "dotted", label = "Analytical soln: Apoptosis")
plt.xlabel('hours')
plt.xlim(0, 150)
plt.ylim(0,1.0)
plt.ylabel("Ratio of cells in each phase to total cells")
plt.title("Ki67 cell cycle with Apoptosis of 0.00105 1/hr")
plt.grid()
plt.legend(loc= "center right")
plt.savefig("figure_ratio_cells_Apop_0.00105.pdf")
plt.show()
