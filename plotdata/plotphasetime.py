import matplotlib.pyplot as plt
import numpy as np
import csv

xdata = []
ydata1 = []  
ydata2 = [] 
ydata3 = []        
with open('./Ki67Cycle_1_20/Timespent.csv','r') as csvfile:
    lines = csv.reader(csvfile, delimiter=',')
    next(lines)
    for row in lines:
        xdata.append(row[0])
        ydata1.append(row[1])
        ydata2.append(row[3])
        ydata3.append(row[5])
new_x = [np.nan if entry == '' or entry == [] else entry for entry in xdata]
new_y = [np.nan if entry == '' or entry == [] else entry for entry in ydata1]
x = np.array(new_x).astype(float)
y1 = np.array(new_y).astype(float)
y2 = np.array(ydata2).astype(float)
y3 = np.array(ydata3).astype(float)

plt.plot(x, y1, color = 'k', linestyle = 'solid',label = "min volume")
plt.plot(x, y2, color = 'r', linestyle = 'solid',label = "median volume")
plt.plot(x, y3, color = 'b', linestyle = 'solid',label = "max volume")
plt.xlabel('Time [min]')
plt.ylabel('Timespent')
plt.title('Timespent', fontsize = 20)
plt.grid()
plt.legend()
plt.show()