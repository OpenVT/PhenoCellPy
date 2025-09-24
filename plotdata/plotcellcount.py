import matplotlib.pyplot as plt
import numpy as np
import csv

xdata = []
ydata = []

with open('./Ki67Cycle_1_20/numcell.csv','r') as csvfile:
    lines = csv.reader(csvfile, delimiter=',')
    next(lines)
    for row in lines:
        xdata.append(row[0])
        ydata.append(row[1])
x = np.array(xdata).astype(float)
y = np.array(ydata).astype(float)

plt.plot(x, y, color = 'g', linestyle = 'dashed',
         marker = 'o',label = "1_20")

xdata = []
ydata = []         
with open('./Ki67Cycle_1_2/numcell.csv','r') as csvfile:
    lines = csv.reader(csvfile, delimiter=',')
    next(lines)
    for row in lines:
        xdata.append(row[0])
        ydata.append(row[1])
x = np.array(xdata).astype(float)
y = np.array(ydata).astype(float)

plt.plot(x, y, color = 'r', linestyle = 'dashed',
         marker = 'o',label = "1_2")
         
xdata = []
ydata = []         
with open('./Ki67Cycle_5_20/numcell.csv','r') as csvfile:
    lines = csv.reader(csvfile, delimiter=',')
    next(lines)
    for row in lines:
        xdata.append(row[0])
        ydata.append(row[1])
x = np.array(xdata).astype(float)*5
y = np.array(ydata).astype(float)

plt.plot(x, y, color = 'k', linestyle = 'dashed',
         marker = 'o',label = "5_20")
         
xdata = []
ydata = []         
with open('./Ki67Cycle_5_2/numcell.csv','r') as csvfile:
    lines = csv.reader(csvfile, delimiter=',')
    next(lines)
    for row in lines:
        xdata.append(row[0])
        ydata.append(row[1])
x = np.array(xdata).astype(float)*5
y = np.array(ydata).astype(float)

plt.plot(x, y, color = 'b', linestyle = 'dashed',
         marker = 'o',label = "5_2")

plt.xlabel('Time [min]')
plt.ylabel('Number of Cells')
plt.title('numcell', fontsize = 20)
plt.grid()
plt.legend()
plt.show()