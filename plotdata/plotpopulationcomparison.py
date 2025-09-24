import matplotlib.pyplot as plt
import numpy as np
import csv

ydata = []         
with open('TF_data.txt','r') as file:
    for line in file:
        processed_line = line.strip().strip('[]')
        str_list = processed_line.split(', ')
        ydata.append(str_list)
x=np.linspace(0,266,16)
y=np.array(ydata[0]).astype(float)

plt.plot(x, y, color = 'b', linestyle = 'dashed',
         marker = 'o',label = "TF")
thresh=y[10]
timestart=x[10] #taking later point when number of cells is higher, we will then try to normalize time around this point

def find_first_passing_index(data_list, threshold):
    for i, value in enumerate(data_list):
        if value == threshold or value > threshold:
            return i
    return -1 # Or raise an error if no value passes

xdata = []
ydata = []

with open('./Ki67Cycle_1_20/numcell.csv','r') as csvfile:
    lines = csv.reader(csvfile, delimiter=',')
    next(lines)
    for row in lines:
        xdata.append(row[0])
        ydata.append(row[1])
x = np.array(xdata).astype(float)/60
y = np.array(ydata).astype(float)

index = find_first_passing_index(y, thresh)
passtime=x[index]
adjust=passtime-timestart
x=x-adjust

plt.plot(x, y, color = 'g', linestyle = 'dashed',
         marker = 'o',label = "CC3D_2D_Vol")
         
xdata = []
ydata = []

with open('./Ki673D.csv','r') as csvfile:
    lines = csv.reader(csvfile, delimiter=',')
    next(lines)
    for row in lines:
        xdata.append(row[0])
        ydata.append(row[1])
x = np.array(xdata).astype(float)/60*5
y = np.array(ydata).astype(float)

index = find_first_passing_index(y, thresh)
passtime=x[index]
adjust=passtime-timestart
x=x-adjust

plt.plot(x, y, color = 'k', linestyle = 'dashed',
         marker = 'o',label = "CC3D_3D_Vol")
         
xdata = []
ydata = []

with open('./Ki67basic2.csv','r') as csvfile:
    lines = csv.reader(csvfile, delimiter=',')
    next(lines)
    for row in lines:
        xdata.append(row[0])
        ydata.append(row[1])
x = np.array(xdata).astype(float)/60*5
y = np.array(ydata).astype(float)

index = find_first_passing_index(y, thresh)
passtime=x[index]
adjust=passtime-timestart
x=x-adjust

plt.plot(x, y, color = 'r', linestyle = 'dashed',
         marker = 'o',label = "CC3D_2D")

plt.xlabel('Time [hr]')
plt.ylabel('Number of Cells')
plt.title('Population', fontsize = 20)
plt.grid()
plt.legend()
plt.show()