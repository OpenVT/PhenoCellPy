import matplotlib.pyplot as plt
import numpy as np
import csv

xdata = []
ydata1 = []  
ydata2 = [] 
ydata3 = []        
with open('./Ki67Cycle_1_20/VolumeMetrics.csv','r') as csvfile:
    lines = csv.reader(csvfile, delimiter=',')
    next(lines)
    for row in lines:
        xdata.append(row[0])
        ydata1.append(row[1])
        ydata2.append(row[3])
        ydata3.append(row[5])
x = np.array(xdata).astype(float)
rx = x.reshape(-1, 10).mean(axis=1) #reducing number of points by 10 taking the mean value
y1 = np.array(ydata1).astype(float)
ry1 = y1.reshape(-1, 10).mean(axis=1)
y2 = np.array(ydata2).astype(float)
ry2 = y2.reshape(-1, 10).mean(axis=1)
y3 = np.array(ydata3).astype(float)
ry3 = y3.reshape(-1, 10).mean(axis=1)

plt.plot(rx/60, ry1, color = 'k', linestyle = 'solid',label = "CC3D median volume")
plt.plot(rx/60, ry2, color = 'r', linestyle = 'solid',label = "CC3D min volume")
plt.plot(rx/60, ry3, color = 'b', linestyle = 'solid',label = "CC3D max volume")

ydata=[]
with open('TF_data.txt','r') as file:
    for line in file:
        processed_line = line.strip().strip('[]')
        str_list = processed_line.split(', ')
        ydata.append(str_list)
x=np.linspace(0,266,16)
y1=np.array(ydata[1]).astype(float)
y2=np.array(ydata[2]).astype(float)
y3=np.array(ydata[3]).astype(float)
y1s=y1/80*200#we scale the value to be normalized to the Phenocellpy doubling values. We know this is double our IC and therefore 80 and 200.
y2s=y2/80*200
y3s=y3/80*200
plt.plot(x, y1s, color = 'k', linestyle = 'dashed',
         marker = 'o',label = "TF")
plt.plot(x, y2s, color = 'r', linestyle = 'dashed',
         marker = 'o',label = "TF")
plt.plot(x, y3s, color = 'b', linestyle = 'dashed',
         marker = 'o',label = "TF")

plt.xlabel('Time [min]')
plt.ylabel('Volmetrics')
plt.title('Volmetrics', fontsize = 20)
plt.grid()
plt.legend()
plt.show()