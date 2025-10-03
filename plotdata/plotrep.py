import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import csv

thresh=5000
#Define theshold check
def find_first_passing_index(data_list, threshold):
    for i, value in enumerate(data_list):
        if value == threshold or value > threshold:
            return i
    return -1 # Or raise an error if no value passes
    
#Define the exponential growth function
def exponential_growth(x, B):
    """
    Exponential growth model: y = A * B^x
    """
    return B**x

#Plot ODE solution from the known cycle time
initial_guess = [1.0351]
x=np.linspace(0,266,16)
plt.plot(x, exponential_growth(x, initial_guess), 
         color='Black', label='ODE Growth Curve')

ydata = []         
with open('./TF_Basic/0/TF_data.txt','r') as file:
    for line in file:
        processed_line = line.strip().strip('[]')
        str_list = processed_line.split(', ')
        ydata.append(str_list)
x=np.array(ydata[0]).astype(float)
y=np.array(ydata[1]).astype(float)

index = find_first_passing_index(y, thresh)
allx=x[0:index]
ally=y[0:index]

#plt.plot(x, y, color = 'b', linestyle = 'dashed',
#         marker = 'o',label = "TF")
         
ydata = []         
with open('./TF_Basic/1/TF_data.txt','r') as file:
    for line in file:
        processed_line = line.strip().strip('[]')
        str_list = processed_line.split(', ')
        ydata.append(str_list)
x=np.array(ydata[0]).astype(float)
y=np.array(ydata[1]).astype(float)

index = find_first_passing_index(y, thresh)

allx=np.concatenate((allx,x[0:index]))
ally=np.concatenate((ally,y[0:index]))

#plt.plot(x, y, color = 'b', linestyle = 'dashed',
#         marker = 'o',label = "TF")

ydata = []
with open('./TF_Basic/2/TF_data.txt','r') as file:
    for line in file:
        processed_line = line.strip().strip('[]')
        str_list = processed_line.split(', ')
        ydata.append(str_list)
x=np.array(ydata[0]).astype(float)
y=np.array(ydata[1]).astype(float)

#plt.plot(x, y, color = 'b', linestyle = 'dashed',
#         marker = 'o',label = "TF")

index = find_first_passing_index(y, thresh)

allx=np.concatenate((allx,x[0:index]))
ally=np.concatenate((ally,y[0:index]))

popt, pcov = curve_fit(exponential_growth, allx, ally, p0=initial_guess)

# Extract the fitted parameters
fitted_b = popt

print(f"Fitted parameters for Tissue_Forge: b={fitted_b}")

#Plot the fitted curve
x=np.linspace(0,266,16)
plt.plot(x, exponential_growth(x, fitted_b), 
         color='blue', label='Tissue Forge Exponential Growth Curve')

xdata = []
ydata = []

with open('./Ki67Cycle_1_20/0/numcell.csv','r') as csvfile:
    lines = csv.reader(csvfile, delimiter=',')
    next(lines)
    for row in lines:
        xdata.append(row[0])
        ydata.append(row[1])
x = np.array(xdata).astype(float)/60
y = np.array(ydata).astype(float)

index = find_first_passing_index(y, thresh)

allx=x[0:index]
ally=y[0:index]
xdata = []
ydata = []

with open('./Ki67Cycle_1_20/1/number_cells.dat','r') as csvfile:
    lines = csv.reader(csvfile, delimiter=',')
    next(lines)
    for row in lines:
        xdata.append(row[0])
        ydata.append(row[1])
x = np.array(xdata).astype(float)/60
y = np.array(ydata).astype(float)

index = find_first_passing_index(y, thresh)

allx=np.concatenate((allx,x[0:index]))
ally=np.concatenate((ally,y[0:index]))

xdata = []
ydata = []

with open('./Ki67Cycle_1_20/2/number_cells.dat','r') as csvfile:
    lines = csv.reader(csvfile, delimiter=',')
    next(lines)
    for row in lines:
        xdata.append(row[0])
        ydata.append(row[1])
x = np.array(xdata).astype(float)/60
y = np.array(ydata).astype(float)

index = find_first_passing_index(y, thresh)

allx=np.concatenate((allx,x[0:index]))
ally=np.concatenate((ally,y[0:index]))

xdata = []
ydata = []
         
with open('./Ki67Cycle_1_20/3/number_cells.dat','r') as csvfile:
    lines = csv.reader(csvfile, delimiter=',')
    next(lines)
    for row in lines:
        xdata.append(row[0])
        ydata.append(row[1])
x = np.array(xdata).astype(float)/60
y = np.array(ydata).astype(float)

index = find_first_passing_index(y, thresh)

allx=np.concatenate((allx,x[0:index]))
ally=np.concatenate((ally,y[0:index]))

plt.scatter(allx, ally, color = 'g', linestyle = 'dashed',
        marker = 'o',label = "CC3D_2D_Vol")
         
popt, pcov = curve_fit(exponential_growth, allx, ally, p0=initial_guess)

# Extract the fitted parameters
fitted_b = popt

print(f"Fitted parameters for 1_20: b={fitted_b}")

#Plot the fitted curve
x=np.linspace(0,266,16)
plt.plot(x, exponential_growth(x, fitted_b), 
         color='green', label='CC3D Exponential Growth Curve; volume args')
         
xdata = []
ydata = []

with open('./Ki67Basic3D/0/Ki673D.csv','r') as csvfile:
    lines = csv.reader(csvfile, delimiter=',')
    next(lines)
    for row in lines:
        xdata.append(row[0])
        ydata.append(row[1])
x = np.array(xdata).astype(float)/60*5
y = np.array(ydata).astype(float)

index = find_first_passing_index(y, thresh)

allx=x[0:index]
ally=y[0:index]

#plt.plot(x, y, color = 'k', linestyle = 'dashed',
#         marker = 'o',label = "CC3D_3D_Vol")
         
xdata = []
ydata = []

with open('./Ki67Basic3D/1/number_cells.dat','r') as csvfile:
    lines = csv.reader(csvfile, delimiter=',')
    next(lines)
    for row in lines:
        xdata.append(row[0])
        ydata.append(row[1])
x = np.array(xdata).astype(float)/60*5
y = np.array(ydata).astype(float)

index = find_first_passing_index(y, thresh)

allx=np.concatenate((allx,x[0:index]))
ally=np.concatenate((ally,y[0:index]))

#plt.plot(x, y, color = 'k', linestyle = 'dashed',
#         marker = 'o',label = "CC3D_3D_Vol")

xdata = []
ydata = []
         
with open('./Ki67Basic3D/2/number_cells.dat','r') as csvfile:
    lines = csv.reader(csvfile, delimiter=',')
    next(lines)
    for row in lines:
        xdata.append(row[0])
        ydata.append(row[1])
x = np.array(xdata).astype(float)/60*5
y = np.array(ydata).astype(float)

index = find_first_passing_index(y, thresh)

allx=np.concatenate((allx,x[0:index]))
ally=np.concatenate((ally,y[0:index]))

#plt.plot(x, y, color = 'k', linestyle = 'dashed',
#         marker = 'o',label = "CC3D_3D_Vol")

popt, pcov = curve_fit(exponential_growth, allx, ally, p0=initial_guess)

# Extract the fitted parameters
fitted_b = popt

print(f"Fitted parameters for 3D: b={fitted_b}")

#Plot the fitted curve
x=np.linspace(0,266,16)
plt.plot(x, exponential_growth(x, fitted_b), 
         color='magenta', label='CC3D Exponential Growth Curve; 3D, volume args')
  
xdata = []
ydata = []

with open('./Ki67Basic_1_20/1/number_cells.dat','r') as csvfile:
    lines = csv.reader(csvfile, delimiter=',')
    next(lines)
    for row in lines:
        xdata.append(row[0])
        ydata.append(row[1])
x = np.array(xdata).astype(float)/60
y = np.array(ydata).astype(float)

index = find_first_passing_index(y, thresh)
allx=x[0:index]
ally=y[0:index]
         
xdata = []
ydata = []

with open('./Ki67Basic_1_20/2/number_cells.dat','r') as csvfile:
    lines = csv.reader(csvfile, delimiter=',')
    next(lines)
    for row in lines:
        xdata.append(row[0])
        ydata.append(row[1])
x = np.array(xdata).astype(float)/60
y = np.array(ydata).astype(float)

index = find_first_passing_index(y, thresh)

allx=np.concatenate((allx,x[0:index]))
ally=np.concatenate((ally,y[0:index]))
         
xdata = []
ydata = []

with open('./Ki67Basic_1_20/3/number_cells.dat','r') as csvfile:
    lines = csv.reader(csvfile, delimiter=',')
    next(lines)
    for row in lines:
        xdata.append(row[0])
        ydata.append(row[1])
x = np.array(xdata).astype(float)/60
y = np.array(ydata).astype(float)

index = find_first_passing_index(y, thresh)

allx=np.concatenate((allx,x[0:index]))
ally=np.concatenate((ally,y[0:index]))

#plt.scatter(allx, ally, color = 'r', linestyle = 'dashed',
#         marker = 'o',label = "CC3D_2D")
         
# Fit the curve
# p0 provides initial guesses for the parameters (a, b)
popt, pcov = curve_fit(exponential_growth, allx, ally, p0=initial_guess)

# Extract the fitted parameters
fitted_b = popt

print(f"Fitted parameters for Basic: b={fitted_b}")

#Plot the fitted curve
x=np.linspace(0,266,16)
plt.plot(x, exponential_growth(x, fitted_b), 
         color='red', label='CC3D Exponential Growth Curve; no args')

plt.yscale('log')
plt.xlabel('Time [hr]')
plt.ylabel('Number of Cells')
plt.title('Population', fontsize = 20)
plt.grid()
plt.legend()
plt.show()