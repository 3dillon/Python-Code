# Import proper libraries for math and file reading below
import numpy as np
from scipy import interpolate
import math

# Flow Properties
R = 0.2
Vinf = 5.0
LiftSpan = 0
Rho = 1.18
Zero = 1*10^-6
pi = np.pi

Lambda = 1
Gamma = LiftSpan / (Rho * Vinf)
#Cl = LiftSpan / (0.5 * Rho * (Vinf^2) * R)

Xmin = -1
Xmax = 1
numX = 101
numX2 = 0 # due to nature of loop, theta variable requires iterated sequence to be reverse from X and Y calculations

Ymin = -1
Ymax = 1
numY = 101
numY2 = 0

dx = (Xmax - Xmin) / (numY - 1)
dy = (Ymax - Ymin) / (numY - 1)

X = []
Y = []
Phi = []
Psi = []

for i in range(numX):

    X.append(Xmin + dx * (numX - 1))   # saves calculated Xvals value into list each iteration of loop
    Y.append(Ymin + dy * (numY - 1))   # saves calculated Yvals value into list each iteration of loop

    # r returns the square-root of the highest value in the list formed from summation of Xvals and Yvals lists
    r = np.sqrt(max(np.add(X, Y)))
    # theta is calculated for each iteration based off of the current numX2 and numY2 values of the loop
    theta = (np.arctan2((Ymin + dy * (numY2 - 1)), (Xmin + dx * (numX2 - 1))))

    Phi.append((Vinf * r * np.cos(theta) * (1+(R/r)**2)))

    numX -= 1   # subtracts 1 from numX so that each iteration of loop has proper value to be used in calculations
    numY -= 1   # subtracts 1 from numY so that each iteration of loop has proper value to be used in calculations
    numX2 += 1   # adds 1 from numX2 so that each iteration of loop has proper value to be used in calculations
    numY2 += 1   # adds 1 from numY2 so that each iteration of loop has proper value to be used in calculations

Xarray = np.flip(np.array([X,]*101).transpose(), axis = 0)
Yarray = np.flip(np.array([Y,]*101), axis = 1)
Phiarray = np.array([Phi,]*101)


print(Xarray, "\n")
#print(Yarray)
print("Theta is: ", theta)
print("r is: ", r)

print(Phiarray)
