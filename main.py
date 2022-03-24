# https://research.wdss.io/planetary-motion/#The-Setup
# https://medium.com/swlh/create-your-own-n-body-simulation-with-python-f417234885e9
# https://courses.physics.ucsd.edu/2018/Winter/physics141/Lectures/Lecture2/volker.pdf

from __future__ import annotations  

import math
from re import A
#from astropy.constants import G as gravitationalConstant
import matplotlib.pyplot as plt
import numpy as np
import Vector
gravitationalConstant = 6.67e-11

class ExoObject:
    def __init__(self, xPos:float, yPos:float, ux:float, uy:float, mass:int, name:str) -> None:
        self.currentPosX:float = xPos
        self.currentPosY:float = yPos
        self.currentVelX:float = ux
        self.currentVelY:float = uy
        #self.currentAccX:float = ax
        #self.currentAccY:float = ay

        self.mass:int = mass
        self.name:str = name

        self.xHistory:list = [xPos]
        self.yHistory:list = [yPos]

    def CalcualteTheta(self, other):
        deltaX = other.currentPosX - self.currentPosX
        deltaY = other.currentPosY - self.currentPosY

        theta = math.atan2(deltaY, deltaX)
        return theta
    
    def CalculateSeperationSquared(self, other):
        deltaX = other.currentPosX - self.currentPosX
        deltaY = other.currentPosY - self.currentPosY

        return (deltaX ** 2) + (deltaY ** 2)

    def CalculateAcceleration(self, other):
        theta = self.CalcualteTheta(other)
        seperationSquared = self.CalculateSeperationSquared(other)
        return (gravitationalConstant * other.mass * math.cos(theta) / seperationSquared,
                gravitationalConstant * other.mass * math.sin(theta) / seperationSquared)

    def CalculateTotalAcceleration(self, bodies:list):
        totalAccelerationX = 0
        totalAccelerationY = 0
        count=0
        for body in bodies:
            if body != self:
                ax, ay = self.CalculateAcceleration(body)
                totalAccelerationX += ax
                totalAccelerationY += ay
                count+=1
        return {'x':totalAccelerationX, 'y':totalAccelerationY}

    def CalculateNewVelocity(self, totalAcceleration:dict, timeStep:float):
        self.currentVelX += (timeStep * totalAcceleration['x'])
        self.currentVelY += (timeStep * totalAcceleration['y'])
    
    def CalculateNewPosition(self, timeStep:float):
        self.currentPosX += self.currentVelX * timeStep
        self.currentPosY += self.currentVelY * timeStep

    def UpdateValues(self, bodies:list, timeStep:float):
        accelerationDict = self.CalculateTotalAcceleration(bodies)
        self.CalculateNewVelocity(accelerationDict, timeStep)
        self.CalculateNewPosition(timeStep)
        self.xHistory.append(self.currentPosX)
        self.yHistory.append(self.currentPosY)


def main():
    sunInitialVelocity = math.sqrt(gravitationalConstant * 5.972e24 / 1.5e11)
    earthInitialVelocity = math.sqrt(gravitationalConstant * 1.989e30 / 1.5e11)
    
    Sun = ExoObject(0,0,0,-sunInitialVelocity,1.989e30,"Sun")
    Earth = ExoObject(1.5e11,0,0,earthInitialVelocity,5.972e24,"Earth")

    solarSystem = [Sun, Earth]

    daysToSeconds = 24*60*60
    runTime = 365 * daysToSeconds
    dt = 0.01 * daysToSeconds

    currentTime = 0

    while currentTime <= runTime:
        for body in solarSystem:
            body.UpdateValues(solarSystem, dt)
        currentTime += dt
  
    
    plt.plot(Earth.xHistory, Earth.yHistory, color='g')
    plt.plot(Sun.xHistory, Sun.yHistory, color='y')
    plt.show()

    #print(Earth.xHistory)

def main2():
    a = Vector.Vector(1,2)

    
if __name__ == "__main__":
    main2()