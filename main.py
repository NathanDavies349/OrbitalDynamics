# https://research.wdss.io/planetary-motion/#The-Setup
# https://medium.com/swlh/create-your-own-n-body-simulation-with-python-f417234885e9
# https://courses.physics.ucsd.edu/2018/Winter/physics141/Lectures/Lecture2/volker.pdf

from __future__ import annotations  

import math
from re import A
#from astropy.constants import G as gravitationalConstant
import matplotlib.pyplot as plt
import numpy as np
from Vector import Vector
gravitationalConstant = 6.67e-11

class Body:
    def __init__(self, initialPosition:Vector, initialVelocity:Vector, mass:int, name:str) -> None:
        self.currentPos:Vector = initialPosition
        self.currentVel:Vector = initialVelocity
        self.mass:int = mass
        self.name:str = name

        self.orbitalHistory:list[Vector] = [initialPosition]

    def CalcualteTheta(self, other:Body):
        deltaX = other.currentPos.x - self.currentPos.x
        deltaY = other.currentPos.y - self.currentPos.y

        theta = math.atan2(deltaY, deltaX)
        return theta
    
    def CalculateSeperationSquared(self, other:Body):
        deltaX = other.currentPos.x - self.currentPos.x
        deltaY = other.currentPos.y - self.currentPos.y
        return (deltaX ** 2) + (deltaY ** 2)

    def CalculateAcceleration(self, other:Body):
        theta = self.CalcualteTheta(other)
        seperationSquared = self.CalculateSeperationSquared(other)
        return (gravitationalConstant * other.mass * math.cos(theta) / seperationSquared,
                gravitationalConstant * other.mass * math.sin(theta) / seperationSquared)

    def CalculateTotalAcceleration(self, bodies:list[Body]):
        totalAccelerationX = 0
        totalAccelerationY = 0
        for body in bodies:
            if body != self:
                ax, ay = self.CalculateAcceleration(body)
                totalAccelerationX += ax
                totalAccelerationY += ay
  
        return Vector(totalAccelerationX, totalAccelerationY)

    def CalculateNewVelocity(self, totalAcceleration:Vector, timeStep:float):
        self.currentVel += (totalAcceleration * timeStep)

    
    def CalculateNewPosition(self, timeStep:float):
        self.currentPos += (self.currentVel * timeStep)

    def UpdateValues(self, bodies:list[Body], timeStep:float):
        acceleration = self.CalculateTotalAcceleration(bodies)
        self.CalculateNewVelocity(acceleration, timeStep)
        self.CalculateNewPosition(timeStep)
        self.orbitalHistory.append(self.currentPos)


def main():
    sunInitialVelocity = math.sqrt(gravitationalConstant * 5.972e24 / 1.5e11)
    earthInitialVelocity = math.sqrt(gravitationalConstant * 1.989e30 / 1.5e11)
    
    Sun = Body(Vector(0,0),Vector(0,-sunInitialVelocity),1.989e30,"Sun")
    Earth = Body(Vector(1.5e11,0),Vector(0,earthInitialVelocity),5.972e24,"Earth")

    solarSystem = [Sun, Earth]

    daysToSeconds = 24*60*60
    runTime = 365 * daysToSeconds
    dt = 0.01 * daysToSeconds

    currentTime = 0

    while currentTime <= runTime:
        for body in solarSystem:
            body.UpdateValues(solarSystem, dt)
        currentTime += dt
  
    
    plt.plot([pos.x for pos in Earth.orbitalHistory], [pos.y for pos in Earth.orbitalHistory], color='g')
    plt.plot([pos.x for pos in Sun.orbitalHistory], [pos.y for pos in Sun.orbitalHistory], color='y')
    plt.show()

    #print(Earth.xHistory)

    
if __name__ == "__main__":
    main()