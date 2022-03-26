# https://research.wdss.io/planetary-motion/#The-Setup
# https://medium.com/swlh/create-your-own-n-body-simulation-with-python-f417234885e9
# https://courses.physics.ucsd.edu/2018/Winter/physics141/Lectures/Lecture2/volker.pdf

from __future__ import annotations  

import math
#from astropy.constants import G as gravitationalConstant
import matplotlib.pyplot as plt
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

    @property
    def xPositionalHistory(self) -> list:
        return [pos.x for pos in self.orbitalHistory]
    
    @property
    def yPositionalHistory(self) -> list:
        return [pos.y for pos in self.orbitalHistory]


class System:
    def __init__(self, bodies:tuple[Body]) -> None:
        self.bodies:tuple[Body] = bodies

    def __iter__(self):
        return (body for body in self.bodies)
    
    @property
    def TotalGpeOfSystem(self):
        # GPE = -GM_1M_2/r
        totalGPE = 0
        for index, body in enumerate(self.bodies):
            for other in self.bodies[index:]: #does not double count bodies in GPE calculation
                if body == other:
                    pass
                else:
                    #print(body.name, other.name)
                    totalGPE += -gravitationalConstant * body.mass * other.mass / math.sqrt(body.CalculateSeperationSquared(other))
        return totalGPE

    @property
    def TotalKeOfSystem(self):
        # KE = mv^2 / 2
        totalKE = 0
        for body in self.bodies:
            totalKE += 0.5 * body.mass * (body.currentVel.Magnitude**2)
        return totalKE

    @property
    def TotalEnergyOfSystem(self):
        return self.TotalGpeOfSystem + self.TotalKeOfSystem


def main():
    astronomicalUnit = 1.5e11
    sunInitialVelocity = math.sqrt(gravitationalConstant * 5.972e24 / 1.5e11)
    earthInitialVelocity = math.sqrt(gravitationalConstant * 1.989e30 / 1.5e11)
    
    Sun = Body(Vector(0,0),Vector(0,-sunInitialVelocity),1.989e30,"Sun")
    Earth = Body(Vector(1.5e11,0),Vector(0,earthInitialVelocity),5.972e24,"Earth")
    
    solarSystem = System((Sun, Earth))

    daysToSeconds = 24*60*60
    runTime = 365 * daysToSeconds
    dt = 0.01 * daysToSeconds

    currentTime = 0

    kineticEnergy = []
    gravitationalPotential = []
    totalEnergy = []
    time = []

    while currentTime <= runTime:
        for body in solarSystem:
            body.UpdateValues(solarSystem, dt)
        kineticEnergy.append(solarSystem.TotalKeOfSystem)
        gravitationalPotential.append(solarSystem.TotalGpeOfSystem)
        totalEnergy.append(solarSystem.TotalEnergyOfSystem)
        currentTime += dt
        time.append(currentTime)
  
    
    plt.plot(Earth.xPositionalHistory, Earth.yPositionalHistory, color='g')
    plt.plot(Sun.xPositionalHistory, Sun.yPositionalHistory, color='y')
    plt.show()

    plt.plot(time, kineticEnergy, color='g')
    plt.plot(time, gravitationalPotential, color='r')
    plt.plot(time, totalEnergy, color='k')
    plt.show()

    
if __name__ == "__main__":
    main()