# https://research.wdss.io/planetary-motion/#The-Setup
# https://medium.com/swlh/create-your-own-n-body-simulation-with-python-f417234885e9
# https://courses.physics.ucsd.edu/2018/Winter/physics141/Lectures/Lecture2/volker.pdf
# https://nctstca.github.io/events/202107-tcassp/lectures/NCTS_TCA_SSP_20210706_Numerical_Simulations_hyschive.pdf

# https://nssdc.gsfc.nasa.gov/planetary/factsheet/

#to Do
#Add scale to plots
#Change plots to animation
#logarithmically scale plots to fit all planets on

from __future__ import annotations
from cProfile import label  

import math
from telnetlib import TM
#from astropy.constants import G as gravitationalConstant
import matplotlib.pyplot as plt
from Vector import Vector
gravitationalConstant = 6.67e-11
planets = {
            "Mercury": {"Mass":0.33e24,
                        "OrbitalVelocity":47.4e3,
                        "DistanceFromSun":57.9e9},
            "Venus": {"Mass":4.87e24,
                      "OrbitalVelocity":35e3,
                      "DistanceFromSun":108.2e9},
            "Earth": {"Mass":5.97e24,
                      "OrbitalVelocity":28.8e3,
                      "DistanceFromSun":149.5e9},
            "Mars": {"Mass":0.642e24,
                     "OrbitalVelocity":24.1e3,
                     "DistanceFromSun":228e9},
            "Jupiter": {"Mass":1898e24,
                        "OrbitalVelocity":13.1e3,
                        "DistanceFromSun":778.5e9,
                        "OrbitalPeriod":4331*24*60*60}
}

class Body:
    def __init__(self, initialPosition:Vector, initialVelocity:Vector, mass:int, name:str) -> None:
        self.currentPos:Vector = initialPosition
        self.currentVel:Vector = initialVelocity
        self.currentAcc:Vector = Vector(0,0)
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
  
        self.currentAcc = Vector(totalAccelerationX, totalAccelerationY)

    def UpdateVelocity(self, timeStep:float):
        self.currentVel += (self.currentAcc * timeStep)

    def UpdatePosition(self, timeStep:float):
        self.currentPos += (self.currentVel * timeStep)

    def UpdateValues(self, timeStep:float):
        self.UpdateVelocity(timeStep)
        self.UpdatePosition(timeStep)
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

def CalculateSolarParameters(planets:dict) -> Body:
    #Sun's motion and barycenter are goverened mostly by Jupiter
    solarMass = 1.989e30
    distanceFromBarycenter = planets["Jupiter"]["DistanceFromSun"]/(1+solarMass/planets["Jupiter"]["Mass"])
    #In the calculation above neglect size of Sun and Jupiter since seperation is much greater
    #Period of the Sun's orbit is equal to that of Jupiter (using the assumption at the start)
    solarVelocity = 2*math.pi*distanceFromBarycenter/planets["Jupiter"]["OrbitalPeriod"]

    return Body(initialPosition=Vector(-distanceFromBarycenter,0),
                initialVelocity=Vector(0,-solarVelocity),
                mass=solarMass,
                name="Sun")
    
def MakeSolarSystem(planets:dict) -> System:
    tempStorage = []
    tempStorage.append(CalculateSolarParameters(planets))
    for planetName, planetAttributes in planets.items():
        tempStorage.append(Body(initialPosition=Vector(planetAttributes["DistanceFromSun"],0),
                                initialVelocity=Vector(0,planetAttributes["OrbitalVelocity"]),
                                mass=planetAttributes["Mass"],
                                name=planetName))
    return System(tuple(tempStorage))

def KickDriftKick(bodies:System, timeStep:float) -> None:
    for body in bodies: #first calculate the acceleration on each body from all other bodies before position is updated
        body.CalculateTotalAcceleration(bodies) #a(t)

    #velocity_half:float
    for body in bodies: #velocity at half time step then update the positions of each body
        body.currentVel += (body.currentAcc * timeStep/2) #v(t + delta_t/2)
        body.currentPos += (body.currentVel * timeStep) #x(t + delta_t)
        body.orbitalHistory.append(body.currentPos)

    for body in bodies: #calculate the acceleration on each body from all other bodies at the updated position
        body.CalculateTotalAcceleration(bodies) #a(t + delta_t)

    for body in bodies: #Now update the velocity for the new position
        body.currentVel += (body.currentAcc * timeStep/2) #v(t + delta_t)
    
def Simplistic(bodies:System, timeStep:float) -> None:
    for body in bodies: #first calculate the acceleration on each body from all other bodies before position is updated
        body.CalculateTotalAcceleration(bodies)
    for body in bodies:
        body.UpdateValues(timeStep)

def ChangeInEnergyPlot(energy:list[float], time:list[float]) -> None:

    deltaOverEnergy = []
    for i in range(len(energy)-1):
        deltaOverEnergy.append((energy[i+1]-energy[i])/energy[i+1])

    plt.scatter(time[1:], deltaOverEnergy)
    plt.show()

def OrbitsPlot(bodies:System) -> None:
    for body in bodies:
        plt.plot(body.xPositionalHistory, body.yPositionalHistory, label=body.name)
    plt.legend()
    plt.show()

def main():
    astronomicalUnit = 1.5e11
    sunInitialVelocity = math.sqrt(gravitationalConstant * 5.972e24 / 1.5e11)
    earthInitialVelocity = math.sqrt(gravitationalConstant * 1.989e30 / 1.5e11)
    
    Sun = Body(Vector(0,0),Vector(0,-sunInitialVelocity),1.989e30,"Sun")
    Earth = Body(Vector(1.5e11,0),Vector(0,earthInitialVelocity),5.972e24,"Earth")
    
    solarSystem = System((Sun, Earth))

    daysToSeconds = 24*60*60
    runTime = 3650 * daysToSeconds
    dt = 0.1 * daysToSeconds

    currentTime = 0
    totalEnergy = []
    time = [currentTime]

    while currentTime <= runTime:
        KickDriftKick(solarSystem, dt)
        #Simplistic(solarSystem, dt)
        

        totalEnergy.append(solarSystem.TotalEnergyOfSystem)
        currentTime += dt
        time.append(currentTime)

    ChangeInEnergyPlot(energy=totalEnergy, time=time[1:])
    OrbitsPlot(solarSystem)

def mainTest():
    solarSystem = MakeSolarSystem(planets)
    daysToSeconds = 24*60*60
    runTime = 36500 * daysToSeconds
    dt = 0.1 * daysToSeconds

    currentTime = 0
    totalEnergy = []
    time = [currentTime]

    while currentTime <= runTime:
        KickDriftKick(solarSystem, dt)
        totalEnergy.append(solarSystem.TotalEnergyOfSystem)
        currentTime += dt
        time.append(currentTime)

    OrbitsPlot(solarSystem)
    
if __name__ == "__main__":
    mainTest()