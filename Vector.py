#Class for vectors to simplify the access of coordinate based values

from __future__ import annotations      

class Vector(object):
    def __init__(self, x:int|float, y:int|float) -> None:
        self.x:int|float = x
        self.y:int|float = y

    def __repr__(self) -> str:
        return f"Vector({self.x},{self.y})"

    def __add__(self, other:Vector) -> Vector:
        return Vector(self.x+other.x, self.y+other.y)

    def __sub__(self, other:Vector) -> Vector:
        return Vector(self.x-other.x, self.y-other.y)

    def __mul__(self, scalar:int|float) -> Vector:
        return Vector(self.x*scalar, self.y*scalar)
    
    def __truediv__(self, scalar:int|float) -> Vector:
        return Vector(self.x/scalar, self.y/scalar)

