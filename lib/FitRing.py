import os, sys
import numpy, math

class FitRing:
    def __init__(self, XArray, YArray):

        #initialize
        self.XArray = numpy.array(XArray)
        self.YArray = numpy.array(YArray)

    def DoFit(self):

        #set variables
        SumX = sum(self.XArray)
        SumY = sum(self.YArray)

        SumX2 = sum([iX**2 for iX in self.XArray])
        SumY2 = sum([iY**2 for iY in self.YArray])
        
        SumXY = sum([iX* iY for (iX, iY) in zip(self.XArray, self.YArray)])

        #calc matrix elements
        F = numpy.array([[SumX2, SumXY, SumX],
                         [SumXY, SumY2, SumY],
                         [SumX, SumY, len(self.XArray)]])

        G = numpy.array([[-sum([iX**3 + iX* iY**2 for (iX, iY) in zip(self.XArray, self.YArray)])],
                         [-sum([iX**2 *iY + iY**3 for (iX, iY) in zip(self.XArray, self.YArray)])],
                         [-sum([iX**2 + iY**2 for (iX, iY) in zip(self.XArray, self.YArray)])]])

        T = numpy.linalg.inv(F).dot(G)

        CenterX = float(T[0]/-2.0)
        CenterY = float(T[1]/-2.0)
        Radius = math.sqrt(CenterX**2 + CenterY**2 - T[2])
        
        return (CenterX, CenterY, Radius)

