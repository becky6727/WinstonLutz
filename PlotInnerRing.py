import os, sys, time
import numpy
import cv2
import ROOT
import argparse
from scipy.signal import argrelmin
import lib.FitRing as FitRing

#options
parser = argparse.ArgumentParser(description = 'options')

parser.add_argument('-i',
                    type = str,
                    dest = 'FileIN',
                    #nargs = '+',
                    help = 'input tiff files')

parser.add_argument('-BG',
                    type = str,
                    dest = 'FileBG',
                    default = None,
                    help = 'input tiff files for background')

parser.add_argument('-dpi',
                    type = int,
                    dest = 'DPI',
                    default = '72',
                    help = 'input tiff files')

parser.add_argument('-p',
                    action = 'store_true',
                    dest = 'isFigOut',
                    help = 'figure file is created or not')

args = parser.parse_args()

FileIN = args.FileIN
FileBG = args.FileBG
DPI = args.DPI
isFigOut = args.isFigOut

if(not(os.path.isfile(FileIN))):
    print 'No such a file!: %s' %(FileIN)
    sys.exit()
    pass

#---- trimming image ----# 
def Trimming(ImgArray, TrimSize = 30):

    #TrimSize = 30
    Width = ImgArray.shape[1]
    Height = ImgArray.shape[0]
    
    Img = ImgArray[TrimSize:(Height - TrimSize), TrimSize:(Width - TrimSize)]

    return Img
#---- end of function ----#

TrimRegion = 30
Img = cv2.imread(FileIN, -1)
Img = Trimming(Img, TrimRegion)

#laplacian filter
kernel = numpy.array([[1, 1, 1],
                      [1, -8, 1],
                      [1, 1, 1]])

ImgLaplace = cv2.filter2D(Img, -1, kernel)
#ImgInnerRing = cv2.cvtColor(ImgLaplace, cv2.COLOR_BGR2GRAY)
ImgInnerRing = cv2.bitwise_not(ImgLaplace)
ImgInnerRing = cv2.cvtColor(ImgInnerRing, cv2.COLOR_BGR2GRAY)

cv2.imwrite('inner.tif', ImgInnerRing)

#set axis(pixel -> mm)
HeightImg = int(Img.shape[0])
WidthImg = int(Img.shape[1])

#SizeImg = numpy.minimum(HeightImg, WidthImg)

inch2mm = 25.4

XArray = []
YArray = []

for x_pixel in range(WidthImg):
    tmp = x_pixel* (inch2mm/DPI)
    XArray.append(tmp)
    pass

for y_pixel in range(HeightImg):
    tmp = y_pixel* (inch2mm/DPI)
    YArray.append(tmp)
    pass

#fit inner ring
XcArray = []
YcArray = []

MinImg = numpy.min(ImgInnerRing[numpy.where(ImgInnerRing > 5.0e4)])

RangePixel = 1.05
RangeStd = 0.5

for i in range(len(YArray)):
    for j in range(len(XArray)):
        if(ImgInnerRing[i][j] > RangePixel* MinImg):
            continue
        XcArray.append(XArray[j])
        YcArray.append(YArray[i])
        pass
    pass

XcArray = numpy.array(XcArray)
YcArray = numpy.array(YcArray)
RcArray = numpy.sqrt(XcArray**2 + YcArray**2)
RatioArray = YcArray/XcArray

#print MinImg
print RcArray
print XcArray
print YcArray
#print RatioArray

RcMean = numpy.min(numpy.array([numpy.mean(RcArray), numpy.median(RcArray)]))
RcStd = RangeStd* numpy.std(RcArray)
#Ratio = numpy.median(RatioArray)

AngleArray =  numpy.arctan2(YcArray - numpy.mean(YcArray), XcArray - numpy.mean(XcArray))

print RcMean
print RcStd
print AngleArray

XcArray = XcArray[numpy.where((RcArray > RcMean - 1.0* RcStd) & (RcArray < RcMean + RcStd))]
YcArray = YcArray[numpy.where((RcArray > RcMean - 1.0* RcStd) & (RcArray < RcMean + RcStd))]

#print XcArray
#print YcArray
#print numpy.sqrt(XcArray**2 + YcArray**2)

ObjInner = FitRing.FitRing(XcArray, YcArray)
XcIn, YcIn, Rin = ObjInner.DoFit()

print XcIn, YcIn, Rin

#plot
c1 = ROOT.TCanvas('c1', 'c1', 0, 0, 1200, 1000)

ROOT.gStyle.SetOptStat(0)

c1.SetFillColor(0)
#c1.SetGridx()
#c1.SetGridy()
#c1.SetLogy()
c1.Divide(2, 2)

#plot 2-D
Ncanvas = 1

c1.cd(Ncanvas)
c1.cd(Ncanvas).SetGridx()
c1.cd(Ncanvas).SetGridy()

#fill histogram
XBinWidth = inch2mm/DPI
MinX = XArray[0]
MaxX = XArray[-1]
#NofXBin = int((MaxX - MinX)/XBinWidth)
NofXBin = WidthImg

YBinWidth = inch2mm/DPI
MinY = YArray[0]
MaxY = YArray[-1]
#NofYBin = int((MaxY - MinY)/YBinWidth)
NofYBin = HeightImg

hMap = ROOT.TH2D('hMap', 'Winston Lutz', NofXBin, MinX, MaxX, NofYBin, MinY, MaxY)

for i in range(len(YArray)):
    for j in range(len(XArray)):
        if(ImgInnerRing[i][j] > RangePixel* MinImg):
            continue
        hMap.Fill(XArray[j], YArray[i], ImgInnerRing[i][j])
        pass
    pass

hMap.Draw('colz')

#plot pixel distribution
Ncanvas = 2

c1.cd(Ncanvas)
c1.cd(Ncanvas).SetGridx()
c1.cd(Ncanvas).SetGridy()
c1.cd(Ncanvas).SetLogy()

PixelBinWidth = 2**7
MinPix = 2**15
MaxPix = 2**16 + PixelBinWidth
NofPixelBin = int((MaxPix - MinPix)/PixelBinWidth)

hPixel = ROOT.TH1D('hPixel', 'Pixel distribution', NofPixelBin, MinPix, MaxPix)

for i in range(len(ImgInnerRing)):
    for j in range(len(ImgInnerRing[i])):
        if(ImgInnerRing[i][j] > RangePixel* MinImg):
            continue
        hPixel.Fill(ImgInnerRing[i][j])
        pass
    pass

hPixel.Draw('hist')

#angle
Ncanvas = 3

c1.cd(Ncanvas)
c1.cd(Ncanvas).SetGridx()
c1.cd(Ncanvas).SetGridy()
#c1.cd(Ncanvas).SetLogy()

AngleBinWidth = 0.155
MinAngle = numpy.min(AngleArray) - AngleBinWidth
MaxAngle = numpy.max(AngleArray) + AngleBinWidth
NofAngleBin = int((MaxAngle - MinAngle)/AngleBinWidth)

hAngle = ROOT.TH1D('hAngle', 'Angle distribution', NofAngleBin, MinAngle, MaxAngle)

for i in range(len(AngleArray)):
    hAngle.Fill(AngleArray[i])
    pass

hAngle.Draw('hist')

c1.Update()
