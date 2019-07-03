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

if(not(os.path.isfile(FileBG))):
    print 'No such a file!: %s' %(FileBG)
    sys.exit()
    pass

Img = cv2.imread(FileIN, -1)
ImgBG = cv2.imread(FileBG, -1)

#set axis(pixel -> mm)
HeightImg = int(Img.shape[0])
WidthImg = int(Img.shape[1])

SizeImg = numpy.minimum(HeightImg, WidthImg)

inch2mm = 25.4

XArray = []
YArray = []

for x_pixel in range(SizeImg):
    tmp = x_pixel* (inch2mm/DPI)
    XArray.append(tmp)
    pass

for y_pixel in range(SizeImg):
    tmp = y_pixel* (inch2mm/DPI)
    YArray.append(tmp)
    pass

#show red-pixel
RedArray = [[] for i in range(SizeImg)]

for i in range(SizeImg):
    for j in range(SizeImg):
        RedArray[i].append(Img[i][j][2])
        pass
    pass

#background
BGArray = []

for i in range(int(ImgBG.shape[0])):
    for j in range(int(ImgBG.shape[1])):
        BGArray.append(ImgBG[i][j][2])
        pass
    pass

BG = numpy.mean(numpy.array(BGArray))

ImgOuterRing = numpy.absolute(RedArray - BG)
#ImgOuterRing = RedArray - BG

#fit outer ring
XcArray = []
YcArray = []

for i in range(len(YArray)):
    for j in range(len(XArray)):
        if((ImgOuterRing[i][j]/numpy.max(ImgOuterRing) < 0.3) or (ImgOuterRing[i][j]/numpy.max(ImgOuterRing) > 0.5)):
            continue
        XcArray.append(XArray[j])
        YcArray.append(YArray[i])
        pass
    pass

XcArray = numpy.array(XcArray)
YcArray = numpy.array(YcArray)
RcArray = numpy.sqrt(XcArray**2 + YcArray**2)

ObjOuter = FitRing.FitRing(XcArray, YcArray)
XcOut, YcOut, Rout = ObjOuter.DoFit()

print XcOut, YcOut, Rout

#plot
c1 = ROOT.TCanvas('c1', 'c1', 0, 0, 1200, 1000)

ROOT.gStyle.SetOptStat(0)

c1.SetFillColor(0)
c1.SetGridx()
c1.SetGridy()
#c1.SetLogy()
#c1.Divide(2, 2)

#plot 2-D
#fill histogram
XBinWidth = inch2mm/DPI
MinX = XArray[0]
MaxX = XArray[-1]
#NofXBin = int((MaxX - MinX)/XBinWidth)
NofXBin = SizeImg

YBinWidth = inch2mm/DPI
MinY = YArray[0]
MaxY = YArray[-1]
#NofYBin = int((MaxY - MinY)/YBinWidth)
NofYBin = SizeImg

hMap = ROOT.TH2D('hMap', 'Winston Lutz', NofXBin, MinX, MaxX, NofYBin, MinY, MaxY)

for i in range(len(YArray)):
    for j in range(len(XArray)):
        if((ImgOuterRing[i][j]/numpy.min(ImgOuterRing) < 0.3) or (ImgOuterRing[i][j]/numpy.min(ImgOuterRing) > 0.5)):
            continue
        hMap.Fill(XArray[j], YArray[i], ImgOuterRing[i][j])
        pass
    pass

hMap.Draw('colz')

c1.Update()
