import os, sys, time
import numpy
import cv2
import ROOT
import argparse
from scipy.signal import argrelmin
import lib.FitRing as FitRing
import re

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

#--- fitting for outer ring ---#
def GetOuterRing(Data, BG, XpixArray, YpixArray):

    #pixels
    NofXpixel = len(XpixArray)
    NofYpixel = len(YpixArray)
    
    #pixel at red-color
    RedArray = [[] for i in range(NofYpixel)]
    
    for i in range(NofYpixel):
        for j in range(NofXpixel):
            RedArray[i].append(Data[i][j][2])
            pass
        pass
    
    #background
    BGArray = []
    
    for i in range(int(BG.shape[0])):
        for j in range(int(BG.shape[1])):
            BGArray.append(BG[i][j][2])
            pass
        pass

    ImgOuterRing = numpy.absolute(RedArray - numpy.mean(numpy.array(BGArray)))

    #for arrays to search center
    XcArray = []
    YcArray = []

    for i in range(NofYpixel):
        for j in range(NofXpixel):

            #FY2019
            if((ImgOuterRing[i][j]/numpy.max(ImgOuterRing) < 0.4) or (ImgOuterRing[i][j]/numpy.max(ImgOuterRing) > 0.5)):
            #if((ImgOuterRing[i][j]/numpy.max(ImgOuterRing) < 0.5) or (ImgOuterRing[i][j]/numpy.max(ImgOuterRing) > 5.0)):
                continue
            XcArray.append(XpixArray[j])
            YcArray.append(YpixArray[i])
            pass
        pass

    XcArray = numpy.array(XcArray)
    YcArray = numpy.array(YcArray)

    ObjFit = FitRing.FitRing(XcArray, YcArray)
    Xc, Yc, Rc = ObjFit.DoFit()

    return Xc, Yc, Rc
#--- end of def ---#

#--- fitting for inner ring ---#
def GetInnerRing(Data, XpixArray, YpixArray):

    #fit parameter
    RangePixel = 1.08
    RangeStd = 0.5

    #pixels
    NofXpixel = len(XpixArray)
    NofYpixel = len(YpixArray)
    
    #laplacian filter(edge filter)
    kernel = numpy.array([[1, 1, 1],
                          [1, -8, 1],
                          [1, 1, 1]])

    ImgLaplace = cv2.filter2D(Data, -1, kernel)

    #invert and convert gray scale
    ImgInnerRing = cv2.bitwise_not(ImgLaplace)
    ImgInnerRing = cv2.cvtColor(ImgInnerRing, cv2.COLOR_BGR2GRAY)

    #for arrays to search center
    XcArray = []
    YcArray = []

    ThPixel = numpy.min(ImgInnerRing[numpy.where(ImgInnerRing > 5.0e4)])
    
    for i in range(NofYpixel):
        for j in range(NofXpixel):
            if(ImgInnerRing[i][j] > RangePixel* ThPixel):
                continue
            XcArray.append(XpixArray[j])
            YcArray.append(YpixArray[i])
            pass
        pass

    XcArray = numpy.array(XcArray)
    YcArray = numpy.array(YcArray)

    #filtering
    XcArray = XcArray[numpy.where((YcArray > numpy.median(YcArray) - 2.5) & (YcArray < numpy.median(YcArray) + 2.5))]
    YcArray = YcArray[numpy.where((YcArray > numpy.median(YcArray) - 2.5) & (YcArray < numpy.median(YcArray) + 2.5))]
    YcArray = YcArray[numpy.where((XcArray > numpy.median(XcArray) - 2.5) & (XcArray < numpy.median(XcArray) + 2.5))]
    XcArray = XcArray[numpy.where((XcArray > numpy.median(XcArray) - 2.5) & (XcArray < numpy.median(XcArray) + 2.5))]

    RcArray = numpy.sqrt(XcArray**2 + YcArray**2)

    RcMean = numpy.min(numpy.array([numpy.mean(RcArray), numpy.median(RcArray)]))
    RcStd = RangeStd* numpy.std(RcArray)
        
    XcArray = XcArray[numpy.where((RcArray > RcMean - 1.0* RcStd) & (RcArray < RcMean + RcStd))]
    YcArray = YcArray[numpy.where((RcArray > RcMean - 1.0* RcStd) & (RcArray < RcMean + RcStd))]

    ObjFit = FitRing.FitRing(XcArray, YcArray)
    Xc, Yc, Rc = ObjFit.DoFit()

    return Xc, Yc, Rc
#--- end of def ---#

#---- trimming image ----# 
def Trimming(ImgArray, TrimSize = 30):

    #TrimSize = 30
    Width = ImgArray.shape[1]
    Height = ImgArray.shape[0]
    
    Img = ImgArray[TrimSize:(Height - TrimSize), TrimSize:(Width - TrimSize)]

    return Img
#---- end of function ----#

#read data
Pattern = './\w+/(\d+)/\w(\d+)_\w(\d+)_\w(\d+).tif'
rePattern = re.compile(Pattern)

FileName = r'%s' %(FileIN)
InfoList = rePattern.findall(FileName)

print 'FY = %s' %(InfoList[0][0])
print '(Enegy = %s MeV, Gantry = %s, Couch = %s)' %(InfoList[0][1], InfoList[0][2], InfoList[0][3])

FY = int(InfoList[0][0])
Energy = int(InfoList[0][1])
Gantry = int(InfoList[0][2])
Couch = int(InfoList[0][3])

Img = cv2.imread(FileIN, -1)
ImgBG = cv2.imread(FileBG, -1)

#triming
if(FY >= 2018):
    TrimRegion = 35
    Img = Trimming(Img, TrimRegion)
    ImgBG = Trimming(ImgBG, TrimRegion)
    pass

#set axis(pixel -> mm)
HeightImg = int(Img.shape[0])
WidthImg = int(Img.shape[1])

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

print 'TIFF size = (%d, %d)' %(Img.shape[0], Img.shape[1])
print 'DPI = %d dpi' %(DPI)

XcOut, YcOut, RcOut = GetOuterRing(Img, ImgBG, XArray, YArray)
XcIn, YcIn, RcIn = GetInnerRing(Img, XArray, YArray)

#diff
diff = numpy.sqrt((XcIn - XcOut)**2 + (YcIn - YcOut)**2)

print 'Outer ring: (%+2.4f, %+2.4f)' %(XcOut, YcOut)
print 'Inner ring: (%+2.4f, %+2.4f)' %(XcIn, YcIn)
print 'difference of radius: %2.4f mm' %(diff)

#pixel at red-color
RedArray = [[] for i in range(HeightImg)]

for i in range(HeightImg):
    for j in range(WidthImg):
        RedArray[i].append(Img[i][j][2])
        pass
    pass

#background
tmp = []

for i in range(int(ImgBG.shape[0])):
    for j in range(int(ImgBG.shape[1])):
        tmp.append(ImgBG[i][j][2])
        pass
    pass

RedArray = numpy.absolute(RedArray - numpy.mean(numpy.array(tmp)))
    
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
NofXBin = WidthImg

YBinWidth = inch2mm/DPI
MinY = YArray[0]
MaxY = YArray[-1]
#NofYBin = int((MaxY - MinY)/YBinWidth)
NofYBin = HeightImg

hMap = ROOT.TH2D('hMap', '%dMeV, G%03d, C%03d' %(Energy, Gantry, Couch),
                 NofXBin, MinX, MaxX,
                 NofYBin, MinY, MaxY)

for i in range(len(YArray)):
    for j in range(len(XArray)):
        hMap.Fill(XArray[j], YArray[i], RedArray[i][j])
        pass
    pass

hMap.Draw('colz')

#make circle(outer)
XOutRingArray1 = []
YOutRingArray1 = []
XOutRingArray2 = []
YOutRingArray2 = []

tmpXArray = numpy.arange(XArray[0], XArray[-1], 0.0005)

for i in range(len(tmpXArray)):
    if(RcOut**2 - (tmpXArray[i] - XcOut)**2 > 0.0):
        tmp = YcOut + numpy.sqrt(RcOut**2 - (tmpXArray[i] - XcOut)**2)
        XOutRingArray1.append(tmpXArray[i])
        YOutRingArray1.append(tmp)
        pass
    pass

for i in range(len(tmpXArray)):
    if(RcOut**2 - (tmpXArray[i] - XcOut)**2 > 0.0):
        tmp = YcOut - numpy.sqrt(RcOut**2 - (tmpXArray[i] - XcOut)**2)
        XOutRingArray2.append(tmpXArray[i])
        YOutRingArray2.append(tmp)
        pass
    pass

gOutRing1 = ROOT.TGraph(len(XOutRingArray1), numpy.array(XOutRingArray1), numpy.array(YOutRingArray1))
gOutRing2 = ROOT.TGraph(len(XOutRingArray2), numpy.array(XOutRingArray2), numpy.array(YOutRingArray2))

gOutRing1.SetLineWidth(3)
gOutRing1.SetLineStyle(2)

gOutRing2.SetLineWidth(3)
gOutRing2.SetLineStyle(2)

gOutRing1.Draw('lsame')
gOutRing2.Draw('lsame')

#make circle(inner)
XInnerRingArray1 = []
YInnerRingArray1 = []
XInnerRingArray2 = []
YInnerRingArray2 = []

for i in range(len(tmpXArray)):
    if(RcIn**2 - (tmpXArray[i] - XcIn)**2 > 0.0):
        tmp = YcIn + numpy.sqrt(RcIn**2 - (tmpXArray[i] - XcIn)**2)
        XInnerRingArray1.append(tmpXArray[i])
        YInnerRingArray1.append(tmp)
        pass
    pass

for i in range(len(tmpXArray)):
    if(RcIn**2 - (tmpXArray[i] - XcIn)**2 > 0.0):
        tmp = YcIn - numpy.sqrt(RcIn**2 - (tmpXArray[i] - XcIn)**2)
        XInnerRingArray2.append(tmpXArray[i])
        YInnerRingArray2.append(tmp)
        pass
    pass

gInnerRing1 = ROOT.TGraph(len(XInnerRingArray1), numpy.array(XInnerRingArray1), numpy.array(YInnerRingArray1))
gInnerRing2 = ROOT.TGraph(len(XInnerRingArray2), numpy.array(XInnerRingArray2), numpy.array(YInnerRingArray2))

gInnerRing1.SetLineWidth(3)
gInnerRing1.SetLineStyle(2)

gInnerRing2.SetLineWidth(3)
gInnerRing2.SetLineStyle(2)

gInnerRing1.Draw('lsame')
gInnerRing2.Draw('lsame')

c1.Update()

#create figure
if(isFigOut):
    if(Energy > 225):
        GIF = './Figures/M%d_G%03d_C%03d_%d.gif' %(Energy, Gantry, Couch, FY)
    else:
        GIF = './Figures/S%d_G%03d_C%03d_%d.gif' %(Energy, Gantry, Couch, FY)
        pass
    c1.Print(GIF)
    pass

