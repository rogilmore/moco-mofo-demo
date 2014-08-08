#!/usr/bin/env python2
""" Motion form (MOFO) stimulus display demo
    Simulates displays in Fesi et al. 2014
"""
#-------------------------------------------------------
#   author: Rick Gilmore, thatrickgilmore@gmail.com
#   based on starField, captureFrames, screensAndWindows demos from PsychoPy2
#
#   Creates an onscreen display and if saveMovie=True,
#   writes it to a series of .jpg images that can be 
#   converted into a movie using MPEG Streamclip or
#   a similar application.

#-------------------------------------------------------
# Version history
# 2014-08-03    v.01    based on moco-demo.py

#-------------------------------------------------------
# Known bugs, desired enhancements
# - Want to generate square dot regions.
# - targetDots does not overwrite bgdDots, so will have to detect target/background dots internally
# - Screen for MOFO was probably 800 x 600. Should change.

#-------------------------------------------------------
# Import dependencies
from psychopy import visual, event, core
from psychopy.tools.coordinatetools import pol2cart
import numpy, math

# Profile of show_figure( )
def show_figure( frames, bgdDotsX, bgdDotsY, targetDotsX, targetDotsY, bgdDir, targetDir, 
        bgdSpd, targetSpd, targetRadius, flowMode ):
            
    for frameN in range(frames): 
        # Calculate background move
        dX = bgdSpd*numpy.cos(bgdDir)
        dY = bgdSpd*numpy.sin(bgdDir)
        bgdDotsX += dX
        bgdDotsY += dY
        
        # If out of range, reposition
        outVertical = ( abs( bgdDotsY ) >= 1 )       # if units = 'norm' then X and Y in [-1, 1]
        outHorizontal = ( abs( bgdDotsX ) >= 1 )
        
        # Wrap out of boundary dots horizontally or vertically
        bgdDotsY[ outVertical ] = -1.0*bgdDotsY[ outVertical ]
        bgdDotsX[ outHorizontal ] = -1.0*bgdDotsX[ outHorizontal ] 
        
        # Calculate figure move
        dX = targetSpd*numpy.cos(targetDir)
        dY = targetSpd*numpy.sin(targetDir)
        targetDotsX += dX
        targetDotsY += dY
    
        # If figure out of range, reposition
        outVertical = ( abs( targetDotsY ) >= targetRadius )       # if units = 'norm' then X and Y in [-1, 1]
        outHorizontal = ( abs( targetDotsX ) >= targetRadius )
        
        # Draw background
        bgdDots.setXYs(numpy.array([bgdDotsX, bgdDotsY]).transpose())
        bgdDots.draw()
        
        # Draw figure
        targetDots.setXYs(numpy.array([targetDotsX, targetDotsY]).transpose())
        if (flowMode == 'flow') :
            targetDots.setFieldPos([ frameN*dX, frameN*dY ])
        targetDots.draw()
        
        # Draw fixation
        fixation.draw()
        
        # If generating movie, save frames
        if saveMovie:
            win.getMovieFrame(buffer='back')
            
        # flip buffer
        win.flip()
    # End for frameN

    # Return
    return bgdDotsX, bgdDotsY, targetDotsX, targetDotsY
    
# end def show_figure


def show_laminar( frames, dotsRadius, dotsTheta, dotsX, dotsY, spd, dir, innerRadius, outerRadius, saveMovie ):
    for frameN in range(frames):
        # Calculate move
        dX = spd*numpy.cos(dir)
        dY = spd*numpy.sin(dir)
        targetDotsX += dX
        targetDotsY += dY

        # Polar coordinates
#        dotsRadius = numpy.sqrt( dotsX**2 + dotsY**2 )
#        dotsTheta = numpy.arctan2( dotsY, dotsX )*degPerRad
#
#        # Determine who's out of range and regenerate values in range
#        outFieldDots = (dotsRadius >= outerRadius)
#        dotsRadius[outFieldDots] = innerRadius + spd
#        
#        inFieldDots = (dotsRadius <= innerRadius)
#        dotsRadius[inFieldDots] = outerRadius - spd
        
        # Convert to Cartestian coordinates and set positions
#        dotsX, dotsY = pol2cart(dotsTheta, dotsRadius)
        dots.setXYs(numpy.array([dotsX, dotsY]).transpose())
        
        # Draw stimuli to window
        dots.draw()
        fixation.draw()
        
        targetDots.setPos([ dX, dY ] )
        targetDots.draw()
        
        if saveMovie:
            win.getMovieFrame(buffer='back')
        win.flip()
    # end for frameN
    return dotsRadius, dotsTheta, dotsX, dotsY
# end def show_laminar

def show_noise( frames, dotsRadius, dotsTheta, dotsX, dotsY, spd, innerRadius, outerRadius, saveMovie ):
  # Incoherent
    for frameN in range(frames):
        # Calculate move
        dotsDtheta = numpy.random.rand(nDots)*360
        dotsX += spd*numpy.cos( dotsDtheta )
        dotsY += spd*numpy.sin( dotsDtheta )
        
        # Polar coordinates
        dotsRadius = numpy.sqrt( dotsX**2 + dotsY**2 )
        dotsTheta = numpy.arctan2( dotsY, dotsX )*degPerRad
        
        # Determine who's out of range and regenerate values in range
        outFieldDots = (dotsRadius >= outerRadius)
        dotsRadius[outFieldDots] = innerRadius + spd
        
        inFieldDots = (dotsRadius <= innerRadius)
        dotsRadius[inFieldDots] = outerRadius - spd

        # Assign dot positions to array, draw, and flip
        dotsX, dotsY = pol2cart(dotsTheta,dotsRadius)
        dots.setXYs(numpy.array([dotsX, dotsY]).transpose())
        dots.draw()
        fixation.draw()
        if saveMovie:
            win.getMovieFrame(buffer='back')
        win.flip()
    # end for frameN
    return dotsRadius, dotsTheta, dotsX, dotsY
# end def show_noise

# Fonts
sans = ['Gill Sans MT', 'Arial','Helvetica','Verdana'] #use the first font found on this list

# Specify constants, parameters for Fesi, Thomas, Gilmore 2014
windowPix = 600     # Window is square 600 x 600
monitorFramesPerSec = 72
updateFramesPerSec = 24
degPerRad = 180.0/math.pi
radPerDeg = 1/degPerRad

nLoops = 5
modulationHz = 1.2
cycleFrames = updateFramesPerSec/modulationHz
halfCycleFrames = int( cycleFrames/2.0 )
dotLifeFrames = 100

stimSpeeds = { '1' : 0.0017,          # 1 deg/s
               '2' : 0.0034,          # 2 deg/s
               '4' : 0.0068,          # 4 deg/s
               '8' : 0.0136,          # 8 deg/s
               '16' : 0.0272          # 16 deg/s
                }

#stimType = {   'radial' : 'radial-movie.jpg',
#                'laminar' : 'laminar-movie.jpg',
#                'rotation' : 'rotation-movie.jpg'
#            }
            
#movieType = 'radial'
#speedType = '16'
#movieName = movieType + '-' + speedType + '-degPerSec-.jpg'

#dotSpeedUnits = stimSpeeds[speedType]
#fieldSizeUnits = 1.5

saveMovie = 1

nDotsTotal = 2600
dotSize = .005     # 7 arcmin or 7/60=0.117 deg. 0.117 deg/dot / 24 deg/field diam

bgdDir = 0.0           # radians
targetDir = math.pi

bgdSpd = .002
targetSpd = .002
flowMode = 'flow'

bgdRadius = 1.0
targetRadius = 0.25
targetAreaProportion = (targetRadius**2)/(bgdRadius**2)

nBgdDots = round( nDotsTotal*(1-targetAreaProportion) )
nTargetDots = nDotsTotal - nBgdDots

# annulus parameters
#innerRadius = 0.19
#outerRadius = 1.0

# Generate initial dot positions
#bgdDotsX = numpy.random.uniform(-1, 1, nBgdDots)
#bgndDotsY = numpy.random.uniform(-1, 1, nBgdDots)
#
#targetDotsX = numpy.random.uniform(-1, 1, nTargetDots)
#targetDotsY = numpy.random.uniform(-1, 1, nTargetDots)

dotsTheta=numpy.random.rand(nBgdDots)*360
dotsRadius=(numpy.random.rand(nBgdDots)**0.5)
bgdDotsX, bgdDotsY = pol2cart(dotsTheta,dotsRadius)

dotsTheta=numpy.random.rand(nTargetDots)*360
dotsRadius=(numpy.random.rand(nTargetDots)**0.5)*targetRadius
targetDotsX, targetDotsY = pol2cart(dotsTheta,dotsRadius)

# Create a window
win =visual.Window( (windowPix,windowPix), allowGUI=False, color=[-1,-1,-1],
    bitsMode=None, units='norm', winType='pyglet')
    
# Create stimuli
bgdDots = visual.ElementArrayStim(win, units='norm', elementTex=None, elementMask='circle', fieldShape='square',
    nElements=nBgdDots, sizes=dotSize)
    
targetDots = visual.ElementArrayStim(win, units='norm', elementTex=None, elementMask='circle', fieldShape='square',
    nElements=nTargetDots, sizes=dotSize)

fixation = visual.TextStim(win,text="+",pos=(0, 0), color=[1,1,1], ori=0, height = 0.1, font=sans)

# Prepare to loop
trialClock = core.Clock()

# Loop
for loopN in range(nLoops):
    bgdDotsX, bgdDotsY, targetDotsX, targetDotsY = show_figure( halfCycleFrames, bgdDotsX, bgdDotsY, targetDotsX, targetDotsY, bgdDir, 
        targetDir, bgdSpd, targetSpd, targetRadius, flowMode )
    targetDir = bgdDir
    bgdDotsX, bgdDotsY, targetDotsX, targetDotsY = show_figure( halfCycleFrames, bgdDotsX, bgdDotsY, targetDotsX, targetDotsY, bgdDir, 
        targetDir, bgdSpd, targetSpd, targetRadius, flowMode )
    targetDir = math.pi 
#end for loopN

# Save movie
#if saveMovie:
#    win.saveMovieFrames(movieName)
#
#
