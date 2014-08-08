  #!/usr/bin/env python2
""" Motion coherence (MOCO) stimulus display demo
    Simulates displays in Fesi et al. 2014
"""
#-------------------------------------------------------------------------
#   author: Rick Gilmore, thatrickgilmore@gmail.com
#   based on starField, captureFrames demos from PsychoPy2
#
#   Creates an onscreen display and if saveMovie=True,
#   writes it to a series of .jpg images that can be 
#   converted into a movie using MPEG Streamclip or
#   a similar application.

#-------------------------------------------------------------------------
# Version history
# 2014-07-31    v.01    first working version
# 2014-08-01    v.02    functions for each half cycle
# 2014-08-02    v.03    added uniform speed radial motion, dictionary for speed, movieName
# 2014-08-02    v.04    made generator functions return dot positions for continuity
# 2014-08-04    v.05    refactored with new move_dots, replot_out_dots functions; added refreshing
# 2014-08-08    v.06    further refactoring, streamlining. Fixed weird coordinate conversion bug.

#-------------------------------------------------------------------------
# Known bugs, desired enhancements
# 2014-08-04    wrapping now encapsulated in separate function, but still not correct
# 2014-08-08    replot-scaled-polar replot mode should work for radial in, but there is
#               density reduction at outer radius. Laminar wrap now broken

#-------------------------------------------------------------------------
# Import dependencies
from psychopy import visual, event, core
from psychopy.tools.coordinatetools import pol2cart, cart2pol
import numpy, math

#-------------------------------------------------------------------------
# define functions

def move_dots( r, th, X, Y, displPerUpdate, dirRads, moveMode ):
    # Lengths of radius, theta, X, Y all equal?
    if not ( len(r)==len(th)==len(X)==len(Y) ):
        print 'Coordinate vectors not of equal length.'
        return radius, thetaRads, X, Y
    else:
        nDots = len(r)
    
    if moveMode == 'radial': 
        r += displPerUpdate
        X, Y = pol2cart(th, r, units='rad')
    elif moveMode == 'rotation-eq-angle':
        th += displPerUpdate
        X, Y = pol2cart(th, r, units='rad')
    elif moveMode == 'rotation-eq-spd':
        X += displPerUpdate*numpy.sin( numpy.pi - th )
        Y += displPerUpdate*numpy.cos( numpy.pi - th )
        th, r = cart2pol( X, Y, units ='rad' )
    elif moveMode == 'laminar': 
        X += displPerUpdate*numpy.cos( dirRads )
        Y += displPerUpdate*numpy.sin( dirRads )
        th, r = cart2pol( X, Y, units ='rad' )
    elif moveMode == 'random':
        dTheta = numpy.random.rand(nDots)*numpy.pi*2.0
        X += displPerUpdate*numpy.cos( dTheta )
        Y += displPerUpdate*numpy.sin( dTheta )
        th, r = cart2pol( X, Y, units ='rad' )
    elif moveMode == 'pol-cart':
        X, Y = pol2cart(th, r, units='rad')
    elif moveMode == 'cart-pol':
        th, r = cart2pol( X, Y, units ='rad' )

    # end if moveMode
    return r, th, X, Y
# end def move_dots
#-------------------------------------------------------------------------

def replot_out_dots( radius, thetaRads, X, Y, spd, dirRads, innerRadius, outerRadius, moveMode, replotMode ):
        
    if not( len(radius)==len(thetaRads)==len(X)==len(Y) ):
        print 'Coordinate vectors not of equal length.'
        return radius, thetaRads, X, Y
    else:
        nDots = len(radius)
    
    # Determine who's out of range
    outField = (radius >= outerRadius)
    inField = (radius <= innerRadius)
        
    # Replot based on replotMode and moveMode
    if replotMode == 'wrap':
        if (moveMode == 'radial') or (moveMode == 'rotation') or (moveMode == 'random'):
            if sum( outField ) :
                radius[ outField ] = innerRadius+.001 
            if sum( inField ) :
                radius[ inField ] = outerRadius-.001
        elif (moveMode == 'laminar'):
            if sum( outField ) :
                thetaRads[outField] += numpy.pi         # Rotate position by pi
                radius[outField] = outerRadius          # Set inside outerRadius
            if sum( inField ) :
                thetaRads[inField] += numpy.pi
                radius[inField] = innerRadius
        X, Y = pol2cart( thetaRads, radius, units='rads' )
    elif replotMode == 'replot-scaled-polar':
        if sum( outField ):
            radius[outField], thetaRads[outField], X[outField], Y[outField] = make_dots( sum( outField ), distributionMode = 'scaled-polar' )
        if sum( inField ):
            radius[inField], thetaRads[inField], X[inField], Y[inField] = make_dots( sum( inField ), distributionMode = 'scaled-polar' )
    elif replotMode == 'replot-radial':
        if sum( outField ):
            radius[ outField ] = new_radius( sum(outField), 0, outerRadius )
            thetaRads[ outField ] = numpy.random.rand(sum(outField))*math.pi*2.0
        if sum( inField ):
            radius[ inField]  = new_radius( sum(inField), 0, outerRadius )
            thetaRads[ inField ] = numpy.random.rand(sum(inField))*math.pi*2.0
        X, Y = pol2cart(thetaRads, radius, units='rad')
    elif replotMode == 'replot-rect':
        if sum( outField ):
            X[ outField ] = numpy.random.rand( sum(outField ) )*2.0*outerRadius - outerRadius # keep in [-1,1]
            Y[ outField ] = numpy.random.rand( sum(outField ) )*2.0*outerRadius - outerRadius # keep in [-1,1]

        # New dot positions outside of outerRadius?
        outAgain = ( numpy.sqrt(X**2 + Y**2) >= outerRadius )
        if sum( outAgain ):
            X[ outAgain ] = numpy.random.rand( sum(outAgain) )*numpy.sqrt(2)*outerRadius*1.25 - numpy.sqrt(2)/2*outerRadius # keep in [-.707, .707]
            Y[ outAgain ] = numpy.random.rand( sum(outAgain) )*numpy.sqrt(2)*outerRadius*1.25 - numpy.sqrt(2)/2*outerRadius # keep in [-1,1]
        
        # Convert changed X,Y to polar
        thetaRads, radius = cart2pol( X, Y, units='rad' )
    
    # end if replotMode
    
    return radius, thetaRads, X, Y
# end def replot_out_dots
#-------------------------------------------------------------------------

def show_dots( frames, frameIndex0, dotsRadius, dotsTheta, dotsX, dotsY, refreshIndex, spd, dirRads, 
                innerRadius, outerRadius, moveMode, replotMode, saveMovie, singleFrameMode ):
    for frameN in range( frames ):
        # refresh dots
        frameIndex = frameN + frameIndex0
        refreshThese = ( refreshIndex == frameIndex )
        if sum( refreshThese ):
            dotsRadius[ refreshThese] = new_radius( sum(refreshThese), innerRadius, outerRadius )
            dotsTheta[ refreshThese ] = numpy.random.rand(sum(refreshThese))*math.pi*2.0
               
        # move dots
        dotsRadius, dotsTheta, dotsX, dotsY = move_dots( dotsRadius, dotsTheta, dotsX, dotsY, spd, dirRads, moveMode )
            
        # reposition if out
        dotsRadius, dotsTheta, dotsX, dotsY = replot_out_dots( dotsRadius, dotsTheta, dotsX, dotsY, spd, dirRads, 
                                                                innerRadius, outerRadius, moveMode, replotMode )
        
        # Check one last time that all dots are in, shouldn't be necessary...
        stillOut = ( dotsRadius >= outerRadius )
        dotsX[ stillOut ] = 0.0
        dotsY[ stillOut ] = 0.0
        
        # Assign dot positions to stimulus
        dots.setXYs( numpy.array( [dotsX, dotsY] ).transpose())
        
        # Draw stimuli to window 
        dots.draw()
        blankCircle.draw()
        fixation.draw() 
        
        # Save movie
        if saveMovie:
            win.getMovieFrame(buffer='back')  
            
        # flip buffer to display
        if singleFrameMode:
            msg =  visual.TextStim(win, text=str(frameN), pos=[0,-0.8], color=[1,1,1])
            msg.draw()
            k = ['']
            while k[0] not in ['space', 'esc']:   
                k = event.waitKeys()
        
        win.flip()
        
    # end for frameN
    return dotsRadius, dotsTheta, dotsX, dotsY, frameN # for next round
# end def show_dots
#-------------------------------------------------------------------------

def make_dots( nDots, min=-1, max=1, distributionMode='uniform-rect' ):
    if distributionMode == 'scaled-polar':      # Adjusts for density change due to outward radial motion
        r = numpy.random.rand(nDots)**(0.5)*(max - min) + min
        th = numpy.random.rand(nDots)*(2*numpy.pi)-numpy.pi
        X, Y = pol2cart( th, r, units='rad' )
    elif distributionMode == 'uniform-polar':
        r = numpy.random.rand(nDots)*(max - min) + min
        th = numpy.random.rand(nDots)*(2*numpy.pi)-numpy.pi
        X, Y = pol2cart( th, r, units='rad' )
    else :
        X, Y = numpy.random.rand( nDots )*(max - min) + min, numpy.random.rand( nDots )*(max - min) + min
        th, r = cart2pol( X, Y, units = 'rads' )
        
    return r, th, X, Y

#-------------------------------------------------------------------------

def new_radius( nDots, innerRadius, outerRadius ):
    dotsRadius = numpy.random.rand(nDots)**(.5)*(outerRadius-innerRadius)+innerRadius
    return dotsRadius


#-------------------------------------------------------------------------

#def make_all (win, halfCycleFrames, dotsRadius0, dotsTheta0, dotsX0, dotsY0, speedRadial, 
#    for m in movieTypes:
#        for s in speedTypes: 
#            movieType = m 
#            speedRadial = stimSpeeds[s]
#            movieName = 'jpg/' + m + '-' + s + '-degPerSec-.jpg'
#            print m, s, speedRadial, movieName
#            for loopN in range(nLoops):
#                if movieType == 'radial':
#                    dotsRadius, dotsTheta, dotsX, dotsY = show_radial_out( halfCycleFrames, dotsRadius0, dotsTheta0, dotsX0, dotsY0, speedRadial, innerRadius, outerRadius, saveMovie )
#                    dotsRadius, dotsTheta, dotsX, dotsY = show_noise( halfCycleFrames, dotsRadius, dotsTheta, dotsX, dotsY, speedRadial, innerRadius, outerRadius, saveMovie )
#                    dotsRadius, dotsTheta, dotsX, dotsY = show_radial_in( halfCycleFrames, dotsRadius, dotsTheta, dotsX, dotsY, speedRadial, innerRadius, outerRadius, saveMovie )
#                    dotsRadius, dotsTheta, dotsX, dotsY = show_noise( halfCycleFrames, dotsRadius, dotsTheta, dotsX, dotsY, speedRadial, innerRadius, outerRadius, saveMovie )
#                elif movieType == 'laminar':
#                    dotsRadius, dotsTheta, dotsX, dotsY = show_laminar( halfCycleFrames, dotsRadius0, dotsTheta0, dotsX0, dotsY0, speedRadial, 0, innerRadius, outerRadius, saveMovie )
#                    dotsRadius, dotsTheta, dotsX, dotsY = show_noise( halfCycleFrames, dotsRadius, dotsTheta, dotsX, dotsY, speedRadial, innerRadius, outerRadius, saveMovie )
#                    dotsRadius, dotsTheta, dotsX, dotsY = show_laminar( halfCycleFrames, dotsRadius, dotsTheta, dotsX, dotsY, speedRadial, math.pi, innerRadius, outerRadius, saveMovie )
#                    dotsRadius, dotsTheta, dotsX, dotsY = show_noise( halfCycleFrames, dotsRadius, dotsTheta, dotsX, dotsY, speedRadial, innerRadius, outerRadius, saveMovie )
#                elif movieType == 'rotation':
#                    dotsRadius, dotsTheta, dotsX, dotsY = show_rotation_fixed_spd( halfCycleFrames, dotsRadius0, dotsTheta0, dotsX0, dotsY0, speedRadial, innerRadius, outerRadius, saveMovie )
#                    dotsRadius, dotsTheta, dotsX, dotsY = show_noise( halfCycleFrames, dotsRadius, dotsTheta, dotsX, dotsY, speedRadial, innerRadius, outerRadius, saveMovie )
#                    dotsRadius, dotsTheta, dotsX, dotsY = show_rotation_fixed_spd( halfCycleFrames, dotsRadius, dotsTheta, dotsX, dotsY, -1.0*speedRadial, innerRadius, outerRadius, saveMovie )
#                    dotsRadius, dotsTheta, dotsX, dotsY = show_noise( halfCycleFrames, dotsRadius, dotsTheta, dotsX, dotsY, speedRadial, innerRadius, outerRadius, saveMovie )
#            # Save movie
#            if saveMovie:
#                win.saveMovieFrames(movieName)
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
# Constants, display parameters

# Modal flags
debugMode = 0
saveMovie = 0
singleFrameMode = 1

# Specify constants, parameters for Fesi, Thomas, Gilmore 2014
windowPix = 600     # Window is square 600 x 600
monitorFramesPerSec = 72
updateFramesPerSec = 24
modulationHz = 1.2
cycleFrames = updateFramesPerSec/modulationHz
halfCycleFrames = int( cycleFrames/2.0 )
stimSpeeds = { '1' : 0.0017,          # 1 deg/s
               '2' : 0.0034,          # 2 deg/s
               '4' : 0.0068,          # 4 deg/s
               '8' : 0.0136,          # 8 deg/s
               '16' : 0.0272          # 16 deg/s
                }
speedType = '4'
dotSpeedUnits = stimSpeeds[speedType]

movieTypes = ['radial', 'laminar', 'rotation']
movieType = movieTypes[1]
movieName = 'jpg/' + movieType + '-' + speedType + '-degPerSec-.jpg'
nMovieLoops = 2

dotLifeFrames = 100
nDots = 2659
dotSize = .005     # 7 arcmin or 7/60=0.117 deg. 0.117 deg/dot / 24 deg/field diam

sans = ['Gill Sans MT', 'Arial','Helvetica','Verdana'] #use the first font found on this list

# annulus parameters
innerRadius = 0.1 # Actually 0.19 units in space, but 0.1 for sqrt replot algorithm
outerRadius = 1.0

# Create a window
win =visual.Window( (windowPix,windowPix), allowGUI=False, color=[-1,-1,-1],
    bitsMode=None, units='norm', winType='pyglet')
    
# Create dot positions, stimuli
refreshIndex = numpy.random.random_integers(0,dotLifeFrames-1,nDots)

dotsRadius0, dotsTheta0, dotsX0, dotsY0 = make_dots( nDots ) # Default is uniform in rectangular region

dots = visual.ElementArrayStim(win, units='norm', elementTex=None, elementMask='circle', fieldShape='circle',
    nElements=nDots, sizes=dotSize, colors=[1,1,1])

fixation = visual.TextStim(win,text="+",pos=(0, 0), color=[1,1,1], ori=0, height = 0.1, font=sans)

blankCircle = visual.GratingStim(win, tex='none',mask='circle',size=0.38,sf=0.0, color=[-1,-1,-1])

# Prepare to loop
trialClock = core.Clock()

#------------------------------------------------------------------------- 
# Testing different stimulus types

# Radial in/out
dotsRadius, dotsTheta, dotsX, dotsY, lastFrame = show_dots(halfCycleFrames, 0, dotsRadius0, dotsTheta0, dotsX0, dotsY0, refreshIndex,
                                        dotSpeedUnits, 0.0, innerRadius, outerRadius, 'radial', 'none', 0, singleFrameMode )
dotsRadius, dotsTheta, dotsX, dotsY, lastFrame = show_dots(halfCycleFrames, lastFrame, dotsRadius, dotsTheta, dotsX, dotsY, refreshIndex,
                                        dotSpeedUnits, 0.0, innerRadius, outerRadius, 'random', 'none', 0, singleFrameMode )
dotsRadius, dotsTheta, dotsX, dotsY, lastFrame = show_dots(halfCycleFrames, lastFrame, dotsRadius0, dotsTheta0, dotsX0, dotsY0, refreshIndex,
                                        -1.0*dotSpeedUnits, 0.0, innerRadius, outerRadius, 'radial', 'none', 0, singleFrameMode )
dotsRadius, dotsTheta, dotsX, dotsY, lastFrame = show_dots(halfCycleFrames, lastFrame, dotsRadius, dotsTheta, dotsX, dotsY, refreshIndex,
                                        dotSpeedUnits, 0.0, innerRadius, outerRadius, 'random', 'none', 0, singleFrameMode )
 

# rotation-eq-spd  left/right  
dotsRadius, dotsTheta, dotsX, dotsY, lastFrame = show_dots(halfCycleFrames, 0, dotsRadius0, dotsTheta0, dotsX0, dotsY0, refreshIndex,
                                        dotSpeedUnits, 0.0, innerRadius, outerRadius, 'rotation-eq-spd', 'none', 0, singleFrameMode )
dotsRadius, dotsTheta, dotsX, dotsY, lastFrame = show_dots(halfCycleFrames, lastFrame, dotsRadius, dotsTheta, dotsX, dotsY, refreshIndex,
                                        dotSpeedUnits, 0.0, innerRadius, outerRadius, 'random', 'none', 0, singleFrameMode )
dotsRadius, dotsTheta, dotsX, dotsY, lastFrame = show_dots(halfCycleFrames, lastFrame, dotsRadius0, dotsTheta0, dotsX0, dotsY0, refreshIndex,
                                        -1.0*dotSpeedUnits, 0.0, innerRadius, outerRadius, 'rotation-eq-spd', 'none', 0, singleFrameMode )
dotsRadius, dotsTheta, dotsX, dotsY, lastFrame = show_dots(halfCycleFrames, lastFrame, dotsRadius, dotsTheta, dotsX, dotsY, refreshIndex,
                                        dotSpeedUnits, 0.0, innerRadius, outerRadius, 'random', 'none', 0, singleFrameMode )

# laminar  left/right  
dotsRadius, dotsTheta, dotsX, dotsY, lastFrame = show_dots(halfCycleFrames, 0, dotsRadius0, dotsTheta0, dotsX0, dotsY0, refreshIndex,
                                        dotSpeedUnits, 0.0, innerRadius, outerRadius, 'laminar', 'wrap', 0, singleFrameMode )
dotsRadius, dotsTheta, dotsX, dotsY, lastFrame = show_dots(halfCycleFrames, lastFrame, dotsRadius, dotsTheta, dotsX, dotsY, refreshIndex,
                                        dotSpeedUnits, 0.0, innerRadius, outerRadius, 'random', 'none', 0, singleFrameMode )
dotsRadius, dotsTheta, dotsX, dotsY, lastFrame = show_dots(halfCycleFrames, lastFrame, dotsRadius0, dotsTheta0, dotsX0, dotsY0, refreshIndex,
                                        -1.0*dotSpeedUnits, 0.0, innerRadius, outerRadius, 'laminar', 'wrap', 0, singleFrameMode )
dotsRadius, dotsTheta, dotsX, dotsY, lastFrame = show_dots(halfCycleFrames, lastFrame, dotsRadius, dotsTheta, dotsX, dotsY, refreshIndex,
                                        dotSpeedUnits, 0.0, innerRadius, outerRadius, 'random', 'none', 0, singleFrameMode )
 