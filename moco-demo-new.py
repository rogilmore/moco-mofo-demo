  #!/usr/bin/env python2
""" Motion coherence (MOCO) stimulus display demo
    Simulates displays in Fesi et al. 2014
"""
#-------------------------------------------------------------------------
#   author: Rick Gilmore, thatrickgilmore@gmail.com
#   based on starField, captureFrames demos from PsychoPy2 http://www.psychopy.org/
#
#   Creates an onscreen display and if saveMovie=True,
#   writes it to a series of .jpg images that can be 
#   converted into a movie using MPEG Streamclip or
#   a similar application.
#   
#   repo @ https://github.com/rogilmore/moco-mofo-demo

#-------------------------------------------------------------------------
# Version history
# 2014-07-31    v.01    first working version
# 2014-08-01    v.02    functions for each half cycle
# 2014-08-02    v.03    added uniform speed radial motion, dictionary for speed, movieName
# 2014-08-02    v.04    made generator functions return dot positions for continuity
# 2014-08-04    v.05    refactored with new move_dots, replot_out_dots functions; added refreshing
# 2014-08-08    v.06    further refactoring, streamlining. Fixed weird coordinate conversion bug.
# 2014-09-02    v.07    refactoring to put parameters in params list
# 2014-09-03    v.08    new function defines, checks dependencies in parameters.
# 2014-09-03    v.09    new hybrid distribution for uniform circles/rings. fixes 2014-09-03.1

#-------------------------------------------------------------------------
# Known bugs, desired enhancements
# 2014-09-03.1  Fixed in v.09. linear motion in ring should have special replot mode to keep uniform density.
# 2014-09-03.2  should generate array of motion types for a stimulus cycle, e.g. radial-out, random, radial-in, random.
# 2014-09-03.3  should print params to separate file for each movie generated.
# 2014-09-03.4  should check params for consistency whenever they are changed.

#-------------------------------------------------------------------------
# To use, manipulate the parameters in define_check_params()
# Laminar flow: 'moveMode' = 'laminar'
# Radial flow:  'moveMode' = 'radial'
# Rotational flow: 'moveMode' = 'rotation-eq-spd' for rotational flows with equal local speeds.
# Change speeds: 'degPerSec' in ['1', '2', '4', '8', '16']
# Laminar direction: 'laminarDirRads', default = 0
#
# Typical cycle consists of motion / random / opposite direction motion / random

#-------------------------------------------------------------------------
# Import dependencies
from psychopy import visual, event, core
from psychopy.tools.coordinatetools import pol2cart, cart2pol
import numpy

#-------------------------------------------------------------------------
# define functions

#-------------------------------------------------------------------------
def move_dots( r, th, X, Y, params ):
    if params['moveMode'] == 'radial': 
        r += params['dotSpeedUnits']
        X, Y = pol2cart(th, r, units='rad')
    elif params['moveMode'] == 'rotation-eq-angle':
        th += params['dotSpeedUnits']
        X, Y = pol2cart(th, r, units='rad')
    elif params['moveMode'] == 'rotation-eq-spd':
        X += params['dotSpeedUnits']*numpy.sin( numpy.pi - th )
        Y += params['dotSpeedUnits']*numpy.cos( numpy.pi - th )
        th, r = cart2pol( X, Y, units ='rad' )
    elif params['moveMode'] == 'laminar': 
        X += params['dotSpeedUnits']*numpy.cos( params['laminarDirRads'] )
        Y += params['dotSpeedUnits']*numpy.sin( params['laminarDirRads'] )
        th, r = cart2pol( X, Y, units ='rad' )
    elif params['moveMode'] == 'random':
        dTheta = numpy.random.rand(params['nDots'])*numpy.pi*2.0
        X += params['dotSpeedUnits']*numpy.cos( dTheta )
        Y += params['dotSpeedUnits']*numpy.sin( dTheta )
        th, r = cart2pol( X, Y, units ='rad' )
    # end if moveMode
    
    return r, th, X, Y
# end def move_dots
#-------------------------------------------------------------------------


#-------------------------------------------------------------------------
def replot_out_dots( radius, thetaRads, X, Y, params ):
        
    # Determine who's out of range
    if ( params['displayRegion'] == 'ring' ):
        outField = (radius >= params['outerRadiusUnits'])
        inField = (radius <= 0)
        out = outField + inField
    elif ( params['displayRegion'] == 'rect' ):
        gtMaxX = (X >= params['maxXunits'])
        ltMinX = (X <= params['minXunits'])
        gtMaxY = (Y >= params['maxYunits'])
        ltMinY = (Y >= params['minYunits'])
        outX = ( gtMaxX | ltMinX )
        outY = ( gtMaxY | gtMaxY )
        out = ( outX | outY )
    else:
        print "Invalid displayRegion"
   
    if params['debugMode']:
        print 'Out= ' + str( sum( outField ) ) + '| In= ' + str( sum( inField ) )+ ' | MaxX = ' + str( numpy.max( X ) )  + ' | MinX = ' + str( numpy.min( X ) ) + ' | MinR = ' + str( numpy.min(radius) ) 
        
    # Replot based on replotMode and moveMode
    if params['replotMode'] == 'wrap':
        if params['displayRegion'] == 'ring':
            if (params['moveMode'] == 'radial') or (params['moveMode'] == 'rotation') or (params['moveMode'] == 'random'):
                if out.any:
                    radius = numpy.mod( radius, params['outerRadiusUnits']  )
                X, Y = pol2cart( thetaRads, radius, units='rads' )
            elif (params['moveMode'] == 'laminar'):
                if out.any: # mirror reverse X, Y coordinates
                    X[ out ] = -1.0*X[out]
                    Y[ out ] = -1.0*Y[out]
                thetaRads, radius = cart2pol( X, Y, units='rad' )
    elif params['replotMode'] == 'replot-scaled-polar':
        if outField.any:
                radius[outField], thetaRads[outField], X[outField], Y[outField] = make_dots( sum( outField ), params, min=0., max = params['innerRadiusUnits'], distributionMode = 'scaled-polar' )
        if inField.any:
            radius[inField], thetaRads[inField], X[inField], Y[inField] = make_dots( sum( inField ), params, min=params['outerRadiusUnits'], max=params['outerRadiusUnits']+params['dotSpeedUnits'], distributionMode = 'scaled-polar' )

    elif params['replotMode'] == 'replot-radial':
        if out.any:
            radius[out], thetaRads[out], X[out], Y[out] = make_dots( sum( out ), params, min=0, max=params['outerRadiusUnits'], distributionMode = 'uniform-polar' )

    elif params['replotMode'] == 'replot-rect':
        if out.any:
            radius[out], thetaRads[out], X[out], Y[out] = make_dots( sum( out ), params, min=params['xMinUnits'], max=params['xMaxUnits'], distributionMode='uniform-rect' )

        # Convert changed X,Y to polar
        thetaRads, radius = cart2pol( X, Y, units='rad' )
    # end if replotMode
    
    return radius, thetaRads, X, Y
# end def replot_out_dots
#-------------------------------------------------------------------------


#-------------------------------------------------------------------------
def show_dots( frameIndex0, dotsRadius, dotsTheta, dotsX, dotsY, refreshIndex, params ):
    # cycle is defined as full on/off or motion/random pair
    for frameN in range( params['halfCycleFrames'] ):
        # refresh dots
        frameIndex = frameN + frameIndex0
        refreshThese = ( refreshIndex == frameIndex )
        if sum( refreshThese ):
            if ( params['displayRegion'] == 'ring' ):
                dotsRadius[ refreshThese ], dotsTheta[ refreshThese ], dotsX[ refreshThese ], dotsY[ refreshThese ] = make_dots( sum( refreshThese ), params, distributionMode='scaled-polar' )
            else:
                dotsRadius[ refreshThese ], dotsTheta[ refreshThese ], dotsX[ refreshThese ], dotsY[ refreshThese ] = make_dots( sum( refreshThese ), params, distributionMode='uniform-rect' )
         
        # move dots
        dotsRadius, dotsTheta, dotsX, dotsY = move_dots( dotsRadius, dotsTheta, dotsX, dotsY, params )

        # reposition if out
        dotsRadius, dotsTheta, dotsX, dotsY = replot_out_dots( dotsRadius, dotsTheta, dotsX, dotsY, params )
                
        # Assign dot positions to stimulus
        dots.setXYs(numpy.array( [dotsX, dotsY] ).transpose())
        
        # Draw stimuli to window 
        dots.draw()
        
        if params['maskCenter']:
            blankCircle.draw()
        
        if params['showFixation']:
            fixation.draw() 
        
        # Save movie
        if params['saveMovie']:
            win.getMovieFrame(buffer='back')  
            
        # flip buffer to display unless single frame mode
        if params['singleFrameMode']:
            msg =  visual.TextStim(win, text=str(frameN), pos=[0.9,-0.9], color=[1,1,1])
            msg.draw()
            k = ['']
            while k[0] not in ['space', 'esc']:   
                k = event.waitKeys()

        win.flip()
        
    # end for frameN
    return dotsRadius, dotsTheta, dotsX, dotsY, frameN # for next round
# end def show_dots
#-------------------------------------------------------------------------


#-------------------------------------------------------------------------
def make_dots( nDots, params, min=-1, max=1, distributionMode='uniform-rect' ):
    if distributionMode == 'scaled-polar':      # Adjusts for density change due to outward radial motion
        r = numpy.random.rand(nDots)**(0.5)*(max - min) + min
        th = numpy.random.rand(nDots)*(2*numpy.pi)-numpy.pi
        X, Y = pol2cart( th, r, units='rad' )
    elif distributionMode == 'uniform-polar':
        r = numpy.random.rand(nDots)*(max - min) + min
        th = numpy.random.rand(nDots)*(2*numpy.pi)-numpy.pi
        X, Y = pol2cart( th, r, units='rad' )
    elif distributionMode == 'fixed-circle':
        th = numpy.random.rand(nDots)*(2*numpy.pi)-numpy.pi
        r = numpy.ones(nDots)*.95
        X, Y = pol2cart( th, r, units='rad' )
    elif distributionMode == 'hybrid-uniform-rect-polar': # dots of uniform density in circular region
        X, Y = numpy.random.rand(nDots)*(max-min)+min, numpy.random.rand(nDots)*(max-min)+min
        th, r = cart2pol( X, Y, units = 'rad' )
        outR = ( r >= params['outerRadiusUnits'] )
        r[ outR ] = numpy.random.rand(sum(outR))*params['outerRadiusUnits']
        X, Y = pol2cart( th, r, units='rad' )
    else : # default is random in [-1,1] and [-1,1] in Cartesian coordinates
        X, Y = numpy.random.rand( nDots )*(max - min) + min, numpy.random.rand( nDots )*(max - min) + min
        th, r = cart2pol( X, Y, units = 'rad' )
        
    return r, th, X, Y
#-------------------------------------------------------------------------


#-------------------------------------------------------------------------
def define_check_params(params):
# Prespecified params
    params = { 'debugMode': False,
           'saveMovie': False,
           'singleFrameMode': False,
           'windowPix': 600,
           'monitorFramesPerSec': 72,
           'updateFramesPerSec': 24,
           'modulationHz': 1.2,
           'nMovieLoops': 1,
           'innerRadiusDeg': 2.4,
           'outerRadiusDeg': 12.,
           'dotLifeFrames': 100,
           'dotDensity': .1,
           'dotAmin' : 7.,
           'displayRegionTypes' : ['ring','rect'],
           'displayRegionIndex' : 0,
           'maskCenter': True,
           'showFixation': True,
           'minTh': -1.0*numpy.pi,
           'maxTh': numpy.pi,
           'minXunits': -1.,
           'maxXunits': 1.,
           'minYunits': -1.,
           'maxYunits': 1.,
           'stimSpeeds' : { '1' : 0.0017,          # 1 deg/s
                            '2' : 0.0034,          # 2 deg/s
                            '4' : 0.0068,          # 4 deg/s
                            '8' : 0.0136,          # 8 deg/s
                            '16' : 0.0272          # 16 deg/s
                          },
            'degPerSec': '8',
            'laminarDirRads': 0,
            'movieTypes': ['radial', 'rotation', 'laminar'],
            'movieTypeIndex': 2,
            'moveModes' : ['radial', 'rotation-eq-spd', 'rotation-eq-angle', 'laminar', 'random'],
            'moveModeIndex': 1,
            'replotModes': ['wrap', 'replot-scaled-polar', 'replot-radial', 'replot-rect'],
            'replotModeIndex': 0, # wrap is default
            'distributionModes': ['uniform-rect', 'scaled-polar', 'uniform-polar', 'fixed-circle', 'hybrid-uniform-rect-polar'],
            'distributionModeIndex': 1
         }

    # Computed parameters
    params['cycleFrames'] = params['updateFramesPerSec']/params['modulationHz']
    params['halfCycleFrames'] = int( params['cycleFrames'] / 2.0 )
    
    params['innerRadiusUnits'] = params['innerRadiusDeg']/params['outerRadiusDeg'] # Normalized to outerRadius
    params['outerRadiusUnits'] = 1.0
    
    params['ringAreaDeg2'] = numpy.pi*params['outerRadiusDeg']**2
    params['ringAreaUnit2'] = numpy.pi*params['outerRadiusUnits']**2
    
    params['rectAreaUnit2'] = (params['maxXunits']-params['minXunits'])*(params['maxYunits']-params['minYunits'])
    rectArea2ringArea = params['rectAreaUnit2'] / params['ringAreaUnit2']
    params['rectAreaDeg2'] = rectArea2ringArea * params['ringAreaDeg2']
    
    params['dotDeg'] = params['dotAmin'] / 60.
    params['dotSizeUnits'] = params['dotDeg']/( params['outerRadiusDeg']/params['outerRadiusUnits'])
    dotAreaDeg2 = (params['dotDeg'])**2
    
    params['maxDotsRing'] = params['ringAreaDeg2'] / dotAreaDeg2
    params['maxDotsRect'] = params['rectAreaDeg2'] / dotAreaDeg2
    
    params['dotSpeedUnits'] = params['stimSpeeds'][params['degPerSec']]
 
    params['replotMode'] = params['replotModes'][ params['replotModeIndex'] ]
    
    params['moveMode'] = params['moveModes'][params['moveModeIndex']]

    params['displayRegion'] = params['displayRegionTypes'][params['displayRegionIndex']]

    # Check dependencies

    if params['moveMode'] == 'radial':
        params['movieTypeIndex'] = 0
    elif params['moveMode'] in ['rotation-eq-spd', 'rotation-eq-angle']:
        params['movieTypeIndex'] = 1
    elif params['moveMode'] == 'laminar':
        params['movieTypeIndex'] = 2
    params['movieName'] = 'jpg/' + params['movieTypes'][params['movieTypeIndex']] + '-' + params['degPerSec'] + '-degPerSec-.jpg'

    if params['displayRegion'] == 'ring':
        params['nDots'] = int( params['dotDensity'] * params['maxDotsRing'] )
        if params['moveMode'] == 'laminar':
            params['distributionModeIndex'] = 4 # 'hybrid-uniform-rect-polar'
            params['dotDistributionMode'] = params['distributionModes'][ params['distributionModeIndex'] ]

            params['replotModeIndex'] = 0 # wrap
            params['replotMode'] = params['replotModes'][ params['replotModeIndex'] ]
        elif params['moveMode'] == 'random':
            params['distributionModeIndex'] = 0 # 'scaled-polar'?
            params['dotDistributionMode'] = params['distributionModes'][ params['distributionModeIndex'] ]
            
            params['replotModeIndex'] = 0 # wrap
            params['replotMode'] = params['replotModes'][ params['replotModeIndex'] ]
        elif params['moveMode'] == 'radial':
            params['distributionModeIndex'] = 1 # 'scaled-polar'
            params['dotDistributionMode'] = params['distributionModes'][ params['distributionModeIndex'] ]
            
            params['replotModeIndex'] = 0 # wrap
            params['replotMode'] = params['replotModes'][ params['replotModeIndex'] ]
        else : # rotation
            params['distributionModeIndex'] = 0 # 'scaled-polar'
            params['dotDistributionMode'] = params['distributionModes'][ params['distributionModeIndex'] ]
            
            params['replotModeIndex'] = 2 # replot-radial
            params['replotMode'] = params['replotModes'][ params['replotModeIndex'] ]

    else : # default to 'rect'
        params['nDots'] = int( params['dotDensity'] * params['maxDotsRect'] )
        if params['moveMode'] == 'linear':
            params['distributionModeIndex'] = 0 # 'uniform-rect'
            params['dotDistributionMode'] = params['distributionModes'][ params['distributionModeIndex'] ]
            
            params['replotModeIndex'] = 0 # wrap
            params['replotMode'] = params['replotModes'][ params['replotModeIndex'] ]
        elif params['moveMode'] == 'random':
            params['distributionModeIndex'] = 0 # 'uniform-rect'
            params['dotDistributionMode'] = params['distributionModes'][ params['distributionModeIndex'] ]
            
            params['replotModeIndex'] = 0 # wrap
            params['replotMode'] = params['replotModes'][ params['replotModeIndex'] ]
        elif params['moveMode'] == 'radial':
            params['distributionModeIndex'] = 0 # 'uniform-rect'
            params['dotDistributionMode'] = params['distributionModes'][ params['distributionModeIndex'] ]
            
            params['replotModeIndex'] = 0 # 'replot-scaled-polar'
            params['replotMode'] = params['replotModes'][ params['replotModeIndex'] ]
        else :
            params['distributionModeIndex'] = 0 # 'uniform-rect'
            params['dotDistributionMode'] = params['distributionModes'][ params['distributionModeIndex'] ]
            
            params['replotModeIndex'] = 2 # replot-radial
            params['replotMode'] = params['replotModes'][ params['replotModeIndex'] ]

    return params
#-------------------------------------------------------------------------


#-------------------------------------------------------------------------
# main program
params = {}
params = define_check_params(params)
   
# fonts for messages
sans = ['Gill Sans MT', 'Arial','Helvetica','Verdana'] #use the first font found on this list

#-------------------------------------------------------------------------
# Begin body of program

# Create a window
win = visual.Window( ( params['windowPix'],params['windowPix'] ), allowGUI=False, color=[-1,-1,-1], bitsMode=None, units='norm', winType='pyglet')
    
# Create dot positions, stimuli
refreshIndex = numpy.random.random_integers( 0, params['dotLifeFrames']-1, params['nDots'] )

dotsRadius0, dotsTheta0, dotsX0, dotsY0 = make_dots( params['nDots'], params, distributionMode=params['dotDistributionMode'] )
dotsRadius0, dotsTheta0, dotsX0, dotsY0 = replot_out_dots( dotsRadius0, dotsTheta0, dotsX0, dotsY0, params )

# Make dot array, fixation, other stims
dots = visual.ElementArrayStim(win, units='norm', elementTex=None, elementMask='circle', fieldShape='square', 
                                nElements=params['nDots'], sizes=params['dotSizeUnits'], colors=[1,1,1])

fixation = visual.TextStim(win,text="+",pos=(0, 0), color=[1,1,1], ori=0, height = 0.1, font=sans)

blankCircle = visual.GratingStim(win, tex='none',mask='circle',size=0.38,sf=0.0, color=[-1,-1,-1])

# Prepare to loop
trialClock = core.Clock()

# Initialize before loop
dotsRadius, dotsTheta, dotsX, dotsY = dotsRadius0, dotsTheta0, dotsX0, dotsY0
lf = 0 # lastFrame index to pass between stimulus segments

for l in range( params['nMovieLoops'] ):

    thisMode = params['moveMode']
    dotsRadius, dotsTheta, dotsX, dotsY, lf = show_dots(lf, dotsRadius, dotsTheta, dotsX, dotsY, refreshIndex,
                                            params )

    params['moveMode'] = 'random'
    dotsRadius, dotsTheta, dotsX, dotsY, lf = show_dots(lf, dotsRadius, dotsTheta, dotsX, dotsY, refreshIndex,
                                            params )
                                            
    params['dotSpeedUnits'] = -1.0*params['dotSpeedUnits']
    params['moveMode'] = thisMode
    dotsRadius, dotsTheta, dotsX, dotsY, lf = show_dots(lf, dotsRadius, dotsTheta, dotsX, dotsY, refreshIndex,
                                            params )

    params['moveMode'] = 'random'
    dotsRadius, dotsTheta, dotsX, dotsY, lf = show_dots(lf, dotsRadius, dotsTheta, dotsX, dotsY, refreshIndex,
                                            params )
# Save movie
if params['saveMovie']:
    print "Saving movie frames to jpg/%s as set of jpg files." % movieName
    win.saveMovieFrames(params['movieName'])
#-------------------------------------------------------------------------
 