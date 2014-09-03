#!/usr/bin/env python2

import numpy

class Dot:
    def __init__(self, r=0, th=0, x=0, y=0, defaultCoord='Cartesian', defaultGenerateMode='randomXY'):
        # Define dots
        if defaultGenerateMode == 'randomXY':
            self.randomXY()
        else:
            self.r = r
            self.th = th
            self.x = x
            self.y = y
        
        self.defaultCoord = defaultCoord
        self.defaultGenerateMode = defaultGenerateMode

        if (r, th != self.cart2pol()):
            if defaultCoord == 'Cartesian':
                 self.r, self.th = self.cart2pol()
            else:
                self.x, self.y = self.pol2cart()

    def __repr__( self ):   
        return 'Dot(r=%s, th=%s, x=%s, y=%s, defaultCoord="%s")' % ( str(self.r), str(self.th), str(self.x), str(self.y), self.defaultCoord )

    def __str__( self ):
        return '[Value: r=%s, th=%s, x=%s, y=%s]' % ( str(self.r), str(self.th), str(self.x), str(self.y) )

    def cart2pol( self ):
        r = numpy.sqrt( self.x**2 + self.y**2 )
        th = numpy.arctan2( self.y, self.x )
        return r, th

    def pol2cart( self ):
        x = self.r*numpy.cos( self.th )
        y = self.r*numpy.sin( self.th )
        return x, y

    def chgR( self, dR=0):
       self.r += dR
       self.x, self.y = self.pol2cart()
   
    def chgTh( self, dTh=0):
       self.th += dTh
       self.x, self.y = self.pol2cart()
       
    def chgXY( self, dX=0, dY=0 ):
       self.x += dX
       self.y += dY
       self.r, self.th = self.cart2pol()

    def randomXY( self, xMin=-1, xMax=1, yMin=-1, yMax=1 ):
        self.x = numpy.random.random()*(xMax-xMin)+xMin
        self.y = numpy.random.random()*(yMax-yMin)+yMin
        self.r, self.th = self.cart2pol()

    def __gt__( self, other ):
        return self.r > other

    def __lt__( self, other ):
        return( self.r < other )

# class Dots:
#     def __init__(self, nDots, defaultCoord='Cartesian', defaultGenerateMode='randomXY'):
#         some_dots = []
#         for d in range( nDots ):
#             some_dots.append( Dot( defaultCoord=defaultCoord, defaultGenerateMode=defaultGenerateMode ) )

#         self = some_dots

    # def __getitem__( self, index ):
    #     return self[index]

    # def __repr__(self):
    #     for d in range( nDots ): # Not correct, need to get len( self ), but let's use it for now
    #         print( self[d] )

# Console tests
if __name__ == '__main__':
    zeroDot = Dot(defaultGenerateMode='fixed')
    cartDot = Dot( x = 3, y = 4, defaultGenerateMode='fixed' )
    polarDot = Dot( r = 1, th = numpy.pi/4, defaultGenerateMode='fixed' )
    incompatibleDot = Dot( r=.5, th=.2, x=1, y=2, defaultGenerateMode='fixed' ) # Cartesian and polar values aren't commensurate
    anotherIncompatibleDot = Dot( r=.5, th=.2, x=1, y=2, defaultCoord='polar', defaultGenerateMode='fixed' )
    
    print( zeroDot )
    print( cartDot )
    print( polarDot )
    print( incompatibleDot )
    print( anotherIncompatibleDot )

    print zeroDot.cart2pol()
    print cartDot.cart2pol()

    print 'Let us change zeroDot by dX=1, dY=2'
    print 'Before change:'
    print( zeroDot )
    zeroDot.chgXY( dX=1, dY=2 )
    print 'After:'
    print( zeroDot )

    print 'Change polarDot by dR=1.2'
    print 'Before change:' 
    print( polarDot )   
    polarDot.chgR( dR = 1.2 )
    print 'After:'
    print( polarDot )

    zeroDot.randomXY()
    print( zeroDot )

    # Testing dot array

    my_dots = []
    nDots = 10
    for d in range( nDots ):
        my_dots.append( Dot() )

    # print a range
    print( my_dots[5:7] )

    # change r for 1
    print( my_dots[2].r )
    my_dots[2].chgR(10)
    print( my_dots[2].r )

    for d in range( len( my_dots ) ):
        my_dots[d].chgR(10)

    print( my_dots[:] )
    
    # change r for a range
    # print( my_dots[4].r, my_dots[5].r )
    # my_dots[4:5].chgR(10) # fails. Can't iterate method over list. Need __iter_ method?
    # print( my_dots[4].r, my_dots[5].r )

    # # Try this
    # my_dots = Dots( 4 )
    # print( my_dots )


    
       
     