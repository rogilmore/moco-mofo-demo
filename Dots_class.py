import numpy

class Dots:
    def __init__( self, r=0, th=0, x=0, y=0, n=1 ):
        if ( type( r ) == type( th ) == type( x ) == type( y ) == type( n ) == int ):
            self.r = r
            self.th = th
            self.x = x
            self.y = y
            self.n = n
        elif ( len( r ) == len( th ) == len( x ) == len( y ) == n ) and ( n > 1 ):
            self.r = r
            self.th = th
            self.x = x
            self.y = y
            self.n = n
        elif ( len( r ) == len( th ) ) and ( r == th == 0 ): # if no r, th specified, generate based on random x, y
            if len( x ) == len( y ) == n:
                self.x, self.y = numpy.rand( n ), numpy.rand( n )
                self.r = numpy.sqrt( self.x**2 + self.y**2 )
                self.th = numpy.arctan2( self.x, self.y )
            else:
                print 'Arguments of unequal length.'
        elif ( len( x ) == len( y ) ) and ( x == y == 0 ): # no x, y specified
            if len( r ) == len( th ) == n:
                self.x, self.y = abs( r )*( numpy.cos( th ), numpy.sin( th ) )
                self.r = r
                self.th = th
            else:
                print 'Arguments of unequal length.'

#    def chgR( self, dR ):
#        return self.r += dR
#    
#    def chgTh( self, dTh):
#        return self.th += dTh
#        
#    def chgXY( self, dX=0, dY=0 ):
#        return self.x += dX, self.y += dY
#                

# Console tests
if __name__ == '__main__':
    oneDot = Dots()
    threeDots = Dots( n=3 )
#    someDots = Dots( r=[.1, -.1, .25], th=[.25, .50, .75] )
    
    print oneDot
    print threeDots
#    print someDots
    
       
     