function [minX,minY] = findGlobalMin( fun, xMin, xMax )

smallX = .0000001;
deltaX = .001;
x = xMin;
y = fun(x+smallX) - fun(x);

minX = x;
minY = y;
%bracket zeos;
while(x < xMax)
        oldX = x;
        oldY = y;
        x = x+ deltaX;
        y = fun(x+smallX) - fun(x);
        if(oldY*y < 0) % we have bracketed a zero
                  [localMinX localMinY] = findZeroOfDeriv(fun,oldX, x);
        end
        if(localMinY < minY)
                minX = localMinX;
                minY = localMinY;
        end
end


function [x,y] = findZeroOfDeriv( fun, xMin, xMax )

smallX = .0000001;
yPrimeMin = fun(xMin+smallX) - fun(xMin);
yPrimeMax = fun(xMax+smallX) - fun(xMax);

while(xMax - xMin > SMALL_NUMBER )
           xTemp = (xMax + xMin)/2;
           yTemp = fun(yTemp+smallX) - fun(yTemp);
           if(yTemp*yMax < 0)
               xMax = xTemp;
           else
               xMin = xTemp;
            end
end
x = xTemp;
y = yTemp;