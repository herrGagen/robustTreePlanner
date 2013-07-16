function [minX,minY] = findGlobalMin( fun, xLeft, xRight )

yLeft = fun(xLeft);
yRight = fun(xRight);

smallX = .0000001;
deltaX = (xRight - xLeft)/10;
x = xLeft;
der = fun(x+smallX) - fun(x);

minX = x;
minDer = der;
localMinX = minX;
localMinDer = minDer;
%bracket zeros;
while(x < xRight)
        oldX = x;
        oldDer = der;
        x = x+ deltaX;
        der = fun(x+smallX) - fun(x);
        if(oldDer*der < 0) % we have bracketed a zero
            [localMinX, localMinDer] = findZeroOfDeriv(fun,oldX, x);                 
        end
        if(localMinDer < minDer)
            minX = localMinX;
            minDer = localMinDer;
        end
end

minY = fun(minX);
if(yLeft < minY)
    minX = xLeft;
    minY = yLeft;
end

if(yRight < minY)
    minX = xRight;
    minY = yRight;
end

function [x der] = findZeroOfDeriv( fun, xMin, xMax )
% Uses dinary search to find zero of the derivative of the function
% handle fun being passed in.

SMALL_NUMBER = abs(xMax-xMin)/2^10;
smallX = .0000001;

yPrimeMax = fun(xMax+smallX) - fun(xMax);

while( abs(xMax - xMin) > SMALL_NUMBER )
           xTemp = (xMax + xMin)/2;
           yPrimeTemp = fun(xTemp+smallX) - fun(xTemp);
           % if derivative in center is different side than derivative on
           % right, make the center the lew left point
           if(yPrimeTemp*yPrimeMax < 0)
               xMin = xTemp;
           else
               xMax = xTemp;
            end
end
x = xTemp;
der = yPrimeTemp;
