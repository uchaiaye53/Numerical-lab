function ret = rk(f, xrange, y0, h)
    
    x = [xrange(1):h:xrange(2)]; 
    y(1) = y0;
    
    for i=2:length(x)
      
        k1 = h*f(x(i-1),y(i-1));
        k2 = h*f( x(i-1)+(h/2), y(i-1)+(k1/2) );
        k3 = h*f( x(i-1)+(h/2), y(i-1)+(k2/2) );
        k4 = h*f( x(i-1)+h, y(i-1)+(k3/2) );
        
        % next y
        y(i) = y(i-1) + (1/6)*(k1+(2*k2)+(2*k3)+k4);
    end
    
    
    ret = [x', y']; 
    return; 
endfunction