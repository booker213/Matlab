function y = reimann(ul,ur)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%for burgers, left hand side has to be quicker that the right hand side for
%a shockwave

s = (ul+ur)/2;

if ul < ur %rarefraction
    
    if ul > 0 
        y = ul;
    elseif ur<0
            y = ur;
            
    else y=0
    end

elseif ul> ur
    if s > 0
        y = ul^2/2;
    elseif s<0
        y =ur^2/2;
    end
elseif ul==ur
y=ul;

end
end

