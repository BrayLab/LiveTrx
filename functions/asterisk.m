function[c] = asterisk(p)
if p < 10^(-10); c = '***'; 
else
    if p < 10^(-5); c = '**';
    else
        if p < 10^(-2); c = '*';
        else
            c = 'ns';
        end
    end
end