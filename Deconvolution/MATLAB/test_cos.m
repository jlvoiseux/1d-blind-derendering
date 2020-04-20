function out = test_cos(input)
    out = (pi/2 -input) - ((pi/2 -input).^3)./ 6 + ((pi/2 -input).^5)./ 120 - ((pi/2 -input).^7)./ 5040 + ((pi/2 -input).^9)./ 362880;    
end