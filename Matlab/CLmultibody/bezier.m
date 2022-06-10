function B = bezier(P,t)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
n = length(P)-1;
B = zeros(length(t),1);
for j = 1:length(t)
    T = t(j);
    y = 0;
    for i = 0:n
        y = y + (nchoosek(n,i)*((1-T)^(n-i))*(T^(i))*P(i+1));
    end
    B(j) = y;
end
