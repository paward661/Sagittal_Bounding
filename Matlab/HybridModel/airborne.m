function dydt = airborne(t,y,u)

g = u(1);

dydt = zeros(6,1);
dydt(1) = y(2);
dydt(2) = -g;
dydt(3) = y(4);
dydt(4) = 0;
dydt(5) = y(6);
dydt(6) = 0;
end

