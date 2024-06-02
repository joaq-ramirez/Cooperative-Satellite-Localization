function dxdt = NLODE(t,x,mu)

Xpos = x(1);
Xvel = x(2);
Ypos = x(3);
Yvel = x(4);

r = norm([Xpos Ypos]);

Xacc = -mu*Xpos/r^3;
Yacc = -mu*Ypos/r^3;

dxdt = [Xvel;Xacc;Yvel;Yacc];

end