function Gravity = GravityModel(Position)


g_para1 = 0.0818191908426^2;
g_para2 = 9.7803267714;
g_para3 = 0.00193185138639;


[Gravity] = -Earth_Gravity(Position,g_para1,g_para2,g_para3);

