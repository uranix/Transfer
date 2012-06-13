algebraic3d
solid box = orthobrick(-1, -1, -1; 1, 1, 1);
solid hot = orthobrick(-1, -1, -1; 0, 1, 1);
solid cold = orthobrick(0, -1, -1; 1, 1, 1);

tlo hot -col=[1,0,0];
tlo cold -col=[0,0,1] -transparent;
