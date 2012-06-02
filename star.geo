algebraic3d
solid box = orthobrick(-1, -1, -1; 1, 1, 1);
solid core = sphere(0, 0, 0; 0.01);

solid air = box and not core;

tlo core -col=[1,0,0];
tlo air -col=[0,0,1] -transparent;
