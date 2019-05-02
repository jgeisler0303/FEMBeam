% Beispiel 6.5 S. 372

data(1).R= 0;
data(1).A= 0.0006;
data(1).Iz= 4.5e-10;
data(1).Iy= 1e-4;
data(1).Ip= 1e-4;
data(1).rho= 8400;
data(1).E= 7e10;
data(1).G= 7e10;


data(2).R= 1;
data(2).A= 0.0006;
data(2).Iz= 4.5e-10;
data(2).Iy= 1e-4;
data(2).Ip= 1e-4;
data(2).rho= 8400;
data(2).E= 7e10;
data(2).G= 7e10;

data(3).R= 2;
data(3).A= 0.0006;
data(3).Iz= 4.5e-10;
data(3).Iy= 1e-4;
data(3).Ip= 1e-4;
data(3).rho= 8400;
data(3).E= 7e10;
data(3).G= 7e10;

sid= FEMBeam2SID(data, 2, [6, 1]);
%  sid= FEMBeam2SID(data, 2);

idx2D= [1 2 6 [1 2 6]+6 [1 2 6]+12];
idx2D_= zeros(1, 18);
idx2D_(idx2D)= 1;
idx2D_([1 2 3 4 18-6+2 18-6+3])= [];
idx2D_= logical(idx2D_);

write_sid(sid, 'test.SID_FEM');
sid= read_sid('test.SID_FEM');
% write_sid(sid, 'test.SID_FEM');
