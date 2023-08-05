ed_blade= FAST2Matlab('/home/jgeisler/Temp/CADynTurb/sim/generated/sim_no_inflow_no_cone/NRELOffshrBsline5MW_Blade.dat', 2);
nnode= 17;

tiledlayout(3, 1)
BlFract= linspace(0, 1, 2*nnode+1);
BlFract= BlFract(2:2:end);

StrcTwst= -1*interp1(ed_blade.BldProp.Table(:, 1), ed_blade.BldProp.Table(:, 3), BlFract);
nexttile
plot(BlFract, StrcTwst)

BldFl1Sh= zeros(6, 1);
BldFl1Amp= zeros(1, nnode);
BldEdgSh= zeros(6, 1);
BldEdgAmp= zeros(1, nnode);
for i= 2:6
    BldFl1Sh(i)= GetFASTPar(ed_blade, sprintf('BldFl1Sh(%d)', i));
    BldFl1Amp= BldFl1Amp + BldFl1Sh(i)*BlFract.^i;
    BldEdgSh(i)= GetFASTPar(ed_blade, sprintf('BldEdgSh(%d)', i));
    BldEdgAmp= BldEdgAmp + BldEdgSh(i)*BlFract.^i;
end

DBldFl1Amp= diff(BldFl1Amp);
DBldEdgAmp= diff(BldEdgAmp);
BldFl1Amp_X= zeros(1, nnode);
BldFl1Amp_Y= zeros(1, nnode);
BldFl1Amp_X(1)= BldFl1Amp(1)*cosd(StrcTwst(1));
BldFl1Amp_Y(1)= BldFl1Amp(1)*sind(StrcTwst(1));
BldEdgAmp_X= zeros(1, nnode);
BldEdgAmp_Y= zeros(1, nnode);
BldEdgAmp_X(1)= -BldEdgAmp(1)*sind(StrcTwst(1));
BldEdgAmp_Y(1)= BldEdgAmp(1)*cosd(StrcTwst(1));
for i= 2:nnode
    BldFl1Amp_X(i)= BldFl1Amp_X(i-1) + DBldFl1Amp(i-1)*cosd(StrcTwst(i));
    BldFl1Amp_Y(i)= BldFl1Amp_Y(i-1) + DBldFl1Amp(i-1)*sind(StrcTwst(i));
    BldEdgAmp_X(i)= BldEdgAmp_X(i-1) - DBldEdgAmp(i-1)*sind(StrcTwst(i));
    BldEdgAmp_Y(i)= BldEdgAmp_Y(i-1) + DBldEdgAmp(i-1)*cosd(StrcTwst(i));
end

nexttile
plot(BlFract, BldFl1Amp, BlFract, BldFl1Amp_X, BlFract, BldFl1Amp_Y, param.R/param.R(end), param.ModalShapes{1})
nexttile
plot(BlFract, BldEdgAmp, BlFract, BldEdgAmp_X, BlFract, BldEdgAmp_Y, param.R/param.R(end), param.ModalShapes{2})


%%
file_F= '/home/jgeisler/Temp/CADynTurb/sim/generated/sim_no_inflow_no_cone/coh_URef-12_maininput.outb';
d_F= loadFAST(file_F);

%%
disp('y')
([d_F.Q_B1F1.Data d_F.Q_B1F2.Data d_F.Q_B1E1.Data]\d_F.TipDyb1.Data)'
([d_F.Q_B1F1.Data d_F.Q_B1E1.Data]\d_F.TipDyb1.Data)'
disp('x')
([d_F.Q_B1F1.Data d_F.Q_B1F2.Data d_F.Q_B1E1.Data]\d_F.TipDxb1.Data)'
([d_F.Q_B1F1.Data d_F.Q_B1E1.Data]\d_F.TipDxb1.Data)'
