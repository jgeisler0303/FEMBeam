function sid= FASTtwr2SID(fst_file, modes)
if ~exist('modes', 'var')
    modes= 1:4;
end

p=  FASTCoeff(fst_file);

n_elem= size(p.TwrFASF, 2)-2;

R= cell(n_elem+2, 1);
Phi= cell(n_elem, 1);
dPhi= cell(n_elem, 1);
ddPhi= cell(n_elem, 1);
IE= cell(n_elem, 1);
for i= 1:n_elem
    R{i+1}= [0 0 p.HNodes(i)]';
    Phi{i}= [p.TwrFASF(1, i+1, 1) p.TwrFASF(2, i+1, 1) 0 0; 0 0 p.TwrSSSF(1, i+1, 1) p.TwrSSSF(2, i+1, 1); 0 0 0 0];
    dPhi{i}= [p.TwrFASF(1, i+1, 2) p.TwrFASF(2, i+1, 2) 0 0; 0 0 p.TwrSSSF(1, i+1, 2) p.TwrSSSF(2, i+1, 2); 0 0 0 0];
    ddPhi{i}= [p.TwrFASF(1, i+1, 3) p.TwrFASF(2, i+1, 3) 0 0; 0 0 p.TwrSSSF(1, i+1, 3) p.TwrSSSF(2, i+1, 3); 0 0 0 0];
    IE{i}= [p.StiffTFA(i) p.StiffTSS(i) 0]';

    Phi{i}= Phi{i}(:, modes);
    dPhi{i}= dPhi{i}(:, modes);
    ddPhi{i}= ddPhi{i}(:, modes);
end

R{1}= [0 0 0]';
R{end}= [0 0 p.TwrFlexL]';

damping= [p.TwrFADmp(1) p.TwrFADmp(2) p.TwrSSDmp(1) p.TwrSSDmp(2)]/100;
damping= damping(modes);

sid= SF2SID(R, Phi, dPhi, ddPhi, p.TElmntMass, IE, damping);

% remove some off-diagonal terms to conform to FAST approach
sid.Me.M0= diag(diag(sid.Me.M0));
for i_elem= 1:n_elem
    for a= 1:3
        sid.frame(i_elem).Phi.M1(a, :, :)= diag(diag(squeeze(sid.frame(i_elem).Phi.M1(a, :, :))));
    end
end
for a= 1:3
    sid.Ct.M1(:, a, :)= diag(diag(squeeze(sid.Ct.M1(:, a, :))));
end
for m= 1:6
    if m<4
        sid.Oe.M1(:, m, :)= squeeze(sid.Oe.M1(:, m, :)) - sid.K0omega{m, m} + diag(diag(sid.K0omega{m, m}));
    else
        c= m-3;
        d= mod(c, 3)+1;
        sid.Oe.M1(:, m, :)= squeeze(sid.Oe.M1(:, m, :)) - sid.K0omega{c, d} - sid.K0omega{c, d}' + 2*diag(diag(sid.K0omega{c, d}));
    end        
end
