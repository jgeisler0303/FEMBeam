function sid= FASTbld2SID(fst_file, modes, twisted_shape)
if ~exist('modes', 'var')
    modes= 1:3;
end
if ~exist('twisted_shape', 'var')
    twisted_shape= true;
end

p=  FASTCoeff(fst_file);

n_elem= size(p.TwistedSF, 4)-2;
n_q= 3;

R= cell(n_elem+2, 1);
Phi= cell(n_elem, 1);
dPhi= cell(n_elem, 1);
ddPhi= cell(n_elem, 1);
IE= cell(n_elem, 1);
for i= 1:n_elem
    R{i+1}= [0 0 p.RNodes(i)+p.HubRad]';

    if twisted_shape
        Phi{i}= zeros(3, n_q);
        dPhi{i}= zeros(3, n_q);
        ddPhi{i}= zeros(3, n_q);
        for i_q= 1:n_q
            Phi{i}(:, i_q)= [squeeze(p.TwistedSF(1, :, i_q, i+1, 1)) 0];
            dPhi{i}(:, i_q)= [squeeze(p.TwistedSF(1, :, i_q, i+1, 2)) 0];
            ddPhi{i}(:, i_q)= [squeeze(p.TwistedSF(1, :, i_q, i+1, 3)) 0];
        end
    else
        Phi{i}= [p.Shape1(i) p.Shape2(i) 0; 0 0 p.Shape(i); 0 0 0];
        dPhi{i}= [p.dShape1(i) p.dShape2(i) 0; 0 0 p.dShape(i); 0 0 0];
        ddPhi{i}= [p.ddShape1(i) p.ddShape2(i) 0; 0 0 p.ddShape(i); 0 0 0];
    end
    IE{i}= [p.StiffBF(1, i) p.StiffBE(1, i) 0]';

    Phi{i}= Phi{i}(:, modes);
    dPhi{i}= dPhi{i}(:, modes);
    ddPhi{i}= ddPhi{i}(:, modes);
end
R{1}= [0 0 p.HubRad]';
R{end}= [0 0 p.TipRad]';

damping= [p.BldFDamp(1, :) p.BldEDamp(1)]/100;
damping= damping(modes);
% twist= p.ThetaS(1, 2:end-1)/180*pi;
twist= p.ThetaS(1, 2:end-1)*0;

sid= SF2SID(R, Phi, dPhi, ddPhi, p.BElmntMass(:, 1), IE, damping, twist);

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
