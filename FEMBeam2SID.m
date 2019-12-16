%  Copyright 2019, Jens Geisler
%  
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%  
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%  
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <http://www.gnu.org/licenses/>.

% Equation and page numbers refer to
% Schwertassek, R & Wallrapp, O. (1999). Dynamik flexibler Mehrk�rpersysteme.

function sid= FEMBeam2SID(data, ref_system, modes, normalize)
ne= length(data)-1;
nk= ne+1;
nqe= 12; % same for all elements
nqk= nqe/2;
nF= nk*nqk;
nq_= nF-nqk;


if ~isfield(data, 'x')
    if isfield(data, 'R')
        for i= 1:nk
            data(i).x= [data(i).R, 0, 0]';
        end
    elseif isfield(data, 'Rx')
        for i= 1:nk
            data(i).x= [data(i).Rx, 0, 0]';
        end            
    elseif isfield(data, 'Ry')
        for i= 1:nk
            data(i).x= [0 data(i).Ry, 0]';
        end            
    elseif isfield(data, 'Rz')
        for i= 1:nk
            data(i).x= [0, 0, data(i).Rz]';
        end            
    end
end
% calculate the tangent (4.24)
if ~isfield(data, 'dx')
    for i= 1:ne
        data(i).dx= data(i+1).x - data(i).x;
    end
end        

for i= 1:ne
    dx= data(i).dx;
    % calculate the tangent (4.24)
    data(i).l= norm(dx);
    data(i).t= dx/data(i).l;
    if i==1
        t_last= data(1).t;
        % start for the normal
        % for t pointing in x or z direction, h will be in y direction
        % for t pointing in y direction, h will be in -x direction
        h_last= null(t_last');
        h_last= h_last(:, 1);
    end

    % calculate normal vector (4.33)
    dt= data(i).t - t_last;
    if norm(dt)<eps
        data(i).h= h_last;
    else
        data(i).h= dt/norm(dt);
    end

    % add twist
    if isfield(data, 'phi')
        % rotation about t
        R= cos(data(i).phi)*eye(3) + sin(data(i).phi)*crossmat(data(i).t) + (1-cos(data(i).phi))*data(i).t*data(i).t';
        data(i).h= R * data(i).h;
    end
    
    % calculate binormal vector (4.34)
    data(i).b= cross(data(i).t, data(i).h);

    % element orientation (4.36)
    data(i).D= [data(i).t data(i).h data(i).b]';

    t_last= data(i).t;
    h_last= data(i).h;
end

[~, extentIdx]= max(sum([data.x], 2));

T= cell(ne, 1);
Gamma= cell(ne, 1);
for e= 1:ne
    T{e}= zeros(nqe, nF);
    Gamma{e}= data(e).D;
    % (5.123) S. 200
    % not sure which is the right version
%     T{e}(1:nqk, (e-1)*nqk+(1:nqk))= blkdiag(data(e).D, eye(3));
%     T{e}(nqk+1:nqe, e*nqk+(1:nqk))= blkdiag(data(e).D, eye(3));
    T{e}(1:nqk, (e-1)*nqk+(1:nqk))= blkdiag(data(e).D, data(e).D);
    T{e}(nqk+1:nqe, e*nqk+(1:nqk))= blkdiag(data(e).D, data(e).D);
end

M= cell(ne, 1);
for e= 1:ne
    le= data(e).l;
    if isfield(data, 'A') && isfield(data, 'rho')
        me01= getElementValue(data, e, 'A') * le * data(e).rho;
    else
        if isfield(data, 'mu')
            me01= getElementValue(data, e, 'mu') * le;
        else
            me01= getElementValue(data, e, 'me');
        end
    end
    me0= me01(1);
    me1= me01(2);
    M{e}= [[(me1+3*me0)/12,0,0,0,0,0,(me1+me0)/12,0,0,0,0,0];[0,(3*me1+10*me0)/35,0,0,0,(7*le*me1+15*le*me0)/420,0,(9*me1+9*me0)/140,0,0,0,-(6*le*me1+7*le*me0)/420];[0,0,(3*me1+10*me0)/35,0,-(7*le*me1+15*le*me0)/420,0,0,0,(9*me1+9*me0)/140,0,(6*le*me1+7*le*me0)/420,0];[0,0,0,(me1+3*me0)/12,0,0,0,0,0,(me1+me0)/12,0,0];[0,0,-(7*le*me1+15*le*me0)/420,0,(3*le^2*me1+5*le^2*me0)/840,0,0,0,-(7*le*me1+6*le*me0)/420,0,-(le^2*me1+le^2*me0)/280,0];[0,(7*le*me1+15*le*me0)/420,0,0,0,(3*le^2*me1+5*le^2*me0)/840,0,(7*le*me1+6*le*me0)/420,0,0,0,-(le^2*me1+le^2*me0)/280];[(me1+me0)/12,0,0,0,0,0,(3*me1+me0)/12,0,0,0,0,0];[0,(9*me1+9*me0)/140,0,0,0,(7*le*me1+6*le*me0)/420,0,(10*me1+3*me0)/35,0,0,0,-(15*le*me1+7*le*me0)/420];[0,0,(9*me1+9*me0)/140,0,-(7*le*me1+6*le*me0)/420,0,0,0,(10*me1+3*me0)/35,0,(15*le*me1+7*le*me0)/420,0];[0,0,0,(me1+me0)/12,0,0,0,0,0,(3*me1+me0)/12,0,0];[0,0,(6*le*me1+7*le*me0)/420,0,-(le^2*me1+le^2*me0)/280,0,0,0,(15*le*me1+7*le*me0)/420,0,(5*le^2*me1+3*le^2*me0)/840,0];[0,-(6*le*me1+7*le*me0)/420,0,0,0,-(le^2*me1+le^2*me0)/280,0,-(15*le*me1+7*le*me0)/420,0,0,0,(5*le^2*me1+3*le^2*me0)/840]];
    
    if isfield(data, 'Mlumped')
        M{e}(1, 1)= M{e}(1, 1) + data(e).Mlumped;
        M{e}(2, 2)= M{e}(2, 2) + data(e).Mlumped;
        M{e}(3, 3)= M{e}(3, 3) + data(e).Mlumped;
        if e==ne
            M{e}(1+6, 1+6)= M{e}(1+6, 1+6) + data(e+1).Mlumped;
            M{e}(2+6, 2+6)= M{e}(2+6, 2+6) + data(e+1).Mlumped;
            M{e}(3+6, 3+6)= M{e}(3+6, 3+6) + data(e+1).Mlumped;
        end
    end
    if isfield(data, 'Ilumped')
        M{e}(4:6, 4:6)= M{e}(4:6, 4:6) + data(e).Ilumped;
        if e==ne
            M{e}(4:6+6, 4:6+6)= M{e}(4:6+6, 4:6+6) + data(e+1).Ilumped;
        end
    end
end

MF= zeros(nF);
for e= 1:ne
    MF= MF + T{e}' * M{e} * T{e};
end

K= cell(ne, 1);
for e= 1:ne
    le= data(e).l;
    if isfield(data, 'E')
        E= data(e).E;
    else
        E= 214e9;
    end
    if isfield(data, 'G')
        G= data(e).G;
    else
        G= E/2/(1+0.3);
    end
    if isfield(data, 'A')
        A01= getElementValue(data, e, 'A');
    else
        if isfield(data, 'EA')
            A01= getElementValue(data, e, 'EA') / E;
        else
            A01= [100, 100]; % very stiff
        end
    end
    A0= A01(1);
    A1= A01(2);

    if isfield(data, 'Iy')
        Iy01= getElementValue(data, e, 'Iy');
    else 
        if isfield(data, 'EIy')
            Iy01= getElementValue(data, e, 'EIy') / E;
        else
            Iy01= [100 100]; % very stiff
        end
    end
    Iy0= Iy01(1);
    Iy1= Iy01(2);

    if isfield(data, 'Iz')
        Iz01= getElementValue(data, e, 'Iz');
    else 
        if isfield(data, 'EIz')
            Iz01= getElementValue(data, e, 'EIz') / E;
        else
            Iz01= [100 100]; % very stiff
        end
    end
    Iz0= Iz01(1);
    Iz1= Iz01(2);

    if isfield(data, 'Ip')
        Ip01= getElementValue(data, e, 'Ip');
    else
        if isfield(data, 'GIp')
            Ip01= getElementValue(data, e, 'GIp') / G;
        else
            Ip01= [100 100]; % very stiff
        end
    end
    Ip0= Ip01(1);
    Ip1= Ip01(2);

    K{e}= [[((A1+A0)*E)/(2*le),0,0,0,0,0,-((A1+A0)*E)/(2*le),0,0,0,0,0];[0,((6*Iz1+6*Iz0)*E)/le^3,0,0,0,((2*Iz1+4*Iz0)*E)/le^2,0,-((6*Iz1+6*Iz0)*E)/le^3,0,0,0,((4*Iz1+2*Iz0)*E)/le^2];[0,0,((6*Iy1+6*Iy0)*E)/le^3,0,-((2*Iy1+4*Iy0)*E)/le^2,0,0,0,-((6*Iy1+6*Iy0)*E)/le^3,0,-((4*Iy1+2*Iy0)*E)/le^2,0];[0,0,0,((Ip1+Ip0)*G)/(2*le),0,0,0,0,0,-((Ip1+Ip0)*G)/(2*le),0,0];[0,0,-((2*Iy1+4*Iy0)*E)/le^2,0,((Iy1+3*Iy0)*E)/le,0,0,0,((2*Iy1+4*Iy0)*E)/le^2,0,((Iy1+Iy0)*E)/le,0];[0,((2*Iz1+4*Iz0)*E)/le^2,0,0,0,((Iz1+3*Iz0)*E)/le,0,-((2*Iz1+4*Iz0)*E)/le^2,0,0,0,((Iz1+Iz0)*E)/le];[-((A1+A0)*E)/(2*le),0,0,0,0,0,((A1+A0)*E)/(2*le),0,0,0,0,0];[0,-((6*Iz1+6*Iz0)*E)/le^3,0,0,0,-((2*Iz1+4*Iz0)*E)/le^2,0,((6*Iz1+6*Iz0)*E)/le^3,0,0,0,-((4*Iz1+2*Iz0)*E)/le^2];[0,0,-((6*Iy1+6*Iy0)*E)/le^3,0,((2*Iy1+4*Iy0)*E)/le^2,0,0,0,((6*Iy1+6*Iy0)*E)/le^3,0,((4*Iy1+2*Iy0)*E)/le^2,0];[0,0,0,-((Ip1+Ip0)*G)/(2*le),0,0,0,0,0,((Ip1+Ip0)*G)/(2*le),0,0];[0,0,-((4*Iy1+2*Iy0)*E)/le^2,0,((Iy1+Iy0)*E)/le,0,0,0,((4*Iy1+2*Iy0)*E)/le^2,0,((3*Iy1+Iy0)*E)/le,0];[0,((4*Iz1+2*Iz0)*E)/le^2,0,0,0,((Iz1+Iz0)*E)/le,0,-((4*Iz1+2*Iz0)*E)/le^2,0,0,0,((3*Iz1+Iz0)*E)/le]];
end

KF= zeros(nF);
for e= 1:ne
    KF= KF + T{e}' * K{e} * T{e};
end


T__= eye(nF);

if ~exist('ref_system', 'var') && ~isempty(ref_system)
    ref_system= 1;
end
switch ref_system
    case 1
        if isfield(data, 'K_foundation')
            if prod(size(data.K_foundation))==length(data.K_foundation)
                data.K_foundation= diag(data.K_foundation);
            end
            for i= 1:6
                if isinf(data.K_foundation(i, i))
                    T__(:, i)= [];
                end
            end
            KF(1:6, 1:6)= KF(1:6, 1:6) + data.K_foundation;
            offsetX= 7 - sum(isinf(diag(data.K_foundation)));
        else    
            % fixed base
            T__(:, 1:nqk)= [];
            offsetX= 1;
        end
    case 2
        % Sehnensystem
        xyz= 1:3;
        free_xyz= setdiff(xyz, extentIdx);
        T__(:, [1 2 3 3+extentIdx nF-6+free_xyz(1) nF-6+free_xyz(2)])= [];
        offsetX= 3;
    otherwise
        error('Reference system number %d not known', ref_system);
end

MF__= T__'*MF*T__;
KF__= T__'*KF*T__;

[V, D]= eig(KF__, MF__);
EF= sqrt(diag(D))/2/pi;

%  C1= cell(ne, 1);
%  for e= 1:ne
%      le= data(e).dx;
%      if isfield(data, 'A') && isfield(data, 'rho')
%          me01= getElementValue(data, e, 'A') * le * data(e).rho;
%      else
%          if isfield(data, 'mu')
%              me01= getElementValue(data, e, 'mu') * le;
%          else
%              me01= getElementValue(data, e, 'me');
%          end
%      end
%      me0= me01(1);
%      me1= me01(2);
%      C1= [[(le*me1+2*le*me0)/(6*le),0,0,0,0,0,(2*le^2*me1+le^2*me0)/(6*le^2),0,0,0,0,0];[0,(3*le*me1+7*le*me0)/(20*le),0,0,0,(2*le*me1+3*le*me0)/60,0,(7*le*me1+3*le*me0)/(20*le),0,0,0,-(3*le*me1+2*le*me0)/60];[0,0,(3*le*me1+7*le*me0)/(20*le),0,-(2*le*me1+3*le*me0)/60,0,0,0,(7*le*me1+3*le*me0)/(20*le),0,(3*le*me1+2*le*me0)/60,0]];

%  end
%  
%  % (5.239) S. 231, (5.166) S. 210
%  CFt= zeros(nF, 3);
%  for e= 1:ne
%      CFt= CFt + T{e}' * C1{e}' * Gamma{e};
%  end
%  CFt= T__'*CFt;

Sr= zeros(nF, 3);
St= zeros(nF, 3);
for i= 1:nk
    % (6.470) S. 365
    St(((i-1)*nqk+1):((i-1)*nqk+3), :)= eye(3);

    % (6.475) S. 366
    R= crossmat(data(i).x);
    Sr(((i-1)*nqk+1):((i-1)*nqk+3), :)= -R;
    Sr(((i-1)*nqk+4):((i-1)*nqk+6), :)= eye(3);
end

% (6.478) S. 366
if exist('modes', 'var') && ~isempty(modes)
    if iscell(modes)
        modes_= [];
        for i= 1:length(modes)
            mode_pair= modes{i};
            switch length(mode_pair)
                case 1
                    modes_= [modes_ mode_pair];
                case 2
                    if abs(EF(mode_pair(1))-EF(mode_pair(2)))/max(EF(mode_pair)) > 1e-5
                        error('Mode pair (%d, %d) does not have identical eigenvalues', mode_pair(1), mode_pair(2));
                    end
                    V= orthogonalizeV(V, mode_pair, offsetX);
                    modes_= [modes_ mode_pair(:)'];                    
                otherwise
                    error('Currently no mode pairs with more than two elements are supported')
            end
        end
        modes= modes_;
    end
    if exist('normalize', 'var')
        switch normalize
            case 1
                V= normalize_to_last(V, modes, offsetX);
            otherwise
                error('Normalization method other than 1 currently not supported');
        end
    end
    
    Se= T__*V(:, modes);
    sid.modes= modes;
else
    Se= T__;
end
%  Se= T__;
nq= size(Se, 2);

% (6.479) S. 367
mE= St'*MF*St;
I0= Sr'*MF*Sr;
Me= Se'*MF*Se;
mc0= Sr'*MF*St;
Ct0= Se'*MF*St; % identical to CFt?
Ct0_= T__'*MF*St; % identical to CFt?
Cr0= Se'*MF*Sr; 

C3= cell(ne, 1);
for e= 1:ne
    le= data(e).l;
    if isfield(data, 'A') && isfield(data, 'rho')
        me01= getElementValue(data, e, 'A') * le * data(e).rho;
    else
        if isfield(data, 'mu')
            me01= getElementValue(data, e, 'mu') * le;
        else
            me01= getElementValue(data, e, 'me');
        end
    end
    me0= me01(1);
    me1= me01(2);
    C3{e}= cell(3, 3);

    C3{e}{1,1}= [[(me1+3*me0)/12,0,0,0,0,0,(me1+me0)/12,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[(me1+me0)/12,0,0,0,0,0,(3*me1+me0)/12,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0]];
    C3{e}{1,2}= [[0,(5*me1+16*me0)/60,0,0,0,(le*me1+2*le*me0)/60,0,(5*me1+4*me0)/60,0,0,0,-(le*me1+le*me0)/60];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,(4*me1+5*me0)/60,0,0,0,(le*me1+le*me0)/60,0,(16*me1+5*me0)/60,0,0,0,-(2*le*me1+le*me0)/60];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0]];
    C3{e}{1,3}= [[0,0,(5*me1+16*me0)/60,0,-(le*me1+2*le*me0)/60,0,0,0,(5*me1+4*me0)/60,0,(le*me1+le*me0)/60,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,(4*me1+5*me0)/60,0,-(le*me1+le*me0)/60,0,0,0,(16*me1+5*me0)/60,0,(2*le*me1+le*me0)/60,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0]];
    C3{e}{2,1}= [[0,0,0,0,0,0,0,0,0,0,0,0];[(5*me1+16*me0)/60,0,0,0,0,0,(4*me1+5*me0)/60,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[(le*me1+2*le*me0)/60,0,0,0,0,0,(le*me1+le*me0)/60,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[(5*me1+4*me0)/60,0,0,0,0,0,(16*me1+5*me0)/60,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[-(le*me1+le*me0)/60,0,0,0,0,0,-(2*le*me1+le*me0)/60,0,0,0,0,0]];
    C3{e}{2,2}= [[0,0,0,0,0,0,0,0,0,0,0,0];[0,(3*me1+10*me0)/35,0,0,0,(7*le*me1+15*le*me0)/420,0,(9*me1+9*me0)/140,0,0,0,-(6*le*me1+7*le*me0)/420];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,(7*le*me1+15*le*me0)/420,0,0,0,(3*le^2*me1+5*le^2*me0)/840,0,(7*le*me1+6*le*me0)/420,0,0,0,-(le^2*me1+le^2*me0)/280];[0,0,0,0,0,0,0,0,0,0,0,0];[0,(9*me1+9*me0)/140,0,0,0,(7*le*me1+6*le*me0)/420,0,(10*me1+3*me0)/35,0,0,0,-(15*le*me1+7*le*me0)/420];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,-(6*le*me1+7*le*me0)/420,0,0,0,-(le^2*me1+le^2*me0)/280,0,-(15*le*me1+7*le*me0)/420,0,0,0,(5*le^2*me1+3*le^2*me0)/840]];
    C3{e}{2,3}= [[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,(3*me1+10*me0)/35,0,-(7*le*me1+15*le*me0)/420,0,0,0,(9*me1+9*me0)/140,0,(6*le*me1+7*le*me0)/420,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,(7*le*me1+15*le*me0)/420,0,-(3*le^2*me1+5*le^2*me0)/840,0,0,0,(7*le*me1+6*le*me0)/420,0,(le^2*me1+le^2*me0)/280,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,(9*me1+9*me0)/140,0,-(7*le*me1+6*le*me0)/420,0,0,0,(10*me1+3*me0)/35,0,(15*le*me1+7*le*me0)/420,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,-(6*le*me1+7*le*me0)/420,0,(le^2*me1+le^2*me0)/280,0,0,0,-(15*le*me1+7*le*me0)/420,0,-(5*le^2*me1+3*le^2*me0)/840,0]];
    C3{e}{3,1}= [[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[(5*me1+16*me0)/60,0,0,0,0,0,(4*me1+5*me0)/60,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[-(le*me1+2*le*me0)/60,0,0,0,0,0,-(le*me1+le*me0)/60,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[(5*me1+4*me0)/60,0,0,0,0,0,(16*me1+5*me0)/60,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[(le*me1+le*me0)/60,0,0,0,0,0,(2*le*me1+le*me0)/60,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0]];
    C3{e}{3,2}= [[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,(3*me1+10*me0)/35,0,0,0,(7*le*me1+15*le*me0)/420,0,(9*me1+9*me0)/140,0,0,0,-(6*le*me1+7*le*me0)/420];[0,0,0,0,0,0,0,0,0,0,0,0];[0,-(7*le*me1+15*le*me0)/420,0,0,0,-(3*le^2*me1+5*le^2*me0)/840,0,-(7*le*me1+6*le*me0)/420,0,0,0,(le^2*me1+le^2*me0)/280];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,(9*me1+9*me0)/140,0,0,0,(7*le*me1+6*le*me0)/420,0,(10*me1+3*me0)/35,0,0,0,-(15*le*me1+7*le*me0)/420];[0,0,0,0,0,0,0,0,0,0,0,0];[0,(6*le*me1+7*le*me0)/420,0,0,0,(le^2*me1+le^2*me0)/280,0,(15*le*me1+7*le*me0)/420,0,0,0,-(5*le^2*me1+3*le^2*me0)/840];[0,0,0,0,0,0,0,0,0,0,0,0]];
    C3{e}{3,3}= [[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,(3*me1+10*me0)/35,0,-(7*le*me1+15*le*me0)/420,0,0,0,(9*me1+9*me0)/140,0,(6*le*me1+7*le*me0)/420,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,-(7*le*me1+15*le*me0)/420,0,(3*le^2*me1+5*le^2*me0)/840,0,0,0,-(7*le*me1+6*le*me0)/420,0,-(le^2*me1+le^2*me0)/280,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,(9*me1+9*me0)/140,0,-(7*le*me1+6*le*me0)/420,0,0,0,(10*me1+3*me0)/35,0,(15*le*me1+7*le*me0)/420,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,(6*le*me1+7*le*me0)/420,0,-(le^2*me1+le^2*me0)/280,0,0,0,(15*le*me1+7*le*me0)/420,0,(5*le^2*me1+3*le^2*me0)/840,0];[0,0,0,0,0,0,0,0,0,0,0,0]];
    
    % not 100% sure if this is correct
    if isfield(data, 'Mlumped')
        C3{e}{1,1}(1, 1)= C3{e}{1,1}(1, 1) + data(e).Mlumped;
        C3{e}{2,2}(2, 2)= C3{e}{2,2}(2, 2) + data(e).Mlumped;
        C3{e}{3,3}(3, 3)= C3{e}{3,3}(3, 3) + data(e).Mlumped;
        if e==ne
            C3{e}{1,1}(1+6, 1+6)= C3{e}{1,1}(1+6, 1+6) + data(e+1).Mlumped;
            C3{e}{2,2}(2+6, 2+6)= C3{e}{2,2}(2+6, 2+6) + data(e+1).Mlumped;
            C3{e}{3,3}(3+6, 3+6)= C3{e}{3,3}(3+6, 3+6) + data(e+1).Mlumped;
        end
    end
end

% (5.252) S. 233, (6.401) S. 338
KFr= cell(3, 1);
Kr= cell(3, 1);
for a= 1:3
    KFr{a}= zeros(nF);
    for e= 1:ne
        lmn= 1:3;
        for l= 1:3
            m= lmn(2);
            n= lmn(3);
            KFr{a}= KFr{a} + T{e}' * (-C3{e}{m, n} + C3{e}{n, m}) * T{e} * Gamma{e}(l, a);
                
            lmn= circshift(lmn, [0 -1]);
        end
    end
    % (6.483) S. 367
    Kr{a}= Se'*KFr{a}*Se;
end

ZF0= zeros(nF, 1);
for k= 1:nk
    ZF0((nqk*(k-1)+1):(nqk*(k-1)+3), 1)= data(k).x;
end
      
% alternative Formulation for Cr0
%  CFr0= [KFr{1}*ZF0 KFr{2}*ZF0 KFr{3}*ZF0];

% (5.268) S. 237
KFom_ab= cell(3, 3); % = C6
for a= 1:3
    for b= 1:3
        KFom_ab{a, b}= zeros(nF);
        for l= 1:3
            for m= 1:3
                for e= 1:ne
                    if l==m
                        m_= l+1;
                        if m_>3, m_= 1; end
                        n_= m_+1;
                        if n_>3, n_= 1; end
                        % (5.266) S. 236
                        Xi= -(C3{e}{m_, m_}+C3{e}{n_, n_});
                    else
                        Xi= C3{e}{m, l};
                    end
                    KFom_ab{a, b}= KFom_ab{a, b} + T{e}' * Gamma{e}(l, a) * Xi * Gamma{e}(m, b) * T{e};
                end
            end
        end
    end
end

% (5.271) S. 237
KFom= cell(6, 1);
Kom= cell(6, 1);
Kom0= zeros(nq, 6);
Kom0_= zeros(nq_, 6);
for i= 1:6
    if i<4
        KFom{i}= KFom_ab{i, i};
    else
        a= i-3;
        b= a+1;
        if b>3, b= 1; end
        KFom{i}= KFom_ab{a, b} + KFom_ab{a, b}';
    end
    Kom{i}= Se'*KFom{i}*Se;

    Kom0(:, i)= Se'*KFom{i}*ZF0;
    Kom0_(:, i)= T__'*KFom{i}*ZF0;
end

% (6.490) S. 368; (6.531) S. 379 or (6.515) S. 375
C4= zeros(3, 3, nq);
for l= 1:nq
    for a= 1:3
        for b= 1:3
            C4(a, b, l)= -Sr(:, a)'*KFr{b}*Se(:, l);
        end
    end
end

Kinv= T__*KF__^-1*T__';

% only axial stiffening
% 6.330 S. 319
Fend_ax= zeros(nF, 1);
Fend_ax(end-nqk+extentIdx, 1)= 1; % tip loads only for now
K0Fend_ax= Se'*geo_stiff(ne, nF, nqk, data, T, -Kinv*Fend_ax)*Se; 
K0t_ax= Se'*geo_stiff(ne, nF, nqk, data, T, -Kinv*T__*Ct0_(:, extentIdx))*Se;

K0omxx= Se'*geo_stiff(ne, nF, nqk, data, T, -Kinv*T__*Kom0_(:, 1))*Se; 
K0omyy= Se'*geo_stiff(ne, nF, nqk, data, T, -Kinv*T__*Kom0_(:, 2))*Se; 
K0omzz= Se'*geo_stiff(ne, nF, nqk, data, T, -Kinv*T__*Kom0_(:, 3))*Se; 
K0omxy= Se'*geo_stiff(ne, nF, nqk, data, T, -Kinv*T__*Kom0_(:, 4))*Se; 
K0omxz= Se'*geo_stiff(ne, nF, nqk, data, T, -Kinv*T__*Kom0_(:, 5))*Se; 
K0omyz= Se'*geo_stiff(ne, nF, nqk, data, T, -Kinv*T__*Kom0_(:, 6))*Se; 


% idx_t= [ones(3, ne); zeros(3, ne)];
% idx_t= logical(idx_t(:));
% idx_r= [zeros(3, ne); ones(3, ne)];
% idx_r= logical(idx_r(:));

if isunix
    [~, user_name] = system('whoami'); % exists on every unix that I know of
    % on my mac, isunix == 1
elseif ispc
    [~, user_name] = system('echo %USERDOMAIN%\%USERNAME%'); % Not as familiar with windows,
    % found it on the net elsewhere, you might want to verify
end
user_name= strrep(user_name, char(10), '');

if exist('modes', 'var')
    sid.comment= sprintf('%d nodes, %d modes, generated by FEMBeam2SID MATLAB script on %s by %s', nk, nq, datestr(now), user_name);
else
    sid.comment= sprintf('%d nodes, %d coordinates, generated by FEMBeam2SID MATLAB script on %s by %s', nk, nq, datestr(now), user_name);
end
% Tabelle 6.9 S. 346
%  sid.refmod % Masse und Name der Koordinaten q_i
sid.refmod.mass= mE(1); % Masse mi des KÃ¶rpers
sid.refmod.nelastq= nq; % Zahl der Koordinaten
coord_name= {'ux' 'uy' 'uz' 'rx' 'ry' 'rz'};
if exist('modes', 'var')
    for i= 1:nq
        sid.refmod.ielastq{i}= sprintf('Eigen Mode %4d : %13f Hz', modes(i), EF(i)); % Name der Koordinaten
    end
else
    for i= 1:nq
        iq= find(T__(:, i), 1);
        in= ceil(iq/6);
        ic= iq - (in-1)*6;
        sid.refmod.ielastq{i}= sprintf('Node %4d, %s', in, coord_name{ic}); % Name der Koordinaten
    end
end
% sid.frame % Daten zur Berechnung der Bewegung von 
for i= 1:nk
    sid.frame(i).node= num2str(i); % Name des Knotens
    sid.frame(i).rframe= 'body ref';  % Name des Bezugssystems
    
    Phi= zeros(3, nF);
    Phi(:, (i-1)*nqk + (1:3))= eye(3);
    Psi= zeros(3, nF);
    Psi(:, (i-1)*nqk + (4:6))= eye(3);

    sid.frame(i).origin= emptyTaylor(1, 3, 1, nq, 0, 3);
    sid.frame(i).origin.M0= data(i).x; % Ort von Ok (6.78)
    for j= 1:nq
        sid.frame(i).origin.M1(:, 1, j)= Phi*Se(:, j);
    end
    sid.frame(i).Phi= emptyTaylor(0, 3, nq, nq, 0, 3);
    sid.frame(i).Phi.M0= Phi*Se; % node translations (6.46),(6.80), (6.337), (5.122), (6.415)
    if i==nk
        sid.frame(i).Phi.order= 1;
        sid.frame(i).Phi.M1(1, :, :)= 0*K0Fend_ax; % stiffening due to force (6.46),(6.80), (6.337), (5.122), (6.415)
        sid.frame(i).Phi.M1(2, :, :)= 0*K0Fend_ax;
        sid.frame(i).Phi.M1(3, :, :)= 0*K0Fend_ax;
        sid.frame(i).Phi.M1(extentIdx, :, :)= K0Fend_ax;        
    end
    
    sid.frame(i).Psi= emptyTaylor(0, 3, nq, nq, 0, 3);
    sid.frame(i).Psi.M0= Psi*Se; % node rotations (6.81), (6.337)
%      sid.frame(i).Psi.M1= 0; % no geometric stiffening due to moments for now

    sid.frame(i).AP= emptyTaylor(1, 3, 3, nq, 0, 3);
    if i==nk
        sid.frame(i).AP.M0= data(i-1).D; % Orientierung von ek (6.79)
    else
        sid.frame(i).AP.M0= data(i).D; % Orientierung von ek (6.79)
    end
    for j= 1:nq
        sid.frame(i).AP.M1(:, :, j)= crossmat(sid.frame(i).Psi.M0(:, j));
    end
   
    sid.frame(i).sigma= emptyTaylor(1, 6, 1, nq, 0, 3);
    sid.frame(i).sigma.M0= zeros(6, 1); % no pretension for now. Spannungen  % (6.413)
    if i<nk
        e= i;
    else
        e= i-1;
    end
    Ke= T{e}' * K{e} * T{e};
    if i<nk
        Ke= -Ke;
    end
    for j= 1:nq
        sid.frame(i).sigma.M1(:, 1, j)= Ke((i-1)*nqk+(1:6), :)*Se(:, j); % left sided nodal force except last node
    end
end
% sid.mdCM
sid.md= emptyTaylor(1, 3, 1, nq, 0, 3);
sid.md.M0= [sum(sum(crossmat([0.5 0 0]).*mc0)) sum(sum(crossmat([0 0.5 0]).*mc0)) sum(sum(crossmat([0 0 0.5]).*mc0))]';
for j= 1:nq
    sid.md.M1(:, 1, j)= Ct0(j, :)';
end
% sid.J
sid.I= emptyTaylor(1, 3, 3, nq, 0, 2);
sid.I.M0= I0;
for i= 1:nq
    sid.I.M1(:, :, i)= -C4(:, :, i) - C4(:, :, i)';
end
% sid.Ct
sid.Ct= emptyTaylor(1, nq, 3, nq, 0, 3);
sid.Ct.M0= Ct0;
sid.Ct.M1(:, 1, :)= 0*K0t_ax;
sid.Ct.M1(:, 2, :)= 0*K0t_ax;
sid.Ct.M1(:, 3, :)= 0*K0t_ax;
sid.Ct.M1(:, extentIdx, :)= K0t_ax;
% sid.Cr
sid.Cr= emptyTaylor(1, nq, 3, nq, 0, 3);
sid.Cr.M0= Cr0;
sid.Cr.M1(:, 1, :)= Kr{1}; % K0r= 0 because only axial stiffening.  Kr= Se*KFr?
sid.Cr.M1(:, 2, :)= Kr{2}; % K0r= 0 because only axial stiffening.  Kr= Se*KFr?
sid.Cr.M1(:, 3, :)= Kr{3}; % K0r= 0 because only axial stiffening.  Kr= Se*KFr?
% sid.Me
sid.Me= emptyTaylor(0, nq, nq, 0, 0, 2);
sid.Me.M0= Me;
% sid.Gr
sid.Gr= emptyTaylor(0, 3, 3*nq, nq, 0, 3);
sid.Gr.M0= -2*reshape(C4, 3, 3*nq); % (6.403) S. 339; or -2*KFr?
%  sid.Gr.M1
% sid.Ge
sid.Ge= emptyTaylor(0, nq, 3*nq, 0, 0, 3);
for i= 1:nq
    sid.Ge.M0(:, 3*(i-1)+(1:3))= 2*[Kr{1}(:, i) Kr{2}(:, i) Kr{3}(:, i)]; % (6.405) S. 340 = 2*C5'
end
%  sid.Oe (6.407)
sid.Oe= emptyTaylor(1, nq, 6, nq, 0, 3);
sid.Oe.M0= Kom0;
sid.Oe.M1(:, 1, :)= Kom{1}+K0omxx;
sid.Oe.M1(:, 2, :)= Kom{2}+K0omyy;
sid.Oe.M1(:, 3, :)= Kom{3}+K0omzz;
sid.Oe.M1(:, 4, :)= Kom{4}+K0omxy;
sid.Oe.M1(:, 5, :)= Kom{5}+K0omxz;
sid.Oe.M1(:, 6, :)= Kom{6}+K0omyz;

%  sid.ksigma
%  sid.ksigma.M0= zeros(size(V, 2), 1); % no pretension for now
%  sid.ksigma.M1= zeros(size(V, 2));

%  sid.Ke
sid.Ke= emptyTaylor(0, nq, nq, 0, 0, 2);
sid.Ke.M0= Se'*KF*Se;
% sid.Ke.M1
sid.ksigma= emptyTaylor(0, nq, 1, nq, 0, 0);
sid.ksigma.M0= [];
sid.De= emptyTaylor(0, nq, nq, 0, 0, 0);
sid.De.M0= [];

sid.EF= EF;

function KFsigma_sharp= geo_stiff(ne, nF, nqk, data, T, zF)
% Spannung: 5.220 S. 226,  Beispiel: 5.287 S. 241
% E*BL : Spannung, (6.311) S. 312 [(6.259) S. 301, (6.300) S. 309]
% K^-1 : Verformung aus (Einheits-)Kraft
% z.B. C4(:, 3) : Einheitskraft
% E bzw. H macht aus Dehnung Spannung
% BL macht aus Knotenkoordinaten Dehnung (eps_xx)
% E bzw. H und BL sind schon in K_G1_ux1 + K_G1_ux2 enthalten

KFsigma_sharp= zeros(nF);
for e= 1:ne
    le= data(e).l;
    if isfield(data, 'E')
        E= data(e).E;
    else
        E= 214e9;
    end
    if isfield(data, 'A')
        A01= getElementValue(data, e, 'A');
    else
        if isfield(data, 'EA')
            A01= getElementValue(data, e, 'EA') / E;
        else
            A01= [100, 100]; % very stiff
        end
    end
    A0= A01(1);
    A1= A01(2);

    K_Ge_ux1_lin= [[0,0,0,0,0,0,0,0,0,0,0,0];[0,-((3*A1+3*A0)*E)/(5*le^2),0,0,0,-(A1*E)/(10*le),0,((3*A1+3*A0)*E)/(5*le^2),0,0,0,-(A0*E)/(10*le)];[0,0,-((3*A1+3*A0)*E)/(5*le^2),0,(A1*E)/(10*le),0,0,0,((3*A1+3*A0)*E)/(5*le^2),0,(A0*E)/(10*le),0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,(A1*E)/(10*le),0,-((A1+3*A0)*E)/30,0,0,0,-(A1*E)/(10*le),0,((A1+A0)*E)/60,0];[0,-(A1*E)/(10*le),0,0,0,-((A1+3*A0)*E)/30,0,(A1*E)/(10*le),0,0,0,((A1+A0)*E)/60];[0,0,0,0,0,0,0,0,0,0,0,0];[0,((3*A1+3*A0)*E)/(5*le^2),0,0,0,(A1*E)/(10*le),0,-((3*A1+3*A0)*E)/(5*le^2),0,0,0,(A0*E)/(10*le)];[0,0,((3*A1+3*A0)*E)/(5*le^2),0,-(A1*E)/(10*le),0,0,0,-((3*A1+3*A0)*E)/(5*le^2),0,-(A0*E)/(10*le),0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,(A0*E)/(10*le),0,((A1+A0)*E)/60,0,0,0,-(A0*E)/(10*le),0,-((3*A1+A0)*E)/30,0];[0,-(A0*E)/(10*le),0,0,0,((A1+A0)*E)/60,0,(A0*E)/(10*le),0,0,0,-((3*A1+A0)*E)/30]];
    K_Ge_ux2_lin= [[0,0,0,0,0,0,0,0,0,0,0,0];[0,((3*A1+3*A0)*E)/(5*le^2),0,0,0,(A1*E)/(10*le),0,-((3*A1+3*A0)*E)/(5*le^2),0,0,0,(A0*E)/(10*le)];[0,0,((3*A1+3*A0)*E)/(5*le^2),0,-(A1*E)/(10*le),0,0,0,-((3*A1+3*A0)*E)/(5*le^2),0,-(A0*E)/(10*le),0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,-(A1*E)/(10*le),0,((A1+3*A0)*E)/30,0,0,0,(A1*E)/(10*le),0,-((A1+A0)*E)/60,0];[0,(A1*E)/(10*le),0,0,0,((A1+3*A0)*E)/30,0,-(A1*E)/(10*le),0,0,0,-((A1+A0)*E)/60];[0,0,0,0,0,0,0,0,0,0,0,0];[0,-((3*A1+3*A0)*E)/(5*le^2),0,0,0,-(A1*E)/(10*le),0,((3*A1+3*A0)*E)/(5*le^2),0,0,0,-(A0*E)/(10*le)];[0,0,-((3*A1+3*A0)*E)/(5*le^2),0,(A1*E)/(10*le),0,0,0,((3*A1+3*A0)*E)/(5*le^2),0,(A0*E)/(10*le),0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,-(A0*E)/(10*le),0,-((A1+A0)*E)/60,0,0,0,(A0*E)/(10*le),0,((3*A1+A0)*E)/30,0];[0,(A0*E)/(10*le),0,0,0,-((A1+A0)*E)/60,0,-(A0*E)/(10*le),0,0,0,((3*A1+A0)*E)/30]];
    ux= T{e}*zF;
    KFsigma_sharp_e= K_Ge_ux1_lin*ux(1) + K_Ge_ux2_lin*ux(1+6);
    % 5.208 S. 222
    KFsigma_sharp= KFsigma_sharp + T{e}'*KFsigma_sharp_e*T{e};
end

function t= emptyTaylor(order, nrow, ncol, nq, nqn, str)
t.order= order;
t.nrow= nrow;
t.ncol= ncol;
t.nq= nq;
t.nqn= nqn;
t.structure= str;

function v= getElementValue(data, i, element)
if ~isfield(data, 'node_interpolation')
    node_interpolation= 'linear';
else
    node_interpolation= data(i).node_interpolation;
end
switch node_interpolation
    case 'mean'
        v= [0.5 0.5]*(data(i).(element)+data(i+1).(element));
    case 'left'
        v= [data(i).(element) data(i).(element)];
    case 'right'
        v= [data(i+1).(element) data(i+1).(element)];
    case 'linear'
        v= [data(i).(element) data(i+1).(element)];
    otherwise
        error('Unknown node interpolation method "%s"', node_interpolation);
end    

function V= orthogonalizeV(V, mode_pair, offsetX)
m1= xyzMagnitude(V(:, mode_pair(1)), offsetX);
[~, idx]= sort(m1, 'descend');
T(1, 1)= sum(V(offsetX+idx(1)-1:6:end, mode_pair(1)));
T(1, 2)= sum(V(offsetX+idx(2)-1:6:end, mode_pair(1)));
T(2, 1)= sum(V(offsetX+idx(1)-1:6:end, mode_pair(2)));
T(2, 2)= sum(V(offsetX+idx(2)-1:6:end, mode_pair(2)));

V1= T(1, 1)*V(:, mode_pair(1)) + T(2, 1)*V(:, mode_pair(2));
V2= T(1, 2)*V(:, mode_pair(1)) + T(2, 2)*V(:, mode_pair(2));

V(:, mode_pair(1))= V1;
V(:, mode_pair(2))= V2;

function V= normalize_to_last(V, modes, offsetX)
for i= 1:length(modes)
    mode= modes(i);
    
    m= xyzMagnitude(V(:, mode), offsetX);
    [~, normtoIdx]= max(m);
    v_= V(offsetX+normtoIdx-1:6:end, mode);
    V(:, mode)= V(:, mode)/v_(end);
end

function m= xyzMagnitude(v, offsetX)
% TODO also consider torsion modes
m(1)= sum(abs(v(offsetX:6:end)));
offsetY= offsetX+1;
m(2)= sum(abs(v(offsetY:6:end)));
offsetZ= offsetY+1;
m(3)= sum(abs(v(offsetZ:6:end)));
