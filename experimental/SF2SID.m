function sid= SF2SID(R, Phi, dPhi, ddPhi, m, IE, damping, twist)

% R is expected to have two entries more than Phi. The first element
% contains the beam root. The last element contains the beam tip position.
% All other elements are the element centers.
% The two extra entries of R are removed in line 24

n_elem= length(m);
n_q= size(Phi{1}, 2);

%% orientation of undeflected elemets
DR= cell(n_elem, 1);
last_node= R{1};
for i_elem= 1:n_elem
    if i_elem==n_elem
        node= R{i_elem+2};
    else
        node= 0.5*(R{i_elem+2} + R{i_elem+1});
    end
    DR{i_elem}= node-last_node;
    last_node= node;
end

R= R(2:end-1);

% E= cell(n_elem, 1);
% for i_elem= 1:n_elem
%     % calculate the tangent
%     l= norm(DR{i_elem});
%     t= DR{i_elem}/l;
%     if i_elem==1
%         t_last= t;
%         t_origin= t;
%         % start for the normal
%         % for t pointing in x or z direction, h will be in y direction
%         % for t pointing in y direction, h will be in -x direction
%         h_last= null(t_last');
%         h_last= h_last(:, 1);
%     end
% 
%     % calculate normal vector (4.33)
%     dt= t - t_last;
%     if norm(dt)<eps
%         h= h_last;
%     else
%         h= dt/norm(dt);
%     end
% 
%     % add twist
%     if exist('twist', 'var')
%         % rotation about t
% %         rot_mat= cos(twist{i_elem})*eye(3) + sin(twist{i_elem})*crossmat(t) + (1-cos(twist{i_elem}))*(t*t');
%         rot_mat= cos(twist(i_elem))*eye(3) + sin(twist(i_elem))*crossmat(t_origin) + (1-cos(twist(i_elem)))*(t_origin*t_origin');
%         h= rot_mat * h;
%     end
%   
%     % calculate binormal vector
%     b= cross(t, h);
% 
%     % element orientation (4.36)
%     E{i_elem}= [t h b]';
% 
%     t_last= t;
%     h_last= h;
% end


%% calculation of basic integrals
mass= sum(m);

mc0= zeros(3, 1);
for i_elem= 1:n_elem
    mc0= mc0 + R{i_elem}*m(i_elem);
end

I0= zeros(3);
for i_elem= 1:n_elem
    I0= I0 + crossmat(R{i_elem})*crossmat(R{i_elem})'*m(i_elem);
end

C1= zeros(3, n_q);
for i_elem= 1:n_elem
    C1= C1 + Phi{i_elem}*m(i_elem);
end

C2= zeros(3, n_q);
for i_elem= 1:n_elem
    C2= C2 + crossmat(R{i_elem})*Phi{i_elem}*m(i_elem);
end

C3= cell(3, 3);
for a= 1:3
    for b= 1:3
        C3{a, b}= zeros(n_q, n_q);
        for i_elem= 1:n_elem
            C3{a, b}= C3{a, b} + Phi{i_elem}(a, :)'*Phi{i_elem}(b, :)*m(i_elem);
        end
    end
end

C4= cell(n_q, 1);
for l= 1:n_q
    C4{l}= zeros(3, 3);
    for i_elem= 1:n_elem
        C4{l}= C4{l} + crossmat(R{i_elem})*crossmat(Phi{i_elem}(:, l))*m(i_elem);
    end
end

C5= cell(n_q, 1);
for l= 1:n_q
    C5{l}= zeros(3, n_q);
    for i_elem= 1:n_elem
        C5{l}= C5{l} + crossmat(Phi{i_elem}(:, l))*Phi{i_elem}*m(i_elem);
    end
end

C6= cell(n_q, n_q);
for l= 1:n_q
    for k= 1:n_q
        C6{k, l}= zeros(3, 3);
        for i_elem= 1:n_elem
            C6{k, l}= C6{k, l} + crossmat(Phi{i_elem}(:, k))*crossmat(Phi{i_elem}(:, l))*m(i_elem);
        end
    end
end

Me= C3{1, 1} + C3{2, 2} + C3{3, 3};

Komega= cell(3);
for a= 1:3
    for b= 1:3
        if a~=b
            Komega{a, b}= C3{a, b};
        else
            c= mod(a, 3)+1;
            d= mod(c, 3)+1;
            Komega{a, b}= -C3{c, c} - C3{d, d};
        end
    end
end

Kr= cell(3, 1);
for a= 1:3
    c= mod(a, 3)+1;
    d= mod(c, 3)+1;
    Kr{a}= -C3{c, d} + C3{c, d}';
end
 
%% stiffness and damping
Ke= zeros(n_q);
for i_elem= 1:n_elem
    stiff= diag(IE{i_elem})*norm(DR{i_elem});
    Ke= Ke + ddPhi{i_elem}' * stiff * ddPhi{i_elem};
end

EF= sqrt(diag(Ke)./diag(Me))/2/pi;

De= Ke*diag(damping(:)./(pi*EF));

%% geometric stiffnesses

% stiffening due to point load
K0F= cell(n_elem, 1);
for i_elem= 1:n_elem
    if i_elem>1
        K0F{i_elem}= K0F{i_elem-1};
    else
        K0F{i_elem}= zeros(3, n_q, n_q);
    end
    for a= 1:3
        e= zeros(1, 3);
        e(a)= 1;
        K0F{i_elem}(a, :, :)= squeeze(K0F{i_elem}(a, :, :)) + e*DR{i_elem}*dPhi{i_elem}'*dPhi{i_elem};
    end
end

% stiffening due to acceleration
% force due to acceleration results from cumulated mass m_cum
m_cum= flip(cumsum(m(end:-1:1))) - 0.5*m;

K0t= cell(3, 1); 
for a= 1:3
    K0t{a}= zeros(n_q);
end
for i_elem= 1:n_elem
    for a= 1:3
        e= zeros(1, 3);
        e(a)= 1;
        K0t{a}= K0t{a} + m_cum(i_elem)*e*DR{i_elem}*dPhi{i_elem}'*dPhi{i_elem};
    end
end

% no stiffening due to rotational acceleration for now
K0r= cell(3, 1);
for a= 1:3
    K0r{a}= zeros(n_q, n_q);
end

% stiffening due to rotational velocity (centrifugal stiffening)
% force due to rotation results from cumulated moment of inertia I_cum
for i_elem= n_elem:-1:1
    I_cum{i_elem}= 0.5*m(i_elem) * (R{i_elem} + 0.5*DR{i_elem});

    if i_elem<n_elem
        I_cum{i_elem}= I_cum{i_elem} + I_cum{i_elem+1} + 0.5*m(i_elem+1)*(R{i_elem+1} - 0.5*DR{i_elem+1});
    end
end

K0omega= cell(3);
for a= 1:3
    for b= 1:3
        K0omega{a, b}= zeros(n_q, n_q);
    end
end
for i_elem= 1:n_elem
    for a= 1:3
        for b= 1:3
            w1= zeros(1, 3);
            w1(a)= 1;
            w2= zeros(1, 3);
            w2(b)= 1;
            K0omega{a, b}= K0omega{a, b} + (-crossmat(w1)*(crossmat(w2)*I_cum{i_elem}))'*DR{i_elem}*dPhi{i_elem}'*dPhi{i_elem};
        end
    end
end
sid.K0omega= K0omega; % save for separate inspection

%% transfert to SID struct
if isunix
    [~, user_name] = system('whoami');
elseif ispc
    [~, user_name] = system('echo %USERDOMAIN%\%USERNAME%');
end
user_name= strrep(user_name, newline, '');

sid.comment= sprintf('%d nodes, %d modes, generated by SF2SID MATLAB script on %s by %s', n_elem, n_q, datestr(now), user_name);

sid.refmod.mass= mass;
sid.refmod.nelastq= n_q;

for k= 1:n_q
    sid.refmod.ielastq{k}= sprintf('Eigen Mode %4d : %13f Hz', k, EF(k));
end

for i_elem= 1:n_elem
    sid.frame(i_elem).node= num2str(i_elem);
    sid.frame(i_elem).rframe= 'body ref';
    
    sid.frame(i_elem).origin= emptyTaylor(1, 3, 1, n_q, 0, 3);
    sid.frame(i_elem).origin.M0= R{i_elem};
    for l= 1:n_q
        sid.frame(i_elem).origin.M1(:, 1, l)= Phi{i_elem}(:, l);
    end
    sid.frame(i_elem).Phi= emptyTaylor(0, 3, n_q, n_q, 0, 3);
    sid.frame(i_elem).Phi.M0= Phi{i_elem};
    sid.frame(i_elem).Phi.M1= K0F{i_elem};
    
    sid.frame(i_elem).Psi= emptyTaylor(0, 3, n_q, n_q, 0, 3);
    sid.frame(i_elem).Psi.M0= dPhi{i_elem};
%     sid.frame(i_elem).Psi.M1= K0L{i_elem};

    sid.frame(i_elem).AP= emptyTaylor(1, 3, 3, n_q, 0, 3);
    sid.frame(i_elem).AP.M0= eye(3); % keep original orientation for now E{i_elem};
    for l= 1:n_q
        sid.frame(i_elem).AP.M1(:, :, l)= crossmat(dPhi{i_elem}(:,l));
    end
   
%     sid.frame(i_elem).sigma= emptyTaylor(1, 6, 1, n_q, 0, 3);
%     sid.frame(i_elem).sigma.M0= zeros(6, 1); % no pretension for now.
%     for l= 1:n_q
%         sid.frame(i_elem).sigma.M1(:, 1, l)= 
%     end
end

sid.md= emptyTaylor(1, 3, 1, n_q, 0, 3);
sid.md.M0= mc0;
for l= 1:n_q
    sid.md.M1(:, 1, l)= C1(:, l);
end

sid.I= emptyTaylor(1, 3, 3, n_q, 0, 2);
sid.I.M0= I0;
for l= 1:n_q
    sid.I.M1(:, :, l)= -C4{l} - C4{l}';
end

sid.Ct= emptyTaylor(1, n_q, 3, n_q, 0, 3);
sid.Ct.M0= C1';
sid.Ct.M1(:, 1, :)= K0t{1};
sid.Ct.M1(:, 2, :)= K0t{2};
sid.Ct.M1(:, 3, :)= K0t{3};

sid.Cr= emptyTaylor(1, n_q, 3, n_q, 0, 3);
sid.Cr.M0= C2';
sid.Cr.M1(:, 1, :)= Kr{1} + K0r{1};
sid.Cr.M1(:, 2, :)= Kr{2} + K0r{2};
sid.Cr.M1(:, 3, :)= Kr{3} + K0r{3};

sid.Me= emptyTaylor(0, n_q, n_q, 0, 0, 2);
sid.Me.M0= Me;

sid.Gr= emptyTaylor(0, 3, 3*n_q, n_q, 0, 3);

for l= 1:n_q
    sid.Gr.M0(:, 3*(l-1)+(1:3))= -2*C4{l};
end
for l= 1:n_q
    for k= 1:n_q
        sid.Gr.M1(:, 3*(l-1)+(1:3), k)= -2*C6{k, l};
    end
end

sid.Ge= emptyTaylor(0, n_q, 3*n_q, 0, 0, 3);
for l= 1:n_q
    sid.Ge.M0(:, 3*(l-1)+(1:3))= 2*[Kr{1}(:, l) Kr{2}(:, l) Kr{3}(:, l)];
end

sid.Oe= emptyTaylor(1, n_q, 6, n_q, 0, 3);
for l= 1:n_q
    sid.Oe.M0(l, :)= [C4{l}(1, 1) C4{l}(2, 2) C4{l}(3, 3) C4{l}(1, 2)+C4{l}(2, 1) C4{l}(2, 3)+C4{l}(3, 2) C4{l}(3, 1)+C4{l}(1, 3)];
end
for m= 1:6
    if m<4
        sid.Oe.M1(:, m, :)= Komega{m, m} + K0omega{m, m};
    else
        c= m-3;
        d= mod(c, 3)+1;
        sid.Oe.M1(:, m, :)= Komega{c, d}+Komega{c, d}' + K0omega{c, d}+K0omega{c, d}';
    end        
end

sid.Ke= emptyTaylor(0, n_q, n_q, 0, 0, 2);
sid.Ke.M0= Ke;
% sid.Ke.M1

% sid.ksigma= emptyTaylor(0, n_q, 1, n_q, 0, 0);
% sid.ksigma.M0= [];

sid.De= emptyTaylor(0, n_q, n_q, 0, 0, 0);
sid.De.M0= De;


function m= crossmat(v)
m= [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
