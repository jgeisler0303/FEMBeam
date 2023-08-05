%  This routine is used to compute rotor (blade and hub) properties
%    KBF(), KBE(), CBF(), CBE(), FreqBF(), FreqBE(), AxRedBld(),
%    TwistedSF(), BldMass(), FirstMom(), SecondMom(), BldCG(),
%    RotMass, RotIner, Hubg1Iner, Hubg2Iner, rSAerCenn1(), and
%    rSAerCenn2(), BElmtMass()
%  tower properties:
%    KTFA(), KTSS(), CTFA(), CTSS(), FreqTFA(), FreqTSS(),
%    AxRedTFA(), AxRedTSS(), TwrFASF(), TwrSSSF(), TwrMass, andTwistedSF
%    TwrTpMass, TElmtMass()
function p= FASTCoeff(fst_file)
fstDataOut = FAST2Matlab(fst_file);
fst_dir= fileparts(fst_file);

EDFile= strrep(GetFASTPar(fstDataOut, 'EDFile'), '"', '');
edDataOut = FAST2Matlab(fullfile(fst_dir, EDFile));

EDTwrFile= strrep(GetFASTPar(edDataOut, 'TwrFile'), '"', '');
edtwrDataOut = FAST2Matlab(fullfile(fst_dir, EDTwrFile));
EDBldFile= strrep(GetFASTPar(edDataOut, 'BldFile(1)'), '"', '');
edbldDataOut = FAST2Matlab(fullfile(fst_dir, EDBldFile));

p.BD4Blades= false; % always use ElastoDyn

p= SetPrimaryParameters( p, edDataOut);
p=SetBladeParameters( p, edDataOut, edbldDataOut);
p= SetTowerParameters( p, edtwrDataOut);

%...............................................................................................................................
% Calculate the distances from point S on a blade to the aerodynamic center in the j1 and j2 directions:
%...............................................................................................................................
% fields= fieldnames(p);

for K = 1:p.NumBl          % Loop through the blades

    for J = 1:p.BldNodes    % Loop through the blade nodes / elements

        TmpDist           = ( 0.25 - p.PitchAxis(K,J) )*p.Chord(J);   % Distance along the chordline from point S (25% chord) to the aerodynamic center of the blade element J--positive towards the trailing edge.
        TmpDistj1         = TmpDist*p.SAeroTwst(J);                   % Distance along the j1-axis   from point S (25% chord) to the aerodynamic center of the blade element J
        TmpDistj2         = TmpDist*p.CAeroTwst(J) ;                  % Distance along the j2-axis   from point S (25% chord) to the aerodynamic center of the blade element J
        p.rSAerCenn1(K,J) = TmpDistj1*p.CThetaS(K,J+1) - TmpDistj2*p.SThetaS(K,J+1);
        p.rSAerCenn2(K,J) = TmpDistj1*p.SThetaS(K,J+1) + TmpDistj2*p.CThetaS(K,J+1);

    end % J - Blade nodes / elements

end    % K - Blades


%...............................................................................................................................
% Calculate the nacelle inertia terms:
%...............................................................................................................................

p.Nacd2Iner = GetFASTPar(edDataOut, 'NacYIner') - p.NacMass*( p.NacCMxn^2 + p.NacCMyn^2 ); % Nacelle inertia about the d2-axis
if ( p.Nacd2Iner < 0.0 )
    error(' NacYIner must not be less than NacMass*( NacCMxn^2 + NacCMyn^2 ).');
end

% Calculate hub inertia about its centerline passing through its c.g..
%   This calculation assumes that the hub for a 2-blader is essentially
%   a uniform cylinder whose centerline is transverse through the cylinder
%   passing through its c.g..  That is, for a 2-blader, Hubg1Iner =
%   Hubg2Iner is the inertia of the hub about both the g1- and g2- axes.  For
%   3-bladers, Hubg1Iner is simply equal to HubIner and Hubg2Iner is zero.
% Also, Initialize RotMass and RotIner to associated hub properties:

if ( p.NumBl == 2 )   % 2-blader
    p.Hubg1Iner = ( InputFileData.HubIner - p.HubMass*( ( p.UndSling - p.HubCM )^2 ) )/( p.CosDel3^2 );
    p.Hubg2Iner = p.Hubg1Iner;
    if ( p.Hubg1Iner < 0.0 )
        error(' HubIner must not be less than HubMass*( UndSling - HubCM )^2 for 2-blader.')
    end
else                    % 3-blader
    p.Hubg1Iner = GetFASTPar(edDataOut, 'HubIner');
    p.Hubg2Iner = 0.0;
end

p.RotMass   = p.HubMass;
p.RotIner   = p.Hubg1Iner;


%...............................................................................................................................

% Initialize several variables to 0.0:

p.TwrMass = 0.0;


for K = 1:p.NumBl          % Loop through the blades


    % Initialize BldMass(), FirstMom(), and SecondMom() using TipMass() effects:

    p.BldMass  (K) = p.TipMass(K);
    p.FirstMom (K) = p.TipMass(K)*p.BldFlexL;
    p.SecondMom(K) = p.TipMass(K)*p.BldFlexL*p.BldFlexL;


    for J = p.BldNodes:-1:1 % Loop through the blade nodes / elements in reverse


        % Calculate the mass of the current element

        p.BElmntMass(J,K) = p.MassB(K,J)*p.DRNodes(J);                        % Mass of blade elementp.FMomAbvNd J


        % Integrate to find some blade properties which will be output in .fsm

        p.BldMass  (K) = p.BldMass  (K) + p.BElmntMass(J,K);
        p.FirstMom (K) = p.FirstMom (K) + p.BElmntMass(J,K)*p.RNodes(J);
        p.SecondMom(K) = p.SecondMom(K) + p.BElmntMass(J,K)*p.RNodes(J)*p.RNodes(J);


        % Integrate to find FMomAbvNd:

        p.FMomAbvNd   (K,J) = ( 0.5*p.BElmntMass(J,K) )*( p.HubRad + p.RNodes(J  ) + 0.5*p.DRNodes(J  ) );

        if ( J == p.BldNodes )   % Outermost blade element
            % Add the TipMass() effects:

            p.FMomAbvNd(K,J) = p.FMomAbvNd(K,J) + p.TipMass(K)*p.TipRad;
        else                       % All other blade elements
            % Add to p.FMomAbvNd(K,J) the effects from the (not yet used) portion of element J+1

            p.FMomAbvNd(K,J) = p.FMomAbvNd(K,J) + p.FMomAbvNd(K,J+1) ...
                + ( 0.5*p.BElmntMass(J+1,K) )*( p.HubRad + p.RNodes(J+1) - 0.5*p.DRNodes(J+1) );
        end


    end % J - Blade nodes / elements in reverse

    if (~p.BD4Blades)
        % Calculate BldCG() using FirstMom() and BldMass(); and calculate
        %   RotMass and RotIner:

        p.BldCG    (K) = p.FirstMom (K) / p.BldMass    (K);
        p.RotMass      = p.RotMass      + p.BldMass    (K);
        p.RotIner      = p.RotIner      + ( p.SecondMom(K) + p.BldMass  (K)*p.HubRad*( 2.0*p.BldCG(K) + p.HubRad ) )*( p.CosPreC(K)^2 );
    end

end % K - Blades


p.MBF= zeros(p.NumBl, 2, 2);
p.MBE= zeros(p.NumBl, 1, 1);
p.KBFCent= zeros(p.NumBl, 2, 2);
p.KBECent= zeros(p.NumBl, 1);
p.KBF= zeros(p.NumBl, 2, 2);
p.KBE= zeros(p.NumBl, 1);
p.TwistedSF= zeros(p.NumBl, 2, 3, p.BldNodes+2,3); % blade no, x/y, BF1/BF2/BE, node, deriv
for K = 1:p.NumBl          % Loop through the blades


    % Initialize the generalized blade masses using tip mass effects:

    p.MBF(K,1,1) = p.TipMass(K);
    p.MBF(K,2,2) = p.TipMass(K);
    p.MBE(K,1,1) = p.TipMass(K);

    for J = 1:p.BldNodes    % Loop through the blade nodes / elements


        % Integrate to find the generalized mass of the blade (including tip mass effects).
        %   Ignore the cross-correlation terms of MBF (i.e. MBF(i,j) where i ~= j) since
        %   these terms will never be used.

        p.Shape1(J) = SHP( p.RNodesNorm(J), p.BldFlexL, p.BldFl1Sh(:,K), 0);
        p.Shape2(J) = SHP( p.RNodesNorm(J), p.BldFlexL, p.BldFl2Sh(:,K), 0);
        p.MBF    (K,1,1) = p.MBF    (K,1,1) + p.BElmntMass(J,K)*p.Shape1(J)*p.Shape1(J);
        p.MBF    (K,2,2) = p.MBF    (K,2,2) + p.BElmntMass(J,K)*p.Shape2(J)*p.Shape2(J);

        p.Shape(J)  = SHP( p.RNodesNorm(J), p.BldFlexL, p.BldEdgSh(:,K), 0);
        p.MBE    (K,1,1) = p.MBE    (K,1,1) + p.BElmntMass(J,K)*p.Shape(J) *p.Shape(J);


        % Integrate to find the generalized stiffness of the blade (not including centrifugal
        %    effects).

        ElmntStff      = p.StiffBF(K,J)*p.DRNodes(J);                       % Flapwise stiffness of blade element J
        p.ddShape1(J) = SHP( p.RNodesNorm(J), p.BldFlexL, p.BldFl1Sh(:,K), 2);
        p.ddShape2(J) = SHP( p.RNodesNorm(J), p.BldFlexL, p.BldFl2Sh(:,K), 2);
        p.KBF    (K,1,1) = p.KBF    (K,1,1) + ElmntStff*p.ddShape1(J)*p.ddShape1(J);
        p.KBF    (K,1,2) = p.KBF    (K,1,2) + ElmntStff*p.ddShape1(J)*p.ddShape2(J);
        p.KBF    (K,2,1) = p.KBF    (K,2,1) + ElmntStff*p.ddShape2(J)*p.ddShape1(J);
        p.KBF    (K,2,2) = p.KBF    (K,2,2) + ElmntStff*p.ddShape2(J)*p.ddShape2(J);

        ElmntStff      = p.StiffBE(K,J)*p.DRNodes(J);                       % Edgewise stiffness of blade element J
        p.ddShape(J)  = SHP( p.RNodesNorm(J), p.BldFlexL, p.BldEdgSh(:,K), 2);
        p.KBE    (K,1,1) = p.KBE    (K,1,1) + ElmntStff*p.ddShape(J) *p.ddShape(J);


        % Integrate to find the centrifugal-term of the generalized flapwise and edgewise
        %   stiffness of the blades.  Ignore the cross-correlation terms of KBFCent (i.e.
        %   KBFCent(i,j) where i ~= j) since these terms will never be used.

        ElmntStff      = p.FMomAbvNd(K,J)*p.DRNodes(J)*p.RotSpeed^2;   % Centrifugal stiffness of blade element J

        p.dShape1(J) = SHP( p.RNodesNorm(J), p.BldFlexL, p.BldFl1Sh(:,K), 1);
        p.dShape2(J) = SHP( p.RNodesNorm(J), p.BldFlexL, p.BldFl2Sh(:,K), 1);
        p.KBFCent(K,1,1) = p.KBFCent(K,1,1) + ElmntStff*p.dShape1(J)*p.dShape1(J);
        p.KBFCent(K,2,2) = p.KBFCent(K,2,2) + ElmntStff*p.dShape2(J)*p.dShape2(J);

        p.dShape(J)  = SHP( p.RNodesNorm(J), p.BldFlexL, p.BldEdgSh(:,K), 1);
        p.KBECent(K,1,1) = p.KBECent(K,1,1) + ElmntStff*p.dShape(J) *p.dShape(J);


        % Calculate the 2nd derivatives of the twisted shape functions:

        Shape  = SHP( p.RNodesNorm(J), p.BldFlexL, p.BldFl1Sh(:,K), 2);
        p.TwistedSF(K,1,1,J+1,3) =  Shape*p.CThetaS(K,J+1);                  % 2nd deriv. of Phi1(J) for blade K
        p.TwistedSF(K,2,1,J+1,3) = -Shape*p.SThetaS(K,J+1);                  % 2nd deriv. of Psi1(J) for blade K

        Shape  = SHP( p.RNodesNorm(J), p.BldFlexL, p.BldFl2Sh(:,K), 2);
        p.TwistedSF(K,1,2,J+1,3) =  Shape*p.CThetaS(K,J+1);                  % 2nd deriv. of Phi2(J) for blade K
        p.TwistedSF(K,2,2,J+1,3) = -Shape*p.SThetaS(K,J+1);                  % 2nd deriv. of Psi2(J) for blade K

        Shape  = SHP( p.RNodesNorm(J), p.BldFlexL, p.BldEdgSh(:,K), 2);
        p.TwistedSF(K,1,3,J+1,3) =  Shape*p.SThetaS(K,J+1);                  % 2nd deriv. of Phi3(J) for blade K
        p.TwistedSF(K,2,3,J+1,3) =  Shape*p.CThetaS(K,J+1);                  % 2nd deriv. of Psi3(J) for blade K


        % Integrate to find the 1st derivatives of the twisted shape functions:
        TwstdSF= zeros(2, 3, 2);
        for I = 1:2     % Loop through Phi and Psi
            for L = 1:3  % Loop through all blade DOFs
                TwstdSF     (  I,L,  2) = p.TwistedSF(K,I,L,J+1,3)*0.5*p.DRNodes(J);
                p.TwistedSF   (K,I,L,J+1,2) = TwstdSF   ( I,L,  2);
            end       % L - All blade DOFs
        end          % I - Phi and Psi

        if ( J ~= 1 )    % All but the innermost blade element
            % Add the effects from the (not yet used) portion of element J-1

            for I = 1:2     % Loop through Phi and Psi
                for L = 1:3  % Loop through all blade DOFs
                    p.TwistedSF(K,I,L,J+1,2) = p.TwistedSF(K,I,L,J+1,2) + p.TwistedSF(K,I,L,J,2) ...
                        + TwstdSFOld( I,L,  2);
                end       % L - All blade DOFs
            end          % I - Phi and Psi
        end


        % Integrate to find the twisted shape functions themselves (i.e., their zeroeth derivative):

        for I = 1:2     % Loop through Phi and Psi
            for L = 1:3  % Loop through all blade DOFs
                TwstdSF     (  I,L, 1 ) = p.TwistedSF(K,I,L,J+1,2)*0.5*p.DRNodes(J);
                p.TwistedSF   (K,I,L,J+1,1) = TwstdSF   ( I,L, 1 );
            end       % L - All blade DOFs
        end          % I - Phi and Psi

        if ( J ~= 1 )    % All but the innermost blade element
            % Add the effects from the (not yet used) portion of element J-1

            for I = 1:2     % Loop through Phi and Psi
                for L = 1:3  % Loop through all blade DOFs
                    p.TwistedSF(K,I,L,J+1,1) = p.TwistedSF(K,I,L,J+1,1) + p.TwistedSF(K,I,L,J,1) ...
                        + TwstdSFOld( I,L, 1 );
                end       % L - All blade DOFs
            end          % I - Phi and Psi
        end


        % Integrate to find the blade axial reduction shape functions:
        AxRdBld= zeros(2, 3);
        for I = 1:3     % Loop through all blade DOFs
            for L = 1:3  % Loop through all blade DOFs
                AxRdBld    (  I,L  ) = 0.5*p.DRNodes(J)*(                          ...
                    p.TwistedSF(K,1,I,J+1,2)*p.TwistedSF(K,1,L,J+1,2) ...
                    + p.TwistedSF(K,2,I,J+1,2)*p.TwistedSF(K,2,L,J+1,2) );
                p.AxRedBld   (K,I,L,J+1) = AxRdBld(I,L);
            end       % L - All blade DOFs
        end          % I - All blade DOFs

        if ( J ~= 1 )    % All but the innermost blade element
            % Add the effects from the (not yet used) portion of element J-1

            for I = 1:3     % Loop through all blade DOFs
                for L = 1:3  % Loop through all blade DOFs
                    p.AxRedBld(K,I,L,J+1) = p.AxRedBld(K,I,L,J+1) + p.AxRedBld(K,I,L,J)   ...
                        + AxRdBldOld(I,L);
                end       % L - All blade DOFs
            end          % I - All blade DOFs
        end


        % Store the TwstdSF and AxRdBld terms of the current element (these will be used for the next element)

        TwstdSFOld = TwstdSF;
        AxRdBldOld = AxRdBld;


    end % J - Blade nodes / elements




    if (p.BD4Blades)

        %p.KBF     ( K,:,:    ) = 0.0

        % the 1st and zeroeth derivatives of the twisted shape functions at the blade root:
        p.TwistedSF(K,:,:,:,2) = 0.0;
        p.TwistedSF(K,:,:,:,1) = 0.0;
        p.AxRedBld( K,:,:,:  ) = 0.0;
    else

        % Apply the flapwise modal stiffness tuners of the blades to KBF():

        for I = 1:2     % Loop through flap DOFs
            for L = 1:2  % Loop through flap DOFs
                p.KBF(K,I,L) = sqrt( p.FStTunr(K,I)*p.FStTunr(K,L) )*p.KBF(K,I,L);
            end       % L - Flap DOFs
        end          % I - Flap DOFs

        % Calculate the blade natural frequencies:

        for I = 1:2     % Loop through flap DOFs
            p.FreqBF(K,I,1) = (1/2/pi)*sqrt(   p.KBF(K,I,I)                   /( p.MBF(K,I,I) - p.TipMass(K) ) );   % Natural blade I-flap frequency w/o centrifugal stiffening nor     tip mass effects
            p.FreqBF(K,I,2) = (1/2/pi)*sqrt(   p.KBF(K,I,I)                   /  p.MBF(K,I,I)                );     % Natural blade I-flap frequency w/o centrifugal stiffening, but w/ tip mass effects
            p.FreqBF(K,I,3) = (1/2/pi)*sqrt( ( p.KBF(K,I,I) + p.KBFCent(K,I,I) )/  p.MBF(K,I,I)                );     % Natural blade I-flap frequency w/  centrifugal stiffening and     tip mass effects
        end          % I - Flap DOFs

        p.FreqBE   (K,1,1) = (1/2/pi)*sqrt(   p.KBE(K,1,1)                   /( p.MBE(K,1,1) - p.TipMass(K) ) );   % Natural blade 1-edge frequency w/o centrifugal stiffening nor      tip mass effects
        p.FreqBE   (K,1,2) = (1/2/pi)*sqrt(   p.KBE(K,1,1)                   /  p.MBE(K,1,1)                );     % Natural Blade 1-edge frequency w/o  centrifugal stiffening, but w/ tip mass effects
        p.FreqBE   (K,1,3) = (1/2/pi)*sqrt( ( p.KBE(K,1,1) + p.KBECent(K,1,1) )/  p.MBE(K,1,1)                );     % Natural Blade 1-edge frequency w/  centrifugal stiffening and      tip mass effects


        % Calculate the generalized damping of the blades:

        for I = 1:2     % Loop through flap DOFs
            for L = 1:2  % Loop through flap DOFs
                p.CBF(K,I,L) = ( 0.01*p.BldFDamp(K,L) )*p.KBF(K,I,L)/( pi*p.FreqBF(K,L,1) );
            end       % L - Flap DOFs
        end          % I - Flap DOFs

        p.CBE      (K,1,1) = ( 0.01*p.BldEDamp(K,1) )*p.KBE(K,1,1)/( pi*p.FreqBE(K,1,1) );


        % Calculate the 2nd derivatives of the twisted shape functions at the blade root:

        Shape  = SHP( 0.0, p.BldFlexL, p.BldFl1Sh(:,K), 2);
        p.TwistedSF(K,1,1,1,3) =  Shape*p.CThetaS(K,1);        % 2nd deriv. of Phi1(0) for blade K
        p.TwistedSF(K,2,1,1,3) = -Shape*p.SThetaS(K,1);        % 2nd deriv. of Psi1(0) for blade K

        Shape  = SHP( 0.0, p.BldFlexL, p.BldFl2Sh(:,K), 2);
        p.TwistedSF(K,1,2,1,3) =  Shape*p.CThetaS(K,1);        % 2nd deriv. of Phi2(0) for blade K
        p.TwistedSF(K,2,2,1,3) = -Shape*p.SThetaS(K,1);        % 2nd deriv. of Psi2(0) for blade K

        Shape  = SHP( 0.0, p.BldFlexL, p.BldEdgSh(:,K), 2);
        p.TwistedSF(K,1,3,1,3) =  Shape*p.SThetaS(K,1);        % 2nd deriv. of Phi3(0) for blade K
        p.TwistedSF(K,2,3,1,3) =  Shape*p.CThetaS(K,1);        % 2nd deriv. of Psi3(0) for blade K


        % Calculate the 2nd derivatives of the twisted shape functions at the tip:

        Shape  = SHP( 1.0, p.BldFlexL, p.BldFl1Sh(:,K), 2);
        p.TwistedSF(K,1,1,p.TipNode+1,3) =  Shape*p.CThetaS(K,p.TipNode+1);        % 2nd deriv. of Phi1(p.TipNode) for blade K
        p.TwistedSF(K,2,1,p.TipNode+1,3) = -Shape*p.SThetaS(K,p.TipNode+1);        % 2nd deriv. of Psi1(p.TipNode) for blade K

        Shape  = SHP( 1.0, p.BldFlexL, p.BldFl2Sh(:,K), 2);
        p.TwistedSF(K,1,2,p.TipNode+1,3) =  Shape*p.CThetaS(K,p.TipNode+1);        % 2nd deriv. of Phi2(p.TipNode) for blade K
        p.TwistedSF(K,2,2,p.TipNode+1,3) = -Shape*p.SThetaS(K,p.TipNode+1);        % 2nd deriv. of Psi2(p.TipNode) for blade K

        Shape  = SHP( 1.0, p.BldFlexL, p.BldEdgSh(:,K), 2);
        p.TwistedSF(K,1,3,p.TipNode+1,3) =  Shape*p.SThetaS(K,p.TipNode+1);        % 2nd deriv. of Phi3(p.TipNode) for blade K
        p.TwistedSF(K,2,3,p.TipNode+1,3) =  Shape*p.CThetaS(K,p.TipNode+1);        % 2nd deriv. of Psi3(p.TipNode) for blade K


        % Integrate to find the 1st and zeroeth derivatives of the twisted shape functions
        %   at the tip:

        for I = 1:2     % Loop through Phi and Psi
            for L = 1:3  % Loop through all blade DOFs
                p.TwistedSF(K,I,L,p.TipNode+1,2) = p.TwistedSF(K,I,L,p.BldNodes+1,2) + TwstdSFOld(I,L,2);
                p.TwistedSF(K,I,L,p.TipNode+1,1) = p.TwistedSF(K,I,L,p.BldNodes+1,1) + TwstdSFOld(I,L,1);
            end       % L - All blade DOFs
        end          % I - Phi and Psi

        % the 1st and zeroeth derivatives of the twisted shape functions at the blade root:
        p.TwistedSF(K,:,:,1,2) = 0.0;
        p.TwistedSF(K,:,:,1,1) = 0.0;
        p.AxRedBld( K,:,:,1  ) = 0.0;

        % Integrate to find the blade axial reduction shape functions at the tip:

        for I = 1:3     % Loop through all blade DOFs
            for L = 1:3  % Loop through all blade DOFs
                p.AxRedBld(K,I,L,p.TipNode+1) = p.AxRedBld(K,I,L,p.BldNodes+1) + AxRdBldOld(I,L);
            end       % L - All blade DOFs
        end          % I - All blade DOFs
    end % p.BD4Blades


end % K - Blades

% setdiff(fieldnames(p), fields)
% fields= fieldnames(p);
% Calculate the tower-top mass:

% p.TwrTpMass = p.RotMass + p.RFrlMass + p.BoomMass + p.TFinMass + p.NacMass + p.YawBrMass;
p.TwrTpMass = p.RotMass + p.NacMass + p.YawBrMass;


for J = p.TwrNodes:-1:1 % Loop through the tower nodes / elements in reverse


    % Calculate the mass of the current element

    p.TElmntMass(J)    = p.MassT(J)*p.DHNodes(J);     % Mass of tower element J


    % Integrate to find the tower mass which will be output in .fsm

    p.TwrMass      = p.TwrMass + p.TElmntMass(J);


    % Integrate to find TMssAbvNd:

    TMssAbvNd   (J) = 0.5*p.TElmntMass(J);

    if ( J == p.TwrNodes )   % Uppermost tower element
        % Add the TwrTpMass effects:

%         TMssAbvNd(J) = TMssAbvNd(J) + p.TwrTpMass;
        TMssAbvNd(J) = TMssAbvNd(J);
    else                       % All other tower elements
        % Add to TMssAbvNd(J) the effects from the (not yet used) portion of element J+1

        TMssAbvNd(J) = 0.5*p.TElmntMass(J+1) + TMssAbvNd(J) + TMssAbvNd(J+1);
    end


end % J - Tower nodes / elements in reverse



% Initialize the generalized tower masses using tower-top mass effects:
p.MTFA= zeros(2, 2);
p.MTSS= zeros(2, 2);
for I = 1:2  % Loop through all tower modes in a single direction
    p.MTFA(I,I) = p.TwrTpMass;
    p.MTSS(I,I) = p.TwrTpMass;
end       % I - All tower modes in a single direction

% set values for tower base (note that we haven't corrctly defined the values for (:,0,2) in the arrays below):
p.TwrFASF(   1:2,1,1:2) = 0.0;
p.TwrSSSF(   1:2,1,1:2) = 0.0;
p.AxRedTFA(1:2,1:2,1)   = 0.0;
p.AxRedTSS(1:2,1:2,1)   = 0.0;

p.KTFA= zeros(2, 2);
p.KTSS= zeros(2, 2);
p.KTFAGrav= zeros(2, 2);
p.KTSSGrav= zeros(2, 2);
p.KTFAGravTT= zeros(2, 2);
p.KTSSGravTT= zeros(2, 2);

for J = 1:p.TwrNodes    % Loop through the tower nodes / elements


    % Calculate the tower shape functions (all derivatives):

    p.TwrFASF(1,J+1,3) = SHP( p.HNodesNorm(J), p.TwrFlexL, p.TwFAM1Sh(:), 2);
    p.TwrFASF(2,J+1,3) = SHP( p.HNodesNorm(J), p.TwrFlexL, p.TwFAM2Sh(:), 2);
    p.TwrFASF(1,J+1,2) = SHP( p.HNodesNorm(J), p.TwrFlexL, p.TwFAM1Sh(:), 1);
    p.TwrFASF(2,J+1,2) = SHP( p.HNodesNorm(J), p.TwrFlexL, p.TwFAM2Sh(:), 1);
    p.TwrFASF(1,J+1,1) = SHP( p.HNodesNorm(J), p.TwrFlexL, p.TwFAM1Sh(:), 0);
    p.TwrFASF(2,J+1,1) = SHP( p.HNodesNorm(J), p.TwrFlexL, p.TwFAM2Sh(:), 0);

    p.TwrSSSF(1,J+1,3) = SHP( p.HNodesNorm(J), p.TwrFlexL, p.TwSSM1Sh(:), 2);
    p.TwrSSSF(2,J+1,3) = SHP( p.HNodesNorm(J), p.TwrFlexL, p.TwSSM2Sh(:), 2);
    p.TwrSSSF(1,J+1,2) = SHP( p.HNodesNorm(J), p.TwrFlexL, p.TwSSM1Sh(:), 1);
    p.TwrSSSF(2,J+1,2) = SHP( p.HNodesNorm(J), p.TwrFlexL, p.TwSSM2Sh(:), 1);
    p.TwrSSSF(1,J+1,1) = SHP( p.HNodesNorm(J), p.TwrFlexL, p.TwSSM1Sh(:), 0);
    p.TwrSSSF(2,J+1,1) = SHP( p.HNodesNorm(J), p.TwrFlexL, p.TwSSM2Sh(:), 0);


    % Integrate to find the generalized mass of the tower (including tower-top mass effects).
    %   Ignore the cross-correlation terms of p.MTFA (i.e. p.MTFA(i,j) where i ~= j) and p.MTSS
    %   since these terms will never be used.


    for I = 1:2     % Loop through all tower DOFs in one direction
        p.MTFA  (I,I) = p.MTFA  (I,I) + p.TElmntMass(J)*p.TwrFASF(I,J+1,1)^2;
        p.MTSS  (I,I) = p.MTSS  (I,I) + p.TElmntMass(J)*p.TwrSSSF(I,J+1,1)^2;
    end          % I - through all tower DOFs in one direction


    % Integrate to find the generalized stiffness of the tower (not including gravitational
    %    effects).

    ElStffFA       = p.StiffTFA(J)*p.DHNodes(J);                        % Fore-aft stiffness of tower element J
    ElStffSS       = p.StiffTSS(J)*p.DHNodes(J);                        % Side-to-side stiffness of tower element J

    for I = 1:2     % Loop through all tower DOFs in one direction
        for L = 1:2  % Loop through all tower DOFs in one direction
            p.KTFA (I,L) = p.KTFA    (I,L) + ElStffFA *p.TwrFASF(I,J+1,3)*p.TwrFASF(L,J+1,3);
            p.KTSS (I,L) = p.KTSS    (I,L) + ElStffSS *p.TwrSSSF(I,J+1,3)*p.TwrSSSF(L,J+1,3);
        end       % L - All tower DOFs in one direction
    end          % I - through all tower DOFs in one direction


    % Integrate to find the gravitational-term of the generalized stiffness of the tower.
    %   Ignore the cross-correlation terms of p.KTFAGrav (i.e. p.KTFAGrav(i,j) where i ~= j)
    %   and p.KTSSGrav since these terms will never be used.

%     ElmntStff      = -TMssAbvNd(J)*p.DHNodes(J)*p.Gravity;              % Gravitational stiffness of tower element J
    ElmntStff      = TMssAbvNd(J)*p.DHNodes(J)*1;              % stiffness of tower element J due to unity acceleration

    for I = 1:2     % Loop through all tower DOFs in one direction
        p.KTFAGrav(I,I) = p.KTFAGrav(I,I) + ElmntStff*p.TwrFASF(I,J+1,2)^2;
        p.KTSSGrav(I,I) = p.KTSSGrav(I,I) + ElmntStff*p.TwrSSSF(I,J+1,2)^2;
    end
    ElmntStff      =  1*p.DHNodes(J)*1;              % due to unitiy mass under unity acceleration stiffness of tower element J due to unity acceleration
    for I= 1:2
        p.KTFAGravTT(I,I) = p.KTFAGravTT(I,I)+ElmntStff*p.TwrFASF(I,J+1,2)^2;
        p.KTSSGravTT(I,I) = p.KTSSGravTT(I,I)+ElmntStff*p.TwrSSSF(I,J+1,2)^2;
    end


    % Integrate to find the tower axial reduction shape functions:
    AxRdTFA= zeros(2, 2);
    AxRdTSS= zeros(2, 2);
    for I = 1:2     % Loop through all tower DOFs in one direction
        for L = 1:2  % Loop through all tower DOFs in one direction
            AxRdTFA (I,L) = 0.5*p.DHNodes(J)*p.TwrFASF(I,J+1,2)*p.TwrFASF(L,J,2);
            AxRdTSS (I,L) = 0.5*p.DHNodes(J)*p.TwrSSSF(I,J+1,2)*p.TwrSSSF(L,J,2);

            p.AxRedTFA(I,L,J+1) = AxRdTFA(I,L);
            p.AxRedTSS(I,L,J+1) = AxRdTSS(I,L);
        end       % L - All tower DOFs in one direction
    end

    if ( J ~= 1 )    % All but the lowermost tower element
        % Add the effects from the (not yet used) portion of element J-1

        for I = 1:2     % Loop through all tower DOFs in one direction
            for L = 1:2  % Loop through all tower DOFs in one direction
                p.AxRedTFA(I,L,J+1) = p.AxRedTFA(I,L,J+1) + p.AxRedTFA(I,L,J)+ AxRdTFAOld(I,L);
                p.AxRedTSS(I,L,J+1) = p.AxRedTSS(I,L,J+1) + p.AxRedTSS(I,L,J)+ AxRdTSSOld(I,L);
            end       % L - All tower DOFs in one direction
        end
    end


    % Store the AxRdTFA and AxRdTSS terms of the current element (these will be used for the next element)

    AxRdTFAOld = AxRdTFA;
    AxRdTSSOld = AxRdTSS;


end % J - Tower nodes / elements


% Apply the modal stiffness tuners of the tower to KTFA() and KTSS():

for I = 1:2     % Loop through all tower DOFs in one direction
    for L = 1:2  % Loop through all tower DOFs in one direction
        p.KTFA(I,L) = sqrt( p.FAStTunr(I)*p.FAStTunr(L) )*p.KTFA(I,L);

        p.KTSS(I,L) = sqrt( p.SSStTunr(I)*p.SSStTunr(L) )*p.KTSS(I,L);
    end       % L - All tower DOFs in one direction
end          % I - through all tower DOFs in one direction


% Calculate the tower natural frequencies:

for I = 1:2     % Loop through all tower DOFs in one direction
    allKTFAGrav= (p.KTFAGrav(I,I) + p.TwrTpMass*p.KTFAGravTT(I,I))*p.Gravity;
    allKTSSGrav= (p.KTSSGrav(I,I) + p.TwrTpMass*p.KTSSGravTT(I,I))*p.Gravity;
    p.FreqTFA(I,1) = (1/2/pi)*sqrt(   p.KTFA(I,I)                  /( p.MTFA(I,I) - p.TwrTpMass ) );  % Natural tower I-fore-aft frequency w/o gravitational destiffening nor tower-top mass effects
    p.FreqTFA(I,2) = (1/2/pi)*sqrt( ( p.KTFA(I,I) + allKTFAGrav )/  p.MTFA(I,I)               );  % Natural tower I-fore-aft frequency w/  gravitational destiffening and tower-top mass effects
    p.FreqTSS(I,1) = (1/2/pi)*sqrt(   p.KTSS(I,I)                  /( p.MTSS(I,I) - p.TwrTpMass ) );  % Natural tower I-side-to-side frequency w/o gravitational destiffening nor tower-top mass effects
    p.FreqTSS(I,2) = (1/2/pi)*sqrt( ( p.KTSS(I,I) + allKTSSGrav )/  p.MTSS(I,I)               );  % Natural tower I-side-to-side frequency w/  gravitational destiffening and tower-top mass effects
end          % I - All tower DOFs in one direction

for I = 1:2  % Loop through all tower modes in a single direction
    p.MTFA(I,I) = p.MTFA(I,I) - p.TwrTpMass;
    p.MTSS(I,I) = p.MTSS(I,I) - p.TwrTpMass;
end       % I - All tower modes in a single direction

% Calculate the generalized damping of the tower:

for I = 1:2     % Loop through all tower DOFs in one direction
    for L = 1:2  % Loop through all tower DOFs in one direction
        p.CTFA(I,L) = ( 0.01*p.TwrFADmp(L) )*p.KTFA(I,L)/( pi*p.FreqTFA(L,1) );

        p.CTSS(I,L) = ( 0.01*p.TwrSSDmp(L) )*p.KTSS(I,L)/( pi*p.FreqTSS(L,1) );
    end       % L - All tower DOFs in one direction
end          % I - All tower DOFs in one direction


% Calculate the tower shape functions (all derivatives) at the tower-top:

p.TwrFASF(1,p.TTopNode+1,3) = SHP( 1.0, p.TwrFlexL, p.TwFAM1Sh(:), 2);
p.TwrFASF(2,p.TTopNode+1,3) = SHP( 1.0, p.TwrFlexL, p.TwFAM2Sh(:), 2);
p.TwrFASF(1,p.TTopNode+1,2) = SHP( 1.0, p.TwrFlexL, p.TwFAM1Sh(:), 1);
p.TwrFASF(2,p.TTopNode+1,2) = SHP( 1.0, p.TwrFlexL, p.TwFAM2Sh(:), 1);
p.TwrFASF(1,p.TTopNode+1,1) = SHP( 1.0, p.TwrFlexL, p.TwFAM1Sh(:), 0);
p.TwrFASF(2,p.TTopNode+1,1) = SHP( 1.0, p.TwrFlexL, p.TwFAM2Sh(:), 0);

p.TwrSSSF(1,p.TTopNode+1,3) = SHP( 1.0, p.TwrFlexL, p.TwSSM1Sh(:), 2);
p.TwrSSSF(2,p.TTopNode+1,3) = SHP( 1.0, p.TwrFlexL, p.TwSSM2Sh(:), 2);
p.TwrSSSF(1,p.TTopNode+1,2) = SHP( 1.0, p.TwrFlexL, p.TwSSM1Sh(:), 1);
p.TwrSSSF(2,p.TTopNode+1,2) = SHP( 1.0, p.TwrFlexL, p.TwSSM2Sh(:), 1);
p.TwrSSSF(1,p.TTopNode+1,1) = SHP( 1.0, p.TwrFlexL, p.TwSSM1Sh(:), 0);
p.TwrSSSF(2,p.TTopNode+1,1) = SHP( 1.0, p.TwrFlexL, p.TwSSM2Sh(:), 0);


% Integrate to find the tower axial reduction shape functions at the tower-top:

for I = 1:2     % Loop through all tower DOFs in one direction
    for L = 1:2  % Loop through all tower DOFs in one direction
        p.AxRedTFA(I,L,p.TTopNode+1) = p.AxRedTFA(I,L,p.TwrNodes+1)+ AxRdTFAOld(I,L);
        p.AxRedTSS(I,L,p.TTopNode+1) = p.AxRedTSS(I,L,p.TwrNodes+1)+ AxRdTSSOld(I,L);
    end       % L - All tower DOFs in one direction
end


% Calculate the turbine mass:

p.TurbMass  = p.TwrTpMass + p.TwrMass;

% setdiff(fieldnames(p), fields)

end

%----------------------------------------------------------------------------------------------------------------------------------
%> SHP calculates the Derive-derivative of the shape function ModShpAry at Fract.
%  NOTE: This function only works for Deriv = 0, 1, or 2.
function shp= SHP(Fract, FlexL, ModShpAry, Deriv)
%..................................................................................................................................

if ( Deriv < 0 || Deriv > 2 )
    error('Function SHP input Deriv=%d is invalid. Deriv must be 0, 1, or 2.', Deriv)
elseif ( Fract < 0.0 || Fract > 1.0 )
    error('Function SHP input Fract=%d does not meet the condition 0<=Fract<=1.', Fract)
end

Swtch        = zeros(3, 1); % Initialize Swtch(:) to 0
Swtch(Deriv+1) = 1;
shp          = 0.0;

for I = 1:size(ModShpAry, 1) % =2,PolyOrd
    J = I + 1;
    CoefTmp = Swtch(1) + Swtch(2)*J + Swtch(3)*I*J;

    if ( (J == 2) && (Deriv == 2) )  %bjj this could be removed as Fract^0 = 1 (0^0 = 1 in Fortran)
        shp =       ModShpAry(I)*CoefTmp                         /( FlexL^Deriv );
    else
        shp = shp + ModShpAry(I)*CoefTmp*( Fract^( J - Deriv ) )/( FlexL^Deriv );
    end
end %I
end

%----------------------------------------------------------------------------------------------------------------------------------
%> This takes the blade input file data and sets the corresponding blade parameters, performing linear interpolation of the
%  input data to the specified blade mesh.
%  This routine assumes p\%HubRad and p\%BldFlexL are already set.
function p=SetBladeParameters( p, edDataOut, edbldDataOut)
%..................................................................................................................................



% ..............................................................................................................................
% Set the blade discretization information here:
% ..............................................................................................................................

for K=1:1 % we're going to assume the discretization is the same for all blades

    if (p.BD4Blades)
        p.BldNodes = 0;
    else
        p.BldNodes = GetFASTPar(edDataOut, 'BldNodes');
    end

    p.TipNode  = p.BldNodes + 1;    % The index for the blade tip and tower top nodes

end

% .......... Allocate arrays for the blade parameters being set in this routine ..........:



if ( ~ p.BD4Blades)

    for K=1:1 % we're going to assume the discretization is the same for all blades

%         if ( allocated( BladeMeshData(K).Chord ) )
% 
%             p.RNodes   = BladeMeshData(K).RNodes - p.HubRad;   % Radius to blade analysis nodes relative to root ( 0 < RNodes(:) < p.BldFlexL ) (Convert RNodes to be relative to the hub)
% 
%             p.DRNodes(1) = 2.0*p.RNodes(1);
%             for J = 2:p.BldNodes
%                 p.DRNodes(J) = 2.0*( p.RNodes(J) - p.RNodes(J-1) ) - p.DRNodes(J-1);
%             end
% 
%             p.Chord     = BladeMeshData(K).Chord;
%             p.AeroTwst  = BladeMeshData(K).AeroTwst;
%             p.CAeroTwst = cosd(p.AeroTwst);
%             p.SAeroTwst = sind(p.AeroTwst);
% 
%         else


            % DRNodes (Let's use constant-spaced nodes for now, but the rest of the code is written to handle variable-spaced nodes--
            %          this will be a future input%):
            p.DRNodes = ones(p.BldNodes, 1)*p.BldFlexL/p.BldNodes; %array

            % RNodes:
            p.RNodes(1) = 0.5*p.DRNodes(1);
            for J=2:p.BldNodes
                p.RNodes(J) = p.RNodes( J - 1 ) + 0.5*( p.DRNodes(J) + p.DRNodes( J - 1 ) );
            end

            % these values aren't used (at least they shouldn't be):
            p.Chord     = zeros(p.BldNodes, 1);
            p.AeroTwst  = zeros(p.BldNodes, 1);
            p.CAeroTwst = ones(p.BldNodes, 1);
            p.SAeroTwst = zeros(p.BldNodes, 1);

%         end


    end


    % ..............................................................................................................................
    % Interpolate the blade properties to this discretization:
    % ..............................................................................................................................

    % Array definitions:

    %    Input      Interp    Description
    %    -----      ------    -----------
    %    p.BlFract    RNodesNorm Fractional radius (0 at root, 1 at tip)
    %    PitchAx    PitchAxis  Pitch axis (0 at LE, 1 at TE)
    %    StrcTwst   ThetaS     Structural twist
    %    BMassDen   MassB      Lineal mass density
    %    FlpStff    StiffBF    Flapwise stiffness
    %    EdgStff    StiffBE    Edgewise stiffness
    %    GJStff     StiffBGJ   Blade torsional stiffness
    %    EAStff     StiffBEA   Blade extensional stiffness
    %    Alpha      BAlpha     Blade flap/twist coupling coefficient
    %    FlpIner    InerBFlp   Blade flap (about local structural yb-axis) mass inertia per unit length
    %    EdgIner    InerBEdg   Blade edge (about local structural xb-axis) mass inertia per unit length
    %    PrecrvRef  RefAxisxb  Blade offset for defining the reference axis from the pitch axis for precurved blades (along xb-axis)
    %    PreswpRef  RefAxisyb  Blade offset for defining the reference axis from the pitch axis for preswept  blades (along yb-axis)
    %    FlpcgOf    cgOffBFlp  Blade flap mass cg offset
    %    EdgcgOf    cgOffBEdg  Blade edge mass cg offset
    %    FlpEAOf    EAOffBFlp  Blade flap elastic axis offset
    %    EdgEAOf    EAOffBEdg  Blade edge elastic axis offset


    % Define RNodesNorm() which is common to all the blades:

    p.RNodesNorm = p.RNodes/p.BldFlexL;  % Normalized radius to analysis nodes relative to hub ( 0 < RNodesNorm(:) < 1 )



    % Perform a linear interpolation of the input data to map to the meshed data for simulation:

    for K=1:p.NumBl
        p.BlFract= getFASTTableColumn(edbldDataOut.BldProp, 'BlFract');

        StrcTwst= getFASTTableColumn(edbldDataOut.BldProp, 'StrcTwst');
        p.ThetaS  (K,1)         = StrcTwst(1);
        p.ThetaS  (K,p.TipNode+1) = StrcTwst(end);

        p.PitchAxis(K, :)= interp1(p.BlFract, getFASTTableColumn(edbldDataOut.BldProp, 'PitchAxis'), p.RNodesNorm);
        p.ThetaS(K, 2:end-1)= interp1(p.BlFract, StrcTwst, p.RNodesNorm);
        p.MassB(K, :)= interp1(p.BlFract, getFASTTableColumn(edbldDataOut.BldProp, 'BMassDen'), p.RNodesNorm);
        p.StiffBF(K, :)= interp1(p.BlFract, getFASTTableColumn(edbldDataOut.BldProp, 'FlpStff'), p.RNodesNorm);
        p.StiffBE(K, :)= interp1(p.BlFract, getFASTTableColumn(edbldDataOut.BldProp, 'EdgStff'), p.RNodesNorm);


        % Set the blade damping and stiffness tuner
        p.BldFDamp(K,1) = GetFASTPar(edbldDataOut, 'BldFlDmp(1)');
        p.BldFDamp(K,2) = GetFASTPar(edbldDataOut, 'BldFlDmp(2)');
        p.BldEDamp(K,1) = GetFASTPar(edbldDataOut, 'BldEdDmp(1)');
        p.FStTunr (K,1) = GetFASTPar(edbldDataOut, 'FlStTunr(1)');
        p.FStTunr (K,2) = GetFASTPar(edbldDataOut, 'FlStTunr(2)');



        % Set the mode shape arrays
        for i= 2:6
            p.BldEdgSh(i-1,K) = GetFASTPar(edbldDataOut, sprintf('BldEdgSh(%d)', i));
            p.BldFl1Sh(i-1,K) = GetFASTPar(edbldDataOut, sprintf('BldFl1Sh(%d)', i));
            p.BldFl2Sh(i-1,K) = GetFASTPar(edbldDataOut, sprintf('BldFl2Sh(%d)', i));
        end


    end % ( Blades )


else

    p.ThetaS  = 0.0;

    % Set the blade damping and stiffness tuner
    p.BldFDamp = 0.0;
    p.BldEDamp = 0.0;
    p.FStTunr  = 0.0;

    % Set the mode shape arrays
    p.BldEdgSh = 0.0;
    p.BldFl1Sh = 0.0;
    p.BldFl2Sh = 0.0;

end

p.CThetaS = cosd(p.ThetaS);
p.SThetaS = sind(p.ThetaS);
end

%----------------------------------------------------------------------------------------------------------------------------------
%> This takes the tower input file data and sets the corresponding tower parameters, performing linear interpolation of the
%  input data to the specified tower mesh.
%  It requires p\%TwrFlexL, and p\%TwrNodes to be set first.
function p= SetTowerParameters( p, edtwrDataOut)
%..................................................................................................................................

%...............................................................................................................................
% Define the tower discretization arrays:
%...............................................................................................................................

% DHNodes (Let's use constant-spaced nodes for now, but the rest of the code is written to handle variable-spaced nodes--
%          this will be a future input%):
p.DHNodes = ones(p.TwrNodes, 1)*p.TwrFlexL/p.TwrNodes;

% HNodes:
p.HNodes(1) = 0.5*p.DHNodes(1);
for J=2:p.TwrNodes
    p.HNodes(J) = p.HNodes( J - 1 ) + 0.5*( p.DHNodes(J) + p.DHNodes( J - 1 ) );
end

% HNodesNorm:
p.HNodesNorm = p.HNodes/p.TwrFlexL;


%...............................................................................................................................
% Interpolate the input data to the tower discretization
%...............................................................................................................................
% Array definitions:

%    Input      Interp    Description
%    -----      ------    -----------
%    HtFract    HNodesNorm Fractional height (0 at top of rigid section, 1 at tower top)
%    TMassDen   MassT      Lineal mass density
%    TwFAStif   StiffTFA   Tower fore-aft stiffness
%    TwSSStif   StiffTSS   Tower side-to-side stiffness
%    TwGJStif   StiffTGJ   Tower torsional stiffness
%    TwEAStif   StiffTEA   Tower extensional stiffness
%    TwFAIner   InerTFA    Tower fore-aft (about yt-axis) mass inertia per unit length
%    TwSSIner   InerTSS    Tower side-to-side (about xt-axis) mass inertia per unit length
%    TwFAcgOf   cgOffTFA   Tower fore-aft mass cg offset
%    TwSScgOf   cgOffTSS   Tower side-to-side mass cg offset

HtFract= getFASTTableColumn(edtwrDataOut.TowProp, 'HtFract');
p.MassT= interp1(HtFract, getFASTTableColumn(edtwrDataOut.TowProp, 'TMassDen'), p.HNodesNorm);
p.StiffTFA= interp1(HtFract, getFASTTableColumn(edtwrDataOut.TowProp, 'TwFAStif'), p.HNodesNorm);
p.StiffTSS= interp1(HtFract, getFASTTableColumn(edtwrDataOut.TowProp, 'TwSSStif'), p.HNodesNorm);


% if ( SetAdmVals )            % An ADAMS model will be created; thus, read in all the cols.
%     for J=1:p.TwrNodes
%         p.StiffTGJ  (J) = InterpStp( p.HNodesNorm(J), InputFileData.HtFract, InputFileData.TwGJStif, InterpInd, InputFileData.NTwInpSt );
%         p.StiffTEA  (J) = InterpStp( p.HNodesNorm(J), InputFileData.HtFract, InputFileData.TwEAStif, InterpInd, InputFileData.NTwInpSt );
%         p.InerTFA   (J) = InterpStp( p.HNodesNorm(J), InputFileData.HtFract, InputFileData.TwFAIner, InterpInd, InputFileData.NTwInpSt );
%         p.InerTSS   (J) = InterpStp( p.HNodesNorm(J), InputFileData.HtFract, InputFileData.TwSSIner, InterpInd, InputFileData.NTwInpSt );
%         p.cgOffTFA  (J) = InterpStp( p.HNodesNorm(J), InputFileData.HtFract, InputFileData.TwFAcgOf, InterpInd, InputFileData.NTwInpSt );
%         p.cgOffTSS  (J) = InterpStp( p.HNodesNorm(J), InputFileData.HtFract, InputFileData.TwSScgOf, InterpInd, InputFileData.NTwInpSt );
%     end % J
% end


%...............................................................................................................................
% Set other tower parameters:
%...............................................................................................................................

p.TTopNode = p.TwrNodes + 1;

%   % these are for HydroDyn ?
%p.DiamT(:) = InputFileData.TwrDiam
%p.CAT(:)   = InputFileData.TwrCA
%p.CDT(:)   = InputFileData.TwrCD
%


for i= 2:6
    p.TwFAM1Sh(i-1) = GetFASTPar(edtwrDataOut, sprintf('TwFAM1Sh(%d)', i));
    p.TwFAM2Sh(i-1) = GetFASTPar(edtwrDataOut, sprintf('TwFAM2Sh(%d)', i));
    p.TwSSM1Sh(i-1) = GetFASTPar(edtwrDataOut, sprintf('TwSSM1Sh(%d)', i));
    p.TwSSM2Sh(i-1) = GetFASTPar(edtwrDataOut, sprintf('TwSSM2Sh(%d)', i));
end

p.FAStTunr(1)= GetFASTPar(edtwrDataOut, 'FAStTunr(1)');
p.FAStTunr(2)= GetFASTPar(edtwrDataOut, 'FAStTunr(2)');
p.SSStTunr(1)= GetFASTPar(edtwrDataOut, 'SSStTunr(1)');
p.SSStTunr(2)= GetFASTPar(edtwrDataOut, 'SSStTunr(2)');

p.TwrFADmp(1)= GetFASTPar(edtwrDataOut, 'TwrFADmp(1)');
p.TwrFADmp(2)= GetFASTPar(edtwrDataOut, 'TwrFADmp(2)');
p.TwrSSDmp(1)= GetFASTPar(edtwrDataOut, 'TwrSSDmp(1)');
p.TwrSSDmp(2)= GetFASTPar(edtwrDataOut, 'TwrSSDmp(2)');

end
%----------------------------------------------------------------------------------------------------------------------------------
%> This takes the primary input file data and sets the corresponding parameters.
function p= SetPrimaryParameters( p, edDataOut)
%..................................................................................................................................


%p.Twr2Shft  = InputFileData.Twr2Shft
%p.HubIner   = InputFileData.HubIner
%p.NacYIner  = InputFileData.NacYIner

%...............................................................................................................................
% Direct copy of variables:
%...............................................................................................................................
p.NumBl     = GetFASTPar(edDataOut, 'NumBl');
p.TipRad    = GetFASTPar(edDataOut, 'TipRad');
p.HubRad    = GetFASTPar(edDataOut, 'HubRad');
p.method    = GetFASTPar(edDataOut, 'method');
p.TwrNodes  = GetFASTPar(edDataOut, 'TwrNodes');

p.PtfmCMxt = GetFASTPar(edDataOut, 'PtfmCMxt');
p.PtfmCMyt = GetFASTPar(edDataOut, 'PtfmCMyt');

p.DT        = GetFASTPar(edDataOut, 'DT');
p.Gravity   = GetFASTPar(edDataOut, 'Gravity');
p.OverHang  = GetFASTPar(edDataOut, 'OverHang');
p.ShftGagL  = GetFASTPar(edDataOut, 'ShftGagL');
p.TowerHt   = GetFASTPar(edDataOut, 'TowerHt');
p.TowerBsHt = GetFASTPar(edDataOut, 'TowerBsHt');
p.PtfmRefzt = GetFASTPar(edDataOut, 'PtfmRefzt');

p.HubMass   = GetFASTPar(edDataOut, 'HubMass');
p.GenIner   = GetFASTPar(edDataOut, 'GenIner');
p.NacMass   = GetFASTPar(edDataOut, 'NacMass');
p.YawBrMass = GetFASTPar(edDataOut, 'YawBrMass');
p.PtfmMass  = GetFASTPar(edDataOut, 'PtfmMass');
p.PtfmRIner = GetFASTPar(edDataOut, 'PtfmRIner');
p.PtfmPIner = GetFASTPar(edDataOut, 'PtfmPIner');
p.PtfmYIner = GetFASTPar(edDataOut, 'PtfmYIner');
p.GBoxEff   = GetFASTPar(edDataOut, 'GBoxEff');
p.GBRatio   = GetFASTPar(edDataOut, 'GBRatio');
p.DTTorSpr  = GetFASTPar(edDataOut, 'DTTorSpr');
p.DTTorDmp  = GetFASTPar(edDataOut, 'DTTorDmp');


p.NTwGages  = GetFASTPar(edDataOut, 'NTwGages');
p.TwrGagNd  = GetFASTPar(edDataOut, 'TwrGagNd');
p.NBlGages  = GetFASTPar(edDataOut, 'NBlGages');
p.BldGagNd  = GetFASTPar(edDataOut, 'BldGagNd');
%p.OutFile   = GetFASTPar(edDataOut, 'OutFile
%p.OutFileFmt= GetFASTPar(edDataOut, 'OutFileFmt %wrbinoutput, wrtxtoutput???
p.OutFmt    = GetFASTPar(edDataOut, 'OutFmt');
p.Tstart    = GetFASTPar(edDataOut, 'Tstart');
%p.DecFact   = GetFASTPar(edDataOut, 'DecFact
% p.NumOuts   = GetFASTPar(edDataOut, 'NumOuts');

if ( p.NumBl == 2 )
    p.UndSling = GetFASTPar(edDataOut, 'UndSling');
    p.TeetMod  = GetFASTPar(edDataOut, 'TeetMod');
    p.TeetDmpP = GetFASTPar(edDataOut, 'TeetDmpP');
    p.TeetDmp  = GetFASTPar(edDataOut, 'TeetDmp');
    p.TeetCDmp = GetFASTPar(edDataOut, 'TeetCDmp');
    p.TeetSStP = GetFASTPar(edDataOut, 'TeetSStP');
    p.TeetHStP = GetFASTPar(edDataOut, 'TeetHStP');
    p.TeetSSSp = GetFASTPar(edDataOut, 'TeetSSSp');
    p.TeetHSSp = GetFASTPar(edDataOut, 'TeetHSSp');
else % Three-bladed turbines don't use these parameters, so set them to zero.
    p.UndSling = 0.0;
    p.TeetMod  = 0;
    p.TeetDmpP = 0.0;
    p.TeetDmp  = 0.0;
    p.TeetCDmp = 0.0;
    p.TeetSStP = 0.0;
    p.TeetHStP = 0.0;
    p.TeetSSSp = 0.0;
    p.TeetHSSp = 0.0;
end


p.TipMass(1)   = GetFASTPar(edDataOut, 'TipMass(1)');
p.TipMass(2)   = GetFASTPar(edDataOut, 'TipMass(2)');
p.TipMass(3)   = GetFASTPar(edDataOut, 'TipMass(3)');


%...............................................................................................................................
% Calculate some indirect inputs:
%...............................................................................................................................
p.TwoPiNB   = (2*pi)/p.NumBl;                                                   % 2*pi/NumBl is used in RtHS().

p.rZT0zt    = p.TowerBsHt - p.PtfmRefzt;                                         % zt-component of position vector rZT0.
p.RefTwrHt  = p.TowerHt   - p.PtfmRefzt;                                         % Vertical distance between ElastoDyn's undisplaced tower height (variable TowerHt) and ElastoDyn's inertia frame reference point (variable PtfmRef).
p.TwrFlexL  = p.TowerHt   - p.TowerBsHt;                                         % Height / length of the flexible portion of the tower.
p.BldFlexL  = p.TipRad    - p.HubRad;                                            % Length of the flexible portion of the blade.
if (p.BD4Blades), p.BldFlexL = 0.0; end

p.rZYzt     = GetFASTPar(edDataOut, 'PtfmCMzt') - p.PtfmRefzt;


for i= 1:p.NumBl
    p.PreCone(i)= GetFASTPar(edDataOut, sprintf('Precone(%d)', i));
end
p.CosPreC  = cosd(p.PreCone);
p.SinPreC  = sind(p.PreCone);

if ( p.NumBl == 2 )
    p.CosDel3  = cosd( GetFASTPar(edDataOut, 'Delta3'));
    p.SinDel3  = sind( GetFASTPar(edDataOut, 'Delta3'));
else
    p.CosDel3  = 1.0;
    p.SinDel3  = 0.0;
end

%...............................................................................................................................

% Calculate the average tip radius normal to the shaft (AvgNrmTpRd)
%   and the swept area of the rotor (ProjArea):

p.AvgNrmTpRd = p.TipRad*sum(p.CosPreC)/p.NumBl;     % Average tip radius normal to the saft.
p.ProjArea   = pi*( p.AvgNrmTpRd^2 );              % Swept area of the rotor projected onto the rotor plane (the plane normal to the low-speed shaft).

p.RotSpeed  = 1;               % unit rotor speed 
p.CShftTilt = cosd( GetFASTPar(edDataOut, 'ShftTilt'));
p.SShftTilt = sind( GetFASTPar(edDataOut, 'ShftTilt'));

p.HubHt     = p.TowerHt + GetFASTPar(edDataOut, 'Twr2Shft') + p.OverHang*p.SShftTilt;


% Direct copy of InputFileData to parameters

%p.FlapforF1  = InputFileData.FlapforF1
%p.FlapforF2  = InputFileData.FlapforF2
%p.EdgeforF   = InputFileData.EdgeforF
%p.TeetforF   = InputFileData.TeetforF
%p.DrTrforF   = InputFileData.DrTrforF
%p.GenforF    = InputFileData.GenforF
%p.YawforF    = InputFileData.YawforF
%p.TwFAforF1  = InputFileData.TwFAforF1
%p.TwFAforF2  = InputFileData.TwFAforF2
%p.TwSSforF1  = InputFileData.TwSSforF1
%p.TwSSforF2  = InputFileData.TwSSforF2
%p.PtfmSgforF = InputFileData.PtfmSgforF
%p.PtfmSwforF = InputFileData.PtfmSwforF
%p.PtfmHvforF = InputFileData.PtfmHvforF
%p.PtfmRforF  = InputFileData.PtfmRforF
%p.PtfmPforF  = InputFileData.PtfmPforF
%p.PtfmYforF  = InputFileData.PtfmYforF
%p.Azimuth   = InputFileData.Azimuth
% p.RotSpeed  = GetFASTPar(edDataOut, 'RotSpeed')/30*pi;
%p.TTDspFA   = InputFileData.TTDspFA
%p.TTDspSS   = InputFileData.TTDspSS
%p.PtfmSurge = InputFileData.PtfmSurge
%p.PtfmSway  = InputFileData.PtfmSway
%p.PtfmHeave = InputFileData.PtfmHeave
%p.PtfmRoll  = InputFileData.PtfmRoll
%p.PtfmPitch = InputFileData.PtfmPitch
%p.PtfmYaw   = InputFileData.PtfmYaw
p.HubCM     = GetFASTPar(edDataOut, 'HubCM');
p.AzimB1Up  = GetFASTPar(edDataOut, 'AzimB1Up');

p.NacCMxn   = GetFASTPar(edDataOut, 'NacCMxn');
p.NacCMyn   = GetFASTPar(edDataOut, 'NacCMyn');
p.NacCMzn   = GetFASTPar(edDataOut, 'NacCMzn');
%p.NcIMUxn   = InputFileData.NcIMUxn
%p.NcIMUyn   = InputFileData.NcIMUyn
%p.NcIMUzn   = InputFileData.NcIMUzn


% plus everything else from FAST_Initialize


end
