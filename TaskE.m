% Course material in EI2410 F�ltteori f�r v�gledare% 2024-02-05 Martin Norgrenclcclear all% Pre-calculated scattering parameters%load('S_Block.mat'); % From a frequency blocking the passage%load('S_Pass.mat'); % From a frequency with almost total passagef=1.6625e+09;c=299792458;% n_max = 70; % Number of modes in the widest of the waveguides w = 2*pi*f;k = w/c;% Sizes of the aperturesa1 = 15.47e-2;  % 1.5cma2 = 75.47*sqrt(1.0001)*1e-2; % 3.14 cma3 = 15e-2;    % 1cm% Location of aperturesxoff1 = 30e-2;      % 1 cmxoff2 = 0;         % 0 cmxoff3 = 7.5e-2;      % 0.5 cm% Length of middle cavitydelta_z = 37.5e-2; % 3cm% Geometrical parameters of the waveguides, used for the precalculated% S-parameters% In each row x_min,x_max,z_min,z_maxXZ = [ xoff1, xoff1+a1, -30e-2, 0; ...       xoff2, xoff2+a2, 0, delta_z; ...       xoff3, xoff3+a3, delta_z, delta_z+37.5e-2];xm = XZ(:,1); xM=XZ(:,2); zm=XZ(:,3); zM=XZ(:,4);a = xM - xm; % Waveguide widthsl1 = 14;l2 = 70;l3 = 14;n1 = (1:l1).';n2 = (1:l2).';n3 = (1:l3).';% Wavenumberskx1 = n1*pi/a(1); Prop1=min(find(k<kx1))-1;kz1 = ( ( k >= kx1 ) - 1i*( k < kx1 ) ) .* sqrt( abs( k^2 - kx1.^2 ) );kx2 = n2*pi/a(2); Prop2=min(find(k<kx2))-1;kz2 = ( ( k >= kx2 ) - 1i*( k < kx2 ) ) .* sqrt( abs( k^2 - kx2.^2 ) );kx3 = n3*pi/a(3); Prop3=min(find(k<kx3))-1;kz3 = ( ( k >= kx3 ) - 1i*( k < kx3 ) ) .* sqrt( abs( k^2 - kx3.^2 ) );% Calculation of S parameters[S11_L, S12_L, S21_L, S22_L] = aperture_S(kx1, kx2, kz1, kz2, l1, l2, k, k, xoff1, a1, a2);[S33_R, S32_R, S23_R, S22_R] = aperture_S(kx3, kx2, kz3, kz2, l3, l2, k, k, xoff3, a3, a2);P = diag(exp(-1i*kz2*delta_z));[S11_LP,S12_LP,S21_LP,S22_LP] = combineLR(S11_L,S12_L,S21_L,S22_L,zeros(l2),P,P,zeros(l2));[S22_PR, S23_PR, S32_PR, S33_PR] = combineLR(zeros(l2),P,P,zeros(l2),S22_R,S23_R,S32_R,S33_R);[S11_13, S12_13, S21_13, S22_13] = combineLR(S11_L,S12_L,S21_L,S22_L,S22_PR,S23_PR,S32_PR,S33_PR);PP_12 = inv(eye(l2) - S22_L*S22_PR) * S21_L;PM_12 = inv(eye(l2) - S22_L*S22_PR) * (S22_L*S23_PR);MP_23 = inv(eye(l2) - S22_R*S22_LP) * (S22_R*S21_LP);MM_23 = inv(eye(l2) - S22_R*S22_LP) * S23_R;% Excitation coefficients of incoming modesa1p = zeros(size(n1)); a1p(1)=1;a3m = zeros(size(n3)); a3m(1)=0; % Scattered and interior modesa1m = S11_13*a1p+S12_13*a3m;a3p = S21_13*a1p+S22_13*a3m;a2p = PP_12*a1p+PM_12*a3m;a2m = MP_23*a1p+MM_23*a3m;Nx_max = 50;Nx1 = round(a(1)/max(a)*Nx_max);Nx2 = round(a(2)/max(a)*Nx_max);Nx3 = round(a(3)/max(a)*Nx_max);Nz1 = round((zM(1)-zm(1))/max(a)*Nx_max);Nz2 = round((zM(2)-zm(2))/max(a)*Nx_max);Nz3 = round((zM(3)-zm(3))/max(a)*Nx_max);xv1 = linspace(xm(1),xM(1),Nx1);z1 = linspace(zm(1),zM(1),Nz1);[X1,Z1] = meshgrid(xv1,z1);EY_1 = zeros(size(X1));for n = 1:l1  moden = 1/sqrt(a(1)*abs(kz1(n)))*sin(kx1(n)*(X1-xm(1)));	  EY_1 = EY_1 + moden .* ( a1p(n)*exp(-1i*kz1(n)*(Z1-zM(1))) + a1m(n)*exp(+1i*kz1(n)*(Z1-zM(1))) );	  end	xv2 = linspace(xm(2),xM(2),Nx2);z2 = linspace(zm(2),zM(2),Nz2);[X2,Z2] = meshgrid(xv2,z2);EY_2 = zeros(size(X2));for n = 1:l2  moden = 1/sqrt(a(2)*abs(kz2(n)))*sin(kx2(n)*(X2-xm(2)));	  EY_2 = EY_2 + moden .* ( a2p(n)*exp(-1i*kz2(n)*(Z2-zm(2))) + a2m(n)*exp(+1i*kz2(n)*(Z2-zM(2))) );	  end	xv3 = linspace(xm(3),xM(3),Nx3);z3 = linspace(zm(3),zM(3),Nz3);[X3,Z3] = meshgrid(xv3,z3);EY_3 = zeros(size(X3));for n = 1:l3  moden = 1/sqrt(a(3)*abs(kz3(n)))*sin(kx3(n)*(X3-xm(3)));	  EY_3 = EY_3 + moden .* ( a3p(n)*exp(-1i*kz3(n)*(Z3-zm(3))) + a3m(n)*exp(+1i*kz3(n)*(Z3-zm(3))) );	  end	E_max = max( [max( max( abs( EY_1 ) ) ), max( max( abs( EY_2 ) ) ), max( max( abs( EY_3 ) ) )] );V = linspace(-E_max,E_max,30);lw = 3;T = 2*pi/w;B = 20;  % Nr of frames - 1figure(6)clfFS=14;film1 = moviein(B+1);for bild = 0:B,  bild  t = bild/B * 2* T;contourf(Z1,X1,imag(EY_1*exp(1i*w*t)),V,'LineStyle','none'); hold on;contour(Z1,X1,real(EY_1*exp(1i*w*t)),V,'-k'); hold on;contourf(Z2,X2,imag(EY_2*exp(1i*w*t)),V,'LineStyle','none'); hold on;contour(Z2,X2,real(EY_2*exp(1i*w*t)),V,'-k'); hold on;contourf(Z3,X3,imag(EY_3*exp(1i*w*t)),V,'LineStyle','none'); hold on;contour(Z3,X3,real(EY_3*exp(1i*w*t)),V,'-k');   axis equal;   zoom on;   tl = title('TE_{m0}-modes in connected waveguides.');    set(tl,'FontSize',FS)  xl=xlabel(strcat('L: ', int2str(min(find(k<kx1))-1),' modes can propagate. M:', ...   int2str(min(find(k<kx2))-1),' modes can propagate. R:', ...   int2str(min(find(k<kx3))-1),' modes can propagate.')); set(xl,'FontSize',FS)   yl=ylabel('Value of E_y and H-lines'); set(yl,'FontSize',FS)  hold off  film1(:,bild+1) = getframe;endmovie( film1, 10, 15 )