r=linspace(0,2, 256);
phi = linspace(0, 2*pi, 256);
[R, Phi] = meshgrid(r, phi);
w0=1;
m=0;
campo = @(r, phi) power((r/w0),m).*exp(-r.^2/w0^2).*exp(1i*m*phi);
campoMatriz=campo(R,Phi);
X=R.*cos(Phi);
Y=R.*sin(Phi);
figure
surf(X, Y, real(campoMatriz), 'EdgeColor', 'none')
axis equal
xlabel('x')
ylabel('y')
zlabel('Re(E)')
title('Real part of the field')
colorbar
view(2)

figure
mesh(X, Y, campoMatriz.*conj(campoMatriz))


E = polarization("V", campoMatriz);
Ex=E(:, :, 1);
Ey=E(:,:, 2);
figure
mesh(X, Y, real(Ex))
title("Re(E_x)")
figure
mesh(X, Y, real(Ey))
title("Re(E_y)")
theta = linspace(0, 2*pi, 100);
I = zeros(size(theta));   
for k = 1:length(theta)
    th = theta(k);
    M = [cos(th)^2,           cos(th)*sin(th);
        sin(th)*cos(th),     sin(th)^2];
    Ex_out = M(1,1).*Ex + M(1,2).*Ey;
    Ey_out = M(2,1).*Ex + M(2,2).*Ey;
    I(k) = sum(abs(Ex_out).^2 + abs(Ey_out).^2, "all");
end
figure
plot(theta, I)
theta = linspace(0, 2*pi, 100);
I = zeros(size(theta));
for k = 1:length(theta)
    th = theta(k);
    M1 = [cos(pi/4)^2 - 1i*sin(pi/4)^2,        (1 + 1i)*cos(pi/4)*sin(pi/4);
        (1 + 1i)*cos(pi/4)*sin(pi/4),        sin(pi/4)^2 - 1i*cos(pi/4)^2];
    M2= [cos(th)^2,           cos(-th)*sin(-th);
        sin(-th)*cos(th),     sin(-th)^2];
    Ex_out1 = M1(1,1).*Ex + M1(1,2).*Ey;
    Ey_out1 = M1(2,1).*Ex + M1(2,2).*Ey;
    Ex_out2 = M2(1,1).*Ex_out1 + M2(1,2).*Ey_out1;
    Ey_out2 = M2(2,1).*Ex_out1 + M2(2,2).*Ey_out1;
    I(k) = sum(abs(Ex_out2).^2 + abs(Ey_out2).^2, "all");
end
figure
plot(theta, I)

analyzers = [1,  0;     
             0,  1;     
             1,  1;     
             1, -1;     
             1,  1i;     
             1, -1i];    

I = zeros(size(analyzers,1),1);

I0 = sum(abs(Ex_out2).^2 + abs(Ey_out2).^2, 'all');

for k = 1:size(analyzers,1)

    a = analyzers(k,1);
    b = analyzers(k,2);

    normFactor = sqrt(abs(a)^2 + abs(b)^2);
    a = a / normFactor;
    b = b / normFactor;

    Eproj = conj(a)*Ex_out2 + conj(b)*Ey_out2;

    I(k) = sum(abs(Eproj).^2, 'all')/I0;

end

IH = I(1);
IV = I(2);
ID = I(3);
IA = I(4);
IR = I(5);
IL = I(6);

S0 = IH + IV;
S1 = IH - IV;
S2 = ID - IA;
S3 = IR - IL;

P = sqrt(S1^2 + S2^2 + S3^2) / S0;

function polarizado = polarization(elected, campoMatriz)
polarizado = zeros(size(campoMatriz,1), size(campoMatriz,2), 2);
if elected == "H"
    jones = [1, 0];
elseif elected == "V"
    jones = [0, 1];
elseif elected == "R"
    jones = (1/sqrt(2)) * [1, 1i];
elseif elected == "L"
    jones = (1/sqrt(2)) * [1, -1i];
elseif elected == "D"
    jones = (1/sqrt(2)) * [1, 1];
elseif elected == "LD"
    jones = (1/sqrt(2)) * [1, -1];
else
    error("invalido")
end
polarizado(:,:,1) = campoMatriz * jones(1); 
polarizado(:,:,2) = campoMatriz * jones(2); 
end



