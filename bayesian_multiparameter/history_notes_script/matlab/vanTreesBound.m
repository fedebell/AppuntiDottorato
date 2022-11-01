Num1 = 30000;
%Q-plates charges
s = [1 2 11 51];

visibilities = [0.900, 0.875, 0.850, 0.765];

%Computation of the lower bound
Itheta = zeros(1, 4);
Ivis = zeros(1, 4);

for i=1:1:4
    Itheta(i) = 2*s(i)^2*(1-sqrt(1-(visibilities(i))^2));
    Ivis(i) = 2*(1-sqrt(1-(visibilities(i))^2))/((visibilities(i))^2*sqrt(1-(visibilities(i))^2));
end

while i<

CRbound = zeros(, 1);

for i=1:1:counterNonZero
    cvx_begin
        variable nu(4)
        minimize G(1, 1)*inv_pos(Itheta(1)*nu(1)+Itheta(2)*nu(2)+Itheta(3)*nu(3)+Itheta(4)*nu(4))+...
        G(2, 2)*inv_pos(Ivis(1)*nu(1))+G(3, 3)*inv_pos(Ivis(2)*nu(2))+G(4, 4)*inv_pos(Ivis(3)*nu(3))+G(5, 5)*inv_pos(Ivis(4)*nu(4));
        subject to
        nu >= 0;
        2*(s(1)*nu(1)+s(2)*nu(2)+s(3)*nu(3)+s(4)*nu(4)) == resources(i, 1);
    cvx_end
    
    CRbound(i) =  G(1, 1)*inv_pos(Itheta(1)*nu(1)+Itheta(2)*nu(2)+Itheta(3)*nu(3)+Itheta(4)*nu(4))+...
        G(2, 2)*inv_pos(Ivis(1)*nu(1))+G(3, 3)*inv_pos(Ivis(2)*nu(2))+G(4, 4)*inv_pos(Ivis(3)*nu(3))+G(5, 5)*inv_pos(Ivis(4)*nu(4));
    
end