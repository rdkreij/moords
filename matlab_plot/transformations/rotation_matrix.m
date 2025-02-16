function M=rotation_matrix(h,p,r,order)

if nargin ==3
    order = [1 2 3]; % as per the ADCP coordinate transformation manual
end

if numel(order) ~= 3 || ~isequal(unique(order), [1 2 3])
    error
end

ch=cosd(h);
sh=sind(h);

cp=cosd(p);
sp=sind(p);

cr=cosd(r);
sr=sind(r);

%%     
H=[ch sh 0;
    -sh ch 0;
    0 0 1]; 


R=[cr 0 sr;
    0 1 0;
    -sr 0 cr];

P=[1 0 0;
    0 cp -sp;
    0 sp cp];
%%
M3(:, :, 1) = H;
M3(:, :, 2) = P;
M3(:, :, 3) = R;
%
M=M3(:, :, order(1))*M3(:, :, order(2))*M3(:, :, order(3)); 