function Mout = quaternion_cluster_interpolation(Min)
% quaternion_cluster_spline function
% Interpolation on the basis of orientiation.
% FUNCTION
%       Mout = quaternion_cluster_spline(Min)
% INPUT
%       Min  = Global positional cluster marker kinematic data (n*9)
%              (collumns should correspond to:
%                     Min = [x1 y1 z1 x2 y2 z2 x3 y3 z3];)
% OUTPUT
%       Mout = Interpolated data
%              (collumns correspond to:
%                    Mout = [x1 y1 z1 x2 y2 z2 x3 y3 z3];)
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Created by N.Kluft (2019)   [nick.kluft@gmail.com]                      %
% Functions adapted from G. Faber & S.M. Bruijn                           %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Make empty matrices
% uses indices to easily link xxxyyyzzz to xyzxyzxyz data
indices = transpose(reshape(1:9,3,3));

% arbitrary marker assignment
mark.m1 = Min(:,1:3);
mark.m2 = Min(:,4:6);
mark.m3 = Min(:,7:9);

%% create orientation matrix
xas= mark.m2-mark.m1;
yas_temp= mark.m3-mark.m1;
% make vector between every marker.
mark.vec12 = xas;
mark.vec21 =-xas;
mark.vec13 = yas_temp;
mark.vec31 =-yas_temp;
mark.vec23 = mark.m3-mark.m2;
mark.vec32 = -mark.vec23;

% create orientation matrix
zas = cross(xas,yas_temp);
xas = xas./norm_col(xas);
yas = yas_temp./norm_col(yas_temp);
zas = zas./norm_col(zas);
R = [xas,yas,zas];

%% transfer R to quaternion
quat = rot2quat(R);
% create old orientation matrix
Rold = quat2rot(quat);
% find the corresponding gaps
igap = any(isnan(Min(:,1:3:end)),2);
igap = find(diff(igap>0));

%% make extrapolation impossible:
% when data at start of traj is missing:
if isnan(quat(1)) && ~isempty(igap)
    igap(1)=[];
    disp('cannot interpolate first part data segment ')
end
% when data at the end of the traj is missing
if isnan(quat(end)) && ~isempty(igap)
    igap(end)=[];
    disp('cannot interpolate first part data segment ')
end

%% interpolation section
if ~isempty(igap)
    % now normalize the two quaternions at either side of the gap
    for ig = 1:2:numel(igap)
        % update mark.m's here when ig>1 (ugly solution)
        if ig>1
            mark.m1 = Min(:,1:3);
            mark.m2 = Min(:,4:6);
            mark.m3 = Min(:,7:9);
        end
        
        % size of gap in samples
        sgap = (igap(ig+1))-igap(ig);
        % repeat and normalize quaternion for size sgap (starting sample)
        pn = repmat(quatnormalize(quat(igap(ig),:)),sgap+2,1);
        % repeat and normalize quaternion for size sgap (ending sample)
        qn = repmat(quatnormalize(quat(igap(ig+1)+1,:)),sgap+2,1);
        % Sanity check here
        if any(isnan([pn;qn]))
            error('Nans used for interpolation (quaterion not assigned)')
        end
        
        % interpolate now over the quaternion
        qi = quatinterp(pn,qn,[0,1/(sgap+1):1/(sgap+1):1-(1/(sgap+1)),1],'slerp');
        % make again rotation matrices from interpolated quaternion
        quat(igap(ig):igap(ig+1)+1,:) = qi;
        
        % return the rotation matrix
        Rnew = quat2rot(quat);
        
        if numel(igap(ig)+1:igap(ig+1))>1
            % find marker(s) missing
            m = find(any(isnan(Min(igap(ig)+1:igap(ig+1),1:3:end)))); % missing marker
            % find marker that is visible during traject
            o = find(~any(isnan(Min(igap(ig)+1:igap(ig+1),1:3:end)))); % observed marker
        else
            % find marker(s) missing
            m = find(isnan(Min(igap(ig)+1:igap(ig+1),1:3:end))); % missing marker
            % find marker that is visible during traject
            o = find(~isnan(Min(igap(ig)+1:igap(ig+1),1:3:end))); % observed marker
        end
        
        if isempty(o)
             disp(['no quaternion spline possible, as none of the clust'...
                 'er markers were visible. Consider a positional spline'...
                 ' first']);
        else
        for im = m
            % reconstruct marker positions
            mark.(['m',num2str(im)]) = mark.(['m',num2str(o(1))])+ ...
                prod_col(Rnew,prod_col(transpose_col(Rold(igap(ig),:)),...
                mark.(['vec',num2str(o(1)),num2str(im)])(igap(ig),:)));
            if numel(o)>1 % take average position if both markers visible
            mark.(['m',num2str(im)]) = (mark.(['m',num2str(im)])+mark.(['m',num2str(o(2))])+ ...
                prod_col(Rnew,prod_col(transpose_col(Rold(igap(ig),:)),...
                mark.(['vec',num2str(o(2)),num2str(im)])(igap(ig),:))))/2;
            end
            Min(igap(ig)+1:igap(ig+1),indices(im,:))...
                = mark.(['m',num2str(im)])(igap(ig)+1:igap(ig+1),:);
        end
        end
    end
end

%% Store newly derived position of markers
Mout = Min;
end

function R = quat2rot( Qrotation )
% qGetR: get a nx9 rotation matrix
% R = quat2rot( Qrotation )
% IN:
%     Qrotation - quaternion describing rotation
%
% OUT:
%     R - rotation matrix
%
% VERSION: 03.03.2012
w = Qrotation( :,1 );
x = Qrotation( :,2 );
y = Qrotation( :,3 );
z = Qrotation( :,4 );
Rxx = 1 - 2*(y.^2 + z.^2);
Rxy = 2*(x.*y - z.*w);
Rxz = 2*(x.*z + y.*w);
Ryx = 2*(x.*y + z.*w);
Ryy = 1 - 2*(x.^2 + z.^2);
Ryz = 2*(y.*z - x.*w );
Rzx = 2*(x.*z - y.*w );
Rzy = 2*(y.*z + x.*w );
Rzz = 1 - 2 *(x.^2 + y.^2);
R = [Rxx, Rxy, Rxz, Ryx, Ryy, Ryz, Rzx, Rzy, Rzz];
end

function Q = rot2quat( R )
% qGetQ: converts 3x3 rotation matrix into equivalent quaternion
% Q = qGetQ( R );
[~,c] = size( R );
if( c ~= 9 )
    fprintf( 'R must be a nx9 matrix\n\r' );
    return;
end
% [ Rxx, Rxy, Rxz ] = R(1,1:3);
% [ Ryx, Ryy, Ryz ] = R(2,1:3);
% [ Rzx, Rzy, Rzz ] = R(3,1:3);
Rxx = R(:,1); Rxy = R(:,2); Rxz = R(:,3);
Ryx = R(:,4); Ryy = R(:,5); Ryz = R(:,6);
Rzx = R(:,7); Rzy = R(:,8); Rzz = R(:,9);
w = sqrt( sum(R(:,[1,5,9]),2) + 1 ) / 2;
% check if w is real. Otherwise, zero it.
if( imag( w ) > 0 )
    w = 0;
end
x = sqrt( 1 + Rxx - Ryy - Rzz ) / 2;
y = sqrt( 1 + Ryy - Rxx - Rzz ) / 2;
z = sqrt( 1 + Rzz - Ryy - Rxx ) / 2;

[~, i ] = max( [w,x,y,z],[],2 );
for n = 1:numel(x)
    switch i(n)
        case 1
            x(n) = ( Rzy(n) - Ryz(n) ) / (4*w(n));
            y(n) = ( Rxz(n) - Rzx(n) ) / (4*w(n));
            z(n) = ( Ryx(n) - Rxy(n) ) / (4*w(n));
            
        case 2
            w(n) = ( Rzy(n) - Ryz(n) ) / (4*x(n));
            y(n) = ( Rxy(n) + Ryx(n) ) / (4*x(n));
            z(n) = ( Rzx(n) + Rxz(n) ) / (4*x(n));
            
        case 3
            w(n) = ( Rxz(n) - Rzx(n) ) / (4*y(n));
            x(n) = ( Rxy(n) + Ryx(n) ) / (4*y(n));
            z(n) = ( Ryz(n) + Rzy(n) ) / (4*y(n));
            
        case 4
            w(n) = ( Ryx(n) - Rxy(n) ) / (4*z(n));
            x(n) = ( Rzx(n) + Rxz(n) ) / (4*z(n));
            y(n) = ( Ryz(n) + Rzy(n) ) / (4*z(n));
    end
end
Q = [ w, x, y, z ];
end

function  [normcol3] = norm_col(data)
% [normcol3] = norm_col(data)
% -------------------------------------------------------------------------
% This function calculates the norm (length) of a vector (V) and copies the result 3
% times. This can be handy when one wants calculate the unit vector (Vunit) at once
% over time vector: Vunit = V./norm_col(V)
% -------------------------------------------------------------------------
% %% Input
% data      : nx3 matrix with xyz data of 1 marker/vector [x y z]
%
% %% Output
% normcol3  : nx3 matrix 3 copies of the norm of the input vectore: [norm norm norm]
% created by Sjoerd M Bruijn & Gert Faber
normcol=sqrt(sum(data.^2,2));
normcol3=[normcol normcol normcol];
end

function  [Rproduct] = prod_col(Rfirst, Rsecond)
% function: [Rproduct] = prod_col(Rfirst, Rsecond)
% -------------------------------------------------------------------------
% This function calculates product of two nx9 matices in the same way a
% multiplication of two 3x3 matrices is done in by matlab (R*R), but than
% over all time samples at once. This is faster and prevents loops over
% time.
% It also calculates the product of a orientation matrix (nx9), and a
% vector (either 1x3 or nx3)
% -------------------------------------------------------------------------
% % INPUT
% Rfirst  : first orientation matrix (n x 9 or 1 x 9)
% Rsecond : second orientation matrix (n x 9 or 1 x 9) or vector (n x 3 or 1 x 3)
%  
% % OUTPUT
% Rproduct: product of Rfirst and Rsecond (n x 9)
%
% (n= number of samples) 
% created by Sjoerd M Bruijn & Gert Faber

[m,n]  =size(Rfirst);
[m1,n1]=size(Rsecond);
if n==3
    error('first argument can not be a vector, only second argument can be a vector!!')
end
if m1==1
    Rsecond=ones(m,1)*Rsecond;
end
if m==1
    Rfirst=ones(m1,1)*Rfirst;
end

no_col=n1/3;

matrix_1=[ 1 4 7
    2 5 8
    3 6 9];
matrix_2=reshape(1:n1,3,no_col);

Rproduct=zeros(max([m m1]),n1);
for i_row=1:3
    row=matrix_1(i_row,:);
    for i_kol=1:no_col
        column=matrix_2(:,i_kol);
        Rproduct(:,matrix_2(i_row,i_kol))=sum    (Rfirst(:,row).*Rsecond(:,column),2);
    end
end
end

function  [transp_R19] = transpose_col(R19)
% [transp_R19] = transpose_col(R19)
% -------------------------------------------------------------------------
% This function calculates the transpose for a nx9 matrix
% equivalent to R' of transpose(R) in case R is a 3x3 matrix. The advantage, 
% however is that no reshaping of the data is necesarry and that the 
% calculation is done over all time samples at once which prevents 
% loops over time samples.
% -------------------------------------------------------------------------
% % INPUT
% R19     : input matrix (n x 9)
%  
% % OUTPUT
% transp_R19: transposed matrix (n x 9)
%
% (n= number of samples) 
% created by Sjoerd M Bruijn & Gert Faber
transp_R19=[R19(:,1:3:end) R19(:,2:3:end) R19(:,3:3:end)];
end
