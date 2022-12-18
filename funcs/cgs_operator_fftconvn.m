function [y, Fk0, siz] = cgs_operator_fftconvn(x, mode, Fk0, k0, siz)
% if isempty(Fk0)
%     siz_x = size(x);
%     siz_k = size(k0);
%     siz_y = siz_x + siz_k-1;
%     siz_valid = siz_y - siz_k+1;
%     siz_filter =  siz_y-siz_valid+1;
%     
%     Fk0 = fftn(k0, siz_y);
%     siz.siz_x = siz_x;
%     siz.siz_k = siz_k;
%     siz.siz_y = siz_y;
%     siz.siz_valid = siz_valid;
%     siz.siz_filter = siz_filter;
% end
% 
% siz_x = siz.siz_x;
% siz_k = siz.siz_k;
% siz_y = siz.siz_y;
% siz_valid = siz.siz_valid;
% siz_filter = siz.siz_filter;
% 
% 
% switch mode
%     case 'A'
%         Fx = fftn(x, siz_y);
%         Fy = Fx.*Fk0;
%         y = real(ifftn(Fy));
%         y = y(floor(siz_filter(1)/2)+(1:siz_valid(1)), floor(siz_filter(2)/2)+(1:siz_valid(2)), floor(siz_filter(3)/2)+(1:siz_valid(3)));
%         
%     case 'Adj'
%         x0 = zeros(siz_y);
%         x0(floor(siz_filter(1)/2)+(1:siz_valid(1)), floor(siz_filter(2)/2)+(1:siz_valid(2)), floor(siz_filter(3)/2)+(1:siz_valid(3))) = x;
%         Fx0 = fftn(x0, siz_y);
%         Fy = Fx0.*conj(Fk0);
%         y = real(ifftn(Fy));
%         y = y(1:siz_x(1), 1:siz_x(2), 1:siz_x(3));

if isempty(Fk0)
    siz_x = size(x);
    siz_k = size(k0);
    siz_y = siz_x + siz_k - 1;
    siz_valid = siz_k - siz_x + 1; % edit
    siz_filter =  siz_y - siz_valid+1;
    
    Fk0 = fftn(k0, siz_y);
    siz.siz_x = siz_x;
    siz.siz_k = siz_k;
    siz.siz_y = siz_y;
    siz.siz_valid = siz_valid;
    siz.siz_filter = siz_filter;
end

siz_x = siz.siz_x;
siz_k = siz.siz_k;
siz_y = siz.siz_y;
siz_valid = siz.siz_valid;
siz_filter = siz.siz_filter;


switch mode
    case 'A'
        Fx = fftn(x, siz_y);
        Fy = Fx.*Fk0;
        y = real(ifftn(Fy));
        y = y(floor(siz_filter(1)/2)+(1:siz_valid(1)), floor(siz_filter(2)/2)+(1:siz_valid(2)), floor(siz_filter(3)/2)+(1:siz_valid(3)));
        
    case 'Adj'
        x0 = zeros(siz_y, 'like', x);
        x0(floor(siz_filter(1)/2)+(1:siz_valid(1)), floor(siz_filter(2)/2)+(1:siz_valid(2)), floor(siz_filter(3)/2)+(1:siz_valid(3))) = x;
        Fx0 = fftn(x0, siz_y);
        Fy = Fx0.*conj(Fk0);
        y = real(ifftn(Fy));
        y = y(1:siz_x(1), 1:siz_x(2), 1:siz_x(3));        
end
end


