function adapt_parameter(obj)
%UNTITLED 이 함수의 요약 설명 위치
%   자세한 설명 위치
num_filter = length(obj.type);
num_sample = length(obj.input);
parameter = obj.init_parameter;
xh = zeros(num_filter, 3);
y1 = zeros(num_filter, 1);
obj.output = zeros(num_sample, 1);
obj.error = zeros(num_sample, 1);
num_sf = 0;
num_peak = 0;
for j = 1:num_filter
    if obj.type(j) == "lsf"
        num_sf = num_sf + 1;
    elseif obj.type(j) == "hsf"
        num_sf = num_sf + 1;
    elseif obj.type(j) == "peak"
        num_peak = num_peak + 1;
    else
        error('Undefined filter type')
    end
end
dy_dx = zeros(num_filter, 1);
dxh_dG = zeros(num_filter,3);
dy_dG = zeros(num_filter,1);
dxh_dfb = zeros(num_filter,3);
dy_dfb = zeros(num_filter,3);
dxh_dfc = zeros(num_filter,3);
dy_dfc = zeros(num_filter,1);

for i = 1:num_sample
    x = obj.input(i);
    num_param = 1;
    for j = 1:num_filter
        if obj.type(j) == "lsf"
            % parameter setting
            G = parameter(num_param);
            fc = parameter(num_param+1);
            V0 = 10^(G/20);
            H0 = V0 - 1;
            if G > 0 % boost
                a = (tan(pi*fc/obj.fs) - 1)/(tan(pi*fc/obj.fs) + 1);
            else % cut
                a = (tan(pi*fc/obj.fs) - V0)/(tan(pi*fc/obj.fs) + V0);
            end
            % filtering
            xh(j,1) = x - a*xh(j,2);
            y1(j,1) = a*xh(j,1) + xh(j,2);
            x = H0/2*(x + y1(j,1)) + x;
            % gradient calculation
            if G > 0
                % dy_dx
                dy_dx(j) = H0/2*(1 + a) + 1;
                % dy_dG
                dy_dG(j) = (x + y1(j,1))/40*10*(G/20)*log(10);
                % dy_dfc
                da_dfc = 2*pi*sec(2*pi*fc/obj.fs)*(sec(2*pi*fc/obj.fs)-tan(2*pi*fc/obj.fs))/obj.fs;
                dxh_dfc(j,1) = -da_dfc*xh(j,2) - a*dxh_dfc(j,2);
                dy1_dfc = da_dfc*xh(j,1) + a*dxh_dfc(j,1) + dxh_dfc(j,2);
                dy_dfc(j) = H0/2*dy1_dfc;
                % circular buffer
                xh(j,2) = xh(j,1);
                dxh_dG(j,2) = dxh_dG(j,1);
                dxh_dfc(j,2) = dxh_dfc(j,1);
            else
                % dy_dx
                dy_dx(j) = H0/2*(1 + a) + 1;
                % dy_dG
                da_dG = (-log(10)*V0*tan(pi*fc/obj.fs))/(10*(tan(pi*fc/obj.fs) + V0)^2);
                dxh_dG(j,1) = -da_dG*xh(j,2) - a*dxh_dG(j,2);
                dy1_dG = da_dG*xh(j,1) + a*dxh_dG(j,1) + dxh_dG(j,2);
                dy_dG(j) = (x + y1(j,1))/40*10^(G/20)*log(10) + H0/2*dy1_dG;
                % dy_dfc
                da_dfc = 2*pi*V0*sec(2*pi*fc/obj.fs)^2/(obj.fs*(tan(pi*fc/obj.fs) + V0)^2);
                dxh_dfc(j,1) = -da_dfc*xh(j,2) - a*dxh_dfc(j,2);
                dy1_dfc = da_dfc*xh(j,1) + a*dxh_dfc(j,1) + dxh_dfc(j,2);
                dy_dfc(j) = H0/2*dy1_dfc;
                % circular buffer
                xh(j,2) = xh(j,1);
                dxh_dG(j,2) = dxh_dG(j,1);
                dxh_dfc(j,2) = dxh_dfc(j,1);
            end
            % calculating parameter index
            num_param = num_param + 2;
        elseif obj.type(j) == "hsf"
            % parameter setting
            G = parameter(num_param);
            fc = parameter(num_param+1);
            V0 = 10^(G/20);
            H0 = V0 - 1;
            if G >= 0 % boost
                a = (tan(pi*fc/obj.fs) - 1)/(tan(pi*fc/obj.fs) + 1);
            else % cut
                a = (V0*tan(pi*fc/obj.fs) - 1)/(V0*tan(pi*fc/obj.fs) + 1);
            end
            % filtering
            xh(j,1) = x - a*xh(j,2);
            y1(j,1) = a*xh(j,1) + xh(j,2);
            x = H0/2*(x - y1(j,1)) + x;
            % gradient calculation
            if G > 0
                % dy_dx
                dy_dx(j) = H0/2*(1 - a) + 1;
                % dy_dG
                dy_dG(j) = (x - y1(j,1))/40*10*(G/20)*log(10);
                % dy_dfc
                da_dfc = 2*pi*sec(2*pi*fc/obj.fs)*(sec(2*pi*fc/obj.fs)-tan(2*pi*fc/obj.fs))/obj.fs;
                dxh_dfc(j,1) = -da_dfc*xh(j,2) - a*dxh_dfc(j,2);
                dy1_dfc = da_dfc*xh(j,1) + a*dxh_dfc(j,1) + dxh_dfc(j,2);
                dy_dfc(j) = -H0/2*dy1_dfc;
                % circular buffer
                xh(j,2) = xh(j,1);
                dxh_dG(j,2) = dxh_dG(j,1);
                dxh_dfc(j,2) = dxh_dfc(j,1);
            else
                % dy_dx
                dy_dx(j) = H0/2*(1 - a) + 1;
                % dy_dG
                da_dG = (log(10)*V0*tan(pi*fc/obj.fs))/(10*(V0*tan(pi*fc/obj.fs) + 1)^2);
                dxh_dG(j,1) = -da_dG*xh(j,2) - a*dxh_dG(j,2);
                dy1_dG = da_dG*xh(j,1) + a*dxh_dG(j,1) + dxh_dG(j,2);
                dy_dG(j) = (x - y1(j,1))/40*10^(G/20)*log(10) - H0/2*dy1_dG;
                % dy_dfc
                da_dfc = 2*pi*V0*sec(2*pi*fc/obj.fs)^2/(obj.fs*(V0*tan(pi*fc/obj.fs) + 1)^2);
                dxh_dfc(j,1) = -da_dfc*xh(j,2) - a*dxh_dfc(j,2);
                dy1_dfc = da_dfc*xh(j,1) + a*dxh_dfc(j,1) + dxh_dfc(j,2);
                dy_dfc(j) = -H0/2*dy1_dfc;
                % circular buffer
                xh(j,2) = xh(j,1);
                dxh_dG(j,2) = dxh_dG(j,1);
                dxh_dfc(j,2) = dxh_dfc(j,1);
            end
            % calculating parameter index
            num_param = num_param + 2;
        elseif obj.type(j) == "peak"
            % parameter setting
            G = parameter(num_param);
            fb = parameter(num_param + 1);
            fc = parameter(num_param + 2);
            V0 = 10^(G/20);
            H0 = V0 - 1;
            d = -cos(2*pi*fc/obj.fs);
            if G > 0 % boost
                a = (tan(pi*fb/obj.fs) - 1)/(tan(pi*fb/obj.fs) + 1);
            else % cut
                a = (tan(pi*fb/obj.fs) - V0)/(tan(pi*fb/obj.fs) + V0);
            end
            % filtering
            xh(j,1) = x - d*(1-a)*xh(j,2) + a*xh(j,3);
            y1(j,1) = -a*xh(j,1) + d*(1 - a)*xh(j,2) + xh(j,3);
            x = H0/2*(x - y1(j,1)) + x;
            % gradient
            if G > 0
                % dy_dx
                dy_dx = H0/2*(1 + a) + 1;
                % dy_dG
                dy_dG(j) = (x - y1(j,1))/40*10^(G/20)*log(10);
                % dy_dfb
                da_dfb = 2*pi*sec(2*pi*fb/obj.fs)*(sec(2*pi*fb/obj.fs)-tan(2*pi*fb/obj.fs))/obj.fs;
                dxh_dfb(j,1) = d*da_dfb*xh(j,2) - d*(1 - a)*dxh_dfb(j,2) + da_dfb*xh(j,3) + a*dxh_dfb(j,3);
                dy1_dfb = -da_dfb*xh(j,1) - a*dxh_dfb(j,1) - d*da_dfb*xh(j,2) + d*(1 - a)*dxh_dfb(j,2) + dxh_dfb(j,3);
                dy_dfb(j) = -H0/2*dy1_dfb;
                % dy_dfc
                dd_dfc = 2*pi*sin(2*pi*fc/obj.fs)/obj.fs;
                dxh_dfc(j,1) = -dd_dfc*(1 - a)*xh(j,2) - d*(1 - a)*dxh_dfc(j,2) + a*dxh_dfc(j,3);
                dy1_dfc = -a*dxh_dfc(j,1) + dd_dfc*(1 - a)*xh(j,2) + d*(1 - a)*dxh_dfc(j,2) + dxh_dfc(j,3);
                dy_dfc(j) = -H0/2*dy1_dfc;
                % circular buffer
                xh(j,3) = xh(j,2);
                xh(j,2) = xh(j,1);
                dxh_dG(j,3) = dxh_dG(j,2);
                dxh_dG(j,2) = dxh_dG(j,1);
                dxh_dfb(j,3) = dxh_dfb(j,2);
                dxh_dfb(j,2) = dxh_dfb(j,1);
                dxh_dfc(j,3) = dxh_dfc(j,2);
                dxh_dfc(j,2) = dxh_dfc(j,1);
            else % cut
                % dy_dx
                dy_dx = H0/2*(1 + a) + 1;
                % dy_dG
                da_dG = -log(10)*V0*tan(pi*fb/obj.fs)/(10*(tan(pi*fb/obj.fs) + V0)^2);
                dxh_dG(j,1) = d*da_dG*xh(j,2) + da_dG*xh(j,3) - d*(1 - a)*dxh_dG(j,2) + a*dxh_dG(j,1);
                dy1_dG = -da_dG*xh(j,1) - a*dxh_dG(j,1) - d*da_dG*xh(j,2) + d*(1 - a)*dxh_dG(j,2) + dxh_dG(j,3);
                dy_dG(j) = (x - y1(j,1))/40*10^(G/20)*log(10) - H0/2*dy1_dG;
                % dy_dfb
                da_dfb = 2*pi*sec(2*pi*fb/obj.fs)*(sec(2*pi*fb/obj.fs)-tan(2*pi*fb/obj.fs))/obj.fs;
                dxh_dfb(j,1) = d*da_dfb*xh(j,2) - d*(1 - a)*dxh_dfb(j,2) + da_dfb*xh(j,3) + a*dxh_dfb(j,3);
                dy1_dfb = -da_dfb*xh(j,1) - a*dxh_dfb(j,1) - d*da_dfb*xh(j,2) + d*(1 - a)*dxh_dfb(j,2) + dxh_dfb(j,3);
                dy_dfb(j) = -H0/2*dy1_dfb;
                % dy_dfc
                dd_dfc = 2*pi*sin(2*pi*fc/obj.fs)/obj.fs;
                dxh_dfc(j,1) = -dd_dfc*(1 - a)*xh(j,2) - d*(1 - a)*dxh_dfc(j,2) + a*dxh_dfc(j,3);
                dy1_dfc = -a*dxh_dfc(j,1) + dd_dfc*(1 - a)*xh(j,2) + d*(1 - a)*dxh_dfc(j,2) + dxh_dfc(j,3);
                dy_dfc(j) = -H0/2*dy1_dfc;
                % circular buffer
                xh(j,3) = xh(j,2);
                xh(j,2) = xh(j,1);
                dxh_dG(j,3) = dxh_dG(j,2);
                dxh_dG(j,2) = dxh_dG(j,1);
                dxh_dfb(j,3) = dxh_dfb(j,2);
                dxh_dfb(j,2) = dxh_dfb(j,1);
                dxh_dfc(j,3) = dxh_dfc(j,2);
                dxh_dfc(j,2) = dxh_dfc(j,1);
            end
            
            % calculating parameter index
            num_param = num_param + 3;
        else
            error('Undefined filter type')
        end
    end
    obj.output(i) = x;
    obj.error(i) = obj.target(i) - obj.output(i);
    dC_dy = -2*obj.error(i);
    num_param = 1;
    for j = 1:num_filter
        if obj.type(j) == "lsf"
            % parameter update
            if mod(i, obj.fs/10) == 0
                parameter(num_param)
            end
            parameter(num_param) = parameter(num_param) - obj.learning_rate(num_param)*sign(dC_dy*dy_dG(j));
            parameter(num_param + 1) = parameter(num_param + 1) - obj.learning_rate(num_param + 1)*sign(dC_dy*dy_dfc(j));
            % calculating parameter index
            num_param = num_param + 2;
        elseif obj.type(j) == "hsf"
            % parameter update
            parameter(num_param) = parameter(num_param) - obj.learning_rate(num_param)*sign(dC_dy*dy_dG(j));
            parameter(num_param + 1) = parameter(num_param + 1) - obj.learning_rate(num_param + 1)*sign(dC_dy*dy_dfc(j));
            % calculating parameter index
            num_param = num_param + 2;
        elseif obj.type(j) == "peak"
            % parameter update
            parameter(num_param) = parameter(num_param) - obj.learning_rate(num_param)*sign(dC_dy*dy_dG(j));
            parameter(num_param + 1) = parameter(num_param + 1) - obj.learning_rate(num_param + 1)*sign(dC_dy*dy_dfb(j));
            parameter(num_param + 2) = parameter(num_param + 2) - obj.learning_rate(num_param + 2)*sign(dC_dy*dy_dfc(j));
            % calculating parameter index
            num_param = num_param + 3;
        else
            error('Undefined filter type')
        end
    end

end
obj.est_parameter = parameter;
end

