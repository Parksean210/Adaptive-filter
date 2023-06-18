function get_est_tf(obj)
%METHOD1 이 메서드의 요약 설명 위치
%   자세한 설명 위치
num_filter = length(obj.type);
num_param = 1;
obj.tf_est_eq = ones(obj.num_fft, 1);
obj.tf_est_ind = zeros(obj.num_fft, num_filter);
obj.sos = zeros(num_filter, 6);
for i = 1:num_filter
    if obj.type(i) == "lsf"
        G = obj.est_parameter(num_param);
        fc = obj.est_parameter(num_param + 1);
        V0 = 10^(G/20);
        H0 = V0 - 1;
        if G >= 0 % boost
            a = (tan(pi*fc/obj.fs) - 1)/(tan(pi*fc/obj.fs) + 1);
        else % cut
            a = (tan(pi*fc/obj.fs) - V0)/(tan(pi*fc/obj.fs) + V0);
        end
        num = [1+(1+a)*H0/2; a+(1+a)*H0/2];
        den = [1; a];
        tf_lsf = obj.z1*num./(obj.z1*den);
        obj.tf_est_ind(:,i) = tf_lsf;
        obj.tf_est_eq = obj.tf_est_eq.*tf_lsf;
        obj.sos(i,:) = [num' 0 den' 0];
        num_param = num_param + 2;
    elseif obj.type(i) == "hsf"
        G = obj.est_parameter(num_param);
        fc = obj.est_parameter(num_param + 1);
        V0 = 10^(G/20);
        H0 = V0 - 1;
        if G >= 0 % boost
            a = (tan(pi*fc/obj.fs) - 1)/(tan(pi*fc/obj.fs) + 1);
        else % cut
            a = (V0*tan(pi*fc/obj.fs) - 1)/(V0*tan(pi*fc/obj.fs) + 1);
        end
        num = [1+(1-a)*H0/2; a+(a-1)*H0/2];
        den = [1; a];
        tf_hsf = obj.z1*num./(obj.z1*den);
        obj.tf_est_ind(:,i) = tf_hsf;
        obj.tf_est_eq = obj.tf_est_eq.*tf_hsf;
        obj.sos(i,:) = [num' 0 den' 0];
        num_param = num_param + 2;
    elseif obj.type(i) == "peak"
        G = obj.est_parameter(num_param);
        fb = obj.est_parameter(num_param + 1);
        fc = obj.est_parameter(num_param + 2);
        V0 = 10^(G/20);
        H0 = V0 - 1;
        d = -cos(2*pi*fc/obj.fs);
        if G >= 0 % boost
            a = (tan(pi*fb/obj.fs) - 1)/(tan(pi*fb/obj.fs) + 1);
        else % cut
            a = (tan(pi*fb/obj.fs) - V0)/(tan(pi*fb/obj.fs) + V0);
        end
        num = [1+(1+a)*H0/2; d*(1-a); -a-(1+a)*H0/2];
        den = [1; d*(1-a); -a];
        tf_pf = obj.z2*num./(obj.z2*den);
        obj.tf_est_ind(:,i) = tf_pf;
        obj.tf_est_eq = obj.tf_est_eq.*tf_pf;
        obj.sos(i,:) = [num' den'];
        num_param = num_param + 3;
    else
        error('Undefined filter type')
    end
end
end