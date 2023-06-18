classdef AdaptiveFilterEqualization < handle
    %UNTITLED 이 클래스의 요약 설명 위치
    %   자세한 설명 위치
    
    properties
        fs
        f
        type
        num_fft
        true_parameter
        init_parameter
        est_parameter
        learning_rate
        tf_eq
        tf_ind
        tf_est_eq
        tf_est_ind
        sos
        input
        target
        output
        error
    end
    
    properties (Access=private)
        z1
        z2
    end
    
    methods
        function obj = AdaptiveFilterEqualization(samplerate, type, num_fft, true_parameter, init_parameter, learning_rate, input)
            %UNTITLED 이 클래스의 인스턴스 생성
            %   자세한 설명 위치
            obj.fs = samplerate;
            obj.type = type;
            obj.num_fft = num_fft;
            obj.f = (0:obj.num_fft-1)'/obj.num_fft*obj.fs;
            obj.true_parameter = true_parameter;
            obj.init_parameter = init_parameter;
            obj.learning_rate = learning_rate;
            obj.input = input;
            obj.z1 = exp(-2*pi*1i/obj.num_fft*(0:obj.num_fft - 1)'*(0:1));
            obj.z2 = exp(-2*pi*1i/obj.num_fft*(0:obj.num_fft - 1)'*(0:2));
        end
        
        % first-order low shelf filter
        [tf_lsf, num, den] = get_tf_lsf(obj, G, fc)
        % first-order high shelf filter
        [tf_hsf, num, den] = get_tf_hsf(obj, G, fc)
        % second-order peak filter
        [tf_pf, num, den] = get_tf_pf(obj, G, fb, fc)
        % total tf
        get_tf(obj)
        % estimated tf
        get_est_tf(obj)
        % filter inputs
        filter_eq(obj)
        % adaptation
        adapt_parameter(obj)
        
        
        
    end
end

