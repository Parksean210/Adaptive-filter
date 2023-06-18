classdef AdaptiveFilterEqualization
    %UNTITLED 이 클래스의 요약 설명 위치
    %   자세한 설명 위치
    
    properties
        fs
        f
        type
        num_fft
    end
    
    properties (Access=private)
        z1
        z2
    end
    
    methods
        function obj = AdaptiveFilterEqualization(samplerate, type, num_fft)
            %UNTITLED 이 클래스의 인스턴스 생성
            %   자세한 설명 위치
            obj.fs = samplerate;
            obj.type = type;
            obj.num_fft = num_fft;
            obj.f = (0:obj.num_fft-1)'/obj.num_fft*obj.fs;
            obj.z1 = exp(-2*pi*1i/obj.num_fft*(0:obj.num_fft - 1)'*(0:1));
            obj.z2 = exp(-2*pi*1i/obj.num_fft*(0:obj.num_fft - 1)'*(0:2));
        end
        
        function tf_lsf = get_tf_lsf(obj, G, fc)
            %METHOD1 이 메서드의 요약 설명 위치
            %   자세한 설명 위치
            V0 = 10^(G/20);
            H0 = V0 - 1;
            
            if G > 0 % boost
                a = (tan(pi*fc/obj.fs) - 1)/(tan(pi*fc/obj.fs) + 1);
            else % cut
                a = (tan(pi*fc/obj.fs) - V0)/(tan(pi*fc/obj.fs) + V0);
            end
            
            A1 = obj.z1*[a; 1]./obj.z1*[1; a];
            
            tf_lsf = 1 + H0/2*(1 + A1);
        end
    end
end

