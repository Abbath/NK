classdef TSingleLayer
    properties
        thickness
        Lambdas
        ref_index
        rx
        tx
        f
        numlam
    end
    methods
        function [f] = TSingleLayer(Ref_index, Thick, NumLam, Waveths)
            if nargin > 0
                f.ref_index = Ref_index;
                f.thickness = Thick;
                f.numlam = NumLam;
                f.Lambdas = zeros(1, WaveSize());
                f.rx = zeros(1, WaveSize());
                f.tx = zeros(1, WaveSize());
                for i = 2:f.numlam+1
                    f.Lambdas(i) = Waveths(i);
                    [f.rx(i), f.tx(i)] = Calc_layer_x(f, f.Lambdas(i));
                end
                f.f = Fresnel(f.ref_index);
            end
        end
        function [rx,tx] = Calc_layer_x(obj, lambda)
            ff = Fresnel([obj.ref_index]);
            fsq = ff * ff;
            fac = 2 * pi * [obj.thickness] / lambda;
            ep = fac * [obj.ref_index];
            ep = exp(ep);
            epsq = ep*ep;
            den = fsq * epsq;
            den = complex(1 - real(den), -imag(den));
            
            num =  epsq * ff;
            num = num - ff;
            rx = num / den;
            
            num = complex(1 - real(fsq), - imag(fsq));
            num = ep*num;
            tx = num / den;
        end
        
    end
end