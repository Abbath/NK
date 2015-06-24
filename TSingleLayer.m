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
            f.ref_index = Ref_index;
            f.thickness = Thick;
            f.numlam = NumLam;
            f.Lambdas = zeros(1, WaveSize());
            f.rx = zeros(1, WaveSize());
            f.tx = zeros(1, WaveSize());
            for i = 1:f.numlam
                Lambdas(i) = Waveths(i);
                [rx(i), tx(i)] = Calc_layer_x(f, Lambdas(i));
            end
            f.f = Fresnel(f.ref_index);
        end
        function [rx,tx] = Calc_layer_x(obj, lambda)
            fac = 0;
            f = Fresnel([obj.ref_index]);
            fsq = f*f;
            fac = 2 * pi * [obj.thickness] / lambda;
            ep = fac * [obj.ref_index];
            ep = exp(ep);
            epsq = ep*ep;
            den = fsq * epsq;
            den = complex(1 - real(den), -imag(den)); 
            
            num =  epsq*f;
            num = num - f;
            rx = num / den;
            
            num = complex(1 - real(fsq), - imag(fsq));
            num = ep*num;
            tx = num / den;
        end
        
    end
end