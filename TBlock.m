classdef TBlock
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        MLayers
        NumItem
        numlam
        d_sav
        lambdas
        n
        Out_F
        In_F
        Base
        alfa
        r_o
        t_o
        r_o_back
        t_o_back
    end
    
    methods
        function [f] = TBlock(NumLam, NumLay, Lambdas, N, D_sav, In, Out)
            f.MLayers(1,LayerSize()+2) = TSingleLayer;
            f.n = zeros(1, LayerSize());
            f.alfa = zeros(1, LayerSize());
            f.d_sav = zeros(1, LayerSize());
            f.lambdas = zeros(1, WaveSize());
            f.r_o = zeros(1, WaveSize());
            f.r_o_back = zeros(1, WaveSize());
            f.t_o_back = zeros(1, WaveSize());
            f.t_o = zeros(1, WaveSize());
            f.NumItem = NumLay;
            f.numlam = NumLam;
            f.Base = Out;
            f.In_F = Fresnel(In);
            f.Out_F = Fresnel(Out);
            for lay = 1:NumLay
                f.n(lay) = N(lay);
                f.d_sav(lay) = D_sav(lay);
            end
            for wav = 1:f.numlam
                f.lambdas(wav) = Lambdas(wav);
            end
        end
        function [obj1] = update(obj, D_sav)
            for lay = 1:obj.NumItem
                obj.d_sav(lay) = D_sav(lay);
            end
            obj1 = obj;
        end
        function [obj1] = AlfaParams(obj)
            obj.alfa(3) = obj.MLayers(3).f;
            for lay = 4:(obj.NumItem-1)
                negF = obj.MLayers(lay-1).f;
                temp_F = obj.MLayers(lay).f;
                obj.alfa(lay) = Compos(temp_F, negF);
            end
            negF = obj.MLayers(obj.NumItem).f;
            obj.alfa(obj.NumItem) = Compos(obj.In_F, negF);
            obj1 = obj;
        end
        function [obj1] = Create_One_Layer(obj, layer)
            obj1 = obj;
            new_layer = TSingleLayer(obj.n(layer), obj.d_sav(layer), obj.numlam, obj.lambdas);
            obj1.MLayers(1,length(obj.MLayers)+1) = TSingleLayer;
            for i = 1:(layer-1)
                obj1.MLayers(i) = obj.MLayers(i);
            end
            obj1.MLayers(layer) = new_layer;
            for i = (layer+1):length(obj1.MLayers)
                obj1.MLayers(i) = obj.MLayers(i-1);
            end
        end
        function [obj1] = substructure(obj, direction, wav)
            pseud_nin = 1.0;
            Ndx = zeros(1, LayerSize());
            theta = zeros(1, LayerSize());
            bottom_level = 2;
            if ~direction
                l_bias = obj.NumItem + bottom_level;
            else
                l_bias = 0;
            end
            thet_fac = pi*4/obj.lambdas(wav);
            Maxlay_plus1 = obj.NumItem + 1;
            Ndx(Maxlay_plus1) = 1;
            Ndx(bottom_level - 1) = 1;
            theta(Maxlay_plus1) = 0;
            if bottom_level <= obj.NumItem
                for lay = bottom_level:obj.NumItem
                    level = l_bias + (2*direction - 1) * lay;
                    Ndx(level) = obj.n(lay);
                    if direction
                        obj1 = Create_One_Layer(obj, lay);
                    end
                    theta(level) = thet_fac*obj.n(lay)*obj.d_sav(lay);
                end
                [f, alpha ] = LayerParameters(obj.NumItem, Ndx);
                if direction
                    [z1, t1] = Zstak(bottom_level, obj.r_o(wav), obj.r_o(wav), alpha, f, theta, obj.NumItem, pseud_nin, 1);
                    obj1.r_o(wav) = z1;
                    obj1.t_o(wav) = t1;
                else
                    [z1, t1] = Zstak(bottom_level, obj.r_o_back(wav), obj.r_o_back(wav), alpha, f, theta, obj.NumItem, pseud_nin, 1);
                    obj1.r_o(wav) = z1;
                    obj1.t_o(wav) = t1;
                end
            end
        end
        function [obj1] = calc_superstructure(obj, wav)
            normal = 1;
            reverse = 0;
            obj1 = substructure(obj, normal, wav);
            obj1 = substructure(obj1, reverse, wav);
        end
        function [obj1] = setup_supers(obj)
            obj1 = obj;
            for wav = 1:obj.numlam;
                obj1 = calc_superstructure(obj1, wav);
            end
        end
        function [obj1] = Add_one(obj, layer)
            Nin = 1.0;
            Ndx = zeros(1,3);
            theta = zeros(1,3);
            Ndx(3) = 1;
            Ndx(2) = obj.n(layer);
            betax = Ndx(2);
            theta(3) = 0;
            [f, alpha] = LayerParameters(1, Ndx);
            obj1 = Create_One_Layer(obj, layer);
            for wav = 1:obj1.numlam
                theta(2) = 4*pi*betax*obj1.d_sav(layer)/obj1.lambdas(wav);
                [z1, t1] = Zstak(1, obj1.r_o(wav), obj1.t_o(wav), alpha, f, theta, 1, Nin, 1);
                obj1.r_o(wav) = z1;
                obj1.t_o(wav) = t1;
            end
        end
        function [obj1] = Back_Off_one(obj, layer)
            Nin = 1.0;
            Ndx = zeros(1,3);
            theta = zeros(1,3);
            if layer > 1 && layer < obj.NumItem 
                N = obj.MLayers(layer).ref_index;
                thick = obj.MLayers(layer).ref_index;
                Ndx(3) = 1;
                Ndx(2) = N;
                theta(3) = 0;
                [f, alpha] = LayerParameters(1, Ndx);
                d_factor = -4*pi*thick;
                for wav = 1:obj.numlam
                    theta_fac = d_factor / obj.lambdas(wav);
                    theta(2) = theta_fac*N;
                    [z1, t1] = Zstak(1, obj.r_o_back(wav), obj.t_o_back(wav), alpha, f, theta, 1, Nin, 1);
                    obj.r_o_back(wav) = z1;
                    obj.t_o_back(wav) = t1;
                    
                    rxx = obj.MLayers(layer).rx(wav);
                    txx = obj.MLayers(layer).tx(wav);
                    
                    r = obj.r_o(wav);
                    t = obj.t_o(wav);
                    
                    c_fac = obj.r_o_back(wav) * rxx;
                    c_fac = complex(1 - real(c_fac), - imag(c_fac));
                    
                    temp = t * c_fac;
                    obj.t_o(wav) = temp / txx;
                    num = obj.t_o(wav) * obj.t_o(wav);
                    temp = num / c_fac;
                    obj.r_o(wav) = r - temp;
                end
                obj1 = obj;
                obj1.MLayers(1, length(obj.MLayers)-1) = TSingleLayer;
                for i = 1:(layer-1)
                    obj1.MLayers(i) = obj.MLayers(i);
                end
                for i = (layer+1):length(obj1.MLayers);
                    obj1.MLayers(i) = obj.MLayers(i+1);
                end
            end
        end
        function [obj1] = setup_subs(obj)
            denom = obj.Base + 1;
            two = 2;
            for wav = 1:obj.numlam
                obj.r_o(wav) = -obj.Out_F;
                temp = two / denom;
                obj.t_o(wav) = temp;
                obj.r_o_back(wav) = 0;
                obj.t_o_back(wav) = 1;
            end
            obj1 = obj;
        end
        
    end
end

