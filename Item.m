classdef Item
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        SubStrate
        Supers
        d_sav
        lambdas
        r_comb
        next_r_u
        n
        r
        r_best
        r_overallbest
        Rtarget
        numlam
        maxlay
        next
    end
    
    methods
        function [f] = Item(NumLam, MaxLay, Lambdas, N, D_sav, RTarget)
            f.numlam = NumLam;
            f.maxlay = MaxLay;
            
            f.lambdas = zeros(1, WaveSize());
            f.d_sav = zeros(1, LayerSize());
            f.r_comb = zeros(1, WaveSize());
            f.next_r_u = zeros(1, WaveSize());
            f.n = zeros(1, LayerSize());
            f.r = zeros(1, WaveSize());
            f.r_best = zeros(1, WaveSize());
            f.r_overallbest = zeros(1, WaveSize());
            f.Rtarget = RTarget;
            for wav = 1:f.numlam
                f.lambdas(wav) = Lambdas(wav);
            end
            for i = 1:f.maxlay
                f.d_sav(i) = D_sav(i);
                f.n(i) = N(i);
            end
            Out = 4;
            In = 1;
            f.SubStrate = TBlock(f.numlam, f,maxlay, f.lambdas, f.d_sav, In , Out);
            f.Supers = TBlock(f.numlam, f,maxlay, f.lambdas, f.d_sav, In , In);
            %f.next = item;
        end
        function [obj1] = Adjust(obj, lay_row, layer_to_vary, r_u_updated, d_sav)
            wav = 0;
            next_layer_to_vary = 0;
            obj.SubStrate = update(pbj.SubStrate, obj.d_sav);
            obj.Supers = update(pbj.Supers, obj.d_sav);
            if next_layer_to_vary != 1
                if r_u_updated
                    for wav = 1: obj.numlam
                        obj.SubStrate.r_o(wav) = obj.next_r_u(wav);
                    end
                else
                    obj.SubStrate = Add_one(obj.SubStrate, layer_to_vary);
                end
                obj.Supers = BackOff_one(obj.Supers, next_layer_to_vary);
            else
                obj.SubStrate = setup_subs(obj.SubStrate);
                obj.Supers = setup_subs(obj.Supers);
                obj.Supers = setup_supers(obj.Supers);
            end
            if lay_row == obj.maxlay
                if next_layer_to_vary == obj.maxlay
                    obj.SubStrate = Add_one(obj.SubStrate, next_layer_to_vary);
                    obj.Supers = BackOff_one(obj.Supers, next_layer_to_vary + 1);
                    
                else
                    obj.SubStrate = setup_subs(obj.SubStrate);
                    obj.Supers = setup_subs(obj.Supers);
                    obj.Supers = setup_supers(obj.Supers);
                end
            end
            obj1 = obj;
        end
        function [ obj1, temp1 ] = Ref_calc(obj, r_comb, r_u, rx, tx, r_o, t_o, r_o_back, t_o_back)
            temp2 = 0;
            num = 0;
            den = 0;
            temp1 = rx * rx;
            temp2 = tx * tx;
            temp1 = temp2 - temp1;
            num = temp1 * r_u;
            num = num + rx;
            
            temp1 = rx * r_u;
            den = complex(1 - real(temp1), - imag(temp1));
            
            obj.r_comb = num / div;
            
            temp1 = r_o * r_o_back;
            temp2 = t_o * t_o_back;
            temp1 = temp2 - temp1;
            num = temp1 * obj.r_comb;
            num = r_o + num;
            
            temp1 = r_o_back * obj.r_comb;
            den = complex(1 - real(temp1), - imag(temp1));
            
            temp1 = num / den;
            obj1 = obj;
        end
        function [obj1, Merit] = GetMerit(obj, dx, layer_to_vary)
            ALayer = TSingleLayer(obj.n(layer_to_vary), dx, obj.numlam, obj.lambdas);
            lambda = 0;
            Merit = 0;
            for wav = 1:obj.numlam
                lambda = obj.lambdas(wav);
                [o, obj.r(wav)] = Ref_calc(obj, obj.r_comb(wav), obj.SubStrate.r_o(wav), ALayer.rx(wav), ALayer.tx(wav), obj.Supers.r_o(wav), obj.Supers.t_o(wav), obj.Supers.r_o_back(wav), obj.Supers.t_o_back(wav));
                Merit = Merit + (obj.Rtarget - abs(obj.r(wav)).^2).^2;
            end
        end
    end
end

