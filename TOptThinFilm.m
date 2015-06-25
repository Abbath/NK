classdef TOptThinFilm
    %UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        WantPrint
        SaveStack
        r_u_updated
        max_film_thickness
        %i
        %iter
        %lay
        next_layer_to_vary
        layer_to_vary
        layer_to_start
        %lay_row
        %lay_col
        %thick_count
        %wav
        numlam
        maxlay
        num_iterations
        numthix
        numthix_plus1
        Merit
        BestMerit
        OverallBestMerit
        maxlam
        minlam
        scale_down
        d_start
        dx
        lambda
        midlam
        r
        temp
        iteration
        d_sav
        d_orig
        opthix_range
        d_inc
        lambdas
        d_best
        n
        n_opt
        Name
        FirstIter
        SecondIter
        ThirdIter
        NumberOfInput
        AcceptableResult
        delta
        Pat
        aver1
        aver2
        aver3
        head
    end
    
    methods
        function f = TOptThinFilm()
            %f = TOptThinFilm;
            f.max_film_thickness = 5;
            f.FirstIter = 59;
            f.SecondIter = 17;
            f.ThirdIter = 7;
            f.AcceptableResult = 0.25;
            f.delta = 0.4;
            f.num_iterations = 3;
            f.layer_to_start = 2;
            f.NumberOfInput = 2;
            f.maxlay = 8;
            f.minlam = 4.4;
            f.maxlam = 5.6;
            f.numlam = 1;
            f.iteration = zeros(1,LayerSize());
            f.d_sav = zeros(1,LayerSize());
            f.d_orig = zeros(1,LayerSize());
            f.opthix_range = zeros(1,LayerSize());
            f.d_inc = zeros(1,LayerSize());
            f.d_best = zeros(1,LayerSize());
            f.n = zeros(1,LayerSize());
            f.lambdas = zeros(1,LayerSize());
            f.n_opt = zeros(1,LayerSize());
            f.Name = 'sample'; %input('Enter file name: ', 's');
            file = fopen(strcat(f.Name,'.stk'), 'r');
            data = fscanf(file, '%f');
            f.maxlay = data(1);
            counter = 3;
            for lay = 1:f.maxlay
                Re = data(4+((lay-1)*3));
                Im = data(5++((lay-1)*3));
                f.n(lay) = complex(Re, Im);
                f.d_orig(lay) = data(6+((lay-1)*3));
                f.n_opt(lay) = f.n(lay);
                counter = counter + 3;
            end
            f.NumberOfInput = data(counter+1);
            f.delta = data(counter+2);
            f.minlam = data(counter+3);
            %LamFac = 0;
            if f.numlam == 1
                LamFac = 0;
            else
                LamFac = (f.maxlam - f.minlam) / (f.numlam - 1);
            end
            for wav = 1:f.numlam
                f.lambdas(wav) = f.minlam + LamFac * (wav - 1);
            end
            for lay = 1:f.maxlay
                f.iteration(lay) = 1;
            end
        end
        function [obj1] = ReadInput(obj)
            file = fopen(strcat(obj.Name, '.fct'), 'r');
            data = fscanf(file, '%f');
            %Rtarget = 0;
            count = 0;
            counter = 0;
            obj.Pat = zeros(1, 30);
            while 1
                if feof(file) || (count == 16)
                    break;
                end
                count = count + 1;
                counter = counter + 1;
                for j = 1:obj.NumberOfInput
                    obj.Pat(j) = data(counter);
                    counter = counter + 1;
                end
                Rtarget = data(counter);
                for j = 1:obj.maxlay
                    obj.n(j) = complex(real(obj.n_opt(j)) + obj.delta * obj.Pat(j),imag(obj.n(j)));
                end
                pt = Item(obj.numlam, obj.maxlay, obj.lambdas, obj.n, obj.d_sav, Rtarget);
                obj = append(obj, pt);
            end
            obj1 = obj;
        end
        function [ obj1, b ] = encode(obj)
            obj.OverallBestMerit = 999.9;
            obj.d_sav = CopyStack(obj.maxlay, obj.d_orig);
            obj.midlam = 0.5*(obj.lambdas(1) + obj.lambdas(obj.numlam));
            for lay = 1:obj.maxlay
                obj.opthix_range(lay) = obj.midlam;
                if obj.opthix_range(lay) < real(obj.n(lay)) * obj.d_sav(lay)
                    obj.opthix_range(lay) = 1.2 * real(obj.n(lay)) * obj.d_sav(lay);
                end
                if obj.opthix_range(lay) > real(obj.n(lay)) * obj.max_film_thickness;
                    obj.opthix_range(lay) = real(obj.n(lay)) * obj.max_film_thickness;
                end
            end
            obj = ReadInput(obj);
            for iter = 1 : obj.num_iterations
                if iter == 1
                    obj.numthix = obj.FirstIter;
                elseif iter == 2
                    obj.numthix = obj.SecondIter;
                else
                    obj.numthix = obj.ThirdIter;
                end
                obj.numthix_plus1 = obj.numthix + 1;
                obj.scale_down = 1.1/obj.numthix;
                obj.d_inc = calc_d_inc(obj, obj.maxlay, obj.numthix, obj.opthix_range, obj.n);
                obj.layer_to_vary = obj.layer_to_start;
                for lay_col = 1:obj.maxlay
                    obj.BestMerit = 999.9;
                    if iter == 1
                        obj.d_sav = CopyStack(obj.maxlay, obj.d_orig);
                    else
                        for lay = 1:obj.maxlay
                            obj.d_sav(lay) = obj.d_orig(lay);
                        end
                    end
                    for i = 1:(length(obj.head)-1)
                        obj.head(i) = Initialize(obj.layer_to_vary, obj.d_sav);
                    end
                    for lay_row = 1:obj.maxlay
                        obj.r_u_updated = 0;
                        if obj.iteration(obj.layer_to_vary) <= iter
                            obj.dx = calc_dstart(obj, obj.d_sav(obj.layer_to_vary), obj.d_inc(obj.layer_to_vary), obj.numthix, obj.max_film_thickness);
                            for thick_count = 1:obj.numthix_plus1
                                obj.Merit = 0.0;
                                obj.Merit = TotalMerit(obj, obj.dx, obj.layer_to_vary);
                                if obj.Merit < obj.BestMerit
                                    obj.r_u_updated = 1;
                                    obj.BestMerit = obj.Merit;
                                    obj.d_sav(obj.layer_to_vary) = obj.dx;
                                    for i = 1:(length(obj.head) - 1)
                                        for j = 1:obj.numlam
                                            tmp = obj.head(i);
                                            tmp.next_r_u(j) = tmp.r_combo(j);
                                            tmp.r_best(j) = tmp.r(j);
                                            obj.head(i) = tmp;
                                        end
                                    end
                                end
                                obj.dx = obj.dx + obj.d_inc(obj.layer_to_vary);
                            end
                        end
                        if obj.BestMerit < obj.OverallBestMerit
                            obj.OverallBestMerit = obj.BestMerit;
                            for lay = 1:obj.maxlay
                                obj.d_best(lay) = obj.d_sav(lay);
                            end
                            for i = 1:(length(obj.head) - 1)
                                for j = 1:obj.numlam
                                    obj.head(i).r_overallbest(j) = obj.head(i).r_best(j);
                                end
                            end
                            
                        end
                        for i = 1:(length(obj.head) - 1)
                            obj.head(i) = Adjust(lay_row, obj.layer_to_vary, obj.r_u_updated, obj.d_sav);
                        end
                        obj.layer_to_vary = mod(obj.layer_to_vary, obj.maxlay) + 1;
                    end
                    obj.layer_to_vary = mod(obj.layer_to_vary, obj.maxlay) + 1;
                end
                disp(obj.OverallBestMerit);
                for lay = 1:obj.maxlay
                    obj.opthix_range(lay) = obj.opthix_range(lay) * obj.scale_down;
                end
            end
            for lay = 1:obj.maxlay
                obj.d_sav(lay) = obj.d_best(lay);
            end
            for i = 1:(length(obj.head) - 1)
                for j = 1:obj.numlam
                    obj.head(i).r(j) = obj.head(i).r_overallbest(j);
                end
            end
            if obj.OverallBestMerit <= obj.AcceptableResult
                b = 1;
            else
                b = 0;
            end
            obj1 = obj;
        end
        function [t, obj] = training(obj)
            t = 1;
        end
        function [d_inc] = calc_d_inc(~, maxlay, numthix, opthix_range, ndx)
            d_inc = zeros(1,maxlay);
            for lay = 1:maxlay
                d_inc(lay) = opthix_range(lay)/(real(ndx(lay))*numthix);
            end
        end
        function [ result ] = calc_dstart(~, d_saved, d_inc, numthix, max_film_thick)
            %thick_count = 0;
            result = d_saved + d_inc * fix(numthix / 2);
            while result > max_film_thick
                result = result - d_inc;
            end
            result = result - d_inc * numthix;
            while result < 0
                result = result + d_inc;
            end
        end
        function [r_u, t_u] = setup_subs(~, NumLam, ndx)
            r_u = zeros(1,NumLam);
            t_u = zeros(1,NumLam);
            f = Fresnel(ndx(1));
            denom = complex(real(ndx(1)) + 1, imag(ndx(1)));
            two = 2;
            for wav = 1:NumLam
                r_u(wav) = complex(-real(f), -imag(f));
                t_u(wav) = two / denom;
            end
        end
        function [r_u_test] = recall(obj,n)
            %r_u = zeros(1,LayerSize());
            %t_u = zeros(1,LayerSize());
            %f = zeros(1,LayerSize());
            %alfa = zeros(1,LayerSize());
            [f, alfa] = LayerParameters(obj.maxlay, n);
            [r_u, t_u] = setup_subs(obj, obj.numlam, obj.n);
            for wav = 1: obj.numlam
                obj.lambda = obj.lambdas(wav);
                theta = zeros(1,LayerSize());
                for lay = 1:obj.maxlay
                    theta(obj.lay) = (4*pi)*real(n(obj.lay))*obj.d_sav(obj.lay)/obj.lambda;
                end
                Zstak(1, r_u(wav), t_u(wav), alfa, f, theta, obj.maxlay, real(n(1)), 1);
            end
            r_u_test = zeros(1,obj.numlam);
            for wav = 1:obj.numlam
                r_u_test(wav) = obj.r_u(wav);
            end
        end
        function [b] = is_empty(obj)
            b = isempty(obj.head);
        end
        function [obj1] = append(obj, item)
            obj1 = obj;
            obj1.head = [obj1.head, item];
        end
        function [obj1] = remove(obj)
            obj1 = obj;
            obj1.head = [];
        end
        function [Merit] = TotalMerit(obj, dx, layer_to_vary)
            count = 0;
            TMerit = 0;
            for i = 1:length(obj.head)
                TMerit = TMerit + GetMerit(obj.head(i), dx, layer_to_vary);
            end
            Merit = TMerit / count;
        end
        function [] = savethick(obj)
            file = fopen(strcat(obj.Name, '.THK'), 'w');
            fprintf(file, '%s\n%d %f %f\n', obj.Name, obj.maxlay, obj.n_opt(obj.maxlay), obj.n_opt(1));
            for lay = 1:obj.maxlay
                fprintf(file, '%f %f %f\n', real(obj.n_opt(lay)), imag(obj.n_opt(lay)), obj.d_sav(lay));
            end
            fclose(file);
        end
        function [] = result(obj)
            file = fopen(strcat(obj.Name, '.RST'), 'w');
            d = clock;
            fprintf(file, '%s\n', obj.Name);
            fprintf(file, 'File date: %d.%d.%d\n', d(1), d(2), d(3));
            fprintf(file, 'File time: %d:%d:%d\n', d(4), d(5), fix(d(6)));
            fprintf(file, 'OverallBestMerit: %f\nFirstIter: %d\nDelta value: %f\nlambdas: ', obj.OverallBestMerit, obj.FirstIter, obj.delta);
            for wav = 1:obj.numlam
                fprintf(file, '%f:', obj.lambdas(wav));
            end
            fprintf(file, '\n');
            count = 0;
            for i = 1:length(obj.head)
                count = count + 1;
                fprintf(file, '%f:', obj.head(i).Rtarget);
                for wav = 1:obj.numlam
                    fprintf(file, '{%f}', (obj.head(i).r(wav))^2);
                end
                fprintf(file, '\n');
            end
            fclose(file);
        end
    end
end


