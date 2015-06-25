function [ z1, t1 ] = Zstak( Layer, z, t, alfa, f, theta, MaxLay, Nin, want_transmission)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
for i = 1:MaxLay
    %if Layer <= MaxLay
        p = Compos(alfa(Layer), z);
        z = exp(theta(Layer)) * p;
        if want_transmission
            if Layer < MaxLay
                half_thet = theta(Layer) / 2;
                p_factor = f(Layer) * p;
                p_factor = complex(1 - real(p_factor), - imag(p_factor));
                z_factor = f(Layer) * p;
                z_factor = complex(1 - real(z_factor), - imag(z_factor));
                temp = t * (p_factor/z_factor);
                t = exp(half_thet) * temp;
            else
                negF = -f(Layer);
                den = Compos(negF, p);
                den = f(Layer) * den;
                den = complex(1 + real(den), imag(den));
                temp = t / den;
                n_fac = 2 * Nin / (Nin+1);
                t = temp * n_fac;
            end
        end
        %[z, t] = Zstak(Layer + 1, z, t, alfa, f , theta, MaxLay, Nin, want_transmission);
    %end
end
z1 = z;
t1 = t;
end

