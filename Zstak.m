function [ z1, t1 ] = Zstak( Layer, z, t, alfa, f, theta, MaxLay, Nin, want_transmission) )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
n_fac = 0;
p = 0;
p_factor = 0;
z_factor = 0;
F = 0;
negF = 0;
half_thet = 0;
den = 0;
temp = 0;
if Layer <= MaxLay+1
    p = Compos(alfa(Layer), z);
    z = exp(theta(Layer)) * p;
    z1 = z;
    if want_transmission
        if Layer < MaxLay+1
            half_thet = theta(Layer) / 2;
            p_factor = f(Layer) * p;
            p_factor = complex(1 - real(p_factor), - imag(p_factor));
            z_factor = f(Layer) * p;
            z_factor = complex(1 - real(z_factor), - imag(z_factor));
            temp = t * (p_factor/z_factor);
            t = exp(half_thet) * temp;
            t1 = t;
        else
            negF = -f(Layer);
            den = Compos(negF. p);
            den = f(Layer) * den;
            den = complex(1 + real(den), imag(den));
            temp = t / den;
            n_fac = 2 * Nin / (Nin+1);
            t = temp * n_fac;
            t1 = t;
        end
    end
    [z1, t1] = Zstak(Layer + 1, z, t, alfa, f , theta, MaxLay, Nin, want_transmission);
end
end

