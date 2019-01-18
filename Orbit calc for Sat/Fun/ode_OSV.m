function ode_OSV_ = ode_OSV(fun, t1, t2, u, Rel)
% tic
% Rel=1e-10;
t_i = t1;
% h = t2;%(t2-t1)*1;
h = 0.25;%(t2-t1)*0.01;
count_y = length(u);
out = 0;
out_out = 0;

meth = 1; % 1 - метод Рунге-Кутта 4го порядка; 2 - метод Рунге-Кутта-Мерсона

h_inside_before(1:count_y) = 0.;
while (t_i < t2)
    out = 0;
    i_out = 0;
    u0=u;
        
    while (out < 1)
        u=u0;
        h_inside(1:count_y) = h;

        if meth == 1
            if (t2 - t_i < h) 
                h = t2 - t_i;
            end;
%             k1 = h * feval(fun,t_i,u); %fun_orbit(t_i, u);
%             k2 = h * feval(fun,t_i + h/2., u + k1/2.); %fun_orbit(t_i + h/2., u + k1/2.);
%             k3 = h * feval(fun,t_i + h/2., u + k2/2.); %fun_orbit(t_i + h/2., u + k2/2.);
%             k4 = h * feval(fun,t_i + h, u + k3); %fun_orbit(t_i + h, u + k3);
%             u = u0 + (k1 + 2.*k2 + 2.*k3 + k4) / 6.;

            h3 = h/3.;
            k1 = h3*feval(fun,t_i,u); %fun_orbit(t_i,u);
            k2 = h3*feval(fun,t_i+h3,u+k1); %fun_orbit(t_i+h3,u+k1);
            k3 = h*feval(fun,t_i+h3,u+(k1+k2)*0.5); %fun_orbit(t_i+h3,u+(k1+k2)*0.5);
            k4 = k1+4.*h3*feval(fun,t_i+h*0.5,u+0.375*(k1+k3)); %fun_orbit(t_i+h*0.5,u+0.375*(k1+k3));
            k5 = h3*feval(fun,t_i+h,u+1.5*(k4-k3)); %fun_orbit(t_i+h,u+1.5*(k4-k3));
            u = u0 + (k4 + k5)*0.5;
            
            out = 1;
        elseif meth == 2
            h3 = h/3.;
            k1 = h3*feval(fun,t_i,u); %fun_orbit(t_i,u);
            k2 = h3*feval(fun,t_i+h3,u+k1); %fun_orbit(t_i+h3,u+k1);
            k3 = h*feval(fun,t_i+h3,u+(k1+k2)*0.5); %fun_orbit(t_i+h3,u+(k1+k2)*0.5);
            k4 = k1+4.*h3*feval(fun,t_i+h*0.5,u+0.375*(k1+k3)); %fun_orbit(t_i+h*0.5,u+0.375*(k1+k3));
            k5 = h3*feval(fun,t_i+h,u+1.5*(k4-k3)); %fun_orbit(t_i+h,u+1.5*(k4-k3));
            u = u0 + (k4 + k5)*0.5;
            Rel_h = abs(0.1*(2.*k4-3.*k3-k5));
            if out_out == 1
                out = 1;
            else
                out_inside = zeros(count_y,1);
                out_out_inside = zeros(count_y,1);
                for q = 1:count_y
                    if Rel_h(q) == 0.
                        out_inside(q) = 1;
                    elseif Rel_h(q) > Rel
%                         h_inside(q) = h*0.5;
                        h_inside(q) = h_inside(q)*0.5;
                    elseif Rel_h(q)*32. < Rel
%                         if h > t2-t_i
                        if h_inside(q) >= t2-t_i
                            h_inside(q) = t2-t_i;
                            out_out_inside(q) = 1;
                        else
%                             h_inside(q) = h*2.;
                            h_inside(q) = h_inside(q)*2.;
                        end
                    else
                        out_inside(q) = 1;
                    end

                    if q == 1
                        h =  h_inside(q);
                    else
                        if h > h_inside(q)
                            h = h_inside(q);
                        end
                    end
                end
                
                if mean(out_inside) == 1
                    out = 1;
                elseif mean(out_out_inside) == 1
                    out_out = 1;
                elseif mean(h_inside_before-h_inside) == 0
                    out = 1;
                end
                h_inside_before = h_inside;
            end
            i_out = i_out + 1;
            if i_out > 1000
                out = 1;
            end
        end
    end
    if out == 1
        t_i = t_i + h;
    end
end
% toc
ode_OSV_ = u;
% end
% 
% % % 
% 
