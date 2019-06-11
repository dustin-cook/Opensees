function[th_CP, th_NC] = parametric_study(L,h,b,s,Av,cov,fc,fy,P)
    
    p_ratio = max(P/(b*h*fc),0.1);
    v_ratio = 1.0;             % typical value
    % min vratio = 0.02
    
    rho_t = min(Av/(b*s),0.0075);
    
    % ASCE 41
    ASCE_a = max(0.042 - 0.043*p_ratio + 0.63*rho_t - 0.023*v_ratio,0);
    if p_ratio <= 0.5
        ASCE_b = max(0.5/(5+p_ratio/0.8*(1/rho_t)*(fc/fy))-0.01,ASCE_a);
    else
        b5 = max(0.5/(5+0.5/0.8*(1/rho_t)*(fc/fy))-0.01,ASCE_a);
        ASCE_b  = max(interp1([0.5,0.7],[b5,0],p_ratio),ASCE_a);
    end
        
    
    ASCE_c = max(0.24-0.4*p_ratio,0);
    
    th_IO = min(0.15*ASCE_a,0.005);
    th_LS = 0.5*ASCE_b;
    th_CP = 0.7*ASCE_b;
    
    
    % EUROCODE
    covm = cov/39.37;   % m
    hm = h/39.37;       % m
    ho = hm - 2*covm;   % m
    bm = b/39.37;       % m
    bo = bm - 2*covm;   % m
    sm = s/39.37;       % m
    Lm = L/39.37;       % m
    fcm = fc/145.04;    % MPa
    fym = fy/145.04;    % MPa
    Pm = P/1000*4.4482; % kN
    
    gam_el = 1.8;
    rho_d = 0;              % diagonal reinforcement 
    sum_bi = 0.006;
    wc = 1;
    wt = 1;
    
    p_ratio_EU = (Pm/1000/(bm*hm*fcm));
%     rho_t_EU   = Av/(b*s);
   
    alpha = (1-sm/(2*bo))*(1-sm/(2*ho))*(1-(sum_bi^2)/(6*ho*bo));
    
    th_NC = (1/gam_el)*0.0145*(0.25^p_ratio_EU)*...
            ((max(0.01,wc)/max(0.01,wt))^0.3)*...
            (fcm^0.2)*((Lm/2/hm)^0.35)*(25^(alpha*rho_t*(fym/fcm)))*...
            (1.275^(100*rho_d));
        

end