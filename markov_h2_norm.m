function norma = markov_h2_norm(Lambda, A, E, C, mu)
    N = numel(A);
    n = size(A{1},1);
    
    P = zeros(1,N);
    
    setlmis([]);
    
    for i=1:N
        P(i) = lmivar(1,[n, 1]);
    end
    
    ct = 1;
    for i=1:N
        lmiterm([-ct 1 1 P(i)],1,1);
        ct = ct + 1;
        lmiterm([ct 1 1 P(i)],A{i}',1,'s');
        lmiterm([ct 1 1 0],C{i}'*C{i});
        for j=1:N
            lmiterm([ct 1 1 P(j)],Lambda(i,j),1);
        end
        ct = ct + 1;
    end

    lmisys = getlmis;

    c = zeros(decnbr(lmisys),1);
    for idx = 1:decnbr(lmisys)
        for i = 1:N
            vPi = defcx(lmisys,idx,P(i));
            c(idx) = c(idx) + mu(i)*trace(E{i}'*vPi*E{i});
        end
    end

    options = [1e-7,2000,0,200,0];
    copt = mincx(lmisys,c,options);
    if(~isempty(copt))
        norma = copt;
    else
        norma = inf;
    end
end
