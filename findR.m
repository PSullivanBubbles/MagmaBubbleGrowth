function Rout = findR(Rin, tin, t)

    if t ==0 
        Rout = Rin(1);
        return ;
    end
    if t >tin(end) 
        Rout = Rin(end);
        return ;
    end
    if t<tin(1)
        R0= 1e-5;
        lever= (tin(1)-t)/(tin(1));
        Rout = R0 + lever*(Rin(1)-R0);
        return;
    end

    for i=1:(size(tin,1)-1)
        if (t>tin(i))&&(t<tin(i+1))
            lever= (t-tin(i))/(tin(i+1)-tin(i));
            Rout = Rin(i) + lever*(Rin(i+1)-Rin(i));

            return;
        end
    end
    Rout = Rin(end);
end