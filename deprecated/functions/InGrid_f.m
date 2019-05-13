function grid = InGrid_f(polysize,polyx,polyy,position,idx)


    nrx = max(polyx)/polysize;
    nry = max(polyy)/polysize;
    n = 1;
    in = inpolygon(position(logical(idx),1),position(logical(idx),2),polyx,polyy);
    if nrx==1 & nry==1
        grid{max(polyx)/polysize,max(polyy)/polysize} = in;
    else
        while n <= 4
            switch n
                case 1
                    polyx_sub = [min(polyx),min(polyx),...
                        min(polyx)+ceil(nrx/2)*polysize,min(polyx)+ceil(nrx/2)*polysize,...
                        min(polyx)];
                    polyy_sub = [min(polyx),min(polyy)+ceil(nry/2)*polysize,...
                        min(polyy)+ceil(nry/2)*polysize,min(polyx),min(polyx)];
                case 2
                    polyx_sub = [min(polyx)+ceil(nrx/2)*polysize,min(polyx)+ceil(nrx/2)*polysize,...
                        max(polyx),max(polyx),...
                        min(polyx)+ceil(nrx/2)*polysize];
                    polyy_sub = [min(polyx),min(polyy)+ceil(nry/2)*polysize,...
                        min(polyy)+ceil(nry/2)*polysize,min(polyx),min(polyx)];
                case 3
                    polyx_sub = [min(polyx),min(polyx),...
                        min(polyx)+ceil(nrx/2)*polysize,min(polyx)+ceil(nrx/2)*polysize,...
                        min(polyx)];
                    polyy_sub = [min(polyy)+ceil(nry/2)*polysiz,max(polyy),...
                        max(polyy),min(polyy)+ceil(nry/2)*polysiz,min(polyy)+ceil(nry/2)*polysiz];
                case 4
            end
            
            InGrid_f(polysize,polyx_sub,polyy_sub,position,idx(in));
            n = n+1;
        end
    
    end
end