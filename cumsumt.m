function out_snap = cumsumt(IC,in_avg,DT)
    [tL,Vl] = size(in_avg);
    out_snap = repmat(IC,[tL+1 1]) + ...
        cat(1,zeros(1,Vl),cumsum(in_avg.*repmat(DT,[1 Vl]),1));
end
