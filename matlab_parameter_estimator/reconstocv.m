% Reconstructs an OCV array from the measured voltage vb
r_edges = find(ib(2:end) & ~ib(1:end-1))-10;
r_edges = r_edges(r_edges>ix(1) & r_edges<ix(end));

ib = flip(ib);
f_edges = find(ib(2:end) & ~ib(1:end-1))-10;
f_edges = length(ib) - f_edges;
f_edges = f_edges(f_edges>ix(1) & f_edges<ix(end));
f_edges = flip(f_edges);
ib = flip(ib);

vb_temp = vb;
vb_temp(1:r_edges(1)) = vb(r_edges(1)); %fill everything until the first rising edge
for j=1:length(f_edges)

    v_start = vb(r_edges(j));
    if j==length(f_edges)
        v_end = vb(ix(end));
    else
        v_end = vb(r_edges(j+1));
        vb_temp(f_edges(j): r_edges(j+1)) = v_end; %fill after the falling edge
    end
    %interpolate inbetween
    vb_temp((r_edges(j)+1) : (f_edges(j)-1)) = linspace(v_start, v_end, f_edges(j)-r_edges(j)-1);
    
end
vb_temp(f_edges(j):end) = v_end; %fill after the last falling edge

 
clear v_end
clear v_start
clear j
clear f_edges
clear r_edges

%OUTPUT: vb_temp