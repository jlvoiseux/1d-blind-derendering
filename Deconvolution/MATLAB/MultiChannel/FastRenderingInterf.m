%% Produces an interferometric render using interferometric values

function l_out = FastRenderingInterf(source_interf, brdf_interf, num_ang, n, T, indexes, temp1)  
    %[~, col] = ind2sub([2*T-1, 2*(2*num_ang-1)-1], indexes);
    temp = zeros(2*T-1, n);
    for j=round(T/2):2*T-1-round(T/2)
        temp5 = zeros(1, n);
        for i=1:sqrt(n)
            for k=1:sqrt(n)
                temp3 = zeros(size(temp1));
                temp3(temp1 > 0.5) = source_interf(:,get_col_num(i,k,sqrt(n))) + source_interf(:,get_col_num(k,i,sqrt(n))); 
                temp4 = conv(temp3, brdf_interf);
                temp5(get_col_num(i,k,sqrt(n))) = temp4(indexes(j));
            end
        end
        temp(j, :) = temp5;
    end
    l_out = zeros(2*T-1, n);
    l_out_temp = zeros(1, 2*T-1);
    for i=1:n
        l_out_temp(end-(1:2*T-1)+1) = temp(:, i);
        l_out_temp = l_out_temp(round(T/2):end-round(T/2));
        l_out_temp = interp(l_out_temp,2);
        if(length(l_out_temp) < 2*T-1)
            diff = 2*T-1-length(l_out_temp);
            l_out_temp = [l_out_temp, zeros(1, diff)];
        elseif(length(l_out_temp) > 2*T-1)
            diff = length(l_out_temp) - (2*T-1);
            l_out_temp = l_out_temp(1:end-diff);
        end
        l_out(:, i) = flip(l_out_temp)';
    end    
end

function out = get_col_num(i,j, obs_num)
    out = (i-1).*obs_num + j;
end
