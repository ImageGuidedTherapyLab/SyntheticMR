% % 7 by 7 matrix: Rl, Rs, Rt, Ro, Rlo, Rso, Rto
nparam = 3;

count = 1;

while count <= nparam*nparam
    for i = 1:nparam
        for j = 1:nparam
            if i == j
                subplot(nparam, nparam, count)
                histogram(post_samples(:,i), 50)
                
%                 subplot(nparam, nparam, count)
%                 scatter(post_samples(:,i), post_samples(:,nparam+1), '.')
            end

            if i < j
                subplot(nparam, nparam, count)
                histogram2(post_samples(:,i), post_samples(:,j), [50 50], 'FaceColor','flat');
            end

            if i > j
                subplot(nparam, nparam, count)
                scatter3(post_samples(:,i), post_samples(:,j), post_samples(:,nparam+1), '.')
            end

            count = count + 1;
        end
    end
end


% subplot(7, 7, 2)
% histogram2(post_samples(:,1), post_samples(:,2), [50 50]) %, 'FaceColor','flat');
% 
% 
% for i = 1:7
% 	for j = i:7
%         scatter3(post_samples(:,i), post_samples(:,j), post_samples(:,8), '.')
% 
%     end
% end
% 
