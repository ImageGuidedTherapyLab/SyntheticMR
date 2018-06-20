
figure;

nparam = size(post_samples,2)-2;
labels = {'M0_{GM}','T1_{GM}','T2_{GM}','M0_{WM}','T1_{WM}','T2_{WM}','M0_{CSF}','T1_{CSF}','T2_{CSF}'};
% count = 1;

% while count <= nparam*nparam
    for i = 1:nparam
        for j = 1:nparam
            if i == j
                subplot(nparam, nparam, nparam*(i-1)+j)%count)
%                 histogram(post_samples(:,i), 50)
                wp = [i];
                posteriors_nofig(post_samples, wp, {prior{wp,1}});
                hold on; yyaxis right; plot(prior{i,3}-1*prior{i,4}:prior{i,4}/100:prior{i,3}+1*prior{i,4},normpdf(prior{i,3}-1*prior{i,4}:prior{i,4}/100:prior{i,3}+1*prior{i,4},prior{i,3},prior{i,4}),'r');
%                 subplot(nparam, nparam, count)
%                 scatter(post_samples(:,i), post_samples(:,nparam+1), '.')
            end

            if i < j
                subplot(nparam, nparam, nparam*(i-1)+j)%count)
                histogram2(post_samples(:,i), post_samples(:,j), [50 50], 'FaceColor','flat');
            end

            if i > j
                subplot(nparam, nparam, nparam*(i-1)+j)%count)
                wp = [i j];
                posteriors_nofig(post_samples, wp, {prior{wp(1),1}, prior{wp(2),1}});
%                 scatter3(post_samples(:,i), post_samples(:,j), post_samples(:,nparam+1), '.')
            end
            if j == 1;
                yyaxis left; ylabel(labels{i});
            end
            if i == 9
                xlabel(labels{j});
            end
%             count = count + 1;
        end
    end
% end

% wp = [1];
% posteriors(post_samples, wp, {prior{wp,1}});
% wp = [2];
% posteriors(post_samples, wp, {prior{wp,1}});
% wp = [3];
% posteriors(post_samples, wp, {prior{wp,1}});
% wp = [1 2];
% posteriors(post_samples, wp, {prior{wp(1),1}, prior{wp(2),1}});

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
