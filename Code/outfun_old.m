function stop = outfun(x,optimValues,state)
stop = false;
phistory.fval=[];
phistory.x=[];
grad=[];
figmarkstr={'bo','gx','r+','c*','ms','yd','kv','w^'};
switch state
    case 'init'
        figure(1);
        hold on
        figure(2);
        hold on
    case 'iter'
        % Concatenate current point and objective function
        % value with history. x must be a row vector.
        phistory.fval = [phistory.fval; optimValues.fval];
        phistory.x = [phistory.x; x];
        % Concatenate current search direction with
        % searchdir.
        grad = [grad;...
            optimValues.gradient'];
        %            plot(x(1),x(2),'o');
        figure(1);
        plot(optimValues.iteration,-optimValues.fval,'o');
        figure(2);
        yyaxis left
        plot(optimValues.iteration,x(1),'o');
        yyaxis right
        for iii=2:length(x)
            plot(optimValues.iteration,x(iii),figmarkstr{iii});
        end
        %            plot(optimValues.iteration,x(2),'o');
        %            plot(optimValues.iteration,x(3),'o');
        %            plot(optimValues.iteration,x(4),'o');
        %            plot(optimValues.iteration,x(5),'o');
        %            plotyy(optimValues.iteration,x(1),optimValues.iteration,x(2),'o');
        %            plot(optimValues.iteration,x(3),optimValues.iteration,x(4),optimValues.iteration,x(5));
        % Label points with iteration number.
        % Add .15 to x(1) to separate label from plotted 'o'
        %            text(x(1)+.15,x(2),num2str(optimValues.iteration));
        fileID=fopen('opt_history.txt','a');
        fprintf(fileID,['%i,%i,',repmat('%f,',[1,length(optimValues.fval)+length(x)+length(optimValues.stepsize)+length(optimValues.gradient)+length(optimValues.trustregionradius)-1]),'%f\n'],...
            [optimValues.iteration,optimValues.funccount,optimValues.fval,x',optimValues.stepsize,optimValues.gradient',optimValues.trustregionradius]);
        fclose(fileID);
    case 'done'
        figure(1);
        hold off
        xlabel('Iteration Number');
        ylabel('Mutual Information');
        figure(2);
        hold off
        xlabel('Iteration Number');
        yyaxis left
        ylabel('Flip Angle (deg)');
        yyaxis right
        ylabel('Delay Times (s)');
        legend('FA','TD2','TD3','TD4','TD5','Location','NorthEast');
    otherwise
end
end