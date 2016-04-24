classdef classSTSVAR < data
    properties
       gamma;
       c;
    end
    properties (SetAccess = private)
       LM = 1.501;
       logL;
       minVal;
       Theta;
       BL;
       B;
       Lambda;
       restrictions_SR
       restrictions_LR
       SR_Matr
       LR_Matr
       StErrors
       B_err
       L_err
       IR_SRun
       IR_LRun
       FEVD1
       FEVD2
       U
       NonLinconstr
       TR_var
       TR_f
       WTest
       restrictions_L
       l_i
       l_j
    end
    methods
        %% constructor method
        function obj = classSTSVAR(gamma, c, r1, r2, r3, R_Matr_SR, R_Matr_LR, l_i, l_j, TR_var)
             if size(TR_var,1) == obj.T(1,1)-obj.lags
                obj.TR_var = TR_var;
             else
                 disp('Error in the length of the transition variable. Cannot proceed.');
                 return;
             end
            obj.gamma = gamma;
            obj.c = c;
            
            for i = 1 :obj.T(1,1)-obj.lags
                 obj.TR_f(i,1) = calc_transition(TR_var(i,1), obj.c, obj.gamma, obj.exp);
            end

            obj.restrictions_SR = r1;
            obj.restrictions_LR = r2;
            if r3 == 1
                obj.restrictions_L = r3;
                P = reshape(1:1:obj.T(1,2)^2, obj.T(1,2), obj.T(1,2));
                obj.l_i = l_i;                      % i, i-i, i-j
                obj.l_i(1,2) = P(l_i, l_i);
                obj.l_i(1,3) = P(l_i, l_j);
                
                obj.l_j = l_j;                      % j, j-i, j-j
                obj.l_j(1,2) = P(l_j, l_i);
                obj.l_j(1,3) = P(l_j, l_j);
            else
                obj.restrictions_L = 0;
            end
            
            obj.SR_Matr = R_Matr_SR;
            obj.LR_Matr = R_Matr_LR;
            %obj.L_R = L_R;
            if obj.restrictions_SR == 1 || obj.restrictions_LR == 1 || obj.restrictions_L == 1
                         obj.options = optimset('Algorithm', 'interior-point', 'Display', 'off', 'MaxFunEvals', 20000, ...
                         'MaxIter', 20000, 'Hessian','bfgs', 'DerivativeCheck','on','Diagnostics','off','GradObj','off','LargeScale','off', 'UseParallel', 'always');
            else                         
                        obj.options = optimset('Algorithm', 'active-set', 'Display', 'off', 'MaxFunEvals', 20000, ...
                        'MaxIter', 20000, 'Hessian','bfgs', 'DerivativeCheck','on','Diagnostics','off','GradObj','off','LargeScale','off', 'UseParallel', 'always');
            end
            
            obj.TR_var = TR_var;
            obj.Lambda0 =  obj.LM *eye(obj.T(1,2)) + diag((0.1:0.1:0.5)) ;
        end

        %% optimization method
        function obj = optimize(obj)
            disp('Estimating...')
            Param0 = [reshape(obj.B0, obj.T(1,2)^2, 1); diag(obj.Lambda0)];
            Theta0 = obj.Theta0;
            lb = lbound(obj.T, obj.s);
            lb = reshape(lb, obj.T(1,2)^2 + obj.T(1,2), 1);
            
            select_matr = zeros(obj.T(1,2)-1, obj.T(1,2));
            
            for k = 1 : (obj.T(1,2)-1)
                select_matr(k, k:k+1) = [1 -1];
            end
 
            val(1,1)=1;
            val(1,2)=1;
            RND = 1;
            
            while val(RND,2) > 0.000000001
                W = eye(obj.T(1,2), obj.T(1,2));
                for w = 1:obj.lags
                    W = W - get_coefficient(obj.Theta0, obj.T, w, obj.trend);
                end  
                
                ineq = @(Theta)[
                select_matr * Theta(obj.T(1,2)^2+1:end, 1)  ;
                ];
            
                if obj.restrictions_SR == 1 && obj.restrictions_LR == 1
                    ceq = @(Theta) [
                        obj.SR_Matr * Theta( 1 :   obj.T(1,2)^2, 1)
                        obj.LR_Matr * kron( eye(obj.T(1,2)), W^-1 ) * Theta( 1 :   obj.T(1,2)^2, 1)
                        ];
                elseif obj.restrictions_SR == 1 && obj.restrictions_LR ~= 1
                        ceq = @(Theta) [
                        obj.SR_Matr * Theta( 1 :   obj.T(1,2)^2, 1)
                        ];
                elseif obj.restrictions_SR ~= 1 && obj.restrictions_LR == 1
                        ceq = @(Theta) [
                        obj.LR_Matr * kron( eye(obj.T(1,2)), W^-1 ) * Theta( 1 :   obj.T(1,2)^2, 1)
                        ];                    

                elseif obj.restrictions_L == 1
                     ceq = [];
                else
                    ceq = @(Theta)[];
                end
            
                if obj.restrictions_SR == 1 || obj.restrictions_LR == 1 
                    obj.NonLinconstr = @(Theta)deal( [] ,ceq(Theta));
                elseif obj.restrictions_L == 1
                    obj.NonLinconstr = @(Theta)deal( [] ,ceq(Theta));
                else
                    %obj.NonLinconstr = @(Theta)deal( [ineq(Theta)], []);
                end
                
                % estimated B and L
                [Param(:, RND), minV] = fmincon( @(Theta) minimize(obj.Y, obj.Z, obj.T, obj.lags, Theta, Theta0, obj.TR_f), Param0 , [], [] , [] , [] , lb , [] , obj.NonLinconstr, obj.options);

                B = reshape(Param(1:obj.T(1,2)^2, RND), obj.T(1,2), obj.T(1,2));
                L = diag(Param(obj.T(1,2)^2+1 : obj.T(1,2)^2 +  obj.T(1,2), RND));
                disp(L);
                Theta0 = GLS(obj.Z, obj.Y, B, L, obj.T, obj.lags, obj.TR_f);
                    
                Param0 = Param(:, RND);
                val(RND+1,1) = minV;
                val(RND+1,2) = abs( ( val(RND+1,1) - val(RND,1) ) / val(RND,1) );
                RND = RND + 1;        
            end
            
            % vectorized structural parameters and VAR parameters
            obj.BL = Param0;
            obj.Theta = Theta0;
            obj.minVal = minV;    
            
            % sorted structural parameters
            if obj.restrictions_SR == 1 || obj.restrictions_LR == 1 || obj.restrictions_L == 1
                obj.B = reshape(Param0(1: obj.T(1,2)^2, 1), obj.T(1,2), obj.T(1,2));
                obj.Lambda = diag(Param0( obj.T(1,2)^2 + 1:  obj.T(1,2)^2 +  obj.T(1,2), 1));
                % normalize signs of B if needed
            else
                [obj.B, obj.Lambda] = sort_BL(Param0, obj.T);
            end
            
            % calculate the full likelihood value for each t.
            obj.logL = log(sqrt(1/((2*pi)^obj.T(1,2)))) *( obj.T(1,1) - obj.lags) - minV ;
            obj.U = residuals2(obj.T, obj.Y, obj.Z, obj.lags, obj.Theta);
            
        end
        %% method to calculate strandard errors
        function obj = std_err(obj)
                param = [reshape(obj.B, obj.T(1,2)^2, 1); diag(obj.Lambda); obj.Theta]; 
                myDelta=1e-4*abs(param)+0.00000001;

                % outer product calculation as in Hamilton page 143
                s = zeros(obj.T(1,1)-obj.lags,numel(param));
                u = residuals2(obj.T, obj.Y, obj.Z, obj.lags, obj.Theta);

                % calculate the initial likelihood value for each t.
                [~, logLikVec] = calculate_L(obj);

                for i=1:numel(param)
                    m=param;
                    m(i)=param(i)+myDelta(i);

                    B_s = reshape(m(1:obj.T(1,2)^2, 1), obj.T(1,2), obj.T(1,2)) ;
                    Lambda_s = diag(m(obj.T(1,2)^2+1:obj.T(1,2)^2 + obj.T(1,2), 1));
                    Theta_s = m(obj.T(1,2)^2 + obj.T(1,2) + 1 : end);
                    u = residuals2(obj.T, obj.Y, obj.Z, obj.lags, Theta_s );

                    % update the likelihood values
                    for t = 1: obj.T(1,1)-obj.lags    
                         Exponent = u(t,:)  * (((1 - obj.TR_f(t,1))*(B_s*B_s') + obj.TR_f(t,1)*B_s*Lambda_s*B_s')^-1 )* u(t,:)';   
                         Det = (det((1 - obj.TR_f(t,1))*(B_s*B_s') + obj.TR_f(t,1) * B_s*Lambda_s*B_s')) ;
                         logLikVec_update(t,i) = log( sqrt(1/(2*pi)^obj.T(1,2))*1/sqrt(Det) * exp(-0.5 * Exponent));
                    end

                    % calculate the approximated derivative for each t and each parameter
                    for t=1:obj.T(1,1) - obj.lags
                        s(t,i)=(logLikVec_update(t,i)-logLikVec(t,1))/(myDelta(i));
                    end

                end

                % form the outer product matrix as in Hamilton  eq 5.8.4
                sum_s_Matrix=zeros(numel(param));
                for idx=1:obj.T(1,1)-obj.lags
                    s_Matrix=s(idx,:)'*s(idx,:);
                    sum_s_Matrix=sum_s_Matrix+s_Matrix;
                end
                OP_Matrix=sum_s_Matrix/obj.T(1,1);
                StErr_all(:, 1) = sqrt(diag(inv(OP_Matrix*obj.T(1,1))));
                obj.StErrors = StErr_all(:, 1);
                obj.B_err = reshape(StErr_all(1:obj.T(1,2)^2, 1), obj.T(1,2), obj.T(1,2));
                obj.L_err = diag(StErr_all(obj.T(1,2)^2+1 : obj.T(1,2)^2 + obj.T(1,2), 1));
                
                 % Wald test for the parameters
                Hess = (OP_Matrix*obj.T(1,1))^-1;
                n = obj.T(1,2)*(obj.T(1,2)-1) - sum(1:(obj.T(1,2)-1)); % number of pairs to be checked
                Test = zeros(n, numel(param));                              
       
                Rest = zeros(n, obj.T(1,2));
                k = 1;
                j = k + 1; 
                for i = 1:n
                    if j <= obj.T(1,2)
                        Rest(i, k) = 1;
                        Rest(i, j) = -1;
                        j = j+1;
                    else
                        k = k + 1;
                        j = k + 1;
                        Rest(i, k) = 1;
                        Rest(i, j) = -1;
                        j = j+1;
                    end
                end
                
                for ii = 1:n
                    Test(ii, obj.T(1,2)^2 + 1: obj.T(1,2)^2 + obj.T(1,2)) = Rest(ii, :);
                    WTest(ii,1)=((Test(ii,:)*param * (Test(ii,:)* Hess * Test(ii,:)')^-1 * Test(ii,:)*param ));
                    WTest(ii,2) = 1-chi2cdf(WTest(ii,1),1);
                end
                obj.WTest = WTest;

        end
        
        %% Impulse responses
        function obj = IR(obj)
            [obj.IR_SRun, obj.IR_LRun] = irs(obj);
        end
        
        %% variance decomposition
        function obj = FVD(obj)
            obj.FEVD1 = FEVD(obj, 1);
            obj.FEVD2 = FEVD(obj, 2);
        end
        
        %% bootstrap
        function obj = mbb_bootstrap(obj)
            B00 = obj.B;
            Lambda00 = obj.Lambda;
            Param00 = [reshape(B00, obj.T(1,2)^2, 1) ; diag(Lambda00)];
            Theta00 = obj.Theta;
            % moving block bootstrap
            
            block_l = 30;                                                    % length of block
            block_nr = round( (obj.T(1,1) - obj.lags) /block_l);
            while (block_l * block_nr < (obj.T(1,1) - obj.lags) )
                block_l = block_l + 1;
                block_nr = round((obj.T(1,1) - obj.lags) /block_l);             
            end
            
            Uhatstar = zeros(1, obj.T(1,2));
            for i0 =    0 : block_nr - 1    % form block_nr blocks
                
                block_id = ( unidrnd(obj.T(1,1) - obj.lags - block_l ) - 1 );          % get random id of the block
                B_il = obj.U(block_id + 1 : block_id  + block_l, :);     % get the corresponding residuals
                Uhatstar = [Uhatstar; B_il];
            end
            Uhatstar = Uhatstar(2 : end , :);
            
            % step 3 centering
            for s = 1: block_l
                for j = 0 : block_nr - 1
                    center = 0;
                    for r = 0 : obj.T(1,1) - obj.lags - block_l
                        center = center + Uhatstar(s+r, :);
                    end
                    Ustar(j*block_l + s, :) = Uhatstar(j*block_l + s,  :) - center / (obj.T(1,1) - obj.lags - block_l + 1);
                end
            end
            Ustar = Ustar(1:obj.T(1,1) - obj.lags, :);
            
            % step 4 
            % generate the observations
            y_star = zeros(obj.T(1,1), obj.T(1,2) );
            for t=obj.lags+1 : obj.T(1,1)
                if obj.trend == 1
                    for i=obj.T(1,2):obj.T(1,2):obj.T(1,2)*obj.lags
                        y_star(t, :)  = y_star(t, :) + (get_coefficient(obj.Theta, obj.T, i/T(1,2), obj.trend) * obj.Y(t-i/obj.T(1,2), :)')' ;
                    end
                    y_star(t, :) = y_star(t, :) + obj.Theta(1:obj.T(1,2))' + t * obj.Theta(obj.T(1,2)+1:2*obj.T(1,2))'  + Ustar(t-obj.lags, :);

                else
                    for i=obj.T(1,2):obj.T(1,2) : obj.T(1,2)*obj.lags
                        y_star(t, :)  = y_star(t, :) + (get_coefficient(obj.Theta, obj.T, i/obj.T(1,2), obj.trend) * y_star(t-i/obj.T(1,2), :)')' ;
                    end
                    y_star(t, :) = y_star(t, :) + obj.Theta(1:obj.T(1,2))' + Ustar(t-obj.lags, :);
                end
            end
            
            %end of generating bootsrap data
            Y_star = y_star(obj.lags+1:end, :);
            for i = 1:obj.T(1,1)-obj.lags;
                Z_star(i,1) = 1;            % constant
                if obj.trend == 1      % with trend
                    Z_star(i,2) = i;
                    for j=1:obj.lags
                        Z_star(i,(j-1)*obj.T(1,2)+2+1 : j*obj.T(1,2)+1+1 )= y_star(obj.lags + i - j , :);  % lagged observations 
                    end
                else
                     for j=1:obj.lags
                        Z_star(i,(j-1)*obj.T(1,2)+2 : j*obj.T(1,2)+1 )= y_star(obj.lags + i - j , :);  % lagged observations 
                    end           
                end
            end
          
            val(1,1) =1;
            val(1,2) =1;
            RND     = 1;            
            while val(RND,2) > 0.000000001  
                lb = lbound(obj.T, obj.s);
                lb = reshape(lb, obj.T(1,2)^2 + obj.T(1,2), 1); 
                [Param(:, RND), minV] = fmincon( @(Theta) minimize(y_star, Z_star, obj.T, obj.lags, Theta, Theta00, obj.TR_f), Param00 , [], [] , [] , [] , lb , [] , obj.NonLinconstr, obj.options);

                B = reshape(Param(1:obj.T(1,2)^2, RND), obj.T(1,2), obj.T(1,2));
                L = diag(Param(obj.T(1,2)^2+1 : obj.T(1,2)^2 +  obj.T(1,2), RND));   
                Theta00 = GLS(Z_star, y_star, B, L, obj.T, obj.lags, obj.TR_f);
                Param00 = Param(:, RND);
                val(RND+1,1) = minV;
                val(RND+1,2) = abs( ( val(RND+1,1) - val(RND,1) ) / val(RND,1) );
                RND = RND + 1;                
            end
            
            % vectorized structural parameters and VAR parameters
            obj.BL = Param00;
            obj.Theta = Theta00;
            obj.minVal = minV;    
            
            % sorted structural parameters
            if obj.restrictions_SR == 1 || obj.restrictions_LR == 1
                obj.B = reshape(Param00(1: obj.T(1,2)^2, 1), obj.T(1,2), obj.T(1,2));
                obj.Lambda = diag(Param00( obj.T(1,2)^2 + 1:  obj.T(1,2)^2 +  obj.T(1,2), 1));
                % normalize signs of B if needed
            else
                [obj.B, obj.Lambda] = sort_BL(Param00, obj.T);
            end
            
            [obj.IR_SRun, obj.IR_LRun] = irs(obj);            
        end
%%%%%%%%%%%%%%%%%%%%%%%        
        function obj = bootstr(obj)
            
            B00 = obj.B;
            Lambda00 = obj.Lambda;
            Param00 = [reshape(B00, obj.T(1,2)^2, 1) ; diag(Lambda00)];
            Theta00 = obj.Theta;
            
             y_star = zeros(obj.T(1,1), obj.T(1,2) );
             y_star(1: obj.lags, :) = obj.Y(1:obj.lags, :); % Pre-sample values

            % Fixed design wild bootstrap series 
            % Produce Rademacher distributed random numbers:
            rade = 2*round(rand(obj.T(1,1)-obj.lags, 1))-1;
    
            %  use the Rademacher distribution
            for t=1:(obj.T(1,1)-obj.lags)
                Ustar(t, :) = (obj.U(t, :)' * rade(t,:) )';
            end

            % Generate wild bootstraped series,
            for t=obj.lags+1 : obj.T(1,1)
                if obj.trend == 1
                    for i=obj.T(1,2):obj.T(1,2):obj.T(1,2)*obj.lags
                        y_star(t, :)  = y_star(t, :) + (get_coefficient(obj.Theta, obj.T, i/T(1,2), obj.trend) * obj.Y(t-i/obj.T(1,2), :)')' ;
                    end
                    y_star(t, :) = y_star(t, :) + obj.Theta(1:obj.T(1,2))' + t * obj.Theta(obj.T(1,2)+1:2*obj.T(1,2))'  + Ustar(t-obj.lags, :);

                else
                    for i=obj.T(1,2):obj.T(1,2) : obj.T(1,2)*obj.lags
                        y_star(t, :)  = y_star(t, :) + (get_coefficient(obj.Theta, obj.T, i/obj.T(1,2), obj.trend) * obj.Y(t-i/obj.T(1,2), :)')' ;
                    end
                    y_star(t, :) = y_star(t, :) + obj.Theta(1:obj.T(1,2))' + Ustar(t-obj.lags, :);
                end
            end
            
            Y_star = y_star(obj.lags+1:end, :);
            for i = 1:obj.T(1,1)-obj.lags;
                Z_star(i,1) = 1;            % constant
                if obj.trend == 1      % with trend
                    Z_star(i,2) = i;
                    for j=1:obj.lags
                        Z_star(i,(j-1)*obj.T(1,2)+2+1 : j*obj.T(1,2)+1+1 )= obj.Y(obj.lags + i - j , :);  % lagged observations 
                    end
                else
                     for j=1:obj.lags
                        Z_star(i,(j-1)*obj.T(1,2)+2 : j*obj.T(1,2)+1 )= obj.Y(obj.lags + i - j , :);  % lagged observations 
                    end           
                end
            end
    
            val(1,1) =1;
            val(1,2) =1;
            RND     = 1;            
            while val(RND,2) > 0.000000001  
                lb = lbound(obj.T, obj.s);
                lb = reshape(lb, obj.T(1,2)^2 + obj.T(1,2), 1); 
                [Param(:, RND), minV] = fmincon( @(Theta) minimize(y_star, Z_star, obj.T, obj.lags, Theta, Theta00, obj.TR_f), Param00 , [], [] , [] , [] , lb , [] , obj.NonLinconstr, obj.options);

                B = reshape(Param(1:obj.T(1,2)^2, RND), obj.T(1,2), obj.T(1,2));
                L = diag(Param(obj.T(1,2)^2+1 : obj.T(1,2)^2 +  obj.T(1,2), RND));   
                Theta00 = GLS(Z_star, y_star, B, L, obj.T, obj.lags, obj.TR_f);
                Param00 = Param(:, RND);
                val(RND+1,1) = minV;
                val(RND+1,2) = abs( ( val(RND+1,1) - val(RND,1) ) / val(RND,1) );
                RND = RND + 1;                
            end
            
            % vectorized structural parameters and VAR parameters
            obj.BL = Param00;
            obj.Theta = Theta00;
            obj.minVal = minV;    
            
            % sorted structural parameters
            if obj.restrictions_SR == 1 || obj.restrictions_LR == 1
                obj.B = reshape(Param00(1: obj.T(1,2)^2, 1), obj.T(1,2), obj.T(1,2));
                obj.Lambda = diag(Param00( obj.T(1,2)^2 + 1:  obj.T(1,2)^2 +  obj.T(1,2), 1));
                % normalize signs of B if needed
            else
                [obj.B, obj.Lambda] = sort_BL(Param00, obj.T);
            end
            
            [obj.IR_SRun, obj.IR_LRun] = irs(obj);
            
        end
        

        %% plot transition
        function plot_tr(obj)
            figure('Name','Optimal transition path','NumberTitle','off')
            plot(obj.TR_f,'k',  'LineWidth', 1.5);
            xlim(gca, [0 obj.T(1,1)-obj.lags]);
            ylim(gca, [0 1.2]);
        end
        %% plot  residuals
        function plot_u(obj)
            periods = 1 : obj.T(1,1) - obj.lags;
            figure('Name','Standardized VAR residuals','NumberTitle','off')
            Cov_M = 1/(obj.T(1,1) - obj.lags*obj.T(1,2) - 1)*obj.U0' * obj.U0;
            for t=1:obj.T(1,1)-obj.lags
              U_VAR(t, :) = (diag(Cov_M).^(-0.5))' .* obj.U0(t, :) ;
            end
            Title_txt = {'u_t^q', 'u_t^{\pi}', 'u_t^c', 'u_t^{\Deltas}', 'u_t^r' };
            for i = 1: obj.T(1,2)
                subplot(1, obj.T(1,2) , i);
                plot(periods, U_VAR(:,i),'k',  'LineWidth', 1);
                xlim(gca, [0 obj.T(1,1)]);
                set(gca, 'XTick', [0; 441])
                set(gca,'XTickLabel', {'1970' , '2007'});
                title(Title_txt(i));
            end
            %Title_txt0 = {'\epsilon_t^1', '\epsilon_t^2', '\epsilon_t^3', '\epsilon_t^4', '\epsilon_t^5' };
            for t=1:obj.T(1,1)-obj.lags
               Omega = ((1- obj.TR_f(t,1))*(obj.B*obj.B') +obj.TR_f(t,1)*(obj.B*obj.Lambda*obj.B'));
                U_st(t, :) = ( diag(Omega).^(-0.5) .* obj.U(t, :)' )' ;
%                 Omega = ((1- calc_transition(t, obj.c, obj.gamma, obj.exp))*eye(obj.T(1,2)) + calc_transition(t, obj.c, obj.gamma, obj.exp) *obj.Lambda);
%                 
%                 U_st(t, :) = (diag(Omega).^(-0.5))' .* (obj.B^(-1) *obj.U(t, :)' )';
            end
            figure('Name','Standardized ST-SVAR residuals','NumberTitle','off')
            for i = 1: obj.T(1,2)
                subplot(1, obj.T(1,2), i);
                plot(periods, U_st(:,i),'k',  'LineWidth', 1);
                xlim(gca, [0 obj.T(1,1)]);
                set(gca, 'XTick', [0; 441])
                set(gca,'XTickLabel', {'1970' , '2007'});
                title(Title_txt(i));
            end      
        end
        
    end
end