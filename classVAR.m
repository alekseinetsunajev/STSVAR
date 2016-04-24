classdef classVAR < data
    properties
       Theta;
       logL;
       minVal;
       B;
       Ksi;
       restrictions_SR
       restrictions_LR
       SR_Matr
       LR_Matr
       IR_SRun
       IR_LRun
       U
        NonLinconstr
%        Bootstrap
%        Bootstrap_LR
    end
    methods
         %% constructor method
        function obj = classVAR(r1, r2, R_Matr_SR, R_Matr_LR)
            obj.restrictions_SR = r1;
            obj.restrictions_LR = r2;
            obj.SR_Matr = R_Matr_SR;
            obj.LR_Matr = R_Matr_LR;
            obj.options = optimset('Algorithm', 'interior-point', 'Display', 'off', 'MaxFunEvals', 20000, ...
                'MaxIter', 20000, 'Hessian','bfgs', 'DerivativeCheck','on','Diagnostics','off','GradObj','off','LargeScale','off' ) ; %, 'UseParallel', 'always');                      

        end
        %%
        function obj = estimateSVAR(obj)
           obj.Theta = obj.Theta0;
           W = eye(obj.T(1,2), obj.T(1,2));
           for w = 1:obj.lags
                 W = W - get_coefficient(obj.Theta0, obj.T, w, obj.trend);
           end  
                
           B0 = chol(1/obj.T(1,1)*(obj.U0'*obj.U0), 'lower');            
           vB = reshape(B0, obj.T(1,2)^2, 1);
              if obj.restrictions_SR == 1 && obj.restrictions_LR == 1
                    ceq = @(vB) [
                        obj.SR_Matr * vB( 1 :   obj.T(1,2)^2, 1)
                        obj.LR_Matr * kron( eye(obj.T(1,2)), W^-1 ) * vB( 1 :   obj.T(1,2)^2, 1)
                        ];
              elseif obj.restrictions_SR == 1 && obj.restrictions_LR ~= 1
                        ceq = @(vB) [
                        obj.SR_Matr * vB( 1 :   obj.T(1,2)^2, 1)
                        ];
              elseif obj.restrictions_SR ~= 1 && obj.restrictions_LR == 1
                        ceq = @(vB) [
                        obj.LR_Matr * kron( eye(obj.T(1,2)), W^-1 ) * vB( 1 :   obj.T(1,2)^2, 1)
                        ];                    
              else
                    ceq = @(vB)[];
              end
                
              if obj.restrictions_SR == 1 || obj.restrictions_LR == 1
                    obj.NonLinconstr = @(Theta)deal( [] ,ceq(Theta));
              else
                    obj.NonLinconstr = @(Theta)deal( [ineq(Theta)], []);
              end
                
            [B, obj.logL] = fmincon( @(vB)LogLikeVAR(vB, obj.T, obj.lags, obj.U0), vB , [], [] , [] , [] , [] , [] , obj.NonLinconstr, obj.options);      
            obj.logL = log(sqrt(1/((2*pi)^obj.T(1,2)))) *( obj.T(1,1) - obj.lags)  - obj.logL;

            obj.B = reshape(B, obj.T(1,2), obj.T(1,2));
            obj.Ksi = W^-1 * obj.B;
            [obj.IR_SRun, obj.IR_LRun] = irs(obj);
            obj.U = residuals2(obj.T, obj.Y, obj.Z, obj.lags, obj.Theta);
        end   
        
        %% bootstrap part
        function obj = mbb_bootstrap(obj)
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
            
               % calculate the parameter vector and var-cov based on bootstrapped
               % series
               Theta_b=(Y_star'*Z_star)/(Z_star'*Z_star);
               U_b = Y_star - (Z_star*Theta_b' ); 
               Sigma_b=U_b'*U_b*1/(obj.T(1,1) - obj.T(1,2));
               % matrix of impact effects
               B0_b= chol(Sigma_b, 'lower');
               obj.Theta =  reshape(Theta_b, obj.lags*(obj.T(1,2)^2) + obj.T(1,2), 1);
               % update the restrctions
               W = eye(obj.T(1,2), obj.T(1,2));
               for w = 1:obj.lags
                     W = W - get_coefficient(obj.Theta, obj.T, w, obj.trend);
               end  
           
               if obj.restrictions_SR == 1 && obj.restrictions_LR == 1
                    ceq = @(vB) [
                        obj.SR_Matr * vB( 1 :   obj.T(1,2)^2, 1)
                        obj.LR_Matr * kron( eye(obj.T(1,2)), W^-1 ) * vB( 1 :   obj.T(1,2)^2, 1)
                        ];
              elseif obj.restrictions_SR == 1 && obj.restrictions_LR ~= 1
                        ceq = @(vB) [
                        obj.SR_Matr * vB( 1 :   obj.T(1,2)^2, 1)
                        ];
              elseif obj.restrictions_SR ~= 1 && obj.restrictions_LR == 1
                        ceq = @(vB) [
                        obj.LR_Matr * kron( eye(obj.T(1,2)), W^-1 ) * vB( 1 :   obj.T(1,2)^2, 1)
                        ];                    
              else
                    ceq = @(vB)[];
              end
                
              if obj.restrictions_SR == 1 || obj.restrictions_LR == 1
                    obj.NonLinconstr = @(vB)deal( [] ,ceq(vB));
              else
                    obj.NonLinconstr = @(vB)deal( [ineq(vB)], []);
              end
               
               vB = fmincon( @(x)LogLikeVAR(x, obj.T, obj.lags, U_b), reshape(B0_b, obj.T(1,2)^2, 1) , [], [] , [] , [] , [] , [] , obj.NonLinconstr, obj.options);
               obj.B = reshape(vB, obj.T(1,2), obj.T(1,2));
               for ii = 1: obj.T(1,2)
                   if obj.B(ii, ii) < 0
                       obj.B(:, ii) = -1 * obj.B(:, ii) ;
                   end
               end
               % produce the IR for the current series
                [obj.IR_SRun, obj.IR_LRun] = irs(obj);
           
        end
%%%%%%%%%%%%%%%%%%%%        
        function obj = bootstr(obj)
            %for n=1:obj.brep
                %disp(n)
                Ustar = obj.U0;
    
                rade = 2*round(rand(obj.T(1,1)-obj.lags, 1))-1;
    
                %  use the Rademacher distribution,
                for t=1:(obj.T(1,1) - obj.lags)
                    Ustar(t, :) = (Ustar(t, :)' * rade(t,:) )';
                end    
    
                % Set initial values:
                 y_star = zeros(obj.T(1,1), obj.T(1,2) );
                 y_star(1: obj.lags, :) = obj.Y(1:obj.lags, :); % Pre-sample values
    
                % Generate Wild bootstrapped series
                for t=obj.lags+1 : obj.T(1,1)
                    if obj.trend == 1
                        for i=obj.T(1,2):obj.T(1,2):obj.T(1,2)*obj.lags
                            y_star(t, :)  = y_star(t, :) + (get_coefficient(obj.Theta0, obj.T, i/T(1,2), obj.trend) * obj.Y(t-i/obj.T(1,2), :)')' ;
                        end
                        y_star(t, :) = y_star(t, :) + obj.Theta0(1:obj.T(1,2))' + t * obj.Theta0(obj.T(1,2)+1:2*obj.T(1,2))'  + Ustar(t-obj.lags, :);

                    else
                        for i=obj.T(1,2):obj.T(1,2) : obj.T(1,2)*obj.lags
                            y_star(t, :)  = y_star(t, :) + (get_coefficient(obj.Theta0, obj.T, i/obj.T(1,2), obj.trend) * obj.Y(t-i/obj.T(1,2), :)')' ;
                        end
                        y_star(t, :) = y_star(t, :) + obj.Theta0(1:obj.T(1,2))' + Ustar(t-obj.lags, :);
                    end
                end

                Y_star = y_star(obj.lags+1:end, :);
                for i = 1:obj.T(1,1)-obj.lags;
                    Z_star(i,1) = 1;            % constant
                    if obj.trend == 1      % with trend
                        Z_star(i,2) = i;
                        for j=1:obj.lags
                            Z_star(i,(j-1)*obj.T(1,2)+2+1 : j*obj.T(1,2)+1+1 ) = obj.Y(obj.lags + i - j , :);
                        end
                    else
                         for j=1:obj.lags
                            Z_star(i,(j-1)*obj.T(1,2)+2 : j*obj.T(1,2)+1 )= obj.Y(obj.lags + i - j , :);  % lagged observations 
                        end           
                    end
                end

               % calculate the parameter vector and var-cov based on bootstrapped
               % series
               Theta_b=(Y_star'*Z_star)/(Z_star'*Z_star);
               U_b = Y_star - (Z_star*Theta_b' ); 
               Sigma_b=U_b'*U_b*1/(obj.T(1,1) - obj.T(1,2));
               % matrix of impact effects
               B0_b= chol(Sigma_b, 'lower');
               obj.Theta =  reshape(Theta_b, obj.lags*(obj.T(1,2)^2) + obj.T(1,2), 1);
               % update the restrctions
               W = eye(obj.T(1,2), obj.T(1,2));
               for w = 1:obj.lags
                     W = W - get_coefficient(obj.Theta, obj.T, w, obj.trend);
               end  
           
               if obj.restrictions_SR == 1 && obj.restrictions_LR == 1
                    ceq = @(vB) [
                        obj.SR_Matr * vB( 1 :   obj.T(1,2)^2, 1)
                        obj.LR_Matr * kron( eye(obj.T(1,2)), W^-1 ) * vB( 1 :   obj.T(1,2)^2, 1)
                        ];
              elseif obj.restrictions_SR == 1 && obj.restrictions_LR ~= 1
                        ceq = @(vB) [
                        obj.SR_Matr * vB( 1 :   obj.T(1,2)^2, 1)
                        ];
              elseif obj.restrictions_SR ~= 1 && obj.restrictions_LR == 1
                        ceq = @(vB) [
                        obj.LR_Matr * kron( eye(obj.T(1,2)), W^-1 ) * vB( 1 :   obj.T(1,2)^2, 1)
                        ];                    
              else
                    ceq = @(vB)[];
              end
                
              if obj.restrictions_SR == 1 || obj.restrictions_LR == 1
                    obj.NonLinconstr = @(vB)deal( [] ,ceq(vB));
              else
                    obj.NonLinconstr = @(vB)deal( [ineq(vB)], []);
              end
               
               vB = fmincon( @(x)LogLikeVAR(x, obj.T, obj.lags, U_b), reshape(B0_b, obj.T(1,2)^2, 1) , [], [] , [] , [] , [] , [] , obj.NonLinconstr, obj.options);
               obj.B = reshape(vB, obj.T(1,2), obj.T(1,2));
               for ii = 1: obj.T(1,2)
                   if obj.B(ii, ii) < 0
                       obj.B(:, ii) = -1 * obj.B(:, ii) ;
                   end
               end
               % produce the IR for the current series
                [obj.IR_SRun, obj.IR_LRun] = irs(obj);
               %[obj.Bootstrap(:, :, n), obj.Bootstrap_LR(:, :, n)] = irs(obj);
            %end
        end

    end
end