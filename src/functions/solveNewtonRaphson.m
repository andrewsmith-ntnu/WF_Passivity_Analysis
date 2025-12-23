function [xs,dxs] = solveNewtonRaphson(fun,ys,xg,cnv,Ns,bta,litmax)
% Solve a system of nonlinear equations using the Newton-Raphson method with backtracking.

x = xg;

% Get the residual.
[y,dy] = fun (x);
Res = y - ys;
Rval = (Res.')*Res;
    

iter  = 0;
while (((real(Rval) > cnv) && (iter < Ns)) || (iter == 0)) 
   iter = iter + 1;

   % Compute the tangent function at the latest point.
   [y,dy] = fun (x);
   dRdx = dy;

%[slap,shp,ifrq] = eigVal (AA);
%return

   % Apply Newton's method on the equations.
   lam = 1;
   dxr = -dRdx\Res;
   lflg = 0;
   liter = 0;
   while ((lflg == 0) && (liter < litmax))
      liter = liter + 1;

      x1 = x + bta(iter)*lam*dxr;

%[bta(iter) lam]
%[x dxr x1]


      % Compute the residual.
      [y1,dy1] = fun (x1);
      Res1 = y1 - ys;
      R1 = (Res1.')*Res1;

%[dxdt dx1dt]
%printf('%5d %5d  %+5.3e  %+5.3e\n',iter,liter,Rval,R1);

      if (real(R1) < Rval)
         % OK!  Prepare for the next iteration.
         lflg = 1;
         Res = Res1;
         Rval = R1;
         x = x1;
      else
         % Backtrack.
         lam = 0.5*lam;
         if (liter == litmax)
            %[iter R1 Rval]
            fprintf('Warning, proceeding without lambda convergence.\n');
%              msgbox('Warning, proceeding without lambda convergence.\n');
            lflg = 1;
            Res = Res1;
            Rval = R1;
            x = x1;
            % added sda            
            xs = 0;
            dxs = 0;

return
         end
      end

   end  % Newton inner.

   if (iter == Ns) && (real(Rval) > cnv)
      %Rval
%       printf('Warning, max iterations, proceeding without Rval convergence.\n');
%       msgbox('Warning, max iterations, proceeding without Rval convergence.\n');
    
     % added sda            
     xs = 0;
     dxs = 0;
return
   end

end  % Newton outer.

xs = x;
dxs = pinv(dy1);