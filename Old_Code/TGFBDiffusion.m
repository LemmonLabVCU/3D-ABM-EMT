function [Ctgfb] = TGFBDiffusion(Param, Ctgfb, TGFB, Dtgfb)
%UNTITLED4 Summary of this function goes here
%%
n           = Param.n;
L           = n*Param.Csize;    % length of the mesh (um)
h           = round(n/2);       % center of the spheroid (pixels)
dt          = 1;                % time step (hr)
dx          = round(L/(n-1));   % distance intervals of the grid (um)
dy = dx; dz = dx;               % distance intervals of the grid (um

Ctmask = TGFB*ones(n+2, n+2, n+2); % make a mask of CTGFB
    
    for i = 2:n
        for j = 2:n
            for k = 2:n
                Ctmask(i,j,k) = Ctgfb(i,j,k);
            end
        end
    end
    
    for i = 2:n           % x-axis
        for j = 2:n       % y-axis
            for k = 2:n   % z-axis
                if i== 2 || i == n || j == 2 || j == n || k == 2 || k == n
                    change(i,j,k) = dt*2*Dtgfb(i-1,j-1,k-1)*...
                    ((Ctmask(i+1,j,k) - Ctmask(i,j,k))/dx/dx + ...
                    (Ctmask(i,j+1,k) - Ctmask(i,j,k))/dy/dy + ...
                    (Ctmask(i,j,k+1) - Ctmask(i,j,k))/dz/dz);
                else
                    change(i,j,k) = dt*Dtgfb(i-1,j-1,k-1)*...
                    ((Ctmask(i+1,j,k) - 2*Ctmask(i,j,k) + Ctmask(i-1,j,k))/dx/dx + ...
                    (Ctmask(i,j+1,k) - 2*Ctmask(i,j,k) + Ctmask(i,j-1,k))/dy/dy + ...
                    (Ctmask(i,j,k+1) - 2*Ctmask(i,j,k) + Ctmask(i,j,k-1))/dz/dz);
                end
            end
        end
    end
    
    change(1,:,:) = 0; change(:,1,:) = 0; change(:,:,1) = 0;
    change(n,:,:) = 0; change(:,n,:) = 0; change(:,:,n) = 0;
    
    Ctgfb = Ctgfb + change;

end

