function [IncHxEz, IncHyEz, IncEzHx, IncEzHy] = get_1D_update_coefficients(dt, dy, dx, epsEnv, muEnv)
    %% calculate update coefficients for 1D update ( PLANE WAVE )
    IncHxEz = (dt / dy) / muEnv;
    IncHyEz = (dt / dx) / muEnv;
    IncEzHx = (dt / dy) / epsEnv;
    IncEzHy = (dt / dx) / epsEnv;
end