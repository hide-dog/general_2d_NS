function common_allocation(cellxmax, cellymax, nval)
    
    Qbase = zeros(cellxmax, cellymax, nval)
    volume = zeros(cellxmax, cellymax)
    
    dx = zeros(cellxmax+1, cellymax)
    dy = zeros(cellxmax, cellymax+1)

    Qcon = zeros(cellxmax, cellymax, nval)
    Qcon_hat = zeros(cellxmax, cellymax, nval)

    mu = zeros(cellxmax, cellymax)
    lambda = zeros(cellxmax, cellymax)

    E_adv_hat = zeros(cellxmax+1,   cellymax, nval)
    F_adv_hat = zeros(  cellxmax, cellymax+1, nval)

    E_vis_hat = zeros(cellxmax+1,   cellymax, nval)
    F_vis_hat = zeros(  cellxmax, cellymax+1, nval)

    RHS = zeros(cellxmax, cellymax, nval)

    return Qbase, volume, dx, dy, Qcon, Qcon_hat, mu, lambda, 
            E_adv_hat, F_adv_hat, E_vis_hat, F_vis_hat, RHS
end 

function allocation_implicit(cellxmax, cellymax, nval)
    Qbasen = zeros(cellxmax, cellymax, nval)
    Qconn  = zeros(cellxmax, cellymax, nval)
    Qconn_hat  = zeros(cellxmax, cellymax, nval)
    Qbasem  = zeros(cellxmax, cellymax, nval)

    dtau = zeros(cellxmax, cellymax)
    lambda_facex = zeros(cellxmax+1, cellymax)
    lambda_facey = zeros(cellxmax, cellymax+1)

    A_adv_hat_p = zeros(cellxmax, cellymax, nval, nval)
    A_adv_hat_m = zeros(cellxmax, cellymax, nval, nval)
    B_adv_hat_p = zeros(cellxmax, cellymax, nval, nval)
    B_adv_hat_m = zeros(cellxmax, cellymax, nval, nval)
    A_beta_shig = zeros(cellxmax, cellymax)
    B_beta_shig = zeros(cellxmax, cellymax)

    jalphaP = zeros(cellxmax, cellymax)
    jbetaP  = zeros(cellxmax, cellymax)

    delta_Q = zeros(cellxmax, cellymax, nval)
    delta_Q_temp = zeros(cellxmax, cellymax, nval)

    norm2 = zeros(nval)

    return Qbasen, Qconn, Qconn_hat, Qbasem, dtau, lambda_facex, lambda_facey,
            A_adv_hat_m, A_adv_hat_p, B_adv_hat_m, B_adv_hat_p, A_beta_shig, B_beta_shig,
            jalphaP, jbetaP, delta_Q, delta_Q_temp, norm2
end