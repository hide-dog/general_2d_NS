using ProgressMeter
using Dates

function main()
    
    start_t = now()
    out_dir  = "result"
    PARAMDAT = "PARAMDAT.json"
    fwrite   = "write"
    
    nval = 4
    Rd   = 287.0
    R    = 8.314

    xmax, ymax, nodes, vecAx, vecAy = read_allgrid()
    out_file_front, out_ext, restartnum, restart_file, init_small, norm_ok,
    time_integ, nt, dt, every_outnum, in_nt, dtau, cfl,
    init_rho, init_u, init_v, init_p, init_T, specific_heat_ratio, Rd, bdcon = input_para(PARAMDAT)
    
    cellxmax = xmax - 1
    cellymax = ymax - 1

    common_allocation(cellxmax, cellymax, nval)

    Qbase, restartnum = set_initQbase(Qbase, cellxmax, cellymax, restart_file, init_rho, init_u, init_v, init_p, init_T,
                                      specific_heat_ratio, out_file_front, out_ext, out_dir, restartnum, Rd, nval)
    
    # init Delta_Qcon_hat
    volume = set_volume(nodes, cellxmax, cellymax, volume)
    dx, dy = set_dx_lts(dx, dy, nodes, cellxmax, cellymax)
    reset_write(fwrite)

    # main loop
    print("threads num : ")
    println(Threads.nthreads())

    # check bd
    check_bd(bdcon)
    
    #throw(UndefVarError(:x))
    # main loop
    prog = Progress(nt,1)
    @time for t in 1:nt
        next!(prog)
                
        evalnum = t+restartnum
        if time_integ == "1"
            # exlicit scheme
            
            Qbase    = set_boundary(Qbase, cellxmax, cellymax, vecAx, vecAy, bdcon, Rd, specific_heat_ratio, nval)
            Qcon     = base_to_conservative(Qbase, Qcon, cellxmax, cellymax, specific_heat_ratio)
            Qcon_hat = setup_Qcon_hat(Qcon, Qcon_hat, cellxmax, cellymax, volume, nval)
            
            # initial_setup
            mu     = set_mu(Qbase, cellxmax, cellymax, specific_heat_ratio, Rd)
            lambda = set_lambda(Qbase, cellxmax, cellymax, mu, specific_heat_ratio, Rd)
            
            # RHS
            # advection_term
            #E_adv_hat, F_adv_hat = AUSM(Qbase,Qcon,cellxmax,cellymax,vecAx,vecAy,specific_heat_ratio)
            E_adv_hat, F_adv_hat = AUSM_plus(Qbase, Qcon, cellxmax, cellymax, vecAx, vecAy, specific_heat_ratio, volume, nval)
                        
            # viscos_term
            E_vis_hat, F_vis_hat = central_diff(Qbase, Qcon, cellxmax, cellymax, mu, lambda,
                                                vecAx, vecAy, specific_heat_ratio, volume, Rd, nval)
            
            RHS = setup_RHS(cellxmax, cellymax, E_adv_hat, F_adv_hat, E_vis_hat, F_vis_hat, nval, volume)
            
            # time integral
            Qcon_hat = time_integration_explicit(dt, Qcon_hat, RHS, cellxmax, cellymax, nval)

            Qcon  = Qhat_to_Q(Qcon, Qcon_hat, cellxmax, cellymax, volume, nval)
            Qbase = conservative_to_base(Qbase, Qcon, cellxmax, cellymax, specific_heat_ratio)

        elseif time_integ == "2"
            output_physicaltime(fwrite, t, dt)
            # Qcn
            Qbasen    = copy(Qbase)
            Qconn     = base_to_conservative(Qbasen, Qconn, cellxmax, cellymax, specific_heat_ratio)
            Qconn_hat = setup_Qcon_hat(Qconn, Qconn_hat, cellxmax, cellymax, volume, nval)

            Qbasem = copy(Qbase)

            for tau in 1:in_nt
                # LHS (A_adv_hat=jacobian)
                Qbasem = set_boundary(Qbasem, cellxmax, cellymax, vecAx, vecAy, bdcon, Rd, specific_heat_ratio, nval)
                
                Qcon     = base_to_conservative(Qbasem, Qcon, cellxmax, cellymax, specific_heat_ratio)
                Qcon_hat = setup_Qcon_hat(Qcon, Qcon_hat, cellxmax, cellymax, volume, nval)

                # initial_setup
                mu     = set_mu(Qbasem, cellxmax, cellymax, specific_heat_ratio, Rd)
                lambda = set_lambda(Qbasem, cellxmax, cellymax, mu, specific_heat_ratio, Rd)
                dtau   = set_lts(Qbasem, cellxmax, cellymax, mu, dx, dy, vecAx, vecAy, volume, specific_heat_ratio, cfl)
                                
                # RHS
                #advection_term
                E_adv_hat, F_adv_hat = AUSM_plus(Qbasem, Qcon, cellxmax, cellymax, vecAx, vecAy, specific_heat_ratio, volume, nval)

                # viscos_term
                E_vis_hat, F_vis_hat = central_diff(Qbasem, Qcon, cellxmax, cellymax, mu, lambda,
                                                vecAx, vecAy, specific_heat_ratio, volume, Rd, nval)
                
                RHS = setup_RHS(cellxmax, cellymax, E_adv_hat, F_adv_hat, E_vis_hat, F_vis_hat, nval, volume)
            
                # lusgs_advection_term
                A_adv_hat_p, A_adv_hat_m, B_adv_hat_p, B_adv_hat_m, A_beta_shig, B_beta_shig = one_wave(Qbasem, Qcon, cellxmax, cellymax, vecAx, vecAy, specific_heat_ratio, volume, nval)
                # lusgs_viscos_term
                jalphaP, jbetaP = central_diff_jacobian(Qbasem, Qcon, cellxmax, cellymax, mu, lambda, vecAx, vecAy, specific_heat_ratio, volume, nval)
                
                # LUSGS
                delta_Q = zeros(cellxmax, cellymax, nval)

                ite = 0
                norm2 = zeros(nval)
                while true
                    delta_Q_temp = copy(delta_Q)
                    
                    delta_Q = lusgs(dt, dtau, Qcon_hat, Qconn_hat, delta_Q, A_adv_hat_p,  A_adv_hat_m,  B_adv_hat_p,  B_adv_hat_m,  A_beta_shig,  B_beta_shig, jalphaP,  jbetaP, RHS, cellxmax, cellymax, volume, nval)
                    
                    res   = set_res(delta_Q, delta_Q_temp, cellxmax, cellymax, nval)
                    norm2 = check_converge(res, RHS, cellxmax, cellymax, init_small, nval)

                    ite += 1
                    if norm2[1] < norm_ok && norm2[2] < norm_ok && norm2[3] < norm_ok && norm2[4] < norm_ok
                        break
                    end
                    if ite % 100 ==0
                        println(" now cal norm2 ")
                        println(norm2)
                    end
                end
                output_innertime(fwrite, tau, norm2, nval)
                
                for l in 1:nval
                    for j in 2:cellymax-1
                        for i in 2:cellxmax-1
                            Qcon_hat[i,j,l] = Qcon_hat[i,j,l] + delta_Q[i,j,l]
                        end
                    end
                end

                Qcon = Qhat_to_Q(Qcon, Qcon_hat, cellxmax, cellymax, volume, nval)
                Qbasem = conservative_to_base(Qbasem, Qcon, cellxmax, cellymax, specific_heat_ratio)

            end            
            Qbase = copy(Qbasem)
        end
        
        if round(evalnum) % every_outnum == 0
            println("\n")
            println("nt_______________________________"*string(round(evalnum)))
            output_result(evalnum, Qbase, cellxmax, cellymax, specific_heat_ratio, out_file_front, out_ext, out_dir, Rd, nval)
        end

        check_divrege(Qbase, cellxmax, cellymax, Rd, fwrite)
    end
    
    end_t = now()
    output_fin(fwrite, start_t, end_t, nt, dt, in_nt, cellxmax, cellymax)
end

# -- main --
main()


