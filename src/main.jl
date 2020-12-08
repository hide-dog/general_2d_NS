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

    Qbase, volume, dx, dy, Qcon, Qcon_hat, mu, lambda, 
    E_adv_hat, F_adv_hat, E_vis_hat, F_vis_hat, RHS    = common_allocation(cellxmax, cellymax, nval)

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
    if time_integ == "1"
        # exlicit scheme
        prog = Progress(nt,1)
        @time for t in 1:nt
            next!(prog)
            
            evalnum = t + restartnum
        
            Qbase    = set_boundary(Qbase, cellxmax, cellymax, vecAx, vecAy, bdcon, Rd, specific_heat_ratio, nval)
            Qcon     = base_to_conservative(Qbase, Qcon, cellxmax, cellymax, specific_heat_ratio)
            Qcon_hat = setup_Qcon_hat(Qcon, Qcon_hat, cellxmax, cellymax, volume, nval)
            
            # initial_setup
            mu     = set_mu(mu, Qbase, cellxmax, cellymax, specific_heat_ratio, Rd)
            lambda = set_lambda(lambda, Qbase, cellxmax, cellymax, mu, specific_heat_ratio, Rd)
            
            # RHS
            # advection_term
            E_adv_hat, F_adv_hat = AUSM_plus(E_adv_hat, F_adv_hat, Qbase, Qcon, cellxmax, cellymax, 
                                            vecAx, vecAy, specific_heat_ratio, volume, nval)
                        
            # viscos_term
            E_vis_hat, F_vis_hat = central_diff(E_vis_hat, F_vis_hat, Qbase, Qcon, cellxmax, cellymax, mu, lambda,
                                                vecAx, vecAy, specific_heat_ratio, volume, Rd, nval)
            
            RHS = setup_RHS(RHS, cellxmax, cellymax, E_adv_hat, F_adv_hat, E_vis_hat, F_vis_hat, nval, volume)
            
            # time integral
            Qcon_hat = time_integration_explicit(dt, Qcon_hat, RHS, cellxmax, cellymax, nval)

            Qcon  = Qhat_to_Q(Qcon, Qcon_hat, cellxmax, cellymax, volume, nval)
            Qbase = conservative_to_base(Qbase, Qcon, cellxmax, cellymax, specific_heat_ratio)

            if round(evalnum) % every_outnum == 0
                println("\n")
                println("nt_______________________________"*string(round(evalnum)))
                output_result(evalnum, Qbase, cellxmax, cellymax, specific_heat_ratio, out_file_front, out_ext, out_dir, Rd, nval)
            end
    
            check_divrege(Qbase, cellxmax, cellymax, Rd, fwrite)
        end
    elseif time_integ == "2"
        Qbasen, Qconn, Qconn_hat, Qbasem, dtau, lambda_facex, lambda_facey,
        A_adv_hat_m, A_adv_hat_p, B_adv_hat_m, B_adv_hat_p, A_beta_shig, B_beta_shig,
        jalphaP, jbetaP, delta_Q, delta_Q_temp, norm2 = allocation_implicit(cellxmax, cellymax, nval)

        prog = Progress(nt,1)
        @time for t in 1:nt
            next!(prog)
            
            evalnum = t + restartnum
        
            output_physicaltime(fwrite, t, dt)
            # Qcn
            for l in 1:nval
                for j in 1:cellymax
                    for i in 1:cellxmax
                        Qbasen[i,j,l] = Qbase[i,j,l]
                        Qbasem[i,j,l] = Qbase[i,j,l]
                    end
                end
            end
            Qconn     = base_to_conservative(Qbasen, Qconn, cellxmax, cellymax, specific_heat_ratio)
            Qconn_hat = setup_Qcon_hat(Qconn, Qconn_hat, cellxmax, cellymax, volume, nval)     

            for tau in 1:in_nt
                # LHS (A_adv_hat=jacobian)
                Qbasem = set_boundary(Qbasem, cellxmax, cellymax, vecAx, vecAy, bdcon, Rd, specific_heat_ratio, nval)
                
                Qcon     = base_to_conservative(Qbasem, Qcon, cellxmax, cellymax, specific_heat_ratio)
                Qcon_hat = setup_Qcon_hat(Qcon, Qcon_hat, cellxmax, cellymax, volume, nval)

                # initial_setup
                mu     = set_mu(mu, Qbasem, cellxmax, cellymax, specific_heat_ratio, Rd)
                lambda = set_lambda(lambda, Qbasem, cellxmax, cellymax, mu, specific_heat_ratio, Rd)
                dtau   = set_lts(dtau, lambda_facex, lambda_facey, Qbase, cellxmax, cellymax, mu, dx, dy,
                                vecAx, vecAy, volume, specific_heat_ratio, cfl)
                                
                # RHS
                #advection_term
                E_adv_hat, F_adv_hat = AUSM_plus(E_adv_hat, F_adv_hat, Qbasem, Qcon, cellxmax, cellymax, vecAx, vecAy, specific_heat_ratio, volume, nval)

                # viscos_term
                E_vis_hat, F_vis_hat = central_diff(E_vis_hat, F_vis_hat, Qbasem, Qcon, cellxmax, cellymax, mu, lambda,
                                                vecAx, vecAy, specific_heat_ratio, volume, Rd, nval)
                
                RHS = setup_RHS(RHS, cellxmax, cellymax, E_adv_hat, F_adv_hat, E_vis_hat, F_vis_hat, nval, volume)
            
                # lusgs_advection_term
                A_adv_hat_p, A_adv_hat_m, B_adv_hat_p, B_adv_hat_m, A_beta_shig, B_beta_shig = one_wave(A_adv_hat_m, A_adv_hat_p, B_adv_hat_m, B_adv_hat_p, A_beta_shig, B_beta_shig,
                                                                                                    Qbase, Qcon, cellxmax, cellymax, vecAx, vecAy, specific_heat_ratio, volume, nval)
                # lusgs_viscos_term
                jalphaP, jbetaP = central_diff_jacobian(jalphaP, jbetaP, Qbase, Qcon, cellxmax, cellymax, mu, lambda,
                                                        vecAx, vecAy, specific_heat_ratio, volume, nval)
                
                # LUSGS
                ite = 0
                while true
                    for l in 1:nval
                        for j in 1:cellymax
                            for i in 1:cellxmax
                                delta_Q_temp[i,j,l] = delta_Q[i,j,l]
                            end
                        end
                    end
                    
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
            for l in 1:nval
                for j in 1:cellymax
                    for i in 1:cellxmax
                        Qbase[i,j,l] = Qbasem[i,j,l]
                    end
                end
            end

            if round(evalnum) % every_outnum == 0
                println("\n")
                println("nt_______________________________"*string(round(evalnum)))
                output_result(evalnum, Qbase, cellxmax, cellymax, specific_heat_ratio, out_file_front, out_ext, out_dir, Rd, nval)
            end
    
            check_divrege(Qbase, cellxmax, cellymax, Rd, fwrite)
        end
    end
        
    end_t = now()
    output_fin(fwrite, start_t, end_t, nt, dt, in_nt, cellxmax, cellymax)
end

# -- main --
main()


