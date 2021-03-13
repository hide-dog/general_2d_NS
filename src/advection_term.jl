# ------------------------------------
# advection term by AUSM
# ------------------------------------
function AUSM_plus(Qbase, Qcon, cellxmax, cellymax, vecAx, vecAy, specific_heat_ratio, volume, nval)
    E_adv_hat = zeros(cellxmax+1,   cellymax, nval)
    F_adv_hat = zeros(  cellxmax, cellymax+1, nval)
    g         = specific_heat_ratio
    temp_vecX = zeros(nval)   # vecAx used for the pressure term
    temp_vecY = zeros(nval)   # vecAx used for the pressure term
    Lpsi = zeros(nval)        
    Rpsi = zeros(nval)

    mdot = 0.0
    ph = 0.0
    
    tvAxx = 0.0
    tvAxy = 0.0
    vAxx = 0.0
    vAxy = 0.0
    tvAyx = 0.0
    tvAyy = 0.0
    vAyx = 0.0
    vAyy = 0.0
    for j in 2:cellymax -1
        for i in 2:cellxmax+1 -1
            volume_av = 0.5*(volume[i-1,j] + volume[i,j])

            # i-1 cell
            # 規格化
            tvAxx = 0.5*(vecAx[i-1,j,1]+vecAx[i,j,1])
            tvAxy = 0.5*(vecAx[i-1,j,2]+vecAx[i,j,2])
            vAxx = tvAxx / (tvAxx^2+tvAxy^2)^0.5
            vAxy = tvAxy / (tvAxx^2+tvAxy^2)^0.5
            
            rhoL = Qbase[i-1,j,1]
            UL = Qbase[i-1,j,2]*vAxx + Qbase[i-1,j,3]*vAxy
            pL = Qbase[i-1,j,4]
            
            # i cell
            # 規格化
            tvAxx = 0.5*(vecAx[i,j,1]+vecAx[i+1,j,1])
            tvAxy = 0.5*(vecAx[i,j,2]+vecAx[i+1,j,2])
            vAxx = tvAxx / (tvAxx^2+tvAxy^2)^0.5
            vAxy = tvAxy / (tvAxx^2+tvAxy^2)^0.5
            
            rhoR = Qbase[i,j,1]
            UR = Qbase[i,j,2]*vAxx + Qbase[i,j,3]*vAxy
            pR = Qbase[i,j,4]
        
            mdot, ph = AUSM_plus_half(rhoL, rhoR, UL, UR, pL, pR, g,i,j)
            
            # flux half
            sqAx = (vecAx[i,j,1]^2 + vecAx[i,j,2]^2)^0.5
            temp_vecX[2] = vecAx[i,j,1] / sqAx
            temp_vecX[3] = vecAx[i,j,2] / sqAx

            for l in 1:nval
                Lpsi[l] = Qcon[i-1,j,l] / Qcon[i-1,j,1]
                Rpsi[l] = Qcon[i,j,l]   / Qcon[i,j,1]
            end
            Lpsi[4] = (Qcon[i-1,j,4] + Qbase[i-1,j,4]) / Qcon[i-1,j,1]
            Rpsi[4] = (Qcon[i,j,4]   + Qbase[i,j,4])   / Qcon[i,j,1]

            if mdot > 0
                for l in 1:nval
                    E_adv_hat[i,j,l] = (mdot * Lpsi[l] + ph * temp_vecX[l]) * sqAx
                end
            else
                for l in 1:nval
                    E_adv_hat[i,j,l] = (mdot * Rpsi[l] + ph * temp_vecX[l]) * sqAx
                end
            end            
        end
    end

    for j in 2:cellymax+1 -1
        for i in 2:cellxmax -1
            volume_av = 0.5*(volume[i,j-1] + volume[i,j])

            # j-1 cell
            # 規格化
            tvAyx = 0.5*(vecAy[i,j-1,1]+vecAy[i,j,1])
            tvAyy = 0.5*(vecAy[i,j-1,2]+vecAy[i,j,2])
            vAyx = tvAyx / (tvAyx^2+tvAyy^2)^0.5
            vAyy = tvAyy / (tvAyx^2+tvAyy^2)^0.5
            
            rhoL = Qbase[i,j-1,1]
            VL = Qbase[i,j-1,2]*vAyx + Qbase[i,j-1,3]*vAyy
            pL = Qbase[i,j-1,4]
            
            # j cell
            tvAyx = 0.5*(vecAy[i,j,1]+vecAy[i,j+1,1])
            tvAyy = 0.5*(vecAy[i,j,2]+vecAy[i,j+1,2])
            vAyx = tvAyx / (tvAyx^2+tvAyy^2)^0.5
            vAyy = tvAyy / (tvAyx^2+tvAyy^2)^0.5
            
            rhoR = Qbase[i,j,1]
            VR = Qbase[i,j,2]*vAyx + Qbase[i,j,3]*vAyy
            pR = Qbase[i,j,4]

            # scheme
            mdot, ph = AUSM_plus_half(rhoL, rhoR, VL, VR, pL, pR, g,i,j)
            
            # flux half
            sqAy = (vecAy[i,j,1]^2 + vecAy[i,j,2]^2)^0.5
            temp_vecY[2] = vecAy[i,j,1] / sqAy
            temp_vecY[3] = vecAy[i,j,2] / sqAy

            for l in 1:nval
                Lpsi[l] = Qcon[i,j-1,l] / Qcon[i,j-1,1]
                Rpsi[l] = Qcon[i,j,l]   / Qcon[i,j,1]
            end
            Lpsi[4] = (Qcon[i,j-1,4] + Qbase[i,j-1,4]) / Qcon[i,j-1,1]
            Rpsi[4] = (Qcon[i,j,4]   + Qbase[i,j,4])   / Qcon[i,j,1]
            
            if mdot > 0
                for l in 1:nval
                    F_adv_hat[i,j,l] = (mdot * Lpsi[l] + ph * temp_vecY[l]) * sqAy
                end
            else
                for l in 1:nval
                    F_adv_hat[i,j,l] = (mdot * Rpsi[l] + ph * temp_vecY[l]) * sqAy
                end
            end
        end
    end
    
    return E_adv_hat, F_adv_hat
end

function AUSM_plus_half(rhoL, rhoR, UL, UR, pL, pR, g,i,j)
    # param
    beta  = 1/8
    alpha = 3/16
    
    # L, R
    aL = (g * pL / rhoL)^0.5
    aR = (g * pR / rhoR)^0.5

    # half
    ah   = 0.5*( aL + aR )
    rhoh = 0.5*( rhoL + rhoR )

    ML = UL/ah
    MR = UR/ah

    M_p4 = 0
    M_m4 = 0
    p_p5 = 0
    p_m5 = 0
    if abs(ML) >= 1
        M_p4 = 0.5*(ML + abs(ML))
        p_p5 = 0.5*(ML + abs(ML)) / ML
    else
        M_p4 = 0.25*(ML+1)^2 + beta*(ML^2-1)^2
        p_p5 = 0.25*(ML+1)^2 * (2-ML) + alpha*ML*(ML^2-1)^2
    end
    
    if abs(MR) >= 1
        M_m4 = 0.5*(MR - abs(MR))
        p_m5 = 0.5*(MR - abs(MR)) / MR
    else
        M_m4 = -0.25*(MR-1)^2 - beta*(MR^2-1)^2
        p_m5 = 0.25*(MR-1)^2 * (2+MR) - alpha*MR*(MR^2-1)^2
    end
    
    Mh = M_p4 + M_m4
    
    # mdot half
    mdot = ah * Mh
    if Mh > 0
        mdot = mdot * rhoL
    else
        mdot = mdot * rhoR
    end
    
    # p half
    ph = p_p5*pL + p_m5*pR

    return mdot, ph
end
