
function base_to_conservative(Qbase, cellxmax, cellymax, specific_heat_ratio, nval)
    """
    Qbase=[rho,u,v,p]
    Qcon=[rho,rhou,rhov,e]
    """    
    Qcon = zeros(cellxmax, cellymax, nval)
    
    for j in 1:cellymax
        for i in 1:cellxmax
            Qcon[i,j,1] = Qbase[i,j,1]
            Qcon[i,j,2] = Qbase[i,j,1]*Qbase[i,j,2]
            Qcon[i,j,3] = Qbase[i,j,1]*Qbase[i,j,3]
            Qcon[i,j,4] = Qbase[i,j,4]/(specific_heat_ratio-1)+Qbase[i,j,1]*(Qbase[i,j,2]^2+Qbase[i,j,3]^2)/2
        end
    end
    return Qcon
end

function conservative_to_base(Qcon, cellxmax, cellymax, specific_heat_ratio, nval)
    """
    Qbase=[rho,u,v,p]
    Qcon=[rho,rhou,rhov,e]
    """
    Qbase = zeros(cellxmax, cellymax, nval)
    for j in 1:cellymax
        for i in 1:cellxmax
            Qbase[i,j,1] = Qcon[i,j,1]
            Qbase[i,j,2] = Qcon[i,j,2]/Qcon[i,j,1]
            Qbase[i,j,3] = Qcon[i,j,3]/Qcon[i,j,1]
            Qbase[i,j,4] = (Qcon[i,j,4]-Qcon[i,j,1]*(Qbase[i,j,2]^2+Qbase[i,j,3]^2)/2)*(specific_heat_ratio-1)
        end
    end
    return Qbase
end

function setup_Qcon_hat(Qcon, cellxmax, cellymax, volume, nval)
    Qcon_hat = zeros(cellxmax, cellymax, nval)
    for l in 1:nval
        for j in 1:cellymax
            for i in 1:cellxmax
                Qcon_hat[i,j,l] = Qcon[i,j,l] * volume[i,j]
            end
        end
    end
    return Qcon_hat
end

function Qhat_to_Q(Qcon_hat, cellxmax, cellymax, volume, nval)
    Qcon = zeros(cellxmax, cellymax, nval)
    for l in 1:nval
        for j in 1:cellymax
            for i in 1:cellxmax
                Qcon[i,j,l] = Qcon_hat[i,j,l] / volume[i,j]
            end
        end
    end
    return Qcon
end