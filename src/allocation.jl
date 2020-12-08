function common_allocation(cellxmax, cellymax, nval)
    
    Qbase = zeros(cellxmax, cellymax, nval)
    volume = zeros(cellxmax, cellymax)
    
    dx = zeros(cellxmax+1, cellymax)
    dy = zeros(cellxmax, cellymax+1)

    Qcon = zeros(cellxmax, cellymax, nval)
    Qcon_hat = zeros(cellxmax, cellymax, nval)

    return Qbase, volume, dx, dy, Qcon, Qcon_hat
end 