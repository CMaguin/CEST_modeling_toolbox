function Z_corrected=correct_B0(Z,xZspec, dw)
    
    Z_corrected=spline(xZspec-dw, Z, xZspec);