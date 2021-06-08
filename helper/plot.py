
def Proj(inv, xx, yy, P, eta):
    """
    Calculates contour values for stereographic projection 

    Parameters:
    inv (function) : Function responsible for returning value of a point
    xx, yy (meshgrid) : Contour meshgrids
    P (Double) : Purity of the stereographic projection
    eta (Double)

    Returns:
    ND Array containing values of Z
    """
    Z = []
    for row in range(len(xx)):
        Z_row = []
        for col in range(len(xx[row])):
            u = (
                2
                * ((2 * P - 1) ** 0.5)
                * xx[row][col]
                / (xx[row][col] ** 2 + yy[row][col] ** 2 + 1)
            )
            v = (
                2
                * ((2 * P - 1) ** 0.5)
                * yy[row][col]
                / (xx[row][col] ** 2 + yy[row][col] ** 2 + 1)
            )
            w = (
                2 * ((2 * P - 1) ** 0.5) / (xx[row][col] ** 2 + yy[row][col] ** 2 + 1)
                - (2 * P - 1) ** 0.5
            )
            Z_row.append(inv(u, v, w, eta))
        Z.append(Z_row)
    return Z
