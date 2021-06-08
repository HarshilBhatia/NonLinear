import numpy as np


def w_next(u, v, w, eta):
    """
    Function for evolution of w

    Parameters:
    u (Double)
    v (Double)
    w (Double)
    eta (Double)

    Returns :
    Next iteration of w
    """
    if(v == 0):
        return eta*(u**2) / (1+w**2) + (1-eta)*(u**2)/(1-w**2)
        
    x1 = (eta * (u ** 2 - v ** 2)) / (1 + w **2)
    if eta == 1:
        return x1
    y1 = (u ** 2 + v ** 2) * (1 - eta) / (1 - w ** 2)

    
    return x1 + y1


def v_next(u, v, w, eta):
    """
    Function for evolution of w

    Parameters:
    u (Double)
    v (Double)
    w (Double)
    eta (Double)

    Returns :
    Next iteration of v
    """
    if( v == 0):
        return 0
    return (-2 * u * v * eta) / (1 + w ** 2)


def u_next(u, v, w, eta):
    """
    Function for evolution of w

    Parameters:
    u (Double)
    v (Double)
    w (Double)
    eta (Double)

    Returns :
    Next iteration of u
    """
    return 2 * w * eta / (1 + w ** 2)


def forward(u, v, w, eta, itr=1, print_itr=False):
    """
    Function for evolution of u,v,w

    Parameters:
    u (Double)
    v (Double)
    w (Double)
    eta (Double)
    itr (int) : Number of iterations, Default : 1
    print_itr (bool) : Prints each iteration, Default : False

    Returns :
    Tuple containing next iteration of u,v,w
    """

    for _ in range(itr):
        U, V, W = u_next(u, v, w, eta), v_next(u, v, w, eta), w_next(u, v, w, eta)
        u, v, w = U, V, W
    return u, v, w

