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
    if v == 0:
        return eta * (u ** 2) / (1 + w ** 2) + (1 - eta) * (u ** 2) / (1 - w ** 2)

    x1 = (eta * (u ** 2 - v ** 2)) / (1 + w ** 2)
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
    if v == 0:
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




def w_inv(Uo, Vo, Wo, eta):
    """
    Function that calculates the inverse of W

    Parameters :
    Uo (double)
    Vo (double)
    Wo (double)
    eta (double)

    Returns:
    Double : value of w
    """
    return (eta - ((eta ** 2 - Uo ** 2) ** (0.5))) / (Uo)


def u_inv(Uo, Vo, Wo, eta, positive=True, invariant=False):
    """
    Function that calculates the inverse of U

    Parameters :
    Uo (double)
    Vo (double)
    Wo (double)
    eta (double)
    positive (bool) : True for the positive inverse map, otherwise false
    Default : True

    invariant (bool) : True for the invariant plane, otherwise false
    Default : False

    Returns:
    Double : value of u
    """

    w = w_inv(Uo, Vo, Wo, eta)
    if invariant:
        q1 = eta / (1 + w ** 2) + (1 - eta) / (1 - w ** 2)
        if positive:
            return +((Wo / q1) ** 0.5)
        return -((Wo / q1) ** 0.5)

    x1 = (1 - eta) / eta
    x2 = (1 + w ** 2) / (1 - w ** 2)
    alp = Vo ** 2 * (x1 ** 2 * x2 ** 2 - 1)
    a = eta / (1 + w ** 2) + (1 - eta) / (1 - w ** 2)
    x = (Wo + (Wo ** 2 - alp) ** 0.5) / (2 * a)

    x = x ** 0.5

    if positive:
        return x
    return -x


def v_inv(Uo, Vo, Wo, eta, positive=True):
    """
    Function that calculates the inverse of V

    Parameters :
    Uo (double)
    Vo (double)
    Wo (double)
    eta (double)
    positive (bool) : True for the positive inverse map, otherwise false
    Default : True


    Returns:
    Double : value of v
    """

    w = w_inv(Uo, Vo, Wo, eta)
    u = u_inv(Uo, Vo, Wo, eta, positive)
    return (Vo * (1 + w ** 2)) / (-2 * eta * u)


def inverse(Uo, Vo, Wo, eta, itr=1, positive=True):
    """
    Function that calculates the inverse of V

    Parameters :
    Uo (Double)
    Vo (Double)
    Wo (Double)
    eta (Double)
    itr (Int): number of iteration of the map
    Default : 1

    positive (bool) : True for the positive inverse map, otherwise false
    Default : True

    Returns:
    Double : Tuple contained iterated values of (u,v,w)
    """
    for _ in range(itr):
        w = w_inv(Uo, Vo, Wo, eta)
        u = u_inv(Uo, Vo, Wo, eta, positive)
        v = v_inv(Uo, Vo, Wo, eta, positive)
        Uo = u
        Vo = v
        Wo = w
    return Uo, Vo, Wo
