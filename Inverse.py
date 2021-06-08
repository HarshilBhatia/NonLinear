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


